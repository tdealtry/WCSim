#!/bin/env python3

import argparse
import sys
from datetime import datetime
import enum

ns_conversion = {'ns':1,
                 'us':1E3,
                 'ms':1E6,
                 's':1E9}

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='KinConverter: convert kin files that have multiple vertices'
                                 ' into a single event (or multiple overlapping events)')
parser.add_argument('--input-filename', required=True, default=argparse.SUPPRESS, type=str,
                    help='Input .kin filename. Output filename(s) will be the same as the input filename & path'
                    ', with [0-9].merge suffix(es) added')
parser.add_argument('--input-time-unit', required=True, default=argparse.SUPPRESS,
                    choices=ns_conversion.keys(), help='The time unit of the input file')
parser.add_argument('--dark-noise-start', required=True, default=argparse.SUPPRESS, type=int,
                    help='When to start the simulation (in ns)')
parser.add_argument('--dark-noise-end', required=True, default=argparse.SUPPRESS, type=int,
                    help='When to end the simulation (in ns)')
parser.add_argument('--event-overlap', required=True, default=argparse.SUPPRESS, type=int,
                    help='How long (in ns) to overlap')
parser.add_argument('--verbose', '--v', type=int, default=0, help='Verbosity level')

subparsers = parser.add_subparsers(title="Actions", dest='mode',
                                   help='The split strategy')
subparsers.required = True
#Use either this
parser_fix = subparsers.add_parser("FIX", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                   help='Use a fixed event duration')
parser.add_argument('--fixed-duration', required=True, default=argparse.SUPPRESS, type=int,
                    help='A fixed duration (in ns) for each event')
#or these
parser_free = subparsers.add_parser("FREE", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    help='Use a variable event duration, based on the amount of hits')
parser_fix.add_argument('--dark-rate-kHz', required=True, default=argparse.SUPPRESS, type=float,
                         help='Dark rate of the PMTs in kHz')
parser_fix.add_argument('--ntubes', required=True, default=argparse.SUPPRESS, type=int,
                         help='Number of PMTs')
parser_fix.add_argument('--nhits-per-MeV', required=True, default=argparse.SUPPRESS, type=float,
                         help='Approximate number of hits per MeV')
parser_fix.add_argument('--max-hits-allowed', required=True, default=argparse.SUPPRESS, type=int,
                         help='Maximum number of hits that can be used fed into WCSim.'
                         'This is generally a memory limit')
parser_fix.add_argument('--step-duration', required=True, default=argparse.SUPPRESS, type=int,
                        help='The duration (in ns) to step through until the --max-hits-allowed is reached')

class Mode(enum.Enum):
    FIX  = 1
    FREE = 2

args = parser.parse_args()
mode = Mode[args.mode]
ToNS = ns_conversion[args.input_time_unit]

def PrintNS(time):
    for x in ['ns', 'us', 'ms', 's']:
        if time < 1000.0:
            return "%f %s" % (time, x)
        time /= 1000.0

def PrintB(mem):
    for x in ['bytes', 'kB', 'MB', 'GB', 'TB']:
        if mem < 1000.0:
            return "%.2f %s" % (mem, x)
        mem /= 1000.0

def PrintMeV(energy):
    for x in ['MeV', 'GeV', 'TeV', 'PeV']:
        if energy < 1000.0:
            return "%.2f %s" % (energy, x)
        energy /= 1000.0

DummyVertex = """$ nuance 0
$ vertex 0 0 0 0
$ track -12 0.00000 0.00000 0.00000 1.00000 -1
$ track 2212 938.27231 0.00000 0.00000 1.00000 -1
$ info 0 0 0
$ track -11 0.511 0 0 0 0
"""

#read a vertex at a time
def GetVertex(seq, group_by, exclude=['']):
    data = []
    for line in seq:
        if isinstance(line, bytes):
            line = line.decode()
        if line.startswith(group_by):
            if data:
                yield data
                data = []
        if line.strip() in exclude:
            continue
        data.append(line)

    if data:
        yield data

#get the time from a vertex
def GetTime(vertex):
    for line in vertex[1:2]:
        if 'vertex' in line:
            return float(line.split()[-1]) * ToNS

#get the energy from a vertex
def GetEnergy(vertex):
    vertex_energy = 0
    #we don't care about the nu/target lines, so slice list
    # until we get past "info"
    for line in vertex[5:]:
        if line.startswith('$ track'):
            s = line.split()
            p = abs(int(s[2]))
            #check for electron/positron
            if p == 11:
                vertex_energy += float(s[3])
            #check for neutron - gives 2.2 MeV gamma
            elif p == 2112:
                vertex_energy += 2.2
            else:
                print(f'Unknown particle {p}')
                sys.exit(1)
    return vertex_energy

#check that the file is time ordered
def IsTimeOrdered(filename):
    print('Checking input file is time ordered...')
    return True
    with open(filename, 'r') as fin:
        for i, vertex in enumerate(GetVertex(fin, "$ begin")):
            #skip the header
            if vertex[0].startswith('#'):
                continue
            try:
                last_time = this_time
            except UnboundLocalError:
                last_time = float(GetTime(vertex))
            this_time = float(GetTime(vertex))
            if this_time < last_time:
                print(PrintNS(this_time), "comes before", PrintNS(last_time))
                return False
    print("Is time ordered")
    return True

#Sort the initial file by time
def SortByTime(filename):
    outfilename = args.input_filename + '.temp'
    with open(filename, 'r') as fin, open(outfilename, 'w') as fout:
        print('TODO SortByTime() not yet implemented')
        sys.exit(-1)

#Get the header
def GetHeader(filename, args):
    with open(filename, 'r') as fin:
        for vertex in GetVertex(fin, "$ begin"):
            #if there's no header, return blank
            if vertex[0].startswith("$ begin"):
                return ''
            #return the header as a single string
            header = ''.join(vertex) + ''\
                     '# Split by kin_converter ' + str(datetime.now()) + '\n'\
                     '# --fixed-duration ' + str(args.fixed_duration) + '\n'\
                     '# --event-overlap ' + str(args.event_overlap) + '\n'
            return header

class CombinedEvent():
    event_header = '$ begin\n'
    event_footer = '$ end\n'\
        '$ stop\n'
    def __init__(self, header, ievent):
        self.header = header
        self.dark_header = None
        self.ievent = ievent
        self.energy = 0
        #treating the contents as a list of strings
        # rather than one long string
        # gives a big performance improvement
        self.contents = []

    def Append(self, vertex):
        self.contents.append(''.join(vertex))
        self.energy += GetEnergy(vertex)

    def SetDarkRange(self, event_start, event_end):
        #write the dark noise range
        self.dark_header = f'# Event {self.ievent}\n'\
            f'# /DarkRate/SetDarkLow  {event_start}\n'\
            f'# /DarkRate/SetDarkHigh {event_end}\n'

    def Write(self):
        if not self.dark_header:
            print('You need to call SetDarkRange() before Write()')
            sys.exit(-1)
        #need to add a dummy vertex, else WCSim/Geant4 will complain
        if not len(self.contents):
            self.contents.append(DummyVertex)
        with open(args.input_filename + '.%09d' % self.ievent, 'w') as fout:
            fout.write(self.header)
            fout.write(self.dark_header)
            fout.write(CombinedEvent.event_header)
            for c in self.contents:
                fout.write(c)
            fout.write(CombinedEvent.event_footer)
        print(PrintB(sys.getsizeof(self.contents)))

    def CalcHitStats(self):
        self.nvertices = len(self.contents)

        #Number of hits
        self.ndarkhits = int(args.dark_rate_kHz * 1E3 * args.ntubes * args.fixed_duration * 1E-9)
        self.nphysicshits = int(args.nhits_per_MeV * self.energy)
        self.nhits = self.ndarkhits + self.nphysicshits
        self.dark_percent = 100. * self.ndarkhits / (self.nhits)

        #Size in memory
        #hit = (13 int + 6 double) * 3 (true + digi + triggered)
        # this assumes each true hit makes exactly 1 digi hit (overestimate)
        # and each digi hit is triggered (overestimate)
        size_hits_mem = (13 * 4 + 6 * 8) * self.nhits * 3
        #true hit output = 3 ints per PMT + 1 int & 1 double per hit
        size_truehits_out = 3 * 4 * args.ntubes + (4 + 8) * self.nhits
        #digi hit output = 1 float + 1 double + 2 int per hit
        size_digihits_out = (4 + 8 + 4 * 2) * self.nhits
        self.size_hits = size_hits_mem + size_truehits_out + size_digihits_out

    def __str__(self):
        try:
            x = self.nhits
        except AttributeError:
            self.CalcHitStats()
        s = f'Event {ievent} corresponds to time range {PrintNS(event_start)} {PrintNS(event_end)} ns\n'\
            f'Contains {self.nvertices} vertices\n'\
            f'Total energy in event {PrintMeV(self.energy)}'
        try:
            s += f' ({PrintMeV(self.energy / self.nvertices)} average)\n'
        except ZeroDivisionError:
            s += '\n'
        f'Total hits in event {self.nhits} ({self.dark_percent:.2f}% dark noise)\n'
        if args.verbose:
            s += f'Physics hits in event {self.nphysicshits}\n'\
                f'Dark hits in event {self.ndarkhits}\n'
        s += f'Estimated memory required for hits in WCSim {PrintB(self.nhits * 40 * 2)}\n'\
            f'Estimated memory required for hits in WCSim {PrintB(self.size_hits)}\n\n'
        return s


##########################################################
#See if the file is time ordered
if not IsTimeOrdered(args.input_filename):
    #If not, sort it
    SortByTime(args.input_filename)

header = GetHeader(args.input_filename, args)
print(header)

#loop over kin file start/stop times
event_start = args.dark_noise_start
last_event_end = args.dark_noise_end
ievent = 0
file_position = 0
while event_start < last_event_end:
    if mode == Mode.FIX:
        event_end = event_start + args.fixed_duration
        next_event_start = event_start + args.fixed_duration - args.event_overlap

    with open(args.input_filename, 'rb') as fin:
        #create the output object
        ce = CombinedEvent(header, ievent)

        #skip forward in the file a bit
        if args.verbose:
            print('Skipping to position in file', file_position)
        fin.seek(file_position)

        #loop over the input file
        for i, vertex in enumerate(GetVertex(fin, "$ begin", ["$ begin", "$ end"])):
            #skip the header and any partial vertices we've found from using seek()
            if not vertex[0].startswith('$ nuance'):
                continue
            #get the event time
            time = GetTime(vertex)
            if args.verbose > 1:
                print(PrintNS(time))
                if args.verbose > 2:
                    print(f"Vertex #{i}")
                    print("".join(vertex))

            if mode == Mode.FIX:
                if time > event_end:
                    break
                if time >= event_start:
                    ce.Append(vertex)
                #save the current file position if it is earlier than required for the next event
                # At most, the position will be the '$ begin' line of the first vertex in the next event
                if time < next_event_start:
                    file_position = fin.tell()

        ce.SetDarkRange(event_start, event_end)
        ce.Write()
        print(ce)
        del ce
    #increment for next event
    event_start = next_event_start
    ievent += 1
