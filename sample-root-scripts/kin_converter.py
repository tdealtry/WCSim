#!/bin/env python3

import argparse
import sys
from datetime import datetime

#All times internally in this code use ns
#Use this to convert between other units (from command line options, or from input files)
ns_conversion = {'ns':1,
                 'us':1E3,
                 'ms':1E6,
                 's':1E9}

class TimeAndUnit(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            value = int(values[0]) * ns_conversion[values[1]]
        except ValueError:
            print(option_string, 'is in form VALUE UNIT')
            print('Require VALUE to be an integer')
            sys.exit(-1)
        except KeyError:
            print(option_string, 'is in form VALUE UNIT')
            print('Require UNIT to be one of:', ' '.join(ns_conversion.keys()))
            sys.exit(-1)
        else:
            setattr(namespace, self.dest, int(value))


parser = argparse.ArgumentParser(description='KinConverter: convert kin files that have multiple vertices into a single event (or multiple overlapping events)')
parser.add_argument('--input-filename', '-i', required=True, type=str,
                    help='Input .kin filename. Output filename(s) will be the same as the input filename & path, with [0-9].merge suffix(es) added')
parser.add_argument('--input-time-unit', required=True, choices=ns_conversion.keys(),
                    help='The time unit of the input file')
parser.add_argument('--dark-noise-start', action=TimeAndUnit, nargs=2, required=True,
                    help='When to start the simulation (in ns)')
parser.add_argument('--dark-noise-end', action=TimeAndUnit, nargs=2, required=True,
                    help='When to end the simulation (in ns)')
parser.add_argument('--event-overlap', action=TimeAndUnit, nargs=2, required=True,
                    help='How long (in ns) to overlap')
parser.add_argument('--verbose','--v',type=int,default=0,help='Verbosity level')

subparsers = parser.add_subparsers(dest='mode', help='Run mode')#, required=True)
#Use either this
parser_fix = subparsers.add_parser('fixed', help='Use a fixed duration for each output event')
parser_fix.add_argument('--fixed-duration', action=TimeAndUnit, nargs=2, required=True,
                        help='A fixed duration (in ns) for each event')
#or these - dark rate, ntubes, NHits per MeV, max allowed hits
parser_free = subparsers.add_parser('free', help='Use a variable duration for each output event, based on dark rate, number of PMTs, NHits per MeV, and limited by the "max number of allowed hits"')
parser_free.add_argument('--dark-rate', type=float, required=True,
                         help='Dark rate (in kHz) per PMT')
parser_free.add_argument('--nPMTs', type=int, required=True,
                         help='Number of PMTs in simulation')
parser_free.add_argument('--nhits-per-MeV', type=float, required=True,
                         help='Number of hits per MeV that are expected. This is used to give an estimate of the number of "physics" hits in the event')
parser_free.add_argument('--max-hits-allowed', type=float, required=True,
                         help='Maximum (expected) number of hits per out .kin file. This is limited by available memory')
#TODO account for case where there are multiple PMT types in the detector
args = parser.parse_args()

ToNS = ns_conversion[args.input_time_unit]

def PrintNS(time):
    for x in ['ns', 'us', 'ms', 's']:
        if time < 1000.0:
            return "%f %s" % (time, x)
        time /= 1000.0

#Dummy vertex to use when there are no true physics events
# in the given time window
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
    for line in vertex:
        if 'vertex' in line:
            return float(line.split()[-1]) * ToNS

#get the total energy from a vertex
def GetEnergy(vertex):
    total_energy = 0
    for line in vertex:
        #get only particles leaving the nucleus (after FSI)
        if line.startswith('$ track') and line.endswith('0'):
            a, b, pdg, energy, x = line.split(None, 4)
        pdg = int(pdg)
        energy = float(energy)
        if abs(pdg) == 11:
            total_energy += energy
        elif pdg == 2112:
            total_energy += 2.2 #assume neutron capture on H giving 2.2 MeV gamma
        else:
            print('Unknown pdg code', pdg)
            sys.exit(1)
    return total_energy

#check that the file is time ordered
def IsTimeOrdered(filename):
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
    header = ''
    with open(filename, 'r') as fin:
        for vertex in GetVertex(fin, "$ begin"):
            #if there's a header in the input file, use that as the start of the header of the output files
            if not vertex[0].startswith("$ begin"):
                header += ''.join(vertex)
            break
    #return the header as a single string
    header += '# Split by kin_converter ' + str(datetime.now()) + '\n' + str(args) + '\n'
    return header

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
    event_end = event_start + args.fixed_duration
    next_event_start = event_start + args.fixed_duration - args.event_overlap
    print("Event", ievent, "corresponds to range", PrintNS(event_start), PrintNS(event_end))
    with open(args.input_filename, 'rb') as fin, open(args.input_filename + '.%09d' % ievent, 'w') as fout:
        #write the original header
        fout.write(header)
        #write the dark noise range
        fout.write('# Event ' + str(ievent) + '\n')
        fout.write('# /DarkRate/SetDarkLow  ' + str(event_start) + '\n')
        fout.write('# /DarkRate/SetDarkHigh ' + str(event_end) + '\n')
        #and the event start
        fout.write('$ begin\n')
        nvertices = 0
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
                    print("Vertex #{}".format(i))
                    print("".join(vertex))
            if time > event_end:
                break
            if time >= event_start:
                fout.write(''.join(vertex))
                nvertices += 1
            #save the current file position if it is earlier than required for the next event
            # At most, the position will be the '$ begin' line of the first vertex in the next event
            if time < next_event_start:
                file_position = fin.tell()
        #need to add a dummy vertex, else WCSim/Geant4 will complain
        if not nvertices:
            fout.write(DummyVertex)
        #and close the event/file
        fout.write('$ end\n')
        fout.write('$ stop\n')
        print('contains', nvertices, 'vertices')
    #increment for next event
    event_start = next_event_start
    ievent += 1
