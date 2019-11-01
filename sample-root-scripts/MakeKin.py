#!/usr/bin/python

from optparse import OptionParser
import random
from math import pi,sin,cos,sqrt
import sys


pid = {"pi0":111, "pi+":211, "k0l":130, "k0s":310, "k+":321,
       "e+":-11, "mu+":-13, "tau+":-15, 
       "nue":12, "nuebar":-12, 
       "numu":14, "numubar":-14, 
       "nutau":16, "nutaubar":-16,
       "p+":2212, "n0":2112,
       "Bi214":1000832140,"Tl208":1000812080,"K40":1000190400}

#holds detector [radius, height] in cm
detectors = {"SuperK":[3368.15/2., 3620.],
             "Cylinder_60x74_20inchBandL_14perCent":[7400./2., 6000.],
             "Cylinder_60x74_20inchBandL_40perCent":[7400./2., 6000.],
             "HyperK":[7080./2., 5480.]}

for pname, no in pid.items():
    if pname.endswith('+'):
        pid[pname.replace('+','-')] = -1*pid[pname]





parser = OptionParser()
optdefault = 1
parser.add_option("-N", "--nfiles", dest="nfiles",
                  help="number of files of particles to produce. Default: %s" \
                      % (optdefault),
                  metavar="#", default=optdefault)
optdefault = 100
parser.add_option("-n", "--npart", dest="npart",
                  help="number of particles to simulate per file. Default: %s" \
                  % (optdefault),
                  metavar="#", default=optdefault)
optdefault = None
parser.add_option("-s", "--seed", dest="seed",
                  help="Random number seed to use. Default: None (use system time)",
                  metavar="SEED", default=optdefault)
optchoices = pid.keys()
optdefault = "mu-"
parser.add_option("-t", "--type", dest="type",
                  help="Particle type to be generated. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  metavar="TYPE",
                  choices=optchoices, default=optdefault)
optdefault = '1000.0'
parser.add_option("-e", "--energy", dest="energy",
                  help="Particle energy to be generated in MeV. Specify range with colon separated pair (e.g. 0:10). Default: %s" \
                      % (optdefault),
                  metavar="ENERGY", default=optdefault)
optchoices = ["center", "randomwater", "randompmt", "minusx", "plusx", "minusz", "plusz"]
optdefault = optchoices[0]
parser.add_option("-v", "--vertex", dest="vertname",
                  help="Type of vertex. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  choices=optchoices, default=optdefault)
optdefault = "./"
parser.add_option("-p", "--pmtlist", dest="pmtlist",
                  help="file with list of PMT positions. End with a '/' for the filename to be constructed automatically as <detector>_pmts.list. Default: %s" \
                      % (optdefault),
                  metavar="FILE", default=optdefault)
optchoices = ["4pi", "towall", "tocap"]
optdefault = optchoices[0]
parser.add_option("-d", "--direction", dest="dirname",
                  help="Type of direction. Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  choices=optchoices, default=optdefault)
optchoices = detectors.keys()
optdefault = "SuperK"
parser.add_option("-w", "--detector", dest="detector",
                  help="Detector water volume to use (for vertex positioning). Choices: %s. Default: %s" \
                      % (optchoices, optdefault),
                  choices=optchoices, default=optdefault)

(options, args) = parser.parse_args()

options.vertname = options.vertname.lower()
options.dirname = options.dirname.lower()
random.seed(options.seed)



def pair_or_single(arg):
    pair = list(set(float(x) for x in arg.split(':')))
    pair.sort()
    if len(pair) > 2 or len(pair) < 1:
        print "Argument %s is not a colon-separted pair of floats, or a single float" % arg
    elif len(pair) == 1:
        return pair, "%.1f" % (pair[0])
    elif len(pair) == 2:
        return pair, "%.1f:%.1f" % (pair[0], pair[1])

nfiles = int(options.nfiles)
npart  = int(options.npart)
energy, energystr = pair_or_single(options.energy)


#Define the particle
particle = {"vertex":(0, 0, 0),
            "time":0,
            "type":pid[options.type],
            "energy":energy[0],
            "direction":(1,0,0)}

def ReadPMTPositions(filename):
    pmts = {}
    with open(filename) as f:
        for l in f.readlines():
            pmtnum,x,y,z = l.split()
            pmts[int(pmtnum)] = [float(x), float(y), float(z)]
    return pmts

randvertwater = False
if options.vertname == "center":
    randvertwater = False
elif options.vertname == "randomwater":
    randvertwater = True
elif options.vertname == "randompmt":
    fname = options.pmtlist
    if options.pmtlist.strip()[-1] == '/':
        fname += options.detector + '_pmts.list'
    pmts = ReadPMTPositions(fname)
    randvertpmt = True
elif options.vertname == "wall":
    print >>sys.stderr, "Wall not implemented yet"
    sys.exit(3)
elif options.vertname == "minusx":
    if options.detector == "SuperK":
        particle["vertex"] = (-1000,0,0)
    else:
        particle["vertex"] = (-1000. * detectors[options.detector][0] / detectors["SuperK"][0], 0, 0)
elif options.vertname == "plusx":
    if options.detector == "SuperK":
        particle["vertex"] = (+1000,0,0)
    else:
        particle["vertex"] = (+1000. * detectors[options.detector][0] / detectors["SuperK"][0], 0, 0)
elif options.vertname == "minusz":
    if options.detector == "SuperK":
        particle["vertex"] = (0, 0, -1000)
    else:
        particle["vertex"] = (0, 0, -1000. * detectors[options.detector][1] / detectors["SuperK"][1])
elif options.vertname == "plusz":
    if options.detector == "SuperK":
        particle["vertex"] = (0, 0, +1000)
    else:
        particle["vertex"] = (0, 0, +1000. * detectors[options.detector][1] / detectors["SuperK"][1])
else:
    print >>sys.stderr, "Don't understand vertex",opttions.vertname
    sys.exit(2)

if options.dirname == "towall":
    randdir = False
    particle["direction"] = (1,0,0)
elif options.dirname == "tocap":
    randdir = False
    particle["direction"] = (0,0,1)
elif options.dirname == "4pi":
    randdir = True
elif options.dirname == "wall":
    print >>sys.stderr, "Wall not implemented yet"
    sys.exit(3)
else:
    print >>sys.stderr, "Don't understand direction",options.dirname
    sys.exit(2)





nu =   {"type":pid["numu"], "energy":energy[0]+1000.0,
        "direction":(1, 0, 0)}
prot = {"type":pid["p+"], "energy":935.9840,
        "direction":(0, 0, 1)}



def partPrint(p, f, recno):
    f.write("$ begin\n")
    f.write("$ nuance 0\n")
    #random energy
    if len(energy) == 2:
        thisenergy  = random.uniform(energy[0], energy[1])
        p ["energy"] = thisenergy
        nu["energy"] = thisenergy + 1000
    if randvertwater:
        rad    = detectors[options.detector][0] - 20. #cm
        height = detectors[options.detector][1] - 20. #cm
        while True:
            x = random.uniform(-rad,rad)
            y = random.uniform(-rad,rad)
            if x**2 + y**2 < rad**2: break
        z = random.uniform(-height/2,height/2)
        f.write("$ vertex %.5f %.5f %.5f %.5f\n" % (x, y, z, p["time"]))
    elif randvertpmt:
        pmt = random.randint(1, len(pmts))
        f.write("$ vertex %.5f %.5f %.5f %.5f\n" % (pmts[pmt][0], pmts[pmt][1], pmts[pmt][2], p["time"]))
    else:
        f.write("$ vertex %.5f %.5f %.5f %.5f\n" % (p["vertex"]+(p["time"],)) )
    printTrack(nu, f, -1)   # "Neutrino" Track
    printTrack(prot, f, -1) # "Target" track
    f.write("$ info 0 0 %i\n" % recno)
    if randdir:
        th = random.random()*2*pi
        u = 1.-2*random.random()
        x = sqrt(1.-u**2)*cos(th)
        y = sqrt(1.-u**2)*sin(th)
        z = u
        p["direction"] = (x, y, z)
        #th = random.random()*pi
        #phi = random.random()*2*pi
        #p["direction"] = (cos(phi)*cos(th), sin(phi)*cos(th), sin(th))
        
    printTrack(p, f)    # Outgoing Particle Track
    f.write("$ end\n")

def printTrack(p, f, code=0):
    f.write("$ track %(type)i %(energy).5f " % p)
    f.write("%.5f %.5f %.5f %i\n" % (p["direction"]+(code,)))


for fileno in range(nfiles):
    typestr = options.type.replace("+","plus").replace("-","minus")
    
    estring = ""
    filename="%s_%sMeV_%s_%s_%s_%03i.kin" % (typestr, energystr, options.vertname, options.dirname, options.detector, fileno)

    outfile = open(filename, 'w')

    print ("Writing %i particles to " % npart) + filename

    for i in range(npart):
        partPrint(particle, outfile, i)
    outfile.write("$ stop")
    outfile.close()
