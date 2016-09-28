from ROOT import gROOT, gStyle, gSystem
from ROOT import TFile, TH1F, TTree, TCanvas
import os.path
import sys

# Simple example of reading a generated Root file

def sample_readfile(filename = "../wcsim.root", verbose=False):
    #Clear global scope
    #gROOT.Reset()
    
    #Load the library with class dictionary info
    #(create with "gmake shared")
    gSystem.Load(os.path.expandvars("${WCSIMDIR}/libWCSimRoot.so"))
    from ROOT import WCSimRootEvent, WCSimRootGeom, WCSimRootTrigger

    #Open the file
    f = TFile(filename, "read");
    if not f.IsOpen():
        print "Error, could not open input file:", filename
        sys.exit(-1)
  
    #Get the a pointer to the tree from the file
    tree = f.wcsimT
  
    #Get the number of events
    nevent = tree.GetEntries();
    if verbose:
        print "nevent", nevent
  
    #Geometry tree
    geotree = f.wcsimGeoT
    if verbose:
        print "Geotree has", geotree.GetEntries(), "entries"
    geotree.GetEntry(0)
    geo = geotree.wcsimrootgeom
    if not geotree.GetEntries() == 1:
        print "Geotree should have exactly 1 entry"
        sys.exit(9)

    h1    = TH1F("PMT Hits", "PMT Hits", 8000, 0, 8000);
    hvtx0 = TH1F("Event VTX0", "Event VTX0", 200, -1500, 1500);
    hvtx1 = TH1F("Event VTX1", "Event VTX1", 200, -1500, 1500);
    hvtx2 = TH1F("Event VTX2", "Event VTX2", 200, -1500, 1500);
  
    num_trig = 0;
  
    print nevent
    #Now loop over events
    for iev, event in enumerate(tree):
        print "Tree entry", iev

        wcsimrootsuperevent = tree.wcsimrootevent
        wcsimrootevent = wcsimrootsuperevent.GetTrigger(0)

        #get the basic event information
        if verbose:
            print "********************************************************"
            print "Evt, date", wcsimrootevent.GetHeader().GetEvtNum(), wcsimrootevent.GetHeader().GetDate()
            print "Mode", wcsimrootevent.GetMode()
            print "Number of subevents", wcsimrootsuperevent.GetNumberOfSubEvents()
            
            print "Vtxvol", wcsimrootevent.GetVtxvol()
            print "Vtx", wcsimrootevent.GetVtx(0), wcsimrootevent.GetVtx(1),wcsimrootevent.GetVtx(2)
            
        hvtx0.Fill(wcsimrootevent.GetVtx(0))
        hvtx1.Fill(wcsimrootevent.GetVtx(1))
        hvtx2.Fill(wcsimrootevent.GetVtx(2))

        ntracks = wcsimrootevent.GetNtrack()
        if verbose:
            print "Jmu", wcsimrootevent.GetJmu()
            print "Nparn", wcsimrootevent.GetNpar()
            print "Ntrack", ntracks

        #get the track information
        for itrack in xrange(ntracks):
            wcsimroottrack = wcsimrootevent.GetTracks().At(itrack)
            if verbose:
                print "Track " + str(itrack) + " ipnu:", wcsimroottrack.GetIpnu()
                print "Track " + str(itrack) + " parent ID:", wcsimroottrack.GetParenttype()
                for idir in xrange(3):
                    print "Track " + str(itrack) + " dir[" + str(idir) + "]:", wcsimroottrack.GetDir(idir)
        #end loop over tracks

    
        #get the Cherenkov hit (i.e. true) information
    
        #Get the number of Cherenkov hits.
        #Note... this is *NOT* the number of photons that hit tubes.
        #It is the number of tubes hit with Cherenkov photons.
        #The number of digitized tubes will be smaller because of the threshold.
        #Each hit "raw" tube has several photon hits.  The times are recorded.
        #See chapter 5 of ../doc/DetectorDocumentation.pdf
        #for more information on the structure of the root file.
        # 
        #The following code prints out the hit times for the first 10 tubes and also
        #adds up the total pe.
        #
        #For digitized info (one time/charge tube after a trigger) use
        #the digitized information.

        ncherenkovhits     = wcsimrootevent.GetNcherenkovhits();
        ncherenkovdigihits = wcsimrootevent.GetNcherenkovdigihits(); 
    
        h1.Fill(ncherenkovdigihits);
        if verbose:
            print "node id:", iev
            print "Ncherenkovhits:",     ncherenkovhits
            print "Ncherenkovdigihits:", ncherenkovdigihits
            print "RAW HITS:"

        #Grab the big arrays of times and parent IDs
        timeArray = wcsimrootevent.GetCherenkovHitTimes()
    
        totalPe = 0
        for ihittube in xrange(ncherenkovhits):
            wcsimrootcherenkovhit = wcsimrootevent.GetCherenkovHits().At(ihittube)

            tubeNumber     = wcsimrootcherenkovhit.GetTubeID()
            timeArrayIndex = wcsimrootcherenkovhit.GetTotalPe(0)
            peForTube      = wcsimrootcherenkovhit.GetTotalPe(1)
            pmt   = geo.GetPMT(tubeNumber-1)

            totalPe += peForTube

            if ihittube < 10 and verbose: #Only print first XX=10 tubes
                print "Total pe:", peForTube, "times (",
                for ihit in xrange(timeArrayIndex, timeArrayIndex + peForTube):
                    HitTime = timeArray.At(ihit)
                    print HitTime.GetTruetime(), ",",
                print ")"

        #End of loop over Cherenkov hits
        if verbose:
            print "Total Pe: ", totalPe
    

        #Look at digitized hit info

        #Get the number of digitized hits
        #These are stored in subevents
   
        if verbose:
            print "DIGITIZED HITS:"

        for itrigger in xrange(wcsimrootsuperevent.GetNumberOfEvents()):

            wcsimrootevent = wcsimrootsuperevent.GetTrigger(itrigger)
            if verbose:
                print "Sub event number =", itrigger
                print "Ncherenkovdigihits", ncherenkovdigihits
     
            if ncherenkovdigihits > 0:
                num_trig += 1

            #Loop through elements in the TClonesArray of WCSimRootCherenkovDigHits
            for idigihit in xrange(ncherenkovdigihits):
                wcsimrootcherenkovdigihit = wcsimrootevent.GetCherenkovDigiHits().At(idigihit)
                if idigihit < 10 and verbose: #Only print first XX=10 tubes
                    print "q, t, tubeid:",wcsimrootcherenkovdigihit.GetQ(), wcsimrootcherenkovdigihit.GetT(), wcsimrootcherenkovdigihit.GetTubeId()
            #End of loop over Cherenkov digihits
        #End of loop over trigger
    
        #reinitialize super event between loops.
        wcsimrootsuperevent.ReInitialize()
    
    #End of loop over events

    #do the drawing
    win_scale = 0.75
    n_wide = 2
    n_high = 2
    c1 = TCanvas("c1", "First canvas", int(500 * n_wide * win_scale), int(500 * n_high * win_scale))
    c1.Draw()
    c1.Divide(2,2)
    c1.cd(1)
    hvtx0.Draw()
    c1.cd(2)
    hvtx1.Draw()
    c1.cd(3)
    hvtx2.Draw()
    c1.cd(4)
    h1.Draw()
  
    print "num_trig", num_trig
#End of sample_readfile()

if __name__ == "__main__":
    sample_readfile(verbose=True)

