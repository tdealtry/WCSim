#include <iostream>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TROOT.h>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimEnumerations.hh"
#endif

using namespace std;

TString create_filename(const char * prefix, TString& filename_string)
{
  //std::cout << "Creating filename from prefix " << prefix << " and filename_string " << filename_string << std::endl;
  TString prefix_string(prefix);
  TString outfilename = prefix_string + filename_string;
  return outfilename;
}

// Simple example of reading a generated Root file
int sample_analysis(char *filename=NULL, bool verbose=false, 
		    long max_nevents = -1, long max_ntriggers = -1, bool write_tree = false,)
{
#if !defined(__MAKECINT__)
  // Load the library with class dictionary info
  // (create with "gmake shared")
  char* wcsimdirenv;
  wcsimdirenv = getenv ("WCSIMDIR");
  if(wcsimdirenv !=  NULL){
    gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
  }else{
    gSystem->Load("../libWCSimRoot.so");
  }
#endif

  // Open the input file
  TFile *file;
  if (filename==NULL){
    file = new TFile("../wcsim.root","read");
  }else{
    file = new TFile(filename,"read");
  }
  if (!file->IsOpen()){
    cout << "Error, could not open input file: " << filename << endl;
    return -1;
  }

  
  //
  // SETUP GEOMETRY TREE FOR READING
  //

  // Geometry tree - only need 1 "event"
  TTree *geotree = (TTree*)file->Get("wcsimGeoT");
  WCSimRootGeom *geo = 0; 
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
  if(geotree->GetEntries() == 0) {
    cout << "geotree not found!" << endl;
    exit(9);
  }
  geotree->GetEntry(0);

  
  //
  // SETUP EVENT TREE FOR READING
  //
  
  // Get the a pointer to the tree from the file
  TTree *tree = (TTree*)file->Get("wcsimT");
  
  // Get the number of events
  long nevents = tree->GetEntries();
  const long nevents_loop = max_nevents > 0 ? TMath::Min(max_nevents, nevents) : nevents;
  printf("nevents - total in file = %d, running over = %d\n",nevents, nevents_loop);
 
  // Create a WCSimRootEvent to put stuff from the tree in
  WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();

  // Set the branch address for reading from the tree
  TBranch *branch = tree->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);

  // Force deletion to prevent memory leak 
  tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

  // start with the main "subevent", as it contains most of the info
  // and always exists.
  WCSimRootTrigger* wcsimrootevent;
  //The 0th trigger contains:
  // - track information
  // - raw hit information
  // - digit information (for trigger 0)
  //Subsequent triggers contain:
  // - digit information (for that trigger) only
  //Therefore if you want to do truth matching on trigger >= 1,
  // you need to keep a copy of trigger 0 to read the track/hit info
  // Store it here
  WCSimRootTrigger* wcsimroottrigger0;

  
  //
  // CREATE OUTPUT FILE(S) (& BOOK HISTOGRAMS, TREES, ETC.)
  //
  
  //create output file for storing histograms and tree
  TString filenameout(filename);
  TFile * fout = new TFile(create_filename(__func__, filenameout).Data(), "RECREATE");
  fout->cd();

  //book 
  TTree * tout = new TTree("digitimes", "digitimes");
  vector<double> tdigitimes;
  tout->Branch("digitimes", &tdigitimes);

  int num_trig = 0;
  
  // Now loop over events
  for (int iev=0; iev<nevents_loop; iev++)
    {
      // Read the event from the tree into the WCSimRootEvent instance
      tree->GetEntry(iev);
      wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);

      //get some of the basic event information
      const int ntriggers = wcsimrootsuperevent->GetNumberOfEvents();
      const int ntriggers_loop = max_ntriggers > 0 ? TMath::Min(max_ntriggers, ntriggers) : ntriggers;
      const int ntracks = wcsimrootevent->GetNtrack();
      const int ncherenkovhits      = wcsimrootevent->GetNcherenkovhits();
      const int ncherenkovhittimes  = wcsimrootevent->GetNcherenkovhittimes();
      const int ntubeshit           = wcsimrootevent->GetNumTubesHit();
      const int ncherenkovdigihits0 = wcsimrootevent->GetNcherenkovdigihits(); 
      const int ntubesdigihit0      = wcsimrootevent->GetNumDigiTubesHit();
      const double vtx0 = wcsimrootevent->GetVtx(0);
      const double vtx1 = wcsimrootevent->GetVtx(1);
      const double vtx2 = wcsimrootevent->GetVtx(2);

      if(verbose || ((iev % 100) == 0)) {
	printf("********************************************************");
	cout << "Event " << iev << " of " << nevents_loop << " has " << ntriggers << " triggers" << endl;
	printf("Evt, date %d %d\n", wcsimrootevent->GetHeader()->GetEvtNum(),
	       wcsimrootevent->GetHeader()->GetDate());
	printf("Mode %d\n", wcsimrootevent->GetMode());
	printf("Number of subevents %d\n",
	       wcsimrootsuperevent->GetNumberOfSubEvents());
	//vertex
	printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
	printf("Vtx %f %f %f\n", vtx0, vtx1, vtx2);
	printf("Jmu %d\n", wcsimrootevent->GetJmu());
	printf("Npar %d\n", wcsimrootevent->GetNpar());
	//tracks
	printf("ntracks=%d\n",ntracks);
	//hits & digits
	printf("Ncherenkovhits (unique PMTs with hits)  %d\n", ncherenkovhits);
	printf("Ncherenkovhittimes (number of raw hits) %d\n", ncherenkovhittimes);
	printf("Ncherenkovdigihits (number of digits) in trigger 0   %d\n", ncherenkovdigihits0);
	printf("NumTubesHit       %d\n", wcsimrootevent->GetNumTubesHit());
	printf("NumDigitizedTubes in trigger 0 %d\n", wcsimrootevent->GetNumDigiTubesHit());
      }

      // Loop through elements in the TClonesArray of WCSimTracks
      for (int itrack = 0; itrack < ntracks; itrack++) {
	TObject *element = (wcsimrootevent->GetTracks())->At(itrack);      
	WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);
	if(verbose) {
	  printf("Track ipnu: %d\n",wcsimroottrack->GetIpnu());
	  printf("Track parent ID: %d\n",wcsimroottrack->GetParenttype());
	  for (int j=0; j<3; j++)
	    printf("Track dir: %d %f\n", j, wcsimroottrack->GetDir(j));
	}
      }  //itrack // End of loop over tracks
    

      //
      // Now look at the raw Cherenkov+noise hits
      //
      if(verbose)
	cout << "RAW HITS:" << endl;

      // Grab the big arrays of times and parent IDs
      TClonesArray *timeArray = wcsimrootevent->GetCherenkovHitTimes();
    
      //calculate total p.e. in event
      int totalPe = 0;
      // Loop through elements in the TClonesArray of WCSimRootCherenkovHits
      for(int ipmt = 0; ipmt < ncherenkovhits; ipmt++) {
	TObject * Hit = (wcsimrootevent->GetCherenkovHits())->At(ipmt);
	WCSimRootCherenkovHit * wcsimrootcherenkovhit =
	  dynamic_cast<WCSimRootCherenkovHit*>(Hit);
	int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
	int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
	int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
	WCSimRootPMT pmt   = geo->GetPMT(tubeNumber-1);
	totalPe += peForTube;
	if(verbose)
	  printf("Total pe for tube %d: %d times( ", tubeNumber, peForTube);
	for(int irawhit = 0; irawhit < peForTube; irawhit++) {
	  TObject * HitTime = (wcsimrootevent->GetCherenkovHitTimes())->At(timeArrayIndex + irawhit);
	  WCSimRootCherenkovHitTime * wcsimrootcherenkovhittime =
	    dynamic_cast<WCSimRootCherenkovHitTime*>(HitTime);
	  double truetime = wcsimrootcherenkovhittime->GetTruetime();
	  if(verbose)
	    printf("%6.2f ", truetime);
	  h1hittime->Fill(truetime);
	  if(wcsimrootcherenkovhittime->GetParentID() == -1) {
	    h1hittime_noise->Fill(truetime);
	    if(hists_per_event)
	      h1event_hittime[ev][1]->Fill(truetime);
	  }
	  else {
	    h1hittime_photon->Fill(truetime);
	    if(hists_per_event)
	      h1event_hittime[ev][0]->Fill(truetime);
	  }
	}//irawhit
	if(verbose)
	  cout << ")" << endl;
      }//ipmt
      if(verbose)
	cout << "Total Pe : " << totalPe << endl;

      //    
      // Now look at digitized hit info
      //
      if(verbose)
	cout << "DIGITIZED HITS:" << endl;

      //
      // Digi hits are arranged in subevents, so loop over these first
      //
      const int ntriggers_loop = max_ntriggers > 0 ? TMath::Min(max_ntriggers, ntriggers) : ntriggers;

      //save a pointer to the 0th WCSimRootTrigger, so can access track/hit information in triggers >= 1
      wcsimroottrigger0 = wcsimrootevent;
      tdigitimes.clear();
      for (int itrigger = 0 ; itrigger < ntriggers_loop; itrigger++) {

	//count the number of noise/photon hits making up all the digits in this trigger
	int n_noise_hits_total = 0, n_photon_hits_total = 0;

	wcsimrootevent = wcsimrootsuperevent->GetTrigger(itrigger);
	if(verbose)
	  cout << "Sub event number = " << itrigger << "\n";

	const int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
	const int ntubeshitdigi      = wcsimrootevent->GetNumDigiTubesHit();
	if(verbose)
	  printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);

	const int            trigger_time = wcsimrootevent->GetHeader()->GetDate();
	const TriggerType_t  trigger_type = wcsimrootevent->GetTriggerType();
	std::vector<Float_t> trigger_info = wcsimrootevent->GetTriggerInfo();

	h1triggertime->Fill(trigger_time);
	if(trigger_info.size() > 0) {
	  if((trigger_type == kTriggerNDigits) || (trigger_type == kTriggerNDigitsTest)) {
	    h1ndigihitstrigger->Fill(trigger_info[0]);
	  }
	}
	h1ndigihits->Fill(ncherenkovdigihits);
	h1ntubeshitdigi->Fill(ntubeshitdigi);

	h1triggertype->Fill(WCSimEnumerations::EnumAsString(trigger_type).c_str(), 1);
	if(verbose) {
	  cout << "Passed trigger "
	       << WCSimEnumerations::EnumAsString(trigger_type)
	       << " with timestamp " << trigger_time
	       << " and " << ncherenkovdigihits
	       << " hits in the saved subevent region";
	  if(trigger_info.size() > 0) {
	    if((trigger_type == kTriggerNDigits) /*|| (trigger_type == kTriggerNHitsSKDETSIM)*/ || (trigger_type == kTriggerNDigitsTest))
	      cout << " (" << trigger_info[0]
		   << " in the 200nsec region)";
	    else if(trigger_type == kTriggerLocalNHits)
	      cout << " (PMT with tubeID " << trigger_info[1]
		   << " fired the trigger with " << trigger_info[0]
		   << " hits)";
	  }
	  cout << endl;
	}

	if(ncherenkovdigihits > 0)
	  num_trig++;


	tout->Fill();
      }//itrigger // End of loop over triggers

      wcsimroottrigger0 = 0;

      // reinitialize super event between loops.
      wcsimrootsuperevent->ReInitialize();

    }//iev // End of loop over events

  cout << "---------------------------------" << endl
       << "Run summary" << endl
       << "nevent (run over) " << nevent << endl
       << "num_trig (run over) " << num_trig << endl;

  //save histograms in .root file
  fout->cd();
  if(write_tree)
    tout->Write();


  return 0;
}
