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

#include "trigger_tools.cxx"
#include "itc_tools.cxx"

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
int itc_analysis(const char *filename=NULL, const bool verbose=false, 
		 const long max_nevents = -1, const int max_ntriggers = -1)
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


  //
  // CREATE TOOL INSTANCE
  //
  vector<itc_tools *> itcconfigs;
  double smallfracs[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9};
  int  largewindows[] = {200, 300, 400, 500, 600, 800, 1000, 1250, 1500, 1750, 2000};
  const int nsmallfracs   = sizeof(smallfracs)   / sizeof(double);
  const int nlargewindows = sizeof(largewindows) / sizeof(int);
  const bool one_time_slice = false;
  for(int ilarge = 0; ilarge < nlargewindows; ilarge++) {
    for(int ismall = 0; ismall < nsmallfracs; ismall++) {
      itcconfigs.push_back(new itc_tools(fout, smallfracs[ismall]*largewindows[ilarge], largewindows[ilarge], 0, !ismall, one_time_slice));
    }//ismall
  }//ilarge
  const int nconfigs = itcconfigs.size();

  //
  // LOOP OVER EVENTS
  //

  //setup event loop counters
  int num_trig = 0;
  
  //loop over events
  for (int iev = 0; iev < nevents_loop; iev++)
    {
      // Read the event from the tree into the WCSimRootEvent instance
      tree->GetEntry(iev);
      wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);

      //
      // GLOBAL EVENT INFO
      //

      const int  ntriggers = wcsimrootsuperevent->GetNumberOfEvents();
      const int  ntriggers_loop = max_ntriggers > 0 ? TMath::Min(max_ntriggers, ntriggers) : ntriggers;
      const int  ntracks = wcsimrootevent->GetNtrack();
      const int  ncherenkovhits      = wcsimrootevent->GetNcherenkovhits();
      const int  ncherenkovhittimes  = wcsimrootevent->GetNcherenkovhittimes();
      const int  ntubeshit           = wcsimrootevent->GetNumTubesHit();
      const int  ncherenkovdigihits0 = wcsimrootevent->GetNcherenkovdigihits(); 
      const int  ntubesdigihit0      = wcsimrootevent->GetNumDigiTubesHit();
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

      //    
      // DIGITS INFO
      //
      if(verbose)
	cout << "DIGITIZED HITS:" << endl;

      // WARNING: just looking in the first trigger
      //save a pointer to the 0th WCSimRootTrigger, so can access track/hit information in triggers >= 1
      wcsimroottrigger0 = wcsimrootevent;

      num_trig++;

      for(int iconfig = 0; iconfig < nconfigs; iconfig++) {
	if(iconfig == 0)
	  itcconfigs[iconfig]->PopulateDigitTimes(wcsimroottrigger0, false, wcsimrootsuperevent);
	else
	  itcconfigs[iconfig]->CopyDigitTimes(itcconfigs[0]);
	itcconfigs[iconfig]->CalcMaxITC();
      }//iconfig

      wcsimroottrigger0 = 0;

      // reinitialize super event between loops.
      wcsimrootsuperevent->ReInitialize();

    }//iev // End of loop over events

  //tools.Write()

  cout << "---------------------------------" << endl
       << "Run summary" << endl
       << "nevent (run over) " << nevents_loop << endl
       << "num_trig (run over) " << num_trig << endl;

  //save histograms in .root file
  fout->cd();
  for(int iconfig = 0; iconfig < nconfigs; iconfig++) {
    itcconfigs[iconfig]->Write();
  }//iconfig  

  return 0;
}
