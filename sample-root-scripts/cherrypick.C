/*********************************************************
 * A macro that reads in a WCSim output file and applies cuts
 *  to produce a copy of the output file with only the events that
 *  pass all cuts
 * This is useful for e.g. using the event display
 * The current cuts implemented are:
 *  - Trigger type
 *  - NDigits in the event
 *********************************************************/

#include <iostream>
#include <algorithm>
#include <vector>
#include <stdio.h>     
#include <stdlib.h>
#include <string>
#include <sstream>
#include <cassert>

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

TString CreateFilename(const char * prefix, TString& filename_string)
{
  TString prefix_string(prefix);
  TString outfilename = prefix_string + filename_string;
  return outfilename;
}

vector<string> Tokenize(string input, string delimiter)
{
  // split a string of 'delimiter' separated values and return each string
  // part as a vector<string> element.
  
  vector<string> string_parts;
  
  while(input.find_first_of(delimiter) < input.length()) {
    string_parts.push_back( input.substr(0, input.find_first_of(delimiter)) );
    input = input.substr(input.find_first_of(delimiter)+1, input.length());
  }
  string_parts.push_back(input);
  return string_parts;
}

template <typename T> T StringToNumber(std::string &Text)
{
  std::istringstream ss(Text);
  T result;
  return ss >> result ? result : 0;
}

// Simple example of reading a generated Root file
int cherrypick(char *filename=NULL, bool verbose=false, string cuts = "", Long64_t max_nevents = 999999999999)
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
  
  // Get the a pointer to the tree from the file
  TTree *tree = (TTree*)file->Get("wcsimT");
  
  // Get the number of events
  Long64_t nevent = tree->GetEntries();
  printf("nevent %d\n",nevent);
  nevent = TMath::Min(nevent, max_nevents); //cut the loop earlier
 
  // Create a WCSimRootEvent to put stuff from the tree in
  WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();

  // Set the branch address for reading from the tree
  TBranch *branch = tree->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);

  // Force deletion to prevent memory leak 
  tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

  // Geometry tree - only need 1 "event"
  TTree *geotree = (TTree*)file->Get("wcsimGeoT");
  WCSimRootGeom *geo = 0; 
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
  if (geotree->GetEntries() == 0) {
    cout << "geotree not found!" << endl;
      exit(9);
  }
  geotree->GetEntry(0);

  // start with the main "subevent", as it contains most of the info
  // and always exists.
  WCSimRootTrigger* wcsimrootevent;

  // make a new file to store the skimmed events
  TString filenameout(filename);
  TFile * fout = new TFile(CreateFilename("skimmed_", filenameout).Data(), "RECREATE");
  fout->cd();
  //and get the trees to put in
  TTree * newgeotree = geotree->CloneTree(-1); //clone all events
  TTree * newtree    = tree->CloneTree(0);     //clone no events (use Fill() later)

  //Parse the options
  vector<string> allcuts = Tokenize(cuts, ",");
  vector<TriggerType_t> allowedTriggers;
  int ndigitsMin = -1, ndigitsMax = 9999999;
  for(vector<string>::iterator it = allcuts.begin(); it != allcuts.end(); it++) {
    if((*it).find("TRIGGERS") == 0) {
      vector<string> triggersStr = Tokenize(*it, ":");
      cout << "Accepting the following triggers:";
      for(vector<string>::iterator it_t = triggersStr.begin() + 1; it_t != triggersStr.end(); it_t++) {
	TriggerType_t trigtype = WCSimEnumerations::TriggerTypeFromString(*it_t);
	if(trigtype == kTriggerUndefined)
	  exit(1);
	allowedTriggers.push_back(trigtype);
	cout << " " << *it_t;
      }
      cout << endl;
      assert(allowedTriggers.size());
    }//TRIGGERS
    else if((*it).find("NDIGITS") == 0) {
      vector<string> ndigitsStr = Tokenize(*it, ":");
      assert(ndigitsStr.size() == 3);
      ndigitsMin = StringToNumber<int>(ndigitsStr[1]);
      ndigitsMax = StringToNumber<int>(ndigitsStr[2]);
    }//NDIGITS
    else {
      cerr << "Unknown cut: " << *it << endl;
      exit(1);
    }
  }//allcuts
  cout << "Accepting NDigits in the range " << ndigitsMin << " " << ndigitsMax << endl;
  
  // Now loop over events
  int num_trig = 0;
  int num_passed = 0;
  for (int ev=0; ev<nevent; ev++) {
    // Read the event from the tree into the WCSimRootEvent instance
    tree->GetEntry(ev);
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);

    if(verbose || ((ev % 100) == 0)) {
      printf("********************************************************\n");
      cout << "Event " << ev << " of " << nevent << " has " << wcsimrootsuperevent->GetNumberOfEvents() << " trigger(s)" << endl;
    }
    num_trig += wcsimrootsuperevent->GetNumberOfEvents();

    if(verbose){
      printf("Evt, date = %d, %d\n", wcsimrootevent->GetHeader()->GetEvtNum(),
	     wcsimrootevent->GetHeader()->GetDate());
      printf("Mode %d\n", wcsimrootevent->GetMode());
      printf("Number of subevents %d\n",
	     wcsimrootsuperevent->GetNumberOfSubEvents());
      
      printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
      printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0),
	     wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
      printf("node id: %i\n", ev);
    }

    for (int itrigger = 0 ; itrigger < wcsimrootsuperevent->GetNumberOfEvents(); itrigger++) {
      if(verbose) {
	printf("Trigger %d in event %d\n", ev, itrigger);
	printf("trigger type: %s\n", WCSimEnumerations::EnumAsString(wcsimrootevent->GetTriggerType()).c_str());
	printf("Ncherenkovhits (unique PMTs with hits)  %d\n", wcsimrootevent->GetNcherenkovhits());
	printf("Ncherenkovhittimes (number of raw hits) %d\n", wcsimrootevent->GetNcherenkovhittimes());
	printf("Ncherenkovdigihits (number of digits)   %d\n", wcsimrootevent->GetNcherenkovdigihits());
	printf("NumTubesHit       %d\n", wcsimrootevent->GetNumTubesHit());
	printf("NumDigitizedTubes %d\n", wcsimrootevent->GetNumDigiTubesHit());
      }
      
      //cut on the trigger type
      bool passedTriggerType = false;
      if(allowedTriggers.size()) {
	TriggerType_t trigger_type = wcsimrootevent->GetTriggerType();
	for(int iat = 0; iat < allowedTriggers.size(); iat++) {
	  if(trigger_type == allowedTriggers[iat])
	    passedTriggerType = true;
	}//iat
      }//allowedTriggers

      //cut on the number of digits
      bool passedNDigits = false;
      const int ndigits = wcsimrootevent->GetNcherenkovdigihits();
      if((ndigits >= ndigitsMin) && (ndigits <= ndigitsMax))
	passedNDigits = true;

      //if all cuts passed, save this event
      if(passedTriggerType && passedNDigits) {
	cout << "PASSED ALL CUTS! Saving event" << endl;
	newtree->Fill();
	num_passed++;
	break; //break on the loop over triggers in the event - only need to save the event once
      }
    }//itrigger
    
    // reinitialize super event between loops.
    wcsimrootsuperevent->ReInitialize();
    
  }//ev // End of loop over events

  
  std::cout << "********************************************************" << endl;
  std::cout << "num_trig " << num_trig << endl;
  std::cout << "num_passed " << num_passed << endl;
  std::cout << "nevent " << nevent << endl;

  //save selected events in .root file
  fout->Write();
  fout->Close();

  
  return 0;
}
