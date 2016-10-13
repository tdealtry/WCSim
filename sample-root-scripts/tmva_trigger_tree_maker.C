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
#include <TMath.h>
#include <TVector3.h>

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

void GetDistributionQuantities(const vector<double> & v, double & av, double & rms, double & skew, double & mode, double & median, vector<int> & histogram)
{
  const size_t n = v.size();
  if(!n) {
    av     = -99999;
    rms    = -99999;
    skew   = -99999;
    mode   = -99999;
    median = -99999;
    return;
  }
  double min = 999999, max = -999999;

  //mean
  av = 0;
  for(size_t i = 0; i < n; i++) {
    double val = v[i];
    av += val;
    if(val < min) min = val;
    if(val > max) max = val;
  }//i
  av /= n;

  //rms
  rms = 0;
  for(size_t i = 0; i < n; i++) {
    double val = v[i];
    rms += (val-av) * (val-av);
  }//i
  rms /= n - 1;
  rms = TMath::Sqrt(rms);

  //skew
  skew = 0;
  for(size_t i = 0; i < n; i++) {
    double val = v[i];
    skew += (val-av) * (val-av) * (val-av);
  }//i
  skew /= rms * rms * rms;
  skew *= n / ((n-1) * (n-2));

  //mode
  const size_t nbins = histogram.size();
  for(size_t i = 0; i < nbins; i++)
    histogram[i] = 0;
  double binsize = (max - min) / nbins;
  for(size_t i = 0; i < n; i++) {
    double val = v[i];
    size_t binnum = floor((val - min) / binsize);
    binnum = TMath::Min(binnum, nbins - 1);
    //cout << val << "\tbinnum " << binnum << "\tmin " << min << "\tmax " << max << "\tbinsize " << binsize << "\tnbins " << nbins << endl;
    histogram.at(binnum)++;
  }//i
  mode = *(std::max_element(histogram.begin(), histogram.end()));
  mode *= binsize;

  //median
  median = v[n / 2];

  return;
}

// Simple example of reading a generated Root file
int tmva_trigger_tree_maker(char *filename=NULL, bool verbose=false, const int tmva_window = 1000, Long64_t max_nevents = 999999999999, int max_ntriggers = -1)
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
  printf("nevent (to actually run over) %d\n",nevent);
 
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
  if(geotree->GetEntries() == 0) {
    cout << "geotree not found!" << endl;
      exit(9);
  }
  geotree->GetEntry(0);

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

  //create output file for storing output tree
  TString filenameout(filename);
  TFile * fout = new TFile(create_filename("analysed_", filenameout).Data(), "RECREATE");
  fout->cd();
  //create output tree & book branches
  TTree * tout = new TTree("low_level_trigger", "Low-level variables to use in a TMVA trigger");
  int    nhits;
  //charge
  double av_charge, rms_charge, skew_charge, mode_charge, median_charge;
  //time
  double av_time, rms_time, skew_time, mode_time, median_time;
  //position
  double av_position_q, rms_position_q, skew_position_q, mode_position_q, median_position_q;
  double av_position_r, rms_position_r, skew_position_r, mode_position_r, median_position_r;
  double av_position_z, rms_position_z, skew_position_z, mode_position_z, median_position_z;
  //angle: take every 2 hit combination & calculate the angle (using 0,0,0)
  double av_angle, rms_angle, skew_angle, mode_angle, median_angle;
  //3hitangle: take every 3 hit combination & calculate the angle
  double av_3hitangle, rms_3hitangle, skew_3hitangle, mode_3hitangle, median_3hitangle;
  //solar_anisotropy: create plane using normal of Sun direction. Solar nu directional -> expect fewer hits closer to Sun
  double av_solar_anisotropy, rms_solar_anisotropy, skew_solar_anisotropy, mode_solar_anisotropy, median_solar_anisotropy;
  double solar_anisotropy_ratio; // calculate NHITS ratio between close & far (tank split in two)
  //truth
  double TRUE_last_physics_hit; //time of last physics hit, relative to start of window (defined as first TRUE physics hit)
  double TRUE_fraction_physics_hit_in_window; //fraction of the physics hits that are used in the histograms
  double TRUE_e_energy;
  TVector3 TRUE_e_direction, TRUE_e_position;
  //branches
  tout->Branch("nhits", &nhits);
  //charge
  tout->Branch("av_charge",     &av_charge);
  tout->Branch("rms_charge",    &rms_charge);
  tout->Branch("skew_charge",   &skew_charge);
  tout->Branch("mode_charge",   &mode_charge);
  tout->Branch("median_charge", &median_charge);
  //time
  tout->Branch("av_time",     &av_time);
  tout->Branch("rms_time",    &rms_time);
  tout->Branch("skew_time",   &skew_time);
  tout->Branch("mode_time",   &mode_time);
  tout->Branch("median_time", &median_time);
  //position
  //q (theta)
  tout->Branch("av_position_q",     &av_position_q);
  tout->Branch("rms_position_q",    &rms_position_q);
  tout->Branch("skew_position_q",   &skew_position_q);
  tout->Branch("mode_position_q",   &mode_position_q);
  tout->Branch("median_position_q", &median_position_q);
  //r
  tout->Branch("av_position_r",     &av_position_r);
  tout->Branch("rms_position_r",    &rms_position_r);
  tout->Branch("skew_position_r",   &skew_position_r);
  tout->Branch("mode_position_r",   &mode_position_r);
  tout->Branch("median_position_r", &median_position_r);
  //z
  tout->Branch("av_position_z",     &av_position_z);
  tout->Branch("rms_position_z",    &rms_position_z);
  tout->Branch("skew_position_z",   &skew_position_z);
  tout->Branch("mode_position_z",   &mode_position_z);
  tout->Branch("median_position_z", &median_position_z);
  //angle
  tout->Branch("av_angle",     &av_angle);
  tout->Branch("rms_angle",    &rms_angle);
  tout->Branch("skew_angle",   &skew_angle);
  tout->Branch("mode_angle",   &mode_angle);
  tout->Branch("median_angle", &median_angle);
  //3hitangle
  tout->Branch("av_3hitangle",     &av_3hitangle);
  tout->Branch("rms_3hitangle",    &rms_3hitangle);
  tout->Branch("skew_3hitangle",   &skew_3hitangle);
  tout->Branch("mode_3hitangle",   &mode_3hitangle);
  tout->Branch("median_3hitangle", &median_3hitangle);
  //solar_anisotropy
  tout->Branch("av_solar_anisotropy",     &av_solar_anisotropy);
  tout->Branch("rms_solar_anisotropy",    &rms_solar_anisotropy);
  tout->Branch("skew_solar_anisotropy",   &skew_solar_anisotropy);
  tout->Branch("mode_solar_anisotropy",   &mode_solar_anisotropy);
  tout->Branch("median_solar_anisotropy", &median_solar_anisotropy);
  tout->Branch("solar_anisotropy_ratio",  &solar_anisotropy_ratio);
  //truth
  tout->Branch("TRUE_last_physics_hit", &TRUE_last_physics_hit);
  tout->Branch("TRUE_fraction_physics_hit_in_window", &TRUE_fraction_physics_hit_in_window);
  tout->Branch("TRUE_e_energy",         &TRUE_e_energy);
  tout->Branch("TRUE_e_direction",      &TRUE_e_direction);
  tout->Branch("TRUE_e_position",       &TRUE_e_position);

  //book vectors to store ALL hit information
  // need to ensure that all vectors are filled at the same time with the same number of entries
  vector<double>   v_digit_charge, v_digit_time, v_digit_time_physics, v_digit_time_offset;
  vector<double>   v_digit_position_q, v_digit_position_r, v_digit_position_z;
  vector<double>   v_digit_angle, v_digit_3hitangle, v_digit_solar_anisotropy;
  vector<int>      v_digit_pmt;
  vector<TVector3> v_digit_pmtvec;
  vector<int>      v_histogram(101, 0);

  WCSimRootPMT * pmtobj = 0;

  int num_trig = 0;  

  // Now loop over events
  for (int ev=0; ev<nevent; ev++)
  {
    v_digit_charge.clear();
    v_digit_time.clear();
    v_digit_time_physics.clear();
    v_digit_time_offset.clear();
    v_digit_pmt.clear();
    v_digit_pmtvec.clear();
    v_digit_position_q.clear();
    v_digit_position_r.clear();
    v_digit_position_z.clear();
    v_digit_angle.clear();
    v_digit_3hitangle.clear();
    v_digit_solar_anisotropy.clear();

    // Read the event from the tree into the WCSimRootEvent instance
    cout << "getting entry " << ev << " of " << tree->GetEntries() << endl;
    tree->GetEntry(ev);
    cout << "gotcha" << endl;
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);

    const int ntriggers = wcsimrootsuperevent->GetNumberOfEvents();

    if(verbose || ((ev % 100) == 0))
      cout << "Event " << ev << " of " << nevent << " has " << ntriggers << " triggers" << endl;

    const double vtx0 = wcsimrootevent->GetVtx(0);
    const double vtx1 = wcsimrootevent->GetVtx(1);
    const double vtx2 = wcsimrootevent->GetVtx(2);
    TRUE_e_position.SetXYZ(vtx0, vtx1, vtx2);
    if(verbose){
      printf("********************************************************");
      printf("Evt, date %d %d\n", wcsimrootevent->GetHeader()->GetEvtNum(),
	     wcsimrootevent->GetHeader()->GetDate());
      printf("Mode %d\n", wcsimrootevent->GetMode());
      printf("Number of subevents %d\n",
	     wcsimrootsuperevent->GetNumberOfSubEvents());
      printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
      printf("Vtx %f %f %f\n", vtx0, vtx1, vtx2);
      printf("Jmu %d\n", wcsimrootevent->GetJmu());
      printf("Npar %d\n", wcsimrootevent->GetNpar());
    }


    // Now read the tracks in the event
    
    // Get the number of tracks
    const int ntracks = wcsimrootevent->GetNtrack();
    if(verbose) printf("ntracks=%d\n",ntracks);
    
    // Loop through elements in the TClonesArray of WCSimTracks
    bool found_electron = false;
    for (int itrack = 0; itrack < ntracks; itrack++)
    {
      TObject *element = (wcsimrootevent->GetTracks())->At(itrack);      
      WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);
      if(verbose){
	printf("Track number: %d\n", itrack);
	printf("Track flag: %d\n", wcsimroottrack->GetFlag());
	printf("Track ipnu: %d\n",wcsimroottrack->GetIpnu());
	printf("Track parent ID: %d\n",wcsimroottrack->GetParenttype());
	for (int j=0; j<3; j++)
	  printf("Track dir: %d %f\n", j, wcsimroottrack->GetDir(j));
      }
      if(itrack == 2 && wcsimroottrack->GetIpnu() == 11 && wcsimroottrack->GetFlag() == 0) {
	if(found_electron) {
	  cerr << "Found multiple true electrons" << endl;
	  exit(-1);
	}
	found_electron = true;
	TRUE_e_direction.SetXYZ(wcsimroottrack->GetDir(0), wcsimroottrack->GetDir(1), wcsimroottrack->GetDir(2));
	TRUE_e_energy = wcsimroottrack->GetE();
      }
    }  //itrack // End of loop over tracks
    if(!found_electron) {
      cerr << "Found no true electrons" << endl;
      exit(-1);
    }
    
    //get number of hits and digits
    const int ncherenkovhits      = wcsimrootevent->GetNcherenkovhits();
    const int ncherenkovhittimes  = wcsimrootevent->GetNcherenkovhittimes();
    const int ntubeshit           = wcsimrootevent->GetNumTubesHit();
    const int ncherenkovdigihits0 = wcsimrootevent->GetNcherenkovdigihits(); 
    const int ntubesdigihit0      = wcsimrootevent->GetNumDigiTubesHit();
    if(verbose){
      printf("node id: %i\n", ev);
      printf("Ncherenkovhits (unique PMTs with hits)  %d\n", ncherenkovhits);
      printf("Ncherenkovhittimes (number of raw hits) %d\n", ncherenkovhittimes);
      printf("Ncherenkovdigihits (number of digits) in trigger 0   %d\n", ncherenkovdigihits0);
      printf("NumTubesHit       %d\n", wcsimrootevent->GetNumTubesHit());
      printf("NumDigitizedTubes in trigger 0 %d\n", wcsimrootevent->GetNumDigiTubesHit());
    }

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
      if(verbose) {
	cout << "Passed trigger "
	     << WCSimEnumerations::EnumAsString(trigger_type)
	     << " with timestamp " << trigger_time
	     << " and " << ncherenkovdigihits
	     << " hits in the saved subevent region";
	if(trigger_info.size() > 0) {
	  if((trigger_type == kTriggerNDigits) || (trigger_type == kTriggerNDigitsTest))
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

      // Loop through elements in the TClonesArray of WCSimRootCherenkovDigHits
      for(int idigipmt = 0; idigipmt < ncherenkovdigihits; idigipmt++) {
	//get the digit
	if(verbose)
	  cout << "Getting digit " << idigipmt << endl;
	TObject *element = (wcsimrootevent->GetCherenkovDigiHits())->At(idigipmt);
    	WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = 
    	  dynamic_cast<WCSimRootCherenkovDigiHit*>(element);
	
	//find out whether this is due to noise, a real photon, or both
	int tube_id = wcsimrootcherenkovdigihit->GetTubeId();
	int timeArrayIndex = -1;
	int peForTube      = -1;
	//loop through the WCSimRootCherenkovHit array to find the entry for this tube ID
	for(int ipmt = 0; ipmt < ncherenkovhits; ipmt++) {
	  if(verbose)
	    cout << "Getting hit " << ipmt << " of " << ncherenkovhits << endl;
	  TObject *Hit = (wcsimroottrigger0->GetCherenkovHits())->At(ipmt);
	  WCSimRootCherenkovHit *wcsimrootcherenkovhit =
	    dynamic_cast<WCSimRootCherenkovHit*>(Hit);
	  int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
	  if(tube_id == tubeNumber) {
	    timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
	    peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
	    break;
	  }
	}//ipmt
	int n_noise_hits = 0, n_photon_hits = 0, n_unknown_hits = 0;
	if(timeArrayIndex == -1) {
	  if(verbose)
	    cout << "No PMT hits found for digit " << idigipmt << " with tube ID " << tube_id << endl;
	}
	else {
	  if(verbose)
	    cout << peForTube << " PMT hits found for digit " << idigipmt << " with tube ID " << tube_id << endl;
	  //loop over the rawhits ids of hits that made up the digit
	  vector<int> rawhit_ids = wcsimrootcherenkovdigihit->GetPhotonIds();
	  if(verbose)
	    cout << rawhit_ids.size() << " rawhits made up this digit" << endl;
	  for(unsigned int irawhit = 0; irawhit < rawhit_ids.size(); irawhit++) {
	    int this_rawhit = rawhit_ids[irawhit];
	    if(verbose)
	      cout << "Attempting to look for rawhit " << this_rawhit+1 << " in WCSimRootCherenkovHitTime array...";
	    if(this_rawhit >= peForTube) {
	      //if(verbose)
	      cerr << " There are only " << peForTube << " rawhits in this PMT" << endl;
	      continue;
	    }
	    //now look in the WCSimRootCherenkovHitTime array to count the number of photon / dark noise hits
	    TObject *Hit = (wcsimroottrigger0->GetCherenkovHitTimes())->At(timeArrayIndex + this_rawhit);
	    WCSimRootCherenkovHitTime *wcsimrootcherenkovhittime =
	      dynamic_cast<WCSimRootCherenkovHitTime*>(Hit);
	    const double hittime  = wcsimrootcherenkovhittime->GetTruetime();
	    const int    parentid = wcsimrootcherenkovhittime->GetParentID();
	    if(verbose)
	      cout << " hit time " << hittime << " " << parentid << endl;
	    if(parentid == -1) {
	      n_noise_hits++;
	    }
	    else if(parentid < 0) {
	      n_unknown_hits++;
	    }
	    else {
	      n_photon_hits++;
	    }
	  }//irawhit
	}//hits in this PMT found

	const double digitime = wcsimrootcherenkovdigihit->GetT();
	const double digipe   = wcsimrootcherenkovdigihit->GetQ();

	//find the PMT position vector
	//GetTubeNo() runs from 1 to NPMT
	//tube_id runs from 0 to NPMT-1
	const int tube_id_to_get = tube_id - 1;
	if(tube_id_to_get < 0 || tube_id_to_get > geo->GetWCNumPMT()) {
	  cerr << "tube_id_to_get " << tube_id_to_get << " GetWCNumPMT " << geo->GetWCNumPMT() << endl;
	  exit(-1);
	}
	pmtobj = geo->GetPMTPointer(tube_id_to_get);
	if(pmtobj->GetTubeNo() != tube_id_to_get + 1) {
	  cerr << "tube_id_to_get + 1 " << tube_id_to_get + 1 << " != WCSimRootPMT->GetTubeNo() " << pmtobj->GetTubeNo() << endl;
	  exit(-1);
	}
	TVector3 tube_vec(pmtobj->GetPosition(0), pmtobj->GetPosition(1), pmtobj->GetPosition(2));

	//fill the vectors
	v_digit_charge.push_back(digipe);
	v_digit_time.push_back(digitime);
	v_digit_pmt.push_back(tube_id);
	v_digit_pmtvec.push_back(tube_vec);
	v_digit_position_q.push_back(tube_vec.Theta());
	v_digit_position_r.push_back(tube_vec.Perp());
	v_digit_position_z.push_back(tube_vec.Z());
	if(n_photon_hits)
	  v_digit_time_physics.push_back(digitime); //includes both physics-only & noise+physics

	if(verbose){
	  printf("q, t, tubeid, nphotonhits, nnoisehits, nunknownhits: %f %f %d  %d %d %d\n",
		 digipe,
		 digitime,
		 wcsimrootcherenkovdigihit->GetTubeId(),
		 n_photon_hits,
		 n_noise_hits,
		 n_unknown_hits);
	}
      }//idigipmt // End of loop over Cherenkov digihits

    }//itrigger // End of loop over triggers

    cout << "Loop over triggers ended" << endl;

    wcsimroottrigger0 = 0;

    // reinitialize super event between loops.
    wcsimrootsuperevent->ReInitialize();

    assert(v_digit_time.size() == v_digit_charge.size());
    assert(v_digit_time.size() == v_digit_pmt.size());
    assert(v_digit_time_physics.size() <= v_digit_time.size());

    // we now have a list of time/charge/pmt for each hit (and time of the physics hits)
    // need to get the min/max physics hit time, which will allow the histograms to be filled
    double min_physics_time = 99999, TRUE_last_physics_hit = -99999;
    if(!v_digit_time_physics.size()) {
      min_physics_time = std::min_element(v_digit_time.begin(), v_digit_time.end()); //dark-noise only events start from the first non-physics hit
      min_physics_time += 10; // + a little (accounts for time smearing, where events are "smeared out" but not "smeared in"
      TRUE_last_physics_hit = min_physics_time;
    }
    else {
      for(size_t idigitime_physics = 0, nphysicshits = v_digit_time_physics.size(); idigitime_physics < nphysicshits; idigitime_physics++) {
	double thistime = v_digit_time_physics[idigitime_physics];
	if(thistime < min_physics_time)
	  min_physics_time = thistime;
	if(thistime > TRUE_last_physics_hit)
	  TRUE_last_physics_hit = thistime;
      }//idigitime_physics
    }

    double max_time = min_physics_time + tmva_window;

    //we have the time window. Now lets check how fraction of physics hits in the window
    int nphysics_in_window = 0;
    for(size_t idigitime_physics = 0, nphysicshits = v_digit_time_physics.size(); idigitime_physics < nphysicshits; idigitime_physics++) {
      double thistime = v_digit_time_physics[idigitime_physics];
      if(thistime <= max_time)
	nphysics_in_window++;
    }//idigitime_physics
    TRUE_fraction_physics_hit_in_window = (double)nphysics_in_window / (double)v_digit_time_physics.size();

    //
    //now it's onto the physics+noise variables (i.e. TMVA inputs)
    //

    //first, remove entries not in the time window (erase later)
    for(int idigit = v_digit_time.size() - 1; idigit >= 0; idigit--) {
      double thistime = v_digit_time[idigit];
      if(thistime > max_time || thistime < min_physics_time) {
	//not in the tmva window
	v_digit_time.erase(v_digit_time.begin() + idigit);
	v_digit_charge.erase(v_digit_charge.begin() + idigit);
	v_digit_pmt.erase(v_digit_pmt.begin() + idigit);
	v_digit_pmtvec.erase(v_digit_pmtvec.begin() + idigit);
	v_digit_position_q.erase(v_digit_position_q.begin() + idigit);
	v_digit_position_r.erase(v_digit_position_r.begin() + idigit);
	v_digit_position_z.erase(v_digit_position_z.begin() + idigit);
      }
      else
	v_digit_time_offset.push_back(thistime - min_physics_time);
    }//idigit

    assert(v_digit_time.size() == v_digit_time_offset.size());

    /*
    cout << "Angle" << endl;

    //next, calculate the angle parameters
    for(size_t idigit = 0, ndigits = v_digit_time.size(); idigit < ndigits; idigit++) {
      for(size_t jdigit = 0; jdigit < ndigits; jdigit++) {
	if(jdigit > idigit)
	  v_digit_angle.push_back(v_digit_pmtvec[idigit].Angle(v_digit_pmtvec[jdigit])); //2pmt angle
	if(idigit == jdigit)
	  continue;
	for(size_t kdigit = 0; kdigit < ndigits; kdigit++) {
	  if(kdigit == idigit || kdigit == jdigit)
	    continue;
	  double l1 = (v_digit_pmtvec[idigit] - v_digit_pmtvec[jdigit]).Mag();
	  double l2 = (v_digit_pmtvec[jdigit] - v_digit_pmtvec[kdigit]).Mag();
	  double l3 = (v_digit_pmtvec[kdigit] - v_digit_pmtvec[idigit]).Mag();
	  double threehitangle = TMath::ACos((l1*l1 + l2*l2 - l3*l3) / (2 * l1 * l2));
	  v_digit_3hitangle.push_back(threehitangle);
	}//kdigit
      }//jdigit
    }//idigit
    */

    //and the solar anisotropy
    int nplane_hi = 0, nplane_lo = 0;
    for(size_t idigit = 0, ndigits = v_digit_time.size(); idigit < ndigits; idigit++) {
      //http://mathworld.wolfram.com/Point-PlaneDistance.html
      double x0 = v_digit_pmtvec[idigit].X();
      double y0 = v_digit_pmtvec[idigit].Y();
      double z0 = v_digit_pmtvec[idigit].Z();
      double a  = TRUE_e_direction.X();
      double b  = TRUE_e_direction.Y();
      double c  = TRUE_e_direction.Z();
      double plane_distance = (a*x0 + b*y0 + c*z0) / TMath::Sqrt(a*a + b*b + c*c);
      (plane_distance < 0) ? nplane_lo++ : nplane_hi++;
      v_digit_solar_anisotropy.push_back(plane_distance);
    }//idigit
    solar_anisotropy_ratio = (double)nplane_hi / (double)nplane_lo;

    //finally, get the mean/rms/median/mode/skew of each distribution
    nhits = v_digit_time.size();
    GetDistributionQuantities(v_digit_charge, av_charge, rms_charge, skew_charge, mode_charge, median_charge, v_histogram);
    GetDistributionQuantities(v_digit_time, av_time, rms_time, skew_time, mode_time, median_time, v_histogram);
    GetDistributionQuantities(v_digit_position_q, av_position_q, rms_position_q, skew_position_q, mode_position_q, median_position_q, v_histogram);
    GetDistributionQuantities(v_digit_position_r, av_position_r, rms_position_r, skew_position_r, mode_position_r, median_position_r, v_histogram);
    GetDistributionQuantities(v_digit_position_z, av_position_z, rms_position_z, skew_position_z, mode_position_z, median_position_z, v_histogram);
    GetDistributionQuantities(v_digit_angle, av_angle, rms_angle, skew_angle, mode_angle, median_angle, v_histogram);
    GetDistributionQuantities(v_digit_3hitangle, av_3hitangle, rms_3hitangle, skew_3hitangle, mode_3hitangle, median_3hitangle, v_histogram);
    GetDistributionQuantities(v_digit_solar_anisotropy, av_solar_anisotropy, rms_solar_anisotropy, skew_solar_anisotropy, mode_solar_anisotropy, median_solar_anisotropy, v_histogram);
    
    // fill the tree
    tout->Fill();

  }//ev // End of loop over events

  cout << "---------------------------------" << endl
       << "Run summary" << endl
       << "nevent (run over) " << nevent << endl
       << "num_trig (run over) " << num_trig << endl;

  //save tree in .root file
  fout->cd();
  tout->Write();
  fout->Close();

  return 0;
}
