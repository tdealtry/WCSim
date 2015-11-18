#include "WCSimWCTrigger.hh"
#include "WCSimWCPMT.hh"
#include "WCSimWCDigi.hh"
#include "WCSimWCHit.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4ios.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "WCSimDetectorConstruction.hh"
#include "WCSimPmtInfo.hh"
#include "WCSimDarkRateMessenger.hh"

#include "TMath.h"

#include <vector>
// for memset
#include <cstring>
#include <iostream>



// *******************************************
// BASE CLASS
// *******************************************

#ifndef WCSIMWCTRIGGER_VERBOSE
//#define WCSIMWCTRIGGER_VERBOSE
#endif
#ifndef WCSIMWCTRIGGER_PMT_NEIGHBOURS_VERBOSE
//#define WCSIMWCTRIGGER_PMT_NEIGHBOURS_VERBOSE
#endif

const double WCSimWCTriggerBase::offset = 950.0; // ns. apply offset to the digit time
const double WCSimWCTriggerBase::LongTime = 1E6; // ns = 1ms. event time


WCSimWCTriggerBase::WCSimWCTriggerBase(G4String name,
				       WCSimDetectorConstruction* inDetector,
				       WCSimWCDAQMessenger* myMessenger)
  :G4VDigitizerModule(name), DAQMessenger(myMessenger), myDetector(inDetector), triggerClassName(""), myPMTs(NULL)
{
  G4String colName = "WCDigitizedCollection";
  collectionName.push_back(colName);

  ReInitialize();

  if(DAQMessenger == NULL) {
    G4cerr << "WCSimWCDAQMessenger pointer is NULL when passed to WCSimWCTriggerBase constructor. Exiting..."
           << G4endl;
    exit(-1);
  }

  digitizeCalled = false;

  nPMTs = myDetector->GetTotalNumPmts();
}

WCSimWCTriggerBase::~WCSimWCTriggerBase()
{
  if(triggerClassName.compare("NHitsThenLocalNHits") == 0)
    delete localNHitsHits;
}

void WCSimWCTriggerBase::GetVariables()
{
  //set the options to class-specific defaults
  multiDigitsPerTrigger    = GetDefaultMultiDigitsPerTrigger();
  ndigitsThreshold         = GetDefaultNDigitsThreshold();
  ndigitsWindow            = GetDefaultNDigitsWindow();
  ndigitsPreTriggerWindow  = GetDefaultNDigitsPreTriggerWindow();
  ndigitsPostTriggerWindow = GetDefaultNDigitsPostTriggerWindow();

  //read the .mac file to override them
  if(DAQMessenger != NULL) {
    DAQMessenger->TellMeAboutTheTrigger(this);
    DAQMessenger->SetTriggerOptions();
  }
  else {
    G4cerr << "WCSimWCDAQMessenger pointer is NULL when used in WCSimWCTriggerBase::GetVariables(). Exiting..." 
	   << G4endl;
    exit(-1);
  }

  G4cout << (multiDigitsPerTrigger ? "Using mutiple digits per PMT per trigger" : "Using a maximum of 1 digit per PMT per trigger" ) << G4endl
	 << "Using NDigits threshold " << ndigitsThreshold
	 << (ndigitsAdjustForNoise ? " (will be adjusted for noise)" : "") << G4endl
	 << "Using NDigits trigger window " << ndigitsWindow << " ns" << G4endl
	 << "Using NDigits event pretrigger window " << ndigitsPreTriggerWindow << " ns" << G4endl
	 << "Using NDigits event posttrigger window " << ndigitsPostTriggerWindow << " ns" << G4endl
	 << "Using LocalNDigits neighbours " << localNHitsNeighbours << G4endl
	 << "Using LocalNDigits threshold " << localNHitsThreshold
	 << (localNHitsAdjustForNoise ? " (will be adjusted for noise)" : "") << G4endl
	 << "Using LocalNDigits window " << localNHitsWindow << " ns" << G4endl
	 << "Using SaveFailures event pretrigger window " << saveFailuresPreTriggerWindow << " ns" << G4endl
	 << "Using SaveFailures event posttrigger window " << saveFailuresPostTriggerWindow << " ns" << G4endl;
}

int WCSimWCTriggerBase::GetPreTriggerWindow(TriggerType_t t)
{
  switch(t) {
  case kTriggerNDigits:
  case kTriggerNDigitsTest:
  case kTriggerNHitsSKDETSIM:
  case kTriggerLocalNHits:
    return ndigitsPreTriggerWindow;
    break;
  case kTriggerFailure:
    return saveFailuresPreTriggerWindow;
    break;
  default:
    G4cerr << "WCSimWCTriggerBase::GetPreTriggerWindow() Unknown trigger type " << t
	   << "\t" << WCSimEnumerations::EnumAsString(t) << G4endl;
    exit(-1);
    break;
  }
}

int WCSimWCTriggerBase::GetPostTriggerWindow(TriggerType_t t)
{
  switch(t) {
  case kTriggerNDigits:
  case kTriggerNDigitsTest:
  case kTriggerNHitsSKDETSIM:
  case kTriggerLocalNHits:
    return ndigitsPostTriggerWindow;
    break;
  case kTriggerFailure:
    return saveFailuresPostTriggerWindow;
    break;
  default:
    G4cerr << "WCSimWCTriggerBase::GetPostTriggerWindow() Unknown trigger type " << t
	   << "\t" << WCSimEnumerations::EnumAsString(t) << G4endl;
    exit(-1);
    break;
  }
}

int WCSimWCTriggerBase::CalculateAverageDarkNoiseOccupancy(int npmts, int window_ns)
{
  double trigger_window_seconds = window_ns * 1E-9;
  double dark_rate_Hz = PMTDarkRate * 1000;
  double average_occupancy = dark_rate_Hz * trigger_window_seconds * npmts;
  
  G4cout << "Average number of PMTs active in a " << window_ns
	 << "ns window with a dark noise rate of " << PMTDarkRate
	 << "kHz is " << average_occupancy
	 << " (" << npmts << " total PMTs)"
	 << G4endl;
  return round(average_occupancy);
}

void WCSimWCTriggerBase::AdjustNDigitsThresholdForNoise()
{
  double average_occupancy = CalculateAverageDarkNoiseOccupancy(nPMTs, ndigitsWindow);
  G4cout << "Updating the NDigits threshold, from " << ndigitsThreshold
	 << " to " << ndigitsThreshold + round(average_occupancy) << G4endl;
  ndigitsThreshold += round(average_occupancy);
}

void WCSimWCTriggerBase::AdjustLocalNDigitsThresholdForNoise()
{
  int npmts = localNHitsNeighbours;
  double average_occupancy = CalculateAverageDarkNoiseOccupancy(npmts, localNHitsWindow);
  G4cout << "Updating the LocalNDigits threshold, from " << localNHitsThreshold
         << " to " << localNHitsThreshold + round(average_occupancy) << G4endl;
  localNHitsThreshold += round(average_occupancy);
}

void WCSimWCTriggerBase::Digitize()
{
  if(!digitizeCalled) {
    if(ndigitsAdjustForNoise)
      AdjustNDigitsThresholdForNoise();
    if(triggerClassName.compare("NHitsThenLocalNHits") == 0) {
      if(localNHitsAdjustForNoise)
	AdjustLocalNDigitsThresholdForNoise();
    }
    digitizeCalled = true;
  }

  //Input is collection of all digitized hits that passed the threshold
  //Output is all digitized hits which pass the trigger
  
  ReInitialize();

  //This is the output digit collection
  DigitsCollection = new WCSimWCTriggeredDigitsCollection ("/WCSim/glassFaceWCPMT",collectionName[0]);

  G4DigiManager* DigiMan = G4DigiManager::GetDMpointer();

  // Get the Digitized hits collection ID
  G4int WCDCID = DigiMan->GetDigiCollectionID("WCDigitizedStoreCollection");
  // Get the PMT Digits Collection
  WCSimWCDigitsCollection* WCDCPMT = 
    (WCSimWCDigitsCollection*)(DigiMan->GetDigiCollection(WCDCID));

  // Do the work  
  if (WCDCPMT) {
    DoTheWork(WCDCPMT);
  }
  
  StoreDigiCollection(DigitsCollection);
}

void WCSimWCTriggerBase::AlgNDigits(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits, bool test)
{

  //if test is true, we run the algorithm with 1/2 the threshold, and kTriggerNDigitsTest
  //for testing multiple trigger algorithms at once
  int this_ndigitsThreshold = ndigitsThreshold;
  TriggerType_t this_triggerType = kTriggerNDigits;
  if(test) {
    this_ndigitsThreshold /= 2;
    this_triggerType = kTriggerNDigitsTest;
  }

  //Now we will try to find triggers
  //loop over PMTs, and Digits in each PMT.  If ndigits > Threshhold in a time window, then we have a trigger

  int ntrig = 0;
  int window_start_time = 0;
  int window_end_time   = WCSimWCTriggerBase::LongTime - ndigitsWindow;
  int window_step_size  = 5; //step the search window along this amount if no trigger is found
  float lasthit;
  std::vector<int> digit_times;
  bool first_loop = true;

  G4cout << "WCSimWCTriggerBase::AlgNDigits. Number of entries in input digit collection: " << WCDCPMT->entries() << G4endl;
#ifdef WCSIMWCTRIGGER_VERBOSE
  int temp_total_pe = 0;
  for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
    temp_total_pe += (*WCDCPMT)[i]->GetTotalPe();
  }
  G4cout << "WCSimWCTriggerBase::AlgNDigits. " << temp_total_pe << " total p.e. input" << G4endl;
#endif

  // the upper time limit is set to the final possible full trigger window
  while(window_start_time <= window_end_time) {
    int n_digits = 0;
    float triggertime; //save each digit time, because the trigger time is the time of the first hit above threshold
    bool triggerfound = false;
    digit_times.clear();
    
    //Loop over each PMT
    for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
      //int tube=(*WCDCPMT)[i]->GetTubeID();
      //Loop over each Digit in this PMT
      for ( G4int ip = 0 ; ip < (*WCDCPMT)[i]->GetTotalPe() ; ip++) {
	int digit_time = (*WCDCPMT)[i]->GetTime(ip);
	//hit in trigger window?
	if(digit_time >= window_start_time && digit_time <= (window_start_time + ndigitsWindow)) {
	  n_digits++;
	  digit_times.push_back(digit_time);
	}
	//G4cout << digit_time << G4endl;
	//get the time of the last hit (to make the loop shorter)
	if(first_loop && (digit_time > lasthit))
	  lasthit = digit_time;
      }//loop over Digits
    }//loop over PMTs

    //if over threshold, issue trigger
    if(n_digits > this_ndigitsThreshold) {
      ntrig++;
      //The trigger time is the time of the first hit above threshold
      std::sort(digit_times.begin(), digit_times.end());
      triggertime = digit_times[this_ndigitsThreshold];
      triggertime -= (int)triggertime % 5;
      TriggerTimes.push_back(triggertime);
      TriggerTypes.push_back(this_triggerType);
      TriggerInfos.push_back(std::vector<Float_t>(1, n_digits));
      triggerfound = true;
    }

#ifdef WCSIMWCTRIGGER_VERBOSE
    if(n_digits)
      G4cout << n_digits << " digits found in 200nsec trigger window ["
	     << window_start_time << ", " << window_start_time + ndigitsWindow
	     << "]. Threshold is: " << this_ndigitsThreshold << G4endl;
#endif

    //move onto the next go through the timing loop
    if(triggerfound) {
      window_start_time = triggertime + GetPostTriggerWindow(TriggerTypes.back());
    }//triggerfound
    else {
      window_start_time += window_step_size;
    }

    //shorten the loop using the time of the last hit
    if(first_loop) {
      int newend = TMath::Max((float)0, lasthit - ((ndigitsWindow - 10)));
#ifdef WCSIMWCTRIGGER_VERBOSE
      G4cout << "Last hit found to be at " << lasthit
	     << ". Changing window_end_time from " << window_end_time
	     << " to " << newend
	     << G4endl;
#endif
      window_end_time = newend;
      first_loop = false;
    }
  }
  
  G4cout << "Found " << ntrig << " NDigit triggers" << G4endl;
  //call FillDigitsCollection() whether any triggers are found or not
  // (what's saved depends on saveFailuresMode)
  FillDigitsCollection(WCDCPMT, remove_hits, this_triggerType);
}

std::vector<int> WCSimWCTriggerBase::FindPMTNearestNeighbours(int ipmt)
{
  if(myPMTs == NULL) {
    myPMTs = myDetector->Get_Pmts();
  }
  WCSimPmtInfo * thisPMT = myPMTs->at(ipmt);
  int thisTubeID = thisPMT->Get_tubeid();
  double thisX   = thisPMT->Get_transx();
  double thisY   = thisPMT->Get_transy();
  double thisZ   = thisPMT->Get_transz();
  //ipmt is the position in the vector. Runs 0->nPMTs-1
  //tubeid is the actual ID of the PMT. Runs 1->nPMTs
  if((ipmt + 1) != thisTubeID)
    G4cerr << "PMT ID is not the expected one!"
	   << " Vector position + 1 " << ipmt+1
	   << " PMT ID " << thisTubeID << G4endl;

  //loop over ALL the other PMTs & save the distance to the current PMT
  std::vector< std::pair<double, int> > distances;
  for(unsigned int i = 0; i < nPMTs; i++) {
    if(i == (unsigned int)ipmt)
      continue;
    WCSimPmtInfo * PMT = myPMTs->at(i);
    double X   = PMT->Get_transx();
    double Y   = PMT->Get_transy();
    double Z   = PMT->Get_transz();
    double distance = sqrt(((thisX - X) * (thisX - X))
			   + ((thisY - Y) * (thisY - Y))
			   + ((thisZ - Z) * (thisZ - Z)));
    distances.push_back(std::pair<double, int>(distance, PMT->Get_tubeid()));
  }//ipmt

  //sorts by pair.first -> have a vector of pairs ordered by distance
  std::sort(distances.begin(), distances.end());

  //get the N nearest neighbours
  std::vector<int> thisNeighbours;
  for(int i = 0; i < localNHitsNeighbours; i++) {
    thisNeighbours.push_back(distances.at(i).second);
  }
  
#ifdef WCSIMWCTRIGGER_PMT_NEIGHBOURS_VERBOSE
  G4cout << "PMT " << ipmt << " has ID " << thisTubeID
	 << " position " 
	 << thisX << " "
	 << thisY << " "
	 << thisZ << G4endl;
#endif
  return thisNeighbours;
}

void WCSimWCTriggerBase::FindAllPMTNearestNeighbours()
{
  if(myPMTs == NULL) {
    myPMTs = myDetector->Get_Pmts();
  }
  for(unsigned int ipmt = 0; ipmt < nPMTs; ipmt++) {
    if(ipmt % (nPMTs/20 + 1) == 0) {
      G4cout << "WCSimWCTriggerBase::FindAllPMTNearestNeighbours at "
	     << ipmt / (float)nPMTs * 100 << "%"
	     << " (" << ipmt << " out of " << nPMTs << ")"
	     << G4endl;
    }
    std::vector<int> neighbours = FindPMTNearestNeighbours(ipmt);
    pmtNeighbours.push_back(neighbours);
  }

#ifdef WCSIMWCTRIGGER_PMT_NEIGHBOURS_VERBOSE
  for(unsigned int i = 0; i < pmtNeighbours.size(); i++) {
    G4cout << "PMT ID " << i+1 << " has neighbours";
    for(unsigned int j = 0; j < pmtNeighbours.at(i).size(); j++) {
      G4cout << " " << pmtNeighbours.at(i).at(j);
    }//j
    G4cout << G4endl;
  }//i
#endif
}

void WCSimWCTriggerBase::AlgNHitsThenLocalNHits(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits)
{
  //Now we will try to find triggers
  //loop over PMTs, and Digits in each PMT.  If nhits > Threshhold in a time window, then we have a trigger

  int ntrig = 0;
  int window_start_time = 0;
  int window_end_time   = WCSimWCTriggerBase::LongTime - ndigitsWindow;
  int window_step_size  = 5; //step the search window along this amount if no trigger is found
  float lasthit;
  std::vector<int> digit_times;
  bool first_loop = true;

  G4cout << "WCSimWCTriggerBase::AlgNHitsThenLocalNHits. Number of entries in input digit collection: " << WCDCPMT->entries() << G4endl;
#ifdef WCSIMWCTRIGGER_VERBOSE
  int temp_total_pe = 0;
  for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
    temp_total_pe += (*WCDCPMT)[i]->GetTotalPe();
  }
  G4cout << "WCSimWCTriggerBase::AlgNHitsThenLocalNHits. " << temp_total_pe << " total p.e. input" << G4endl;
#endif

  // the upper time limit is set to the final possible full trigger window
  while(window_start_time <= window_end_time) {
    int n_digits = 0, local_n_digits = 0;
    float triggertime; //save each digit time, because the trigger time is the time of the first hit above threshold
    bool triggerfound = false;
    digit_times.clear();
    memset(localNHitsHits, 0, nPMTs * sizeof(int));
    /*
    for(unsigned int i = 0; i < nPMTs; i++)
      if(localNHitsHits[i] != 0)
	G4cerr << "memset didn't work! " << i << " " << localNHitsHits[i] << G4endl;
    */
    std::vector<int> pmts_hit;
    
    //Loop over each PMT
    for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
      int tubeid = (*WCDCPMT)[i]->GetTubeID();
      bool tubeidsaved = false;
      //Loop over each Digit in this PMT
      for ( G4int ip = 0 ; ip < (*WCDCPMT)[i]->GetTotalPe() ; ip++) {
	int digit_time = (*WCDCPMT)[i]->GetTime(ip);
	//hit in trigger window?
	if(digit_time >= window_start_time) {
	  //count the digits for the NHits
	  if(digit_time <= (window_start_time + ndigitsWindow)) {
	    n_digits++;
	    digit_times.push_back(digit_time);
	  }
	  //save the digits & tubes for the local NHits
	  if(digit_time <= (window_start_time + localNHitsWindow)) {
	    local_n_digits++;
	    localNHitsHits[tubeid-1]++; //tubeID goes from 1->nPMTs
	    if(!tubeidsaved) {
	      pmts_hit.push_back(tubeid);
	      tubeidsaved = true;
	    }
	  }
	}
	//G4cout << digit_time << G4endl;
	//get the time of the last hit (to make the loop shorter)
	if(first_loop && (digit_time > lasthit))
	  lasthit = digit_time;
      }//loop over Digits
    }//loop over PMTs

    //if over threshold, issue trigger
    if(n_digits > ndigitsThreshold) {
      ntrig++;
      //The trigger time is the time of the first hit above threshold
      std::sort(digit_times.begin(), digit_times.end());
      triggertime = digit_times[ndigitsThreshold];
      triggertime -= (int)triggertime % 5;
      TriggerTimes.push_back(triggertime);
      TriggerTypes.push_back(kTriggerNDigits);
      TriggerInfos.push_back(std::vector<Float_t>(1, n_digits));
      triggerfound = true;
      G4cout << "Found NDigits trigger with " << n_digits
	     << " digits in the " << ndigitsWindow
	     << " ns trigger decision window" << G4endl;
    }//NHits trigger

    //NHits failed. Try local NHits
    //no point looking at all PMTs if there aren't enough hits in the whole detector in the time window
    else if(local_n_digits > localNHitsThreshold) {
      std::vector<Float_t> local_info;

      //loop over all PMTs with hits
      for(unsigned int ip = 0; ip < pmts_hit.size(); ip++) {
	int this_pmtid = pmts_hit[ip];
	int nlocal = localNHitsHits[this_pmtid - 1];
#ifdef WCSIMWCTRIGGER_VERBOSE
	G4cout << "PMT ID " << this_pmtid << " has " << nlocal << " hits and neighbours";
#endif
	//add the neighbours hits
	std::vector<int> thisNeighbours = pmtNeighbours.at(this_pmtid - 1);
	for(int in = 0; in < localNHitsNeighbours; in++) {
	  nlocal += localNHitsHits[thisNeighbours[in] - 1];
#ifdef WCSIMWCTRIGGER_VERBOSE
	  G4cout << " " << thisNeighbours[in] << "(" << localNHitsHits[thisNeighbours[in] - 1] << ")";
#endif
	}//in
#ifdef WCSIMWCTRIGGER_VERBOSE
	G4cout << " total = " << nlocal  << G4endl;
#endif
	//if over threshold, issue trigger
	if(nlocal > localNHitsThreshold) {
	  local_info.push_back(nlocal); //the number of hits locally to the tube
	  local_info.push_back(ip+1);   //the tube ID of the tube that fired the trigger
	  triggerfound = true;
	  G4cout << "Found LocalNDigits trigger with " << nlocal
		 << " neighbours of PMT " << this_pmtid
		 << " with digits in the " << localNHitsWindow
		 << " ns trigger decision window" << G4endl;
	}//trigger found on a specific PMT
      }//ip

      if(local_info.size()) {
	ntrig++;
	//The trigger time is the time in the middle of the local NHits trigger window                                                                                                                                                         
	//TODO potentially make this something dependant on digit times                                                                                                                                                                        
	triggertime = window_start_time + ((double)localNHitsWindow / 2);
	triggertime -= (int)triggertime % 5;
	TriggerTimes.push_back(triggertime);
	TriggerTypes.push_back(kTriggerLocalNHits);
	TriggerInfos.push_back(local_info);
	triggerfound = true;
      }//trigger(s) found on any PMT
    }//local NHits trigger

    //Cheat sheet
    // arrays/vectors that run from 0 to npmts-1
    //pmts_hit, localNHitsHits, pmtNeighbours
    // arrays/vectors that run from 1 to npmts
    //thisNeighbours (an element of pmtNeighbours)

#ifdef WCSIMWCTRIGGER_VERBOSE
    if(n_digits)
      G4cout << n_digits << " digits found in 200nsec trigger window ["
	     << window_start_time << ", " << window_start_time + ndigitsWindow
	     << "]. Threshold is: " << ndigitsThreshold << G4endl;
#endif

    //move onto the next go through the timing loop
    if(triggerfound) {
      window_start_time = triggertime + GetPostTriggerWindow(TriggerTypes.back());
    }//triggerfound
    else {
      window_start_time += window_step_size;
    }

    //shorten the loop using the time of the last hit
    if(first_loop) {
      int newend = TMath::Max((float)0, lasthit - ((ndigitsWindow - 10)));
#ifdef WCSIMWCTRIGGER_VERBOSE
      G4cout << "Last hit found to be at " << lasthit
	     << ". Changing window_end_time from " << window_end_time
	     << " to " << newend
	     << G4endl;
#endif
      window_end_time = newend;
      first_loop = false;
    }
  }
  
  G4cout << "Found " << ntrig << " NHit or LocalNHit triggers" << G4endl;
  //call FillDigitsCollection() whether any triggers are found or not
  // (what's saved depends on saveFailuresMode)
  FillDigitsCollection(WCDCPMT, remove_hits, kTriggerUndefined);
}

void WCSimWCTriggerBase::FillDigitsCollection(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits, TriggerType_t save_triggerType)
{
  //Adds the digits within the trigger window to the output WCSimWCDigitsCollection
  // optionally removes digits from the input digits collection (when running different Alg* methods concurently) 
  // so they are not used in subsequent trigger decisions or saved twice
  //Also, only save digits of a specific type (again for when running different Alg* methods concurently)

  // Add dummy triggers / exit without saving triggers as required
  //
  //saveFailuresMode = 0 - save only triggered events
  //saveFailuresMode = 1 - save both triggered & not triggered events
  //saveFailuresMode = 2 - save only not triggered events
  if(TriggerTimes.size()) {
    if(saveFailuresMode == 2)
      return;
  }
  else {
    if(saveFailuresMode == 0)
      return;
    TriggerTypes.push_back(kTriggerFailure);
    TriggerTimes.push_back(saveFailuresTime);
    TriggerInfos.push_back(std::vector<Float_t>(1, -1));
    save_triggerType = kTriggerFailure;
  }

  //Get the PMT info for hit time smearing
  G4String WCIDCollectionName = myDetector->GetIDCollectionName();
  WCSimPMTObject * PMT = myDetector->GetPMTPointer(WCIDCollectionName);

  //Loop over trigger times
  for(unsigned int itrigger = 0; itrigger < TriggerTimes.size(); itrigger++) {
    TriggerType_t triggertype = TriggerTypes[itrigger];
    //check if we've already saved this trigger
    if(triggertype != save_triggerType && save_triggerType != kTriggerUndefined)
      continue;
    float         triggertime = TriggerTimes[itrigger];
    std::vector<Float_t> triggerinfo = TriggerInfos[itrigger];

    //these are the boundary of the trigger gate: we want to add all digits within these bounds to the output collection
    float lowerbound = triggertime + GetPreTriggerWindow(triggertype);
    float upperbound = triggertime + GetPostTriggerWindow(triggertype);

#ifdef WCSIMWCTRIGGER_VERBOSE
    G4cout << "Saving trigger " << itrigger << " of type " << WCSimEnumerations::EnumAsString(triggertype)
	   << " in time range [" << lowerbound << ", " << upperbound << "]"
	   << " with trigger time " << triggertime
	   << " and additional trigger info";
    for(std::vector<Float_t>::iterator it = triggerinfo.begin(); it != triggerinfo.end(); ++it)
      G4cout << " " << *it;
    G4cout << G4endl;
#endif

    //loop over PMTs
    for (G4int i = 0; i < WCDCPMT->entries(); i++) {
      int tube=(*WCDCPMT)[i]->GetTubeID();
      //loop over digits in this PMT
      for ( G4int ip = 0; ip < (*WCDCPMT)[i]->GetTotalPe(); ip++){
	int digit_time  = (*WCDCPMT)[i]->GetTime(ip);
	if(digit_time >= lowerbound && digit_time <= upperbound) {
	  //hit in event window
	  //add it to DigitsCollection

	  //first smear the charge & time
	  float peSmeared = (*WCDCPMT)[i]->GetPe(ip);
	  float Q = (peSmeared > 0.5) ? peSmeared : 0.5;
	  G4double digihittime = -triggertime
	    + WCSimWCTriggerBase::offset
	    + digit_time
	    + PMT->HitTimeSmearing(Q);
	  if(digihittime < 0)
	    continue;

	  //get the composition information for the triggered digit
	  //WCDCPMT stores this information in pairs of (digit id, photon id)
	  //need to loop to ensure we get all the photons associated with the current digit (digit ip)
	  std::vector< std::pair<int,int> > digitized_composition = (*WCDCPMT)[i]->GetDigiCompositionInfo();
	  std::vector<int> triggered_composition;
	  for(std::vector< std::pair<int,int> >::iterator it = digitized_composition.begin(); it != digitized_composition.end(); ++it) {
	    if((*it).first == ip) {
	      triggered_composition.push_back((*it).second);
	    }
	    else if ((*it).first > ip)
	      break;
	  }//loop over digitized_composition

#ifdef WCSIMWCTRIGGER_VERBOSE
	  G4cout << "Saving digit on PMT " << tube
		 << " time " << digihittime
		 << " pe "   << peSmeared
		 << " digicomp";
	  for(unsigned int iv = 0; iv < triggered_composition.size(); iv++)
	    G4cout << " " << triggered_composition[iv];
	  G4cout << G4endl;
#endif
	  assert(triggered_composition.size());

	  //add hit
	  if ( DigiHitMap[tube] == 0) {
	    //this PMT has no digits saved yet; create a new WCSimWCDigiTrigger
	    WCSimWCDigiTrigger* Digi = new WCSimWCDigiTrigger();
	    Digi->SetTubeID(tube);
	    Digi->AddGate  (itrigger);
	    Digi->SetTime  (itrigger,digihittime);
	    Digi->SetPe    (itrigger,peSmeared);
	    Digi->AddPe    ();
	    Digi->AddDigiCompositionInfo(itrigger,triggered_composition);
	    DigiHitMap[tube] = DigitsCollection->insert(Digi);
	  }
	  else {
	    //this PMT has digits saved already; add information to the WCSimWCDigiTrigger
	    (*DigitsCollection)[DigiHitMap[tube]-1]->AddGate(itrigger);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->SetTime(itrigger, digihittime);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->SetPe  (itrigger, peSmeared);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->AddPe  ();
	    (*DigitsCollection)[DigiHitMap[tube]-1]->AddDigiCompositionInfo(itrigger,triggered_composition);
	  }
	  if(remove_hits)
	    (*WCDCPMT)[i]->RemoveDigitizedGate(ip);

	  //we've found a digit on this PMT. If we're restricting to just 1 digit per trigger window (e.g. SKI)
	  // then ignore later digits and break. This takes us to the next PMT
	  if(!multiDigitsPerTrigger)
	    break;
	}//digits within trigger window
      }//loop over Digits
    }//loop over PMTs
  }//loop over Triggers
  G4cout << "WCSimWCTriggerBase::FillDigitsCollection. Number of entries in output digit collection: " << DigitsCollection->entries() << G4endl;

}




// *******************************************
// CONTAINER CLASS
// *******************************************

G4Allocator<WCSimWCDigiTrigger> WCSimWCDigiTriggerAllocator;

WCSimWCDigiTrigger::WCSimWCDigiTrigger()
{
  Gates.clear();
  tubeID = 0;
  pe.clear();
  time.clear();
  fDigiComp.clear();
  totalPe = 0;
}

WCSimWCDigiTrigger::~WCSimWCDigiTrigger(){;}

WCSimWCDigiTrigger::WCSimWCDigiTrigger(const WCSimWCDigiTrigger& right)
  :G4VDigi()
{
  // in principle assignment = is defined for containers...
  Gates = right.Gates;
  tubeID = right.tubeID;
  pe     = right.pe;
  time   = right.time;
  fDigiComp = right.fDigiComp;
  totalPe = right.totalPe;
}

const WCSimWCDigiTrigger& WCSimWCDigiTrigger::operator=(const WCSimWCDigiTrigger& right)
{
  Gates = right.Gates;
  tubeID = right.tubeID;
  pe     = right.pe;
  time   = right.time;
  fDigiComp = right.fDigiComp;
  totalPe = right.totalPe;
  return *this;
}

void WCSimWCDigiTrigger::Print()
{
  G4cout << "TubeID: " << tubeID
         << ", Number of Gates: " << NumberOfGates()
	 << G4endl;
  std::multimap<int,float>::iterator it_pe   = pe.begin();
  std::multimap<int,float>::iterator it_time = time.begin();
  for( ; it_pe != pe.end(), it_time != time.end(); ++it_pe, ++it_time) {
    if(it_pe->first != it_time->first) {
      G4cerr << "WCSimWCDigiTrigger::Print() pe and time gate counters disagree!" << G4endl;
      exit(-1);
    }
    G4cout  << "Gate = " << it_pe->first
            << " PE: "   << it_pe->second
            << " Time: " << it_time->second
	    << G4endl;
  }
}



// *******************************************
// DERIVED CLASS
// *******************************************

WCSimWCTriggerNDigits::WCSimWCTriggerNDigits(G4String name,
					 WCSimDetectorConstruction* myDetector,
					 WCSimWCDAQMessenger* myMessenger)
  :WCSimWCTriggerBase(name, myDetector, myMessenger)
{
  triggerClassName = "NDigits";
  GetVariables();
}

WCSimWCTriggerNDigits::~WCSimWCTriggerNDigits()
{
}

void WCSimWCTriggerNDigits::DoTheWork(WCSimWCDigitsCollection* WCDCPMT) {
  //Apply an NDigits trigger
  bool remove_hits = false;
  AlgNDigits(WCDCPMT, remove_hits);
}

// *******************************************
// DERIVED CLASS
// *******************************************

WCSimWCTriggerNDigits2::WCSimWCTriggerNDigits2(G4String name,
					 WCSimDetectorConstruction* myDetector,
					 WCSimWCDAQMessenger* myMessenger)
  :WCSimWCTriggerBase(name, myDetector, myMessenger)
{
  triggerClassName = "NDigits2";
  GetVariables();
}

WCSimWCTriggerNDigits2::~WCSimWCTriggerNDigits2(){
}


void WCSimWCTriggerNDigits2::DoTheWork(WCSimWCDigitsCollection* WCDCPMT) {
  //This calls 2 trigger algorithms; the second algorithm is called on hits that failed the first algorithm
  //(for a second trigger working on hits that passed a pretrigger, FillDigitsCollection() needs to have a new option)

  //Make a copy of the input DigitsCollection, so we can remove hits from the copy
  WCSimWCDigitsCollection* WCDCPMTCopy = new WCSimWCDigitsCollection(*WCDCPMT);
  
  //Apply an NDigits trigger
  bool remove_hits = true;
  AlgNDigits(WCDCPMTCopy, remove_hits);

  //Apply an NDigits trigger with a lower threshold & different saved trigger type
  remove_hits = false;
  bool ndigits_test = true;
  AlgNDigits(WCDCPMTCopy, remove_hits, ndigits_test);
}



// *******************************************
// DERIVED CLASS
// *******************************************

WCSimWCTriggerNHitsThenLocalNHits::WCSimWCTriggerNHitsThenLocalNHits(G4String name,
								     WCSimDetectorConstruction* myDetector,
								     WCSimWCDAQMessenger* myMessenger)
  :WCSimWCTriggerBase(name, myDetector, myMessenger)
{
  triggerClassName = "NHitsThenLocalNHits";
  GetVariables();

  FindAllPMTNearestNeighbours();
  //reserve an array to store the number of digits on each PMT
  localNHitsHits = new int[nPMTs];
  //and fill it with 0s
  memset(localNHitsHits, 0, nPMTs * sizeof(int));
}

WCSimWCTriggerNHitsThenLocalNHits::~WCSimWCTriggerNHitsThenLocalNHits()
{
}

void WCSimWCTriggerNHitsThenLocalNHits::DoTheWork(WCSimWCDigitsCollection* WCDCPMT)
{
  //Apply an NHitsThenLocalNHits trigger
  bool remove_hits = false;
  AlgNHitsThenLocalNHits(WCDCPMT, remove_hits);
}
