#include "WCSimWCTriggerBase.hh"
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

#include <vector>
// for memset
#include <cstring>
#include <iostream>

#ifndef WCSIMWCTRIGGERBASE_VERBOSE
//#define WCSIMWCTRIGGERBASE_VERBOSE
#endif

const double WCSimWCTriggerBase::offset = 950.0 ; // ns. apply offset to the digit time
const double WCSimWCTriggerBase::eventgateup = 950.0 ; // ns. save eventgateup ns after the trigger time
const double WCSimWCTriggerBase::eventgatedown = -400.0 ; // ns. save eventgateup ns before the trigger time
const double WCSimWCTriggerBase::LongTime = 100000.0 ; // ns = 0.1ms. event time


WCSimWCTriggerBase::WCSimWCTriggerBase(G4String name,
				       WCSimDetectorConstruction* myDetector,
				       WCSimWCDAQMessenger* myMessenger)
  :G4VDigitizerModule(name)
{
  G4String colName = "WCDigitizedCollection";
  this->myDetector = myDetector;
  collectionName.push_back(colName);
  DigiHitMap.clear();
  TriggerTimes.clear();
  TriggerTypes.clear(); 
  TriggerInfos.clear(); 

  if(myMessenger) {
    DAQMessenger = myMessenger;
    DAQMessenger->TellMeAboutTheTrigger(this);
    DAQMessenger->TellTrigger();
  }

  digitizeCalled = false;
}
  
WCSimWCTriggerBase::~WCSimWCTriggerBase(){
}

void WCSimWCTriggerBase::AdjustNHitsThresholdForNoise()
{
  int npmts = this->myDetector->GetTotalNumPmts();
  double trigger_window_seconds = nhitsWindow * 1E-9;
  double dark_rate_Hz = PMTDarkRate * 1000;
  double average_occupancy = dark_rate_Hz * trigger_window_seconds * npmts;
  
  G4cout << "Average number of PMTs in detector active in a " << nhitsWindow
	 << "ns window with a dark noise rate of " << PMTDarkRate
	 << "kHz is " << average_occupancy
	 << " (" << npmts << " total PMTs)"
	 << G4endl
	 << "Updating the NHits threshold, from " << nhitsThreshold
	 << " to " << nhitsThreshold + round(average_occupancy) << G4endl;
  nhitsThreshold += round(average_occupancy);
}

void WCSimWCTriggerBase::Digitize()
{
  if(nhitsAdjustForNoise && !digitizeCalled) {
    AdjustNHitsThresholdForNoise();
    digitizeCalled = true;
  }

  //Input is collection of all digitized hits that passed the threshold
  //Output is all digitized hits which pass the trigger
  
  DigiHitMap.clear();
  TriggerTimes.clear();
  TriggerTypes.clear(); 
  TriggerInfos.clear(); 

  //This is the output digit collection
  DigitsCollection = new WCSimWCDigitsCollection ("/WCSim/glassFaceWCPMT",collectionName[0]);

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

void WCSimWCTriggerBase::AlgNHits(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits, bool test)
{

  //if test is true, we run the algorithm with 1/2 the threshold, and kTriggerNHitsTest
  //for testing multiple trigger algorithms at once
  int this_nhitsThreshold = nhitsThreshold;
  TriggerType_t this_triggerType = kTriggerNHits;
  if(test) {
    this_nhitsThreshold /= 2;
    this_triggerType = kTriggerNHitsTest;
  }

  //Now we will try to find triggers
  //loop over PMTs, and Digits in each PMT.  If nhits > Threshhold in a time window, then we have a trigger

  int ntrig = 0;
  int window_start_time = 0;
  int window_end_time   = WCSimWCTriggerBase::LongTime - nhitsWindow;
  int window_step_size  = 5; //step the search window along this amount if no trigger is found
  float lasthit;
  std::vector<int> digit_times;
  bool first_loop = true;

  G4cout << "WCSimWCTriggerBase::AlgNHits. Number of entries in input digit collection: " << WCDCPMT->entries() << G4endl;
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
  int temp_total_pe = 0;
  for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
    temp_total_pe += (*WCDCPMT)[i]->GetTotalPe();
  }
  G4cout << "WCSimWCTriggerBase::AlgNHits. " << temp_total_pe << " total p.e. input" << G4endl;
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
	if(digit_time >= window_start_time && digit_time <= (window_start_time + nhitsWindow)) {
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
    if(n_digits > this_nhitsThreshold) {
      ntrig++;
      //The trigger time is the time of the first hit above threshold
      std::sort(digit_times.begin(), digit_times.end());
      triggertime = digit_times[this_nhitsThreshold];
      triggertime -= (int)triggertime % 5;
      TriggerTimes.push_back(triggertime);
      TriggerTypes.push_back(this_triggerType);
      TriggerInfos.push_back(std::vector<Float_t>(1, n_digits));
      triggerfound = true;
    }

#ifdef WCSIMWCTRIGGERBASE_VERBOSE
    if(n_digits)
      G4cout << n_digits << " digits found in 200nsec trigger window ["
	     << window_start_time << ", " << window_start_time + nhitsWindow
	     << "]. Threshold is: " << this_nhitsThreshold << G4endl;
#endif

    //move onto the next go through the timing loop
    if(triggerfound) {
      window_start_time = triggertime + WCSimWCTriggerBase::eventgateup;
    }//triggerfound
    else {
      window_start_time += window_step_size;
    }

    //shorten the loop using the time of the last hit
    if(first_loop) {
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
      G4cout << "Last hit found to be at " << lasthit
	     << ". Changing window_end_time from " << window_end_time
	     << " to " << lasthit - (nhitsWindow - 10)
	     << G4endl;
#endif
      window_end_time = lasthit - (nhitsWindow - 10);
      first_loop = false;
    }
  }
  
  //call FillDigitsCollection() if at least one trigger was issued
  G4cout << "Found " << ntrig << " NHit triggers" << G4endl;
  FillDigitsCollection(WCDCPMT, remove_hits, this_triggerType);
}

void WCSimWCTriggerBase::AlgNHitsThenITC(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits)
{

  //Now we will try to find triggers
  //loop over PMTs, and Digits in each PMT. 
  // If nhits > Threshhold in a time window, then we have a trigger
  // If this cut fails, attempt an ITC ratio trigger

  int ntrig = 0;
  int window_start_time = 0;
  int window_end_time   = WCSimWCTriggerBase::LongTime - nhitsWindow;
  int window_step_size  = 5; //step the search window along this amount if no trigger is found
  float lasthit;
  std::vector<int> digit_times;
  std::vector<int> digit_times_itc_small, digit_times_itc_large;
  bool first_loop = true;

  G4cout << "WCSimWCTriggerBase::AlgNHits. Number of entries in input digit collection: " << WCDCPMT->entries() << G4endl;
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
  int temp_total_pe = 0;
  for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
    temp_total_pe += (*WCDCPMT)[i]->GetTotalPe();
  }
  G4cout << "WCSimWCTriggerBase::AlgNHits. " << temp_total_pe << " total p.e. input" << G4endl;
#endif

  // the upper time limit is set to the final possible full trigger window
  while(window_start_time <= window_end_time) {
    int n_digits = 0;
    int n_digits_itc_small = 0, n_digits_itc_large = 0;
    float triggertime; //save each digit time, because the trigger time is the time of the first hit above threshold
    bool triggerfound = false;
    digit_times.clear();
    
    //Loop over each PMT & count NDigits in window [window_start_time, window_start_time + nhitsWindow]
    //Also count in two extra windows for the ITC cut
    // [window_start_time, window_start_time + itcSmallWindow]  
    // [window_start_time - itcLargeWindowLow, window_start_time + itcLargeWindowHigh]  
    for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
      //int tube=(*WCDCPMT)[i]->GetTubeID();
      //Loop over each Digit in this PMT
      for ( G4int ip = 0 ; ip < (*WCDCPMT)[i]->GetTotalPe() ; ip++) {
	int digit_time = (*WCDCPMT)[i]->GetTime(ip);
	//hit in trigger window?
	if(digit_time >= window_start_time && digit_time <= (window_start_time + nhitsWindow)) {
	  n_digits++;
	  digit_times.push_back(digit_time);
	}
        //hit in the small ITC window?
        if((digit_time >= window_start_time) && (digit_time <= (window_start_time + itcSmallWindow))) {
          n_digits_itc_small++;
          digit_times_itc_small.push_back(digit_time);
        }
        //hit in the large ITC window?
        if((digit_time >= (window_start_time - itcLargeWindowLow)) && digit_time <= (window_start_time + itcLargeWindowHigh)) {
          n_digits_itc_large++;
          digit_times_itc_large.push_back(digit_time);
        }
	//G4cout << digit_time << G4endl;
	//get the time of the last hit (to make the loop shorter)
	if(first_loop && (digit_time > lasthit))
	  lasthit = digit_time;
      }//loop over Digits
    }//loop over PMTs

    //if over threshold, issue trigger
    if(n_digits > nhitsThreshold) {
      ntrig++;
      //The trigger time is the time of the first hit above threshold
      std::sort(digit_times.begin(), digit_times.end());
      triggertime = digit_times[nhitsThreshold];
      triggertime -= (int)triggertime % 5;
      //save the trigger information
      TriggerTimes.push_back(triggertime);
      TriggerTypes.push_back(kTriggerNHits);
      TriggerInfos.push_back(std::vector<Float_t>(1, n_digits));
      triggerfound = true;
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
      G4cout << EnumAsString(kTriggerNHits) << " trigger passed with time " << triggertime << G4endl;
#endif
    }//NDigits trigger passed

    //The simple NHits trigger hasn't been passed. See if the ITC ratio trigger can be passed
    double itc_ratio = (double)n_digits_itc_small / (double)n_digits_itc_large;
    if(!triggerfound && itc_ratio > itcRatioThreshold) {
      ntrig++;
      //The trigger time is the time of the first hit TEMPORARY
      std::sort(digit_times.begin(), digit_times.end());
      triggertime = digit_times[0];
      triggertime -= (int)triggertime % 5;
      std::vector<Float_t> triggerinfo;
      triggerinfo.push_back(itc_ratio);
      triggerinfo.push_back(n_digits_itc_small);
      triggerinfo.push_back(n_digits_itc_large);
      //save the trigger information
      TriggerTimes.push_back(triggertime);
      TriggerTypes.push_back(kTriggerITCRatio);
      TriggerInfos.push_back(triggerinfo);
      triggerfound = true;
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
      G4cout << EnumAsString(kTriggerITCRatio) << " trigger passed with time " << triggertime << G4endl;
#endif
    }//ITC trigger passed
    

#ifdef WCSIMWCTRIGGERBASE_VERBOSE
    if(n_digits)
      G4cout << n_digits << " digits found in 200nsec trigger window ["
	     << window_start_time << ", " << window_start_time + nhitsWindow
	     << "]. Threshold is: " << nhitsThreshold << G4endl;
#endif

    //move onto the next go through the timing loop
    if(triggerfound) {
      window_start_time = triggertime + WCSimWCTriggerBase::eventgateup;
    }//triggerfound
    else {
      window_start_time += window_step_size;
    }

    //shorten the loop using the time of the last hit
    if(first_loop) {
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
      G4cout << "Last hit found to be at " << lasthit
	     << ". Changing window_end_time from " << window_end_time
	     << " to " << lasthit - (nhitsWindow - 10)
	     << G4endl;
#endif
      window_end_time = lasthit - (nhitsWindow - 10);
      first_loop = false;
    }
  }
  
  //call FillDigitsCollection() if at least one trigger was issued
  G4cout << "Found " << ntrig << " NHitThenITC triggers" << G4endl;
  FillDigitsCollection(WCDCPMT, remove_hits, kTriggerUndefined);
}

void WCSimWCTriggerBase::FillDigitsCollection(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits, TriggerType_t save_triggerType)
{

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

  WCSimPMTObject * PMT = myDetector->GetPMTPointer("glassFaceWCPMT"); //for hit time smearing

  //Loop over trigger times
  for(unsigned int itrigger = 0; itrigger < TriggerTimes.size(); itrigger++) {
    TriggerType_t triggertype = TriggerTypes[itrigger];
    //check if we've already saved this trigger
    if(triggertype != save_triggerType && save_triggerType != kTriggerUndefined)
      continue;
    float         triggertime = TriggerTimes[itrigger];
    std::vector<Float_t> triggerinfo = TriggerInfos[itrigger];
    float lowerbound = triggertime + WCSimWCTriggerBase::eventgatedown;
    float upperbound = triggertime + WCSimWCTriggerBase::eventgateup;

#ifdef WCSIMWCTRIGGERBASE_VERBOSE
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
      //loop over digits
      for ( G4int ip = 0; ip < (*WCDCPMT)[i]->GetTotalPe(); ip++){
	int digit_time  = (*WCDCPMT)[i]->GetTime(ip);
	if(digit_time >= lowerbound && digit_time <= upperbound) {
	  //hit in event window
	  //add it to DigitsCollection
	  float peSmeared = (*WCDCPMT)[i]->GetPe(ip);
	  float Q = (peSmeared > 0.5) ? peSmeared : 0.5;
	  G4double digihittime = -triggertime
	    + WCSimWCTriggerBase::offset
	    + digit_time
	    + PMT->HitTimeSmearing(Q);
	  if(digihittime < 0)
	    continue;

	  //int parentID    = (*WCDCPMT)[i]->GetParentID(ip);
	  std::vector< std::pair<int,int> > digitized_composition = (*WCDCPMT)[i]->GetDigiCompositionInfo();
	  std::vector< std::pair<int,int> > triggered_composition;
	  for(std::vector< std::pair<int,int> >::iterator it = digitized_composition.begin(); it != digitized_composition.end(); ++it) {
	    if((*it).first == ip) {
	      triggered_composition.push_back(std::make_pair(itrigger, (*it).second));
	    }
	    else if ((*it).first > ip)
	      break;
	  }//loop over digitized_composition
	  //add hit
	  if ( DigiHitMap[tube] == 0) {
	    WCSimWCDigi* Digi = new WCSimWCDigi();
	    Digi->SetTubeID(tube);
	    //Digi->AddParentID(parentID);
	    Digi->AddGate  (itrigger,triggertime);
	    Digi->SetTime  (itrigger,digihittime);
	    Digi->SetPe    (itrigger,peSmeared);
	    Digi->AddPe    (digihittime);
	    Digi->AddDigiCompositionInfo(triggered_composition);
	    DigiHitMap[tube] = DigitsCollection->insert(Digi);
	  }
	  else {
	    //(*DigitsCollection)[DigiHitMap[tube]-1]->AddParentID(parentID);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->AddGate(itrigger, triggertime);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->SetTime(itrigger, digihittime);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->SetPe  (itrigger, peSmeared);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->AddPe  (digihittime);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->AddDigiCompositionInfo(triggered_composition);
	  }
	  if(remove_hits)
	    (*WCDCPMT)[i]->RemoveDigitizedGate(ip);
	}//digits within trigger window
      }//loop over Digits
    }//loop over PMTs
  }//loop over Triggers
  G4cout << "WCSimWCTriggerBase::FillDigitsCollection. Number of entries in output digit collection: " << DigitsCollection->entries() << G4endl;

}
