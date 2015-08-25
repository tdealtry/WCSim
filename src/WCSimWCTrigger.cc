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
#ifndef WCSIMWCTRIGGERBASE_PMT_NEIGHBOURS_VERBOSE
//#define WCSIMWCTRIGGERBASE_PMT_NEIGHBOURS_VERBOSE
#endif

const double WCSimWCTriggerBase::offset = 950.0 ; // ns. apply offset to the digit time
const double WCSimWCTriggerBase::eventgateup = 950.0 ; // ns. save eventgateup ns after the trigger time
const double WCSimWCTriggerBase::eventgatedown = -400.0 ; // ns. save eventgateup ns before the trigger time
const double WCSimWCTriggerBase::LongTime = 100000.0 ; // ns = 0.1ms. event time


WCSimWCTriggerBase::WCSimWCTriggerBase(G4String name,
				       WCSimDetectorConstruction* inDetector,
				       WCSimWCDAQMessenger* myMessenger)
  :G4VDigitizerModule(name), myDetector(inDetector), triggerClassName("")
{
  G4String colName = "WCDigitizedCollection";
  collectionName.push_back(colName);

  ReInitialize();

  if(myMessenger != NULL) {
    DAQMessenger = myMessenger;
    DAQMessenger->TellMeAboutTheTrigger(this);
    DAQMessenger->SetTriggerOptions();
  }
  else {
    G4cerr << "WCSimWCDAQMessenger pointer is NULL when passed to WCSimWCTriggerBase constructor. Exiting..."
           << G4endl;
    exit(-1);
  }
  digitizeCalled = false;
}

WCSimWCTriggerBase::~WCSimWCTriggerBase()
{
  if(triggerClassName.compare("NHitsThenLocalNHits") == 0)
    delete localNHitsHits;
}

int WCSimWCTriggerBase::CalculateAverageDarkNoiseOccupancy(int npmts, int window)
{
  double trigger_window_seconds = window * 1E-9;
  double dark_rate_Hz = PMTDarkRate * 1000;
  double average_occupancy = dark_rate_Hz * trigger_window_seconds * npmts;
  
  G4cout << "Average number of PMTs active in a " << window
	 << "ns window with a dark noise rate of " << PMTDarkRate
	 << "kHz is " << average_occupancy
	 << " (" << npmts << " total PMTs)"
	 << G4endl;
  return round(average_occupancy);
}

void WCSimWCTriggerBase::Digitize()
{
  if(!digitizeCalled) {
    nPMTs = this->myDetector->GetTotalNumPmts();
    if(nhitsAdjustForNoise) {
      int nnoisehits = CalculateAverageDarkNoiseOccupancy(nPMTs, nhitsWindow);
      G4cout << "Updating the NHits threshold, from " << nhitsThreshold
	     << " to " << nhitsThreshold + nnoisehits << G4endl;
      nhitsThreshold += nnoisehits;
    }
    if(localNHitsAdjustForNoise) {
      int nnoisehits = CalculateAverageDarkNoiseOccupancy(localNHitsNeighbours + 1, localNHitsWindow);
      G4cout << "Updating the NHits threshold, from " << localNHitsThreshold
	     << " to " << localNHitsThreshold + nnoisehits << G4endl;
      localNHitsThreshold += nnoisehits;
    }
    if(triggerClassName.compare("NHitsThenLocalNHits") == 0) {
      FindAllPMTNearestNeighbours();
      //reserve an array to store the number of digits on each PMT
      localNHitsHits = new int[nPMTs];
      //and fill it with 0s
      memset(localNHitsHits, 0, nPMTs);
    }
    digitizeCalled = true;
  }
    

  //Input is collection of all digitized hits that passed the threshold
  //Output is all digitized hits which pass the trigger
  
  ReInitialize();

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
#ifdef WCSIMWCTRIGGER_VERBOSE
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

#ifdef WCSIMWCTRIGGER_VERBOSE
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
#ifdef WCSIMWCTRIGGER_VERBOSE
      G4cout << "Last hit found to be at " << lasthit
	     << ". Changing window_end_time from " << window_end_time
	     << " to " << lasthit - (nhitsWindow - 10)
	     << G4endl;
#endif
      window_end_time = lasthit - (nhitsWindow - 10);
      first_loop = false;
    }
  }
  
  G4cout << "Found " << ntrig << " NHit triggers" << G4endl;
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
    if(i == (unsigned int)ipmt) continue;
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
  
#ifdef WCSIMWCTRIGGERBASE_PMT_NEIGHBOURS_VERBOSE
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

#ifdef WCSIMWCTRIGGERBASE_PMT_NEIGHBOURS_VERBOSE
  for(unsigned int i = 0; i < pmtNeighbours.size(); i++) {
    G4cout << "PMT ID " << i << " has neighbours";
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
  int window_end_time   = WCSimWCTriggerBase::LongTime - nhitsWindow;
  int window_step_size  = 5; //step the search window along this amount if no trigger is found
  float lasthit;
  std::vector<int> digit_times;
  bool first_loop = true;

  G4cout << "WCSimWCTriggerBase::AlgNHitsThenLocalNHits. Number of entries in input digit collection: " << WCDCPMT->entries() << G4endl;
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
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
    memset(localNHitsHits, 0, nPMTs);
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
	  if(digit_time <= (window_start_time + nhitsWindow)) {
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
    if(n_digits > nhitsThreshold) {
      ntrig++;
      //The trigger time is the time of the first hit above threshold
      std::sort(digit_times.begin(), digit_times.end());
      triggertime = digit_times[nhitsThreshold];
      triggertime -= (int)triggertime % 5;
      TriggerTimes.push_back(triggertime);
      TriggerTypes.push_back(kTriggerNHits);
      TriggerInfos.push_back(std::vector<Float_t>(1, n_digits));
      triggerfound = true;
    }//NHits trigger

    //NHits failed. Try local NHits
    //no point looking at all PMTs if there aren't enough hits in the time window
    if(local_n_digits > localNHitsThreshold) {
      //loop over all PMTs with hits
      for(unsigned int ip = 0; ip < pmts_hit.size(); ip++) {
	int this_pmtid = pmts_hit[ip];
	int nlocal = localNHitsHits[this_pmtid - 1];
	//add the neighbours hits
	std::vector<int> thisNeighbours = pmtNeighbours.at(this_pmtid - 1);
	for(int in = 0; in < localNHitsNeighbours; in++) {
	  nlocal += localNHitsHits[thisNeighbours[in - 1]];
	}//in
	//if over threshold, issue trigger
	if(nlocal > localNHitsThreshold) {
	  ntrig++;
	  //The trigger time is the time in the middle of the local NHits trigger window
	  //TODO potentially make this something dependant on digit times
	  triggertime = window_start_time + ((double)localNHitsWindow / 2);
	  triggertime -= (int)triggertime % 5;
	  TriggerTimes.push_back(triggertime);
	  TriggerTypes.push_back(kTriggerLocalNHits);
	  std::vector<Float_t> local_info;
	  local_info.push_back(nlocal); //the number of hits locally to the tube
	  local_info.push_back(ip+1);   //the tube ID of the tube that fired the trigger
	  TriggerInfos.push_back(local_info);
	  triggerfound = true;
	  break; //may want to continue loop & find all local triggers instead of breaking
	}//trigger found
      }//ip
    }//local NHits trigger

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
  
  G4cout << "Found " << ntrig << " NHit triggers" << G4endl;
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
    float lowerbound = triggertime + WCSimWCTriggerBase::eventgatedown;
    float upperbound = triggertime + WCSimWCTriggerBase::eventgateup;

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

	  //int parentID    = (*WCDCPMT)[i]->GetParentID(ip);

	  //get the composition information for the triggered digit
	  //WCDCPMT stores this information in pairs of (digit id, photon id)
	  //need to loop to ensure we get all the photons associated with the current digit (digit ip)
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
	    //this PMT has no digits saved yet; create a new WCSimWCDigi
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
	    //this PMT has digits saved already; add information to the WCSimWCDigi
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




// *******************************************
// DERIVED CLASS
// *******************************************

WCSimWCTriggerNHits::WCSimWCTriggerNHits(G4String name,
					 WCSimDetectorConstruction* myDetector,
					 WCSimWCDAQMessenger* myMessenger)
  :WCSimWCTriggerBase(name, myDetector, myMessenger)
{
  triggerClassName = "NHits";
}

WCSimWCTriggerNHits::~WCSimWCTriggerNHits()
{
}

void WCSimWCTriggerNHits::DoTheWork(WCSimWCDigitsCollection* WCDCPMT)
{
  //Apply an NHits trigger
  bool remove_hits = false;
  AlgNHits(WCDCPMT, remove_hits);
}



// *******************************************
// DERIVED CLASS
// *******************************************

WCSimWCTriggerNHits2::WCSimWCTriggerNHits2(G4String name,
					 WCSimDetectorConstruction* myDetector,
					 WCSimWCDAQMessenger* myMessenger)
  :WCSimWCTriggerBase(name, myDetector, myMessenger)
{
  triggerClassName = "NHits2";
}

WCSimWCTriggerNHits2::~WCSimWCTriggerNHits2()
{
}


void WCSimWCTriggerNHits2::DoTheWork(WCSimWCDigitsCollection* WCDCPMT)
{
  //This calls 2 trigger algorithms; the second algorithm is called on hits that failed the first algorithm
  //(for a second trigger working on hits that passed a pretrigger, FillDigitsCollection() needs to have a new option)

  //Make a copy of the input DigitsCollection, so we can remove hits from the copy
  WCSimWCDigitsCollection* WCDCPMTCopy = new WCSimWCDigitsCollection(*WCDCPMT);
  
  //Apply an NHits trigger
  bool remove_hits = true;
  AlgNHits(WCDCPMTCopy, remove_hits);

  //Apply an NHits trigger with a lower threshold & different saved trigger type
  remove_hits = false;
  bool nhits_test = true;
  AlgNHits(WCDCPMTCopy, remove_hits, nhits_test);
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
