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
#include "TVector3.h"
#include "TAxis.h"

#include <vector>
// for memset
#include <cstring>
#include <fstream>

#define WCSIMWCTRIGGER_DEBUG_GEOM_INFO 1
#if WCSIMWCTRIGGER_DEBUG_GEOM_INFO >= 1
#include "TFile.h"
#include "TTree.h"
#endif


// *******************************************
// BASE CLASS
// *******************************************

#ifndef WCSIMWCTRIGGER_VERBOSE
//#define WCSIMWCTRIGGER_VERBOSE
#endif
#ifndef WCSIMWCTRIGGER_PMT_NEIGHBOURS_VERBOSE
//#define WCSIMWCTRIGGER_PMT_NEIGHBOURS_VERBOSE
#endif
#ifndef WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE
#define WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE 1
#endif





const double WCSimWCTriggerBase::offset = 950.0; // ns. apply offset to the digit time
const double WCSimWCTriggerBase::LongTime = 1E6; // ns = 1ms. event time


WCSimWCTriggerBase::WCSimWCTriggerBase(G4String name,
				       WCSimDetectorConstruction* inDetector,
				       WCSimWCDAQMessenger* myMessenger)
  :G4VDigitizerModule(name), myPMTs(NULL), DAQMessenger(myMessenger), myDetector(inDetector), triggerClassName("")
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
  if(triggerClassName.compare("NHitsThenRegions") == 0) {
    delete localNHitsHits;
    delete regionsHits;
  }
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
	 << "Using side PMT regions with " << regionsNBinsP << " phi bins and " << regionsNBinsZ << " Z bins" << G4endl
	 << "Using top/bottom PMT regions with " << regionsNCentralSectors << " central sectors and " 
	 <<  regionsNRings << " rings with " << regionsNRingSectors << " times more sectors per ring " << G4endl
	 << "Using ITC ratio threshold " << itcRatioThreshold << G4endl
	 << "Using ITC ratio small window " << itcSmallWindow << G4endl
	 << "Using ITC ratio large window (lo) " << itcLargeWindowLow << G4endl
	 << "Using ITC ratio large window (hi) " << itcLargeWindowHigh << G4endl
	 << "Using SaveFailures event pretrigger window " << saveFailuresPreTriggerWindow << " ns" << G4endl
	 << "Using SaveFailures event posttrigger window " << saveFailuresPostTriggerWindow << " ns" << G4endl
	 << (writeGeom ? "Will write out geometry information for triggering to a file, then exit" : "") << G4endl;

  if(writeGeom) {
    WriteGeomInfo();
    exit(0);
  }
  else {
    ReadGeomInfo();
  }
  if(saveFailuresMode == 0)
    G4cout << "Saving only triggered digits" << G4endl;
  else if(saveFailuresMode == 1)
    G4cout << "Saving both triggered and not-triggered digits" << G4endl;
  else if(saveFailuresMode == 2)
    G4cout << "Saving only not-triggered digits" << G4endl;
  if(saveFailuresMode > 0)
    G4cout << "Using SaveFailures trigger time" << saveFailuresTime << " ns" << G4endl
	   << "Using SaveFailures event pretrigger window " << saveFailuresPreTriggerWindow << " ns" << G4endl
	   << "Using SaveFailures event posttrigger window " << saveFailuresPostTriggerWindow << " ns" << G4endl;
}

int WCSimWCTriggerBase::GetPreTriggerWindow(TriggerType_t t)
{
  switch(t) {
  case kTriggerNDigits:
  case kTriggerNDigitsTest:
  case kTriggerLocalNHits:
  case kTriggerITCRatio:
  case kTriggerRegions:
  case kTriggerAnisotropy:
    return ndigitsPreTriggerWindow;
    break;
  case kTriggerFailure:
    return saveFailuresPreTriggerWindow;
    break;
  case kTriggerNoTrig:
    return - 2147483647 / 2;
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
  case kTriggerLocalNHits:
  case kTriggerITCRatio:
  case kTriggerRegions:
  case kTriggerAnisotropy:
    return ndigitsPostTriggerWindow;
    break;
  case kTriggerFailure:
    return saveFailuresPostTriggerWindow;
    break;
  case kTriggerNoTrig:
    return 2147483647 / 2;
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
      G4cout << n_digits << " digits found in " << ndigitsWindow << " nsec trigger window ["
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
  }//while <= window_end_time
  
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

/*
std::vector<int> WCSimWCTriggerBase::FindRegionsNearestNeighbours(int ir)
{
  int thisTubeID = pmtBlocks[ir].first;
  TVector3 thisV = pmtBlocks[ir].second;
    
  //loop over ALL the other regions & save the distance to the current region
  std::vector< std::pair<double, int> > distances;
  for(unsigned int i = 0; i < pmtBlocks.size(); i++) {
    if(i == (unsigned int)ir)
      continue;
    TVector3 V = pmtBlocks[i].second;
    double distance = sqrt(((thisV.X() - V.X()) * (thisV.X() - V.X())) +
			   ((thisV.Y() - V.Y()) * (thisV.Y() - V.Y())) +
			   ((thisZ.Y() - V.Z()) * (thisZ.Y() - V.Z())));
    distances.push_back(std::pair<double, int>(distance, i));
  }//ipmt

  //sorts by pair.first -> have a vector of pairs ordered by distance
  std::sort(distances.begin(), distances.end());

  /*
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
*/

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
  }//ipmt

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

/*
void WCSimWCTriggerBase::FindAllRegionsNearestNeighbours()
{
  for(int ir = 0; ir < pmtBlocks.size(); ir++) {
    std::vector<int> neighbours = FindRegionNearestNeighbours(ir);
    regionsNeighbours.push_back(neighbours);
  }//ir

#ifdef WCSIMWCTRIGGER_PMT_NEIGHBOURS_VERBOSE
  for(unsigned int i = 0; i < regionsNeighbours.size(); i++) {
    G4cout << "REGIONS ID " << i << " has neighbours";
    for(unsigned int j = 0; j < regionsNeighbours.at(i).size(); j++) {
      G4cout << " " << regionsNeighbours.at(i).at(j);
    }//j
    G4cout << G4endl;
  }//i
#endif
}
*/

void WCSimWCTriggerBase::PopulatePMTAreas()
{
  const double overlap = 0;

  pmtBlocks.clear();
  if(myPMTs == NULL) {
    myPMTs = myDetector->Get_Pmts();
  }

  //Get the position of all PMTs in the detector
  std::vector<double> vR, vP, vZ;
  std::vector<TVector3> allPMTs;
  double RMax = 0;
  for(unsigned int ipmt = 0; ipmt < nPMTs; ipmt++) {
    if(ipmt % (nPMTs/10 + 1) == 0) {
      G4cout << "WCSimWCTriggerBase::PopulatePMTAreas() at "
             << ipmt / (float)nPMTs * 100 << "%"
             << " (" << ipmt << " out of " << nPMTs << ")"
             << G4endl;
    }
    WCSimPmtInfo * thisPMT = myPMTs->at(ipmt);
    int thisTubeID = thisPMT->Get_tubeid();
    double thisX   = thisPMT->Get_transx();
    double thisY   = thisPMT->Get_transy();
    double thisZ   = thisPMT->Get_transz();
    if(TMath::Abs(thisZ) < 1E-10)
      thisZ = 0;
    TVector3 thisV3(thisX, thisY, thisZ);
    double thisR   = thisV3.Perp();
    double thisP   = thisV3.Phi();
    //ipmt is the position in the vector. Runs 0->nPMTs-1
    //tubeid is the actual ID of the PMT. Runs 1->nPMTs
    if((int)(ipmt + 1) != thisTubeID)
      G4cerr << "PMT ID is not the expected one!"
	     << " Vector position + 1 " << ipmt+1
	     << " PMT ID " << thisTubeID << G4endl;

#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 4
    G4cout << "PMT " << ipmt << " has ID " << thisTubeID
	   << " position "
	   << thisX << " "
	   << thisY << " "
	   << thisZ << "\t"
	   << thisR << " "
	   << thisP << " "
	   << thisZ << G4endl;
#endif

    if(thisR > RMax)
      RMax = thisR;
    vR.push_back(thisR);
    vP.push_back(thisP);
    vZ.push_back(thisZ);
    allPMTs.push_back(thisV3);
  }//ipmt
  
#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 3
  G4cout << "R positions" << G4endl;
  PrintVectorCount(vR);
  G4cout << "phi positions" << G4endl;
  PrintVectorCount(vP);
  G4cout << "Z positions" << G4endl;
  PrintVectorCount(vZ);
#endif

  //make a definition of top/bottom (use Z position)
  std::set<double> sZ(vZ.begin(), vZ.end());
  int nzmax = 0, nzmin = 0;
  double zmin, zmax;
  for(std::set<double>::iterator it = sZ.begin(); it != sZ.end(); ++it) {
    double thisZ = *it;
    double thiscount = std::count(vZ.begin(), vZ.end(), thisZ);
    if(thisZ < 0 && thiscount > nzmin) {
      zmin  = thisZ;
      nzmin = thiscount;
    }
    else if(thisZ > 0 && thiscount > nzmax) {
      zmax  = thisZ;
      nzmax = thiscount;
    }
  }//it
#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 2
  G4cout << "top    defined as z = " << zmax << " with " << nzmax << " PMTs" << G4endl
	 << "bottom defined as z = " << zmin << " with " << nzmin << " PMTs" << G4endl;
#endif

  //separate PMTs into top/bottom/side
  std::vector<std::pair<int, TVector3> > topPMTs, bottomPMTs, sidePMTs;
  for(unsigned int iv = 0; iv < allPMTs.size(); iv++) {
    int thispmtid = iv + 1;
    double thisZ = allPMTs[iv].Z();
    if(TMath::Abs(thisZ - zmax) < 1E-6)
      topPMTs.push_back(std::make_pair(thispmtid, allPMTs[iv]));
    else if(TMath::Abs(thisZ - zmin) < 1E-6)
      bottomPMTs.push_back(std::make_pair(thispmtid, allPMTs[iv]));
    else
      sidePMTs.push_back(std::make_pair(thispmtid, allPMTs[iv]));
  }//iv
#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 1
  G4cout << topPMTs.size() << " top PMTs + "
	 << bottomPMTs.size() << " bottom PMTs + "
	 << sidePMTs.size() << " side PMTs = "
	 << allPMTs.size() << " total PMTs " << G4endl;
#endif

#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 3
  G4cout << "SIDE" << G4endl;
  PrintVector3Count(sidePMTs,   false, true, true);
  G4cout << "TOP" << G4endl;
  PrintVector3Count(topPMTs,    true,  true, false);
  G4cout << "BOTTOM" << G4endl;
  PrintVector3Count(bottomPMTs, true,  true, false);
#endif

  //split the side PMTS into rectangles
  // note that the side PMTs are arranged in rows of constant z, and columns of constant phi
  TAxis aZ(regionsNBinsZ, zmin, zmax);
  TAxis aP(regionsNBinsP, -TMath::Pi(), +TMath::Pi());
  std::vector<int> thesePMTs;
  TVector3 avposition;
  int nsides = 0;
  for(int iz = 1; iz <= regionsNBinsZ; iz++) {
    double thiszlo = aZ.GetBinLowEdge(iz) - overlap;
    double thiszhi = aZ.GetBinUpEdge(iz)  - overlap;
    for(int ip = 1; ip <= regionsNBinsP; ip++) {
      thesePMTs.clear();
      double thisplo = aP.GetBinLowEdge(ip) - overlap;
      double thisphi = aP.GetBinUpEdge(ip)  - overlap;
      for(unsigned int ipmt = 0; ipmt < sidePMTs.size(); ipmt++) {
	double thisZ = sidePMTs[ipmt].second.Z();
	double thisP = sidePMTs[ipmt].second.Phi();
	if(thisZ >= thiszlo && thisZ < thiszhi &&
	   thisP >= thisplo && thisP < thisphi)
	  thesePMTs.push_back(sidePMTs[ipmt].first);
	avposition += sidePMTs[ipmt].second;
      }//ipmt
      int thispmts = thesePMTs.size();
      avposition.SetX(avposition.X() / thispmts);
      avposition.SetY(avposition.Y() / thispmts);
      avposition.SetZ(avposition.Z() / thispmts);
      pmtBlocks.push_back(std::make_pair(thesePMTs, avposition));
#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 2
      G4cout << "Side PMT block " << pmtBlocks.size() << " has " << thispmts << " entries. "
	     << thiszlo << " <= Z < " << thiszhi << " && "
	     << thisplo << " <= P < " << thisphi
	     << " Av position " << avposition.X() << "," << avposition.Y() << "," << avposition.Z()
	     << G4endl;
#endif
      nsides += thispmts;
    }//ip
  }//iz
#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 1
  G4cout << nsides << " total entries in side PMTs (including overlaps)" << G4endl;
#endif

  //split the top & bottom PMTs into groups
  // note all have identical Z
  // essentially no pattern in phi (majority are unique)
  // essentially no pattern in R (some rings exist, but very few)
  const int R = RMax + 1E-6;
  std::vector<double> radii; //0th entry is a dummy index. 1st entry (central circle) not counted as 'ring'
  radii.push_back(0);
  const double Z = 1 + regionsNRingSectors;
  
  //calculate the values of radii (i.e. the radius of each segment)
  if(regionsNRings == 1) {
    radii.push_back(TMath::Sqrt(R*R / Z));
  }
  else if(regionsNRings == 2) {
    radii.push_back(TMath::Sqrt(R*R / (Z*Z - regionsNRingSectors)));
    radii.push_back(TMath::Sqrt(radii[1]*radii[1]*Z));
  }
  else if(regionsNRings == 3) {
    radii.push_back(TMath::Sqrt(R*R / (Z*Z*Z - 2*Z*regionsNRingSectors)));
    radii.push_back(TMath::Sqrt(radii[1]*radii[1]*Z));
    radii.push_back(TMath::Sqrt(radii[2]*radii[2]*Z - radii[1]*radii[1]*regionsNRingSectors));
  }
  else if(regionsNRings == 4) {
    radii.push_back(TMath::Sqrt(R*R / (Z*Z*Z*Z - 3*Z*Z*regionsNRingSectors + regionsNRingSectors*regionsNRingSectors)));
    radii.push_back(TMath::Sqrt(radii[1]*radii[1]*Z));
    radii.push_back(TMath::Sqrt(radii[2]*radii[2]*Z - radii[1]*radii[1]*regionsNRingSectors));
    radii.push_back(TMath::Sqrt(radii[3]*radii[3]*Z - radii[2]*radii[2]*regionsNRingSectors));
  }
  else if(regionsNRings == 5) {
    radii.push_back(TMath::Sqrt(R*R / (Z*Z*Z*Z*Z - 4*Z*Z*Z*regionsNRingSectors + 3*Z*regionsNRingSectors*regionsNRingSectors)));
    radii.push_back(TMath::Sqrt(radii[1]*radii[1]*Z));
    radii.push_back(TMath::Sqrt(radii[2]*radii[2]*Z - radii[1]*radii[1]*regionsNRingSectors));
    radii.push_back(TMath::Sqrt(radii[3]*radii[3]*Z - radii[2]*radii[2]*regionsNRingSectors));
    radii.push_back(TMath::Sqrt(radii[4]*radii[4]*Z - radii[3]*radii[3]*regionsNRingSectors));
  }
  else {
    G4cerr << "regionsNRings value of " << regionsNRings << " not implemented" << G4endl;
    exit(-1);
  }
  radii.push_back(R);
  assert((int)radii.size() == regionsNRings + 2);

#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 2
  for(int i = 0; i <= regionsNRings; i++)
    G4cout << "Ring " << i << " radius: " << radii[i] << G4endl;
#endif

  //make the segments
  std::vector<int> theseTopPMTs, theseBottomPMTs;
  TVector3 avtopposition, avbottomposition;
  double total_area = 0;
  int ntops = 0;
  int nbottoms = 0;
  for(int ir = regionsNRings; ir >= 0; ir--) {
    const int nSectors = regionsNCentralSectors * TMath::Power(regionsNRingSectors, ir);
    double thisrlo = radii[ir];
    double thisrhi = radii[ir+1];
    for(int is = 0; is < nSectors; is++) {
      double thisplo = ((is  )*2*TMath::Pi() / (double)nSectors) - TMath::Pi();
      double thisphi = ((is+1)*2*TMath::Pi() / (double)nSectors) - TMath::Pi();
      //top
      theseTopPMTs.clear();
      for(unsigned int ipmt = 0; ipmt < topPMTs.size(); ipmt++) {
	double thisR = topPMTs[ipmt].second.Perp();
	double thisP = topPMTs[ipmt].second.Phi();
	if(thisR >= thisrlo && thisR < thisrhi &&
	   thisP >= thisplo && thisP < thisphi)
	  theseTopPMTs.push_back(topPMTs[ipmt].first);
        avtopposition += topPMTs[ipmt].second;
      }//ipmt
      int ntoppmts = theseTopPMTs.size();
      avtopposition.SetX(avtopposition.X() / ntoppmts);
      avtopposition.SetY(avtopposition.Y() / ntoppmts);
      avtopposition.SetZ(avtopposition.Z() / ntoppmts);
      pmtBlocks.push_back(std::make_pair(theseTopPMTs, avtopposition));
      //bottom
      theseBottomPMTs.clear();
      for(unsigned int ipmt = 0; ipmt < bottomPMTs.size(); ipmt++) {
	double thisR = bottomPMTs[ipmt].second.Perp();
	double thisP = bottomPMTs[ipmt].second.Phi();
	if(thisR >= thisrlo && thisR < thisrhi &&
	   thisP >= thisplo && thisP < thisphi)
	  theseBottomPMTs.push_back(bottomPMTs[ipmt].first);
        avbottomposition += sidePMTs[ipmt].second;
      }//ipmt
      int nbottompmts = theseBottomPMTs.size();
      avbottomposition.SetX(avbottomposition.X() / nbottompmts);
      avbottomposition.SetY(avbottomposition.Y() / nbottompmts);
      avbottomposition.SetZ(avbottomposition.Z() / nbottompmts);
      pmtBlocks.push_back(std::make_pair(theseBottomPMTs, avbottomposition));

      double area = TMath::Pi() * (TMath::Power(thisrhi, 2) - TMath::Power(thisrlo, 2)) / (double)nSectors;
#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 2
      G4cout << "Ring " << ir << " Sector " << is << ":"
	     << thisrlo << " <= R < " << thisrhi << " && "
	     << thisplo << " <= P < " << thisphi
	     << " (area: " << area
	     << ") TopPMTs: " << ntoppmts
	     << " BottomPMTs: " << nbottompmts
	     << G4endl;
#endif
      ntops    += ntoppmts;
      nbottoms += nbottompmts;
      total_area += area;
    }//is
  }//ir
#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 2
  G4cout << "Total area " << total_area << " = area from detector radius " 
	 << TMath::Pi() * TMath::Power(R, 2) << G4endl;
#endif
#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 1
  G4cout << ntops << " total top PMTs (including overlaps)" << G4endl;
  G4cout << nbottoms << " total bottom PMTs (including overlaps)" << G4endl;
#endif

  G4cout << pmtBlocks.size() << " total PMT blocks" << G4endl;
  std::vector<int> pmts_per_block;
  unsigned int blocklo = 9999, blockhi = 0;
  for(unsigned int ib = 0; ib < pmtBlocks.size(); ib++) {
    const unsigned int npmtsinblock = pmtBlocks[ib].first.size();
    if(npmtsinblock > blockhi)
      blockhi = npmtsinblock;
    if(npmtsinblock < blocklo)
      blocklo = npmtsinblock;
    pmts_per_block.push_back(npmtsinblock);
#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 2
    G4cout << "Block " << ib << "\t" << npmtsinblock << " PMTs:";
    for(unsigned int ip = 0; ip < npmtsinblock; ip++) {
      G4cout << " " << pmtBlocks[ib][ip];
    }//ip
    G4cout << G4endl;
#endif
  }//ib
  G4cout << "PMTs per block ranges from " << blocklo << " to " << blockhi << G4endl;
#if WCSIMWCTRIGGER_POPULATE_PMT_AREAS_VERBOSE >= 1
  PrintVectorCount(pmts_per_block);
#endif
}

void WCSimWCTriggerBase::PrintVector3Count(std::vector<std::pair<int, TVector3> > & v, bool r, bool p, bool z)
{
  std::vector<double> vz, vp, vr;
  for(unsigned int iv = 0; iv < v.size(); iv++) {
    if(r)
      vr.push_back(v[iv].second.Perp());
    if(z)
      vz.push_back(v[iv].second.Z());
    if(p)
      vp.push_back(v[iv].second.Phi());
  }//iv
  if(z) {
    G4cout << "Z positions" << G4endl;
    PrintVectorCount(vz);
  }
  if(r) {
    G4cout << "R positions" << G4endl;
    PrintVectorCount(vr);
  }
  if(p) {
    G4cout << "phi positions" << G4endl;
    PrintVectorCount(vp);
  }
}

template<typename T> void WCSimWCTriggerBase::PrintVectorCount(std::vector<T> & v)
{
  std::set<T> s(v.begin(), v.end());
  for(typename std::set<T>::iterator it = s.begin(); it != s.end(); ++it) {
    G4cout << *it << "\t" << std::count(v.begin(), v.end(), *it) << "\t" << s.count(*it) << G4endl;
  }
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

std::vector<int> WCSimWCTriggerBase::FindRegion(const int tubeid)
{
  std::vector<int> regions;
  std::vector<int>::iterator it;
  for(unsigned int ir = 0; ir < pmtBlocks.size(); ir++) {
    it = std::find(pmtBlocks[ir].first.begin(), pmtBlocks[ir].first.end(), tubeid);
    if(it != pmtBlocks[ir].first.end())
      regions.push_back(ir);
  }//ir
  if(!regions.size()) {
    G4cerr << "WCSimWCTriggerBase::FindRegion() >> no region found for PMT with ID " << tubeid << G4endl;
    exit(-1);
  }
  return regions;
}

void WCSimWCTriggerBase::AlgNHitsThenRegions(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits)
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

  G4cout << "WCSimWCTriggerBase::AlgNHitsThenRegions. Number of entries in input digit collection: " << WCDCPMT->entries() << G4endl;
#ifdef WCSIMWCTRIGGER_VERBOSE
  int temp_total_pe = 0;
  for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
    temp_total_pe += (*WCDCPMT)[i]->GetTotalPe();
  }
  G4cout << "WCSimWCTriggerBase::AlgNHitsThenRegions. " << temp_total_pe << " total p.e. input" << G4endl;
#endif

  // the upper time limit is set to the final possible full trigger window
  while(window_start_time <= window_end_time) {
    int n_digits = 0, local_n_digits = 0;
    float triggertime; //save each digit time, because the trigger time is the time of the first hit above threshold
    bool triggerfound = false;
    digit_times.clear();
    memset(localNHitsHits, 0, nPMTs            * sizeof(int));
    memset(regionsHits,    0, pmtBlocks.size() * sizeof(int));
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
      }//ip //loop over Digits
    }//i //loop over PMTs

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

    //NHits failed. Try region-based local NHits
    else {
      std::vector<Float_t> local_info;

      //get a unique list of regions hit
      std::vector<int> regions_hit;
      for(size_t ip = 0; ip < pmts_hit.size(); ip++) {
	int tubeid = pmts_hit[ip];
	std::vector<int> theseregions = FindRegion(tubeid);
	//#ifdef WCSIMWCTRIGGER_VERBOSE
	G4cout << "Tube ID " << tubeid << " has " << localNHitsHits[tubeid-1] << " in region(s)";
	//#endif
	for(unsigned int ir = 0; ir < theseregions.size(); ir++) {
	  int thisregion = theseregions[ir];
	  regions_hit.push_back(thisregion);
	  regionsHits[thisregion] += localNHitsHits[tubeid-1]; //tubeID goes from 1->nPMTs
	  //#ifdef WCSIMWCTRIGGER_VERBOSE
	  G4cout << " " << thisregion;
	  //#endif
	}//ir
	//#ifdef WCSIMWCTRIGGER_VERBOSE
	G4cout << G4endl;
	//#endif
      }//ip //loop over PMTs with hits
      std::set<int> unique_regions_hit(regions_hit.begin(), regions_hit.end());
      
      //loop over all regions with hits
      for(std::set<int>::iterator it = unique_regions_hit.begin(); it != unique_regions_hit.end(); ++it) {
	int this_region = *it;
	int nlocal = regionsHits[this_region];
	//#ifdef WCSIMWCTRIGGER_VERBOSE
	G4cout << "Region ID " << this_region << " has " << nlocal << " hits" << G4endl;
	//#endif
      }//it //end loop over regions with hits
    }//region trigger

    /*
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
	TriggerTypes.push_back(kTriggerRegions);
	TriggerInfos.push_back(local_info);
	triggerfound = true;
      }//trigger(s) found on any PMT
    }//local NHits trigger
    */



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
  
  G4cout << "Found " << ntrig << " NHit or NHitRegions triggers" << G4endl;

  //call FillDigitsCollection() whether any triggers are found or not
  // (what's saved depends on saveFailuresMode)
  FillDigitsCollection(WCDCPMT, remove_hits, kTriggerUndefined);
}

void WCSimWCTriggerBase::AlgNHitsThenITC(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits)
{
  //Now we will try to find triggers
  //loop over PMTs, and Digits in each PMT.
  // If nhits > Threshhold in a time window, then we have a trigger
  // If this cut fails, attempt an ITC ratio trigger

  int ntrig = 0;
  int window_start_time = 0;
  int window_end_time   = WCSimWCTriggerBase::LongTime - ndigitsWindow;
  int window_step_size  = 5; //step the search window along this amount if no trigger is found
  float lasthit;
  std::vector<int> digit_times;
  std::vector<int> digit_times_itc_small, digit_times_itc_large;
  bool first_loop = true;

  G4cout << "WCSimWCTriggerBase::AlgNDigitsThenITC. Number of entries in input digit collection: " << WCDCPMT->entries() << G4endl;
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
  int temp_total_pe = 0;
  for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
    temp_total_pe += (*WCDCPMT)[i]->GetTotalPe();
  }
  G4cout << "WCSimWCTriggerBase::AlgNDigitsThenITC. " << temp_total_pe << " total p.e. input" << G4endl;
#endif

  // the upper time limit is set to the final possible full trigger window
  while(window_start_time <= window_end_time) {
    int n_digits = 0;
    int n_digits_itc_small = 0, n_digits_itc_large = 0;
    float triggertime; //save each digit time, because the trigger time is the time of the first hit above threshold
    bool triggerfound = false;
    digit_times.clear();

    //Loop over each PMT & count NDigits in window [window_start_time, window_start_time + ndigitsWindow]
    //Also count in two extra windows for the ITC cut
    // [window_start_time, window_start_time + itcSmallWindow]
    // [window_start_time - itcLargeWindowLow, window_start_time + itcLargeWindowHigh]
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
    if(n_digits > ndigitsThreshold) {
      ntrig++;
      //The trigger time is the time of the first digit above threshold
      std::sort(digit_times.begin(), digit_times.end());
      triggertime = digit_times[ndigitsThreshold];
      triggertime -= (int)triggertime % 5;
      //save the trigger information
      TriggerTimes.push_back(triggertime);
      TriggerTypes.push_back(kTriggerNDigits);
      TriggerInfos.push_back(std::vector<Float_t>(1, n_digits));
      triggerfound = true;
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
      G4cout << EnumAsString(kTriggerNDigits) << " trigger passed with time " << triggertime << G4endl;
#endif
    }//NDigits trigger passed

    //The simple NHits trigger hasn't been passed. See if the ITC ratio trigger can be passed
    double itc_ratio = (double)n_digits_itc_small / (double)n_digits_itc_large;
    if(!triggerfound && itc_ratio > itcRatioThreshold) {
      ntrig++;
      //The trigger time is the time of the last hit, in order to be close to the NDigits trigger time (i.e. the first hit above the threshold)
      std::sort(digit_times.begin(), digit_times.end());
      triggertime = digit_times.back();
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
      G4cout << n_digits << " digits found in " << ndigitsWindow << "nsec trigger window ["
             << window_start_time << ", " << window_start_time + ndigitsWindow
             << "]. Threshold is: " << nhitsThreshold
             << G4endl
             << n_digits_itc_small << "(" << n_digits_itc_large << ") digits found in "
             << itcSmallWindow << "(" << itcLargeWindowHigh - itcLargeWindowLow << ") window giving ITC ratio"
             << itc_ratio << ". Threshold is: " << itcRatioThreshold
             << G4endl;
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
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
      G4cout << "Last hit found to be at " << lasthit
             << ". Changing window_end_time from " << window_end_time
             << " to " << lasthit - (ndigitsWindow - 10)
             << G4endl;
#endif
      window_end_time = lasthit - (ndigitsWindow - 10);
      first_loop = false;
    }
  }//while <= window_end_time

  G4cout << "Found " << ntrig << " NHitThenITC triggers" << G4endl;
  //call FillDigitsCollection() whether any triggers are found or not
  // (what's saved depends on saveFailuresMode)
  FillDigitsCollection(WCDCPMT, remove_hits, kTriggerUndefined);
}

void WCSimWCTriggerBase::AlgNoTrigger(WCSimWCDigitsCollection* WCDCPMT)
{
  //Does not doanything, just writes out all hits
  std::vector<Float_t> triggerinfo;
  Int_t Ndigits=0;
  for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
    for ( G4int ip = 0 ; ip < (*WCDCPMT)[i]->GetTotalPe() ; ip++) {
      Ndigits++;
    }
  }
  triggerinfo.push_back(Ndigits);
  TriggerTypes.push_back(kTriggerNoTrig);
  TriggerInfos.push_back(triggerinfo);
  TriggerTimes.push_back(0.);

  FillDigitsCollection(WCDCPMT, false, kTriggerNoTrig);
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

	  //first apply time offsets
	  float peSmeared = (*WCDCPMT)[i]->GetPe(ip);
	  G4double digihittime = -triggertime
	    + WCSimWCTriggerBase::offset
	    + digit_time;

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

void WCSimWCTriggerNHitsThenLocalNHits::WriteGeomInfo()
{
  localNHitsNeighbours = nPMTs - 1; //don't include self!
  //fill the arrays with neighbours
  FindAllPMTNearestNeighbours();

#if WCSIMWCTRIGGER_DEBUG_GEOM_INFO >= 1
  TFile f("det.root","RECREATE");
  TTree t("det","");
  std::vector<int> neighbours;
  int pmt;
  t.Branch("pmt", &pmt);
  t.Branch("neighbours",&neighbours); 
#endif

  std::ofstream ofs("det.bin", std::ofstream::out | std::ofstream::binary | std::ofstream::trunc);
  for(int ipmt = 0, npmt = pmtNeighbours.size(); ipmt < npmt; ipmt++) {
    if(ipmt % (nPMTs/20 + 1) == 0) {
      G4cout << "WCSimWCTriggerNHitsThenLocalNHits::WriteGeomInfo at "
             << ipmt / (float)nPMTs * 100 << "%"
             << " (" << ipmt << " out of " << nPMTs << ")"
             << G4endl;
    }
    int pmtid = ipmt + 1;
    ofs.write(reinterpret_cast<const char *>(&pmtid), sizeof(pmtid));
    ofs.write(reinterpret_cast<const char *>(pmtNeighbours[ipmt].data()), sizeof(int) * pmtNeighbours[ipmt].size());
#if WCSIMWCTRIGGER_DEBUG_GEOM_INFO >= 2
    neighbours = pmtNeighbours[ipmt];
    pmt = ipmt;
    t.Fill();
    G4cout << pmtid << G4endl;
    for(int in = 0, nn = pmtNeighbours[ipmt].size(); in < nn; in++) {
      G4cout << " " << pmtNeighbours[ipmt][in];
    }//in
    G4cout << G4endl;
#endif
  }//ipmt

  ofs.close();
#if WCSIMWCTRIGGER_DEBUG_GEOM_INFO >= 1
  t.Write();
#endif


}

void WCSimWCTriggerNHitsThenLocalNHits::ReadGeomInfo()
{
  std::ifstream ifs("det.bin", std::ifstream::in | std::ifstream::binary);
  int pmtid;
  for(unsigned int ipmt = 0; ipmt < nPMTs; ipmt++) {
    if(ipmt % (nPMTs/20 + 1) == 0) {
      G4cout << "WCSimWCTriggerNHitsThenLocalNHits::ReadGeomInfo at "
             << ipmt / (float)nPMTs * 100 << "%"
             << " (" << ipmt << " out of " << nPMTs << ")"
             << G4endl;
    }
    std::vector<int> neighbours;
    ifs.read((char*)&pmtid, sizeof(int));
#if WCSIMWCTRIGGER_DEBUG_GEOM_INFO >= 2
    G4cout << pmtid << G4endl;
#endif
    if(pmtid != (int)(ipmt + 1)) {
      G4cerr << "PMT ID mistmatch. Expected " << ipmt + 1 << " got " << pmtid << G4endl;
      exit(-1);
    }
    for(unsigned int in = 0; in < nPMTs - 1; in++) {
      ifs.read((char*)&pmtid, sizeof(int));
#if WCSIMWCTRIGGER_DEBUG_GEOM_INFO >= 2
      G4cout << " " << pmtid;
#endif
      if((int)in < localNHitsNeighbours)
	neighbours.push_back(pmtid);
    }//in
#if WCSIMWCTRIGGER_DEBUG_GEOM_INFO >= 2
    G4cout << G4endl;
#endif
    pmtNeighbours.push_back(neighbours);
  }//ipmt
  ifs.close();

#if WCSIMWCTRIGGER_DEBUG_GEOM_INFO >= 1
  TFile f("detread.root","RECREATE");
  TTree t("det","");
  std::vector<int> neighbours;
  int pmt;
  t.Branch("pmt", &pmt);
  t.Branch("neighbours",&neighbours);
  for(int ipmt = 0, npmt = pmtNeighbours.size(); ipmt < npmt; ipmt++) {
    neighbours = pmtNeighbours[ipmt];
    pmt = ipmt;
    t.Fill();
  }//ipmt
  t.Write();
#endif
}


// *******************************************
// DERIVED CLASS
// *******************************************

WCSimWCTriggerNHitsThenAnisotropy::WCSimWCTriggerNHitsThenAnisotropy(G4String name,
								     WCSimDetectorConstruction* myDetector,
								     WCSimWCDAQMessenger* myMessenger)
  :WCSimWCTriggerNHitsThenLocalNHits(name, myDetector, myMessenger)
{
  triggerClassName = "NHitsThenAnisotropy";
}

WCSimWCTriggerNHitsThenAnisotropy::~WCSimWCTriggerNHitsThenAnisotropy()
{
}

void WCSimWCTriggerNHitsThenAnisotropy::DoTheWork(WCSimWCDigitsCollection* WCDCPMT)
{
  //Apply an NHitsThenAnisotropy trigger
  bool remove_hits = false;
  //AlgNHitsThenAnisotropy(WCDCPMT, remove_hits);
}

void WCSimWCTriggerNHitsThenAnisotropy::ReadGeomInfo()
{
  localNHitsNeighbours = nPMTs / 2;
  WCSimWCTriggerNHitsThenLocalNHits::ReadGeomInfo();
}


// *******************************************
// DERIVED CLASS
// *******************************************

WCSimWCTriggerNHitsThenRegions::WCSimWCTriggerNHitsThenRegions(G4String name,
							       WCSimDetectorConstruction* myDetector,
							       WCSimWCDAQMessenger* myMessenger)
  :WCSimWCTriggerBase(name, myDetector, myMessenger)
{
  triggerClassName = "NHitsThenRegions";
  GetVariables();
  PopulatePMTAreas();

  //reserve an array to store the number of digits on each PMT
  localNHitsHits = new int[nPMTs];
  //and fill it with 0s
  memset(localNHitsHits, 0, nPMTs * sizeof(int));

  //reserve an array to store the number of digits in each PMT region
  regionsHits = new int[pmtBlocks.size()];
  //and fill it with 0s
  memset(regionsHits, 0, pmtBlocks.size() * sizeof(int));
}

WCSimWCTriggerNHitsThenRegions::~WCSimWCTriggerNHitsThenRegions()
{
}

void WCSimWCTriggerNHitsThenRegions::DoTheWork(WCSimWCDigitsCollection* WCDCPMT)
{
  //Apply an NHitsThenRegions trigger
  bool remove_hits = false;
  AlgNHitsThenRegions(WCDCPMT, remove_hits);
}

// *******************************************
// DERIVED CLASS
// *******************************************
WCSimWCTriggerNHitsThenITC::WCSimWCTriggerNHitsThenITC(G4String name,
					 WCSimDetectorConstruction* myDetector,
					 WCSimWCDAQMessenger* myMessenger)
  :WCSimWCTriggerBase(name, myDetector, myMessenger)
{
  triggerClassName = "NHitsThenITC";
  GetVariables();
}

WCSimWCTriggerNHitsThenITC::~WCSimWCTriggerNHitsThenITC()
{
}

void WCSimWCTriggerNHitsThenITC::DoTheWork(WCSimWCDigitsCollection* WCDCPMT) {
  //Apply an NHitsThenITC trigger
  bool remove_hits = false;
  AlgNHitsThenITC(WCDCPMT, remove_hits);
}

// *******************************************
// DERIVED CLASS
// *******************************************
WCSimWCTriggerNoTrigger::WCSimWCTriggerNoTrigger(G4String name,
						 WCSimDetectorConstruction* myDetector,
						 WCSimWCDAQMessenger* myMessenger)
  :WCSimWCTriggerBase(name, myDetector, myMessenger)
{
  triggerClassName = "NoTrigger";
  GetVariables();
}

WCSimWCTriggerNoTrigger::~WCSimWCTriggerNoTrigger()
{
}

void WCSimWCTriggerNoTrigger::DoTheWork(WCSimWCDigitsCollection* WCDCPMT) {
  //Apply a dummy 'pass everything' trigger
  bool remove_hits = false;
  AlgNoTrigger(WCDCPMT);
}
