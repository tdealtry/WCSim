#include "trigger_tools.hxx"

#include <algorithm>
#include <iostream>

#ifndef __DIGIT_TIME_VERBOSE__
//#define __DIGIT_TIME_VERBOSE__
#endif

trigger_tools::trigger_tools()
{
  digit_times         = NULL;
  digit_times_physics = NULL;
  digit_times_noise   = NULL;
  digit_times_mix     = NULL;
}

trigger_tools::~trigger_tools()
{
  CleanupDigitTimes();
}

void trigger_tools::CleanupDigitTimes()
{
  if(digit_times != NULL) {
    delete digit_times;
    digit_times = NULL;
  }
  if(digit_times_physics != NULL) {
    delete digit_times_physics;
    digit_times_physics = NULL;
  }
  if(digit_times_noise != NULL) {
    delete digit_times_noise;
    digit_times_noise = NULL;
  }
  if(digit_times_mix != NULL) {
    delete digit_times_mix;
    digit_times_mix = NULL;
  }
}

void trigger_tools::CreateDigitTimes()
{
  if(digit_times == NULL) {
    digit_times = new vector<double>;
  }
  if(digit_times_physics == NULL) {
    digit_times_physics = new vector<double>;
  }
  if(digit_times_noise == NULL) {
    digit_times_noise = new vector<double>;
  }
  if(digit_times_mix == NULL) {
    digit_times_mix = new vector<double>;
  }
}

void trigger_tools::SortDigitTimes()
{
  if(digit_times != NULL) {
    std::sort(digit_times->begin(), digit_times->end());
  }
  if(digit_times_physics != NULL) {
    std::sort(digit_times_physics->begin(), digit_times_physics->end());
  }
  if(digit_times_noise != NULL) {
    std::sort(digit_times_noise->begin(), digit_times_noise->end());
  }
  if(digit_times_mix != NULL) {
    std::sort(digit_times_mix->begin(), digit_times_mix->end());
  }
}

void trigger_tools::PrintDigitTimes(DigiType_t digitype)
{
  vector<double> * v;
  switch(digitype) {
  case kDigiTypePhysics:
    v = digit_times_physics;
    break;
  case kDigiTypeNoise:
    v = digit_times_noise;
    break;
  case kDigiTypeMix:
    v = digit_times_mix;
    break;
  case kDigiTypeUndefined:
    v = digit_times;
    break;
  case kDigiTypeError:
    exit(-1);
    break;
  }
  cout << "trigger_tools::PrintDigitTimes() Printing digit type " << digitype << endl
       << "vector has " << v->size() << " entries" << endl;
  for(size_t i = 0; i < v->size(); i++) {
    cout << v->at(i);
    if((i % 10) == 9)
      cout << endl;
    else
      cout << "\t";
  }//i
}

void trigger_tools::FillDigitTimes(double digitime, DigiType_t digitype)
{
  digit_times->push_back(digitime);
  switch(digitype) {
  case kDigiTypePhysics:
    digit_times_physics->push_back(digitime);
    break;
  case kDigiTypeNoise:
    digit_times_noise->push_back(digitime);
    break;
  case kDigiTypeMix:
    digit_times_mix->push_back(digitime);
    break;
  case kDigiTypeUndefined:
    break;
  case kDigiTypeError:
    exit(-1);
    break;
  }
}

void trigger_tools::PopulateDigitTimes(WCSimRootTrigger * trigger, bool append, WCSimRootEvent * event)
{
  if(!append)
    CleanupDigitTimes();
  CreateDigitTimes();
  
  // Loop through elements in the TClonesArray of WCSimRootCherenkovDigHits
  const long ncherenkovdigihits = trigger->GetNcherenkovdigihits();
  for(long idigipmt = 0; idigipmt < ncherenkovdigihits; idigipmt++) {
    //get the digit
#ifdef __DIGIT_TIME_VERBOSE__
      cout << "Getting digit " << idigipmt << endl;
#endif
    TObject * Digit = (trigger->GetCherenkovDigiHits())->At(idigipmt);
    WCSimRootCherenkovDigiHit * wcsimrootcherenkovdigihit = 
      dynamic_cast<WCSimRootCherenkovDigiHit *>(Digit);
    //get the charge, time, PMT
    const double digitime   = wcsimrootcherenkovdigihit->GetT();
    //const double digipe     = wcsimrootcherenkovdigihit->GetQ();
    //const int    digitubeid = wcsimrootcherenkovdigihit->GetTubeId();

    DigiType_t digitype = GetDigitType(wcsimrootcherenkovdigihit, event);
    FillDigitTimes(digitime, digitype);
  }//idigipmt
  SortDigitTimes();
}

trigger_tools::DigiType_t trigger_tools::GetDigitType(WCSimRootCherenkovDigiHit * wcsimrootcherenkovdigihit, WCSimRootEvent * event)
{
  if(event == NULL)
    return kDigiTypeUndefined;

  const int digitubeid = wcsimrootcherenkovdigihit->GetTubeId();

  TClonesArray * hits     = event->GetTrigger(0)->GetCherenkovHits();
  TClonesArray * hittimes = event->GetTrigger(0)->GetCherenkovHitTimes();
  
  //first, find the correct PMT in the WCSimRootCherenkovHit array
  int timeArrayIndex = -1; //position in the WCSimRootCherenkovHitTime for this PMT
  int peForTube      = -1; //number of pe (i.e. entries) for this tube in the WCSimRootCherenkovHitTime
  const long ncherenkovhits = event->GetTrigger(0)->GetNcherenkovhits();
  for(int ipmt = 0; ipmt < ncherenkovhits; ipmt++) {
#ifdef __DIGIT_TIME_VERBOSE__
    cout << "Getting hit " << ipmt << " of " << ncherenkovhits << endl;
#endif
    TObject * Hit = hits->At(ipmt);
    WCSimRootCherenkovHit *wcsimrootcherenkovhit =
      dynamic_cast<WCSimRootCherenkovHit *>(Hit);
    int hittubeid = wcsimrootcherenkovhit->GetTubeID();
    if(digitubeid == hittubeid) {
      timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
      peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
      break;
    }
  }//ipmt

  //check that a true hit has been found
  if(timeArrayIndex == -1) {
    cerr << "No PMT hits found for digit with tube ID " << digitubeid << endl;
    return kDigiTypeError;
  }

#ifdef __DIGIT_TIME_VERBOSE__
    cout << peForTube << " PMT hits found for digit " << idigipmt << " with tube ID " << tube_id << endl;
#endif

  //loop over the rawhits ids of hits that made up the digit,
  // and use them to find the hits in the WCSimRootCherenkovHitTime
  int n_noise_hits = 0, n_photon_hits = 0;
  vector<int> rawhit_ids = wcsimrootcherenkovdigihit->GetPhotonIds();
#ifdef __DIGIT_TIME_VERBOSE__
    cout << rawhit_ids.size() << " rawhits made up this digit" << endl;
#endif
    for(unsigned int irawhit = 0; irawhit < rawhit_ids.size(); irawhit++) {
    int this_rawhit = rawhit_ids[irawhit];
#ifdef __DIGIT_TIME_VERBOSE__
      cout << "Attempting to look for rawhit " << this_rawhit+1 << " in WCSimRootCherenkovHitTime array...";
#endif
    if(this_rawhit >= peForTube) {
      cerr << " There are only " << peForTube << " rawhits in this PMT" << endl;
      return kDigiTypeError;
    }
    //now look in the WCSimRootCherenkovHitTime array to count the number of photon / dark noise hits
    TObject * Hit = (hittimes)->At(timeArrayIndex + this_rawhit);
    WCSimRootCherenkovHitTime * wcsimrootcherenkovhittime =
      dynamic_cast<WCSimRootCherenkovHitTime *>(Hit);
    //const double hittime  = wcsimrootcherenkovhittime->GetTruetime();
    const int    parentid = wcsimrootcherenkovhittime->GetParentID();
#ifdef __DIGIT_TIME_VERBOSE__
      cout << " hit time " << hittime << " " << parentid << endl;
#endif
    if(parentid == -1) {
      n_noise_hits++;
    }
    else if(parentid < -1) {
      cerr << "Unknown hit parent id " << parentid << endl;
      return kDigiTypeError;
    }
    else {
      n_photon_hits++;
    }
  }//irawhit

  //work out what to return
  if(n_noise_hits && !n_photon_hits) {
    return kDigiTypeNoise;
  }
  else if(!n_noise_hits && n_photon_hits) {
    return kDigiTypePhysics;
  }
  else {
    return kDigiTypeMix;
  }
}

/*
vector<double> * trigger_tools::GetDigitTimes(DigiType_t digitype)
{
  switch(digitype) {
  case kDigiTypePhysics:
    return digit_times_physics;
    break;
  case kDigiTypeNoise:
    return digit_times_noise;
    break;
  case kDigiTypeMix:
    return digit_times_mix;
    break;
  case kDigiTypeUndefined:
    return digit_times;
    break;
  case kDigiTypeError:
    exit(-1);
    return 0;
    break;
  }
}
*/

void trigger_tools::CalcMaxITC(itc_tools * itc)
{
  itc->CalcMaxITC(digit_times, digit_times_physics, digit_times_mix);
}
