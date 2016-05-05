#include trigger_tools.hxx

#ifndef __DIGIT_TIME_VERBOSE__
//#define __DIGIT_TIME_VERBOSE__
#endif

trigger_tools::trigger_tools()
{
  digit_times         = null;
  digit_times_physics = null;
  digit_times_noise   = null;
  digit_times_mix     = null;
}

trigger_tools::~trigger_tools()
{
  CleanupDigitTimes();
}

void trigger_tools::CleanupDigitTimes()
{
  if(digit_times != null) {
    delete digit_times;
    digit_times = null;
  }
  if(digit_times_physics != null) {
    delete digit_times_physics;
    digit_times_physics = null;
  }
  if(digit_times_noise != null) {
    delete digit_times_noise;
    digit_times_noise = null;
  }
  if(digit_times_mix != null) {
    delete digit_times_mix;
    digit_times_mix = null;
  }
}

void trigger_tools::CreateDigitTimes()
{
  if(digit_times == null) {
    digit_times = new vector<double>;
  }
  if(digit_times_physics == null) {
    digit_times_physics = new vector<double>;
  }
  if(digit_times_noise == null) {
    digit_times_noise = new vector<double>;
  }
  if(digit_times_mix == null) {
    digit_times_mix = new vector<double>;
  }
}

void trigger_tools::FillDigitTimes(double digitime, DigiType_t digitype)
{
  digit_times.push_back(digitime);
  switch(filltype) {
  case kDigiTypePhysics:
    digi_times_physics.push_back(digitime);
    break;
  case kDigiTypeNoise:
    digi_times_noise.push_back(digitime);
    break;
  case kDigiTypeMix:
    digi_times_mix.push_back(digitime);
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

    DigiType_t filltype = GetDigitType(wcsimrootcherenkovdigihit, event);
    FillDigitTimes(digitime, filltype);
  }//idigipmt
}

DigiType_t trigger_tools::GetDigitType(WCSimRootCherenkovDigiHit * wcsimrootcherenkovdigihit, WCSimRootEvent * event)
{
  if(event == null)
    return kDigiTypeUndefined;

  const int digitubeid = wcsimrootcherenkovdigihit->GetTubeId();

  TClonesArray * hits     = event->GetTrigger(0)->GetCherenkovHits();
  TClonesArray * hittimes = event->GetTrigger(0)->GetCherenkovHitTimes();
  
  //first, find the correct PMT in the WCSimRootCherenkovHit array
  int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0); //position in the WCSimRootCherenkovHitTime for this PMT
  int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1); //number of pe (i.e. entries) for this tube in the WCSimRootCherenkovHitTime
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

  if(verbose)
    cout << peForTube << " PMT hits found for digit " << idigipmt << " with tube ID " << tube_id << endl;
  
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
    if(verbose)
      cout << " hit time " << hittime << " " << parentid << endl;
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
