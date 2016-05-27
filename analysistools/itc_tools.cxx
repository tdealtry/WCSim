#include "itc_tools.hxx"

#include <algorithm>
#include <iostream>

using std::sort;
using std::cout;
using std::endl;
using std::cerr;

itc_tools::itc_tools(TFile * f, int smallwindow, int largewindow, int offset, bool calculatendigits, bool onetimeslice, int verbosity) :
  smallwindow(smallwindow), largewindow(largewindow), offset(offset),
  calcndigits(calculatendigits), trigger_tools(f, onetimeslice, verbosity)
{
  f->cd();

  if(calcndigits) {
    h1_max_ndigits        = new TH1D(TString::Format("h1_max_ndigits%d", largewindow), ";Max NDigits;Number of entries", 1000, 0, 1000);
    h1_max_ndigits_1slice = new TH1D(TString::Format("h1_max_ndigits%d_1slice", largewindow), ";NDigits (one slice);Number of entries", 1000, 0, 1000);
    h1_max_ndigits_time   = new TH1D(TString::Format("h1_max_ndigits%d_time", largewindow), ";Time of max NDigits;Number of entries", 2000, -10000, +10000);
  }
  h1_max_itc        = new TH1D(TString::Format("h1_max_itc%d_%d_%d", smallwindow, largewindow, offset), ";Max ITC ratio;Number of entries", 100, 0, 1);
  h1_max_itc_1slice = new TH1D(TString::Format("h1_max_itc%d_%d_%d_1slice", smallwindow, largewindow, offset), ";ITC ratio (one slice);Number of entries", 100, 0, 1);
  h1_max_itc_time   = new TH1D(TString::Format("h1_max_itc%d_%d_%d_time", smallwindow, largewindow, offset), ";Time of max ITC ratio;Number of entries", 12000, -10000, +110000);

  tree = new TTree(TString::Format("t_itc%d_%d_%d", smallwindow, largewindow, offset).Data(), TString::Format("ITC tree small %d large %d offset %d", smallwindow, largewindow, offset).Data());
  tree->Branch("var_max_ndigits", &var_max_ndigits);
  tree->Branch("var_max_ndigits_time", &var_max_ndigits_time);
  tree->Branch("var_max_ndigits_1slice", &var_max_ndigits_1slice);
  tree->Branch("var_max_itc", &var_max_itc);
  tree->Branch("var_max_itc_time", &var_max_itc_time);
  tree->Branch("var_max_itc_1slice", &var_max_itc_1slice);
}

itc_tools::~itc_tools()
{
  if(calcndigits) {
    delete h1_max_ndigits;
    delete h1_max_ndigits_1slice;
    delete h1_max_ndigits_time;
  }
  delete h1_max_itc;
  delete h1_max_itc_1slice;
  delete h1_max_itc_time;
  delete tree;
}

void itc_tools::CalcMaxITC()
{
  if(verbosity > 2)
    cout << "Calling itc_tools::CalcMaxITC() with " << smallwindow << "\t" << largewindow << endl;

  //add a little bit - the start is digits smeared early - there no earlier digits smeared into this gap
  vector<double>::iterator mintime_it = std::lower_bound(digit_times->begin(), digit_times->end(), digit_times->at(0) + 10);
  int mintime_pos = mintime_it - digit_times->begin();
  double mintime = *mintime_it;
  vector<double>::iterator maxtime_it = std::lower_bound(digit_times->begin(), digit_times->end(), digit_times->back() - largewindow - 10);
  maxtime_it--;
  int maxtime_pos = maxtime_it - digit_times->begin();
  double maxtime = *maxtime_it;

  /*
  //check
  cout << "mintime is " << mintime << " at position " << mintime_pos << endl;
  for(int i = 0; i < mintime_pos + 2; i++) {
    cout << i << "\t" << digit_times->at(i);
    if(i == mintime_pos)
      cout << "\t*" << endl;
    cout << endl;
  }
  */

  double first_physics = mintime;
  if(digit_times_physics->size() && digit_times_mix->size())
    first_physics = TMath::Min(digit_times_physics->at(0), digit_times_mix->at(0));
  else if( digit_times_physics->size() && !digit_times_mix->size())
    first_physics = digit_times_physics->at(0);
  else if(!digit_times_physics->size() &&  digit_times_mix->size())
    first_physics = digit_times_mix->at(0);
  if(mintime < first_physics) {
    mintime = first_physics;
    mintime_it = std::lower_bound(digit_times->begin(), digit_times->end(), mintime);
    mintime_pos = mintime_it - digit_times->begin();
  }
  vector<double>::iterator first_physics_it = mintime_it;
  if(verbosity > 3)
    cout << "itc_tools::CalcMaxITC() Physics starts at " << first_physics << endl;

  double maxndigitstime = -9999, maxitctime = -9999, maxitc = -1;
  int maxndigits = -1;
  for(vector<double>::iterator it = mintime_it; it != maxtime_it && it != digit_times->end(); it++) {
    double time = *it;
    if(verbosity > 4)
      cout << "itc_tools::CalcMaxITC() vector position " << it - digit_times->begin() << " of " << digit_times->size() << " with time " << time << endl;
    int ndigits_small = std::lower_bound(it, digit_times->end(), time + smallwindow) - it;
    int ndigits_large = std::lower_bound(it, digit_times->end(), time + largewindow) - it;
    double thisitc = (double)ndigits_small / (double)ndigits_large;
    if(thisitc > maxitc) {
      maxitc = thisitc;
      maxitctime = time;
    }
    if(verbosity > 4)
      cout << "itc_tools::CalcMaxITC() " << ndigits_small << "\t" << ndigits_large << "\t" << thisitc << "\t" << maxitc << endl;
    if(ndigits_large > maxndigits) {
      maxndigits = ndigits_large;
      maxndigitstime = time;
    }
    if(one_time_slice || (it == first_physics_it)) {
      if(calcndigits)
	h1_max_ndigits_1slice->Fill(ndigits_large);
      h1_max_itc_1slice->Fill(thisitc);
      var_max_ndigits_1slice = ndigits_large;
      var_max_itc_1slice = thisitc;
    }
    if(one_time_slice)
      break;
  }//it
  if(calcndigits) {
    h1_max_ndigits->Fill(maxndigits);
    h1_max_ndigits_time->Fill(maxndigitstime);
  }
  if(verbosity > 3)
    cout << "itc_tools::CalcMaxITC() filling itc with " << maxitc << "\t" << maxitctime << endl;
  h1_max_itc->Fill(maxitc);
  h1_max_itc_time->Fill(maxitctime);
  var_max_ndigits = maxndigits;
  var_max_ndigits_time = maxndigitstime;
  var_max_itc = maxitc;
  var_max_itc_time = maxitctime;
  tree->Fill();
}

void itc_tools::Write()
{
  f->cd();
  if(calcndigits) {
    h1_max_ndigits->Write();
    h1_max_ndigits_1slice->Write();
    h1_max_ndigits_time->Write();
  }
  h1_max_itc->Write();
  h1_max_itc_1slice->Write();
  h1_max_itc_time->Write();
  tree->Write();
}

/*
void itc_tools::TestCalcITC()
{
  double times[] = {1,5,5,5,5,5,15,15,15,15,15,24.5,25,25,25,25,25,35,35,35,35,35,51,52,52,53,53,500,502,502,799,800,801,1000};
  vector<double> v1(times, times + sizeof(times) / sizeof(double));
  vector<double> v2(times, times + sizeof(times) / sizeof(double));
  vector<double> v3(times, times + sizeof(times) / sizeof(double));
  CalcMaxITC(&v1, &v2, &v3);
}
*/
