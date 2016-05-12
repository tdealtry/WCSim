#include "tmva_tools.hxx"

#include <cassert>

#include "TMath.h"

#ifndef __DIGIT_TIME_VERBOSE__
//#define __DIGIT_TIME_VERBOSE__
#endif

tmva_tools::tmva_tools(TFile * f, bool onetimeslice) :
  f(f), one_time_slice(onetimeslice), trigger_tools()
{
  distnames[0] = "charge";
  distnames[1] = "time";
  distnames[2] = "position_q";
  distnames[3] = "position_r";
  distnames[4] = "position_z";
  distnames[5] = "angle";
  distnames[6] = "3hitangle";
  distnames[7] = "solar_anisotropy";

  quantnames[0] = "mean";
  quantnames[1] = "rms";
  quantnames[2] = "skew";
  quantnames[3] = "mode";
  quantnames[4] = "median" ;
  
  f->cd();

  digit_angle = NULL;
  digit_3hitangle = NULL;
  digit_solar_anisotropy = NULL;

  tree = new TTree(TString::Format("t_tmva"), TString::Format("TMVA input tree"));
  tree->Branch("nhits", &var_nhits);
  for(int id = 0; id < tmva_tools::ndistributions; id++) {
    for(int iq = 0; iq < tmva_tools::quantities_per_distribution; iq++) {
      tree->Branch(TString::Format("%s_%s", tmva_tools::distnames[id].c_str(), tmva_tools::quantnames[iq].c_str()), &(var_distributions[id][iq]));
    }//iq
  }//id
  tree->Branch("solar_anisotropy_ratio", &var_solar_anisotropy_ratio);
  tree->Branch("TRUE_last_physics_hit",  &var_TRUE_last_physics_hit);
  tree->Branch("TRUE_fraction_physics_hit_in_window", &var_TRUE_fraction_physics_hit_in_window);
  tree->Branch("TRUE_e_energy",    &TRUE_e_energy);
  tree->Branch("TRUE_e_direction", &TRUE_e_direction);
  tree->Branch("TRUE_e_vertex",    &TRUE_e_vertex);
}

tmva_tools::~tmva_tools()
{
  CleanupVectors();
  delete tree;
}

void tmva_tools::CleanupVectors()
{
  if(digit_angle != NULL) {
    delete digit_angle;
    digit_angle = NULL;
  };
  if(digit_3hitangle != NULL) {
    delete digit_3hitangle;
    digit_3hitangle = NULL;
  };
  if(digit_solar_anisotropy != NULL) {
    delete digit_solar_anisotropy;
    digit_solar_anisotropy = NULL;
  };
}

void tmva_tools::CreateVectors()
{
  if(digit_angle == NULL) {
    digit_angle = new vector<double>;
  };
  if(digit_3hitangle == NULL) {
    digit_3hitangle = new vector<double>;
  };
  if(digit_solar_anisotropy == NULL) {
    digit_solar_anisotropy = new vector<double>;
  };
}

void tmva_tools::FillTree(double duration)
{
  CleanupVectors();
  CreateVectors();

  //Find the start of the window to make distributions from
  double initial_time = digit_times->at(0) + 10; //add a little, noise "smeared out" of window, but not "smeared in"
  //if it's a physics event, initial_time is the time of first physics
  if(digit_times_physics->size() && digit_times_mix->size())
    initial_time = TMath::Min(digit_times_physics->at(0), digit_times_mix->at(0));
  else if( digit_times_physics->size() && !digit_times_mix->size())
    initial_time = digit_times_physics->at(0);
  else if(!digit_times_physics->size() &&  digit_times_mix->size())
    initial_time = digit_times_mix->at(0);

  //cout << "Initial time is " << initial_time << endl;

  //first, remove entries not in the time window
  size_t counter = 0;
  for(int idigit = digit_times_nosort->size() - 1; idigit >= 0; idigit--) {
    double thistime = digit_times_nosort->at(idigit);
    if(thistime > (initial_time + duration) || thistime < initial_time) {
      //not in the tmva window
      //cout << "erasing at position " << idigit << " with time " << thistime << endl;
      digit_times_nosort->erase(digit_times_nosort->begin() + idigit);
      digit_charges->erase(digit_charges->begin() + idigit);
      digit_ipmt->erase(digit_ipmt->begin() + idigit);
      digit_pmtvec->erase(digit_pmtvec->begin() + idigit);
      digit_position_q->erase(digit_position_q->begin() + idigit);
      digit_position_r->erase(digit_position_r->begin() + idigit);
      digit_position_z->erase(digit_position_z->begin() + idigit);
    }
    else
      counter++;
  }//idigit
  assert(counter == digit_times_nosort->size());
  var_nhits = counter;

  //cout << "NHITS (in window) is " << var_nhits << endl;

  FillTruthInfo(initial_time, duration);
  FillAngleVectors();
  FillAnisotropyVector();

  GetDistributionQuantities(digit_charges,      var_distributions[0], false, 101);
  GetDistributionQuantities(digit_times_nosort, var_distributions[1], true,  101);
  GetDistributionQuantities(digit_position_q,   var_distributions[2], false, 101);
  GetDistributionQuantities(digit_position_r,   var_distributions[3], false, 101);
  GetDistributionQuantities(digit_position_z,   var_distributions[4], false, 101);
  GetDistributionQuantities(digit_angle,        var_distributions[5], false, 101);
  GetDistributionQuantities(digit_3hitangle,    var_distributions[6], false, 101);
  GetDistributionQuantities(digit_solar_anisotropy, var_distributions[7], false, 101);

  tree->Fill();
}

void tmva_tools::FillTruthInfo(double initial_time, double duration)
{
  double last = -9999999;
  int nphysics_window = 0;
  const int nphysics_total = digit_times_physics->size() + digit_times_mix->size();
  for(vector<double>::iterator it = digit_times_physics->begin(); it != digit_times_physics->end(); ++it) {
    double thistime = *it;
    if(thistime < (initial_time + duration) && thistime > initial_time) {
      nphysics_window++;
      if(thistime > last)
	last = thistime;
    }
  }
  for(vector<double>::iterator it = digit_times_mix->begin(); it != digit_times_mix->end(); ++it) {
    double thistime = *it;
    if(thistime < (initial_time + duration) && thistime > initial_time) {
      nphysics_window++;
      if(thistime > last)
	last = thistime;
    }
  }
  var_TRUE_fraction_physics_hit_in_window = (double)nphysics_window / (double)nphysics_total;
  var_TRUE_last_physics_hit = last;
}

void tmva_tools::FillAngleVectors()
{
  for(int idigit = 0; idigit < var_nhits; idigit++) {
    for(int jdigit = 0; jdigit < var_nhits; jdigit++) {
      if(jdigit > idigit)
	digit_angle->push_back(digit_pmtvec->at(idigit).Angle(digit_pmtvec->at(jdigit))); //2pmt angle
      else if(idigit == jdigit)
	continue;
      for(int kdigit = 0; kdigit < var_nhits; kdigit++) {
	if(kdigit == idigit || kdigit == jdigit)
	  continue;
	double l1 = (digit_pmtvec->at(idigit) - digit_pmtvec->at(jdigit)).Mag();
	double l2 = (digit_pmtvec->at(jdigit) - digit_pmtvec->at(kdigit)).Mag();
	double l3 = (digit_pmtvec->at(kdigit) - digit_pmtvec->at(idigit)).Mag();
	double threehitangle = TMath::ACos((l1*l1 + l2*l2 - l3*l3) / (2 * l1 * l2));
	digit_3hitangle->push_back(threehitangle);
      }//kdigit
    }//jdigit
  }//idigit
}

void tmva_tools::FillAnisotropyVector()
{
  int nplane_hi = 0, nplane_lo = 0;
  for(int idigit = 0; idigit < var_nhits; idigit++) {
    //http://mathworld.wolfram.com/Point-PlaneDistance.html
    double x0 = digit_pmtvec->at(idigit).X();
    double y0 = digit_pmtvec->at(idigit).Y();
    double z0 = digit_pmtvec->at(idigit).Z();
    double a  = TRUE_e_direction.X();
    double b  = TRUE_e_direction.Y();
    double c  = TRUE_e_direction.Z();
    double plane_distance = (a*x0 + b*y0 + c*z0) / TMath::Sqrt(a*a + b*b + c*c);
    (plane_distance < 0) ? nplane_lo++ : nplane_hi++;
    digit_solar_anisotropy->push_back(plane_distance);
  }//idigit
  var_solar_anisotropy_ratio = (double)nplane_hi / (double)nplane_lo;
}

void tmva_tools::Write()
{
  f->cd();
  tree->Write();
}

void tmva_tools::GetDistributionQuantities(vector<double> * v, double * quantities, const bool mean_subtraction, const size_t nbins)
{
  vector<int> histogram(nbins, 0);
  const size_t n = v->size();
  for(int i = 0; i < tmva_tools::quantities_per_distribution; i++)
    quantities[i] = -99999;
  if(!n) {
    return;
  }
  double min = 99999999, max = -99999999;

  //mean
  double mean = 0;
  for(size_t i = 0; i < n; i++) {
    double val = v->at(i);
    mean += val;
    if(val < min) min = val;
    if(val > max) max = val;
  }//i
  mean /= n;

  if(mean_subtraction) {
    min =  99999999;
    max = -99999999;
    for(size_t i = 0; i < n; i++) {
      v->at(i) = v->at(i) - mean;
      if(v->at(i) < min) min = v->at(i);
      if(v->at(i) > max) max = v->at(i);
    }
    mean = 0;
  }

  //rms
  double rms = 0;
  for(size_t i = 0; i < n; i++) {
    double val = v->at(i);
    rms += (val-mean) * (val-mean);
  }//i
  rms /= n - 1;
  rms = TMath::Sqrt(rms);

  //skew
  double skew = 0;
  for(size_t i = 0; i < n; i++) {
    double val = v->at(i);
    skew += (val-mean) * (val-mean) * (val-mean);
  }//i
  skew /= rms * rms * rms;
  skew *= n / ((n-1) * (n-2));

  //mode
  double binsize = (max - min) / nbins;
  for(size_t i = 0; i < n; i++) {
    double val = v->at(i);
    size_t binnum = floor((val - min) / binsize);
    binnum = TMath::Min(binnum, nbins - 1);
    //cout << val << "\tbinnum " << binnum << "\tmin " << min << "\tmax " << max << "\tbinsize " << binsize << "\tnbins " << nbins << endl;
    histogram.at(binnum)++;
  }//i
  double mode = *(std::max_element(histogram.begin(), histogram.end()));
  mode *= binsize;

  //median
  std::sort(v->begin(), v->end());
  double median = v->at(n / 2);

  quantities[0] = mean;
  quantities[1] = rms;
  quantities[2] = skew;
  quantities[3] = mode;
  quantities[4] = median;

  return;
}
