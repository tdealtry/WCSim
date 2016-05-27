#include "tmva_tools.hxx"

#include <cassert>
#include <algorithm>
#include <iostream>

#include "TMath.h"

using std::sort;
using std::cout;
using std::endl;
using std::cerr;

tmva_tools::tmva_tools(TFile * f, bool onetimeslice, int verbosity) :
  trigger_tools(f, onetimeslice, verbosity), filledfirstevent(false)
{
  distnames[0]  = "charge";
  distnames[1]  = "time";
  distnames[2]  = "position_q";
  distnames[3]  = "position_r";
  distnames[4]  = "position_z";
  distnames[5]  = "angle";
  distnames[6]  = "cosangle";
  distnames[7]  = "beta2";
  distnames[8]  = "beta3";
  distnames[9]  = "beta4";
  distnames[10] = "beta5";
  distnames[11] = "3hitangle";
  distnames[12] = "solar_anisotropy";

  quantnames[0] = "mean";
  quantnames[1] = "rms";
  quantnames[2] = "skew";
  quantnames[3] = "mode";
  quantnames[4] = "median" ;

  digit_angle = NULL;
  digit_cosangle = NULL;
  digit_beta2 = NULL;
  digit_beta3 = NULL;
  digit_beta4 = NULL;
  digit_beta5 = NULL;
  digit_3hitangle = NULL;
  digit_solar_anisotropy = NULL;

  CreateVectors();
  CreateTreeVectors();

  tree = new TTree(TString::Format("t_tmva"), TString::Format("TMVA input tree"));
  tree->Branch("nhits", &var_nhits);
  for(int id = 0; id < tmva_tools::ndistributions; id++) {
    for(int iq = 0; iq < tmva_tools::quantities_per_distribution; iq++) {
      tree->Branch(TString::Format("%s_%s", tmva_tools::distnames[id].c_str(), tmva_tools::quantnames[iq].c_str()), &(var_distributions[id][iq]));
    }//iq
  }//id
  tree->Branch("beta_14", &var_beta14);
  tree->Branch("solar_anisotropy_ratio", &var_solar_anisotropy_ratio);

  tree->Branch("TRUE_last_physics_hit",  &var_TRUE_last_physics_hit);
  tree->Branch("TRUE_fraction_physics_hit_in_window", &var_TRUE_fraction_physics_hit_in_window);
  tree->Branch("TRUE_e_energy",    &TRUE_e_energy);
  tree->Branch("TRUE_e_direction", &TRUE_e_direction);
  tree->Branch("TRUE_e_vertex",    &TRUE_e_vertex);

  tree->Branch("digit_charges", &tree_digit_charges);
  tree->Branch("digit_times", &tree_digit_times_nosort);
  tree->Branch("digit_position_q", &tree_digit_position_q);
  tree->Branch("digit_position_r", &tree_digit_position_r);
  tree->Branch("digit_position_z", &tree_digit_position_z);

  tree->Branch("digit_angle", &tree_digit_angle);
  tree->Branch("digit_cosangle", &tree_digit_cosangle);
  tree->Branch("digit_beta2", &tree_digit_beta2);
  tree->Branch("digit_beta3", &tree_digit_beta3);
  tree->Branch("digit_beta4", &tree_digit_beta4);
  tree->Branch("digit_beta5", &tree_digit_beta5);
  tree->Branch("digit_3hitangle", &tree_digit_3hitangle);
  tree->Branch("digit_solar_anisotropy", &tree_digit_solar_anisotropy);
}

tmva_tools::~tmva_tools()
{
  CleanupVectors();
  CleanupTreeVectors();
  delete tree;
  f = NULL;
}

void tmva_tools::CreateTreeVectors()
{
  tree_digit_angle = new vector<double>;
  tree_digit_cosangle = new vector<double>;
  tree_digit_beta2 = new vector<double>;
  tree_digit_beta3 = new vector<double>;
  tree_digit_beta4 = new vector<double>;
  tree_digit_beta5 = new vector<double>;
  tree_digit_3hitangle = new vector<double>;
  tree_digit_solar_anisotropy = new vector<double>;
  tree_digit_charges = new vector<double>;
  tree_digit_times_nosort = new vector<double>;
  tree_digit_position_q = new vector<double>;
  tree_digit_position_r = new vector<double>;
  tree_digit_position_z = new vector<double>;
}

void tmva_tools::CleanupTreeVectors()
{
  delete tree_digit_angle;
  delete tree_digit_cosangle;
  delete tree_digit_beta2;
  delete tree_digit_beta3;
  delete tree_digit_beta4;
  delete tree_digit_beta5;
  delete tree_digit_3hitangle;
  delete tree_digit_solar_anisotropy;
  delete tree_digit_charges;
  delete tree_digit_times_nosort;
  delete tree_digit_position_q;
  delete tree_digit_position_r;
  delete tree_digit_position_z;
}

void tmva_tools::ClearTreeVectors()
{
  tree_digit_angle->clear();
  tree_digit_cosangle->clear();
  tree_digit_beta2->clear();
  tree_digit_beta3->clear();
  tree_digit_beta4->clear();
  tree_digit_beta5->clear();
  tree_digit_3hitangle->clear();
  tree_digit_solar_anisotropy->clear();
  tree_digit_charges->clear();
  tree_digit_times_nosort->clear();
  tree_digit_position_q->clear();
  tree_digit_position_r->clear();
  tree_digit_position_z->clear();
}

void tmva_tools::FillTreeVectors()
{
  (*tree_digit_angle) = (*digit_angle);
  (*tree_digit_cosangle) = (*digit_cosangle);
  (*tree_digit_beta2) = (*digit_beta2);
  (*tree_digit_beta3) = (*digit_beta3);
  (*tree_digit_beta4) = (*digit_beta4);
  (*tree_digit_beta5) = (*digit_beta5);
  (*tree_digit_3hitangle) = (*digit_3hitangle);
  (*tree_digit_solar_anisotropy) = (*digit_solar_anisotropy);
  (*tree_digit_charges) = (*digit_charges);
  (*tree_digit_times_nosort) = (*digit_times_nosort);
  (*tree_digit_position_q) = (*digit_position_q);
  (*tree_digit_position_r) = (*digit_position_r);
  (*tree_digit_position_z) = (*digit_position_z);
}

void tmva_tools::CleanupVectors()
{
  if(digit_angle != NULL) {
    delete digit_angle;
    digit_angle = NULL;
  };
  if(digit_cosangle != NULL) {
    delete digit_cosangle;
    digit_cosangle = NULL;
  };
  if(digit_beta2 != NULL) {
    delete digit_beta2;
    digit_beta2 = NULL;
  };
  if(digit_beta3 != NULL) {
    delete digit_beta3;
    digit_beta3 = NULL;
  };
  if(digit_beta4 != NULL) {
    delete digit_beta4;
    digit_beta4 = NULL;
  };
  if(digit_beta5 != NULL) {
    delete digit_beta5;
    digit_beta5 = NULL;
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
  if(digit_cosangle == NULL) {
    digit_cosangle = new vector<double>;
  };
  if(digit_beta2 == NULL) {
    digit_beta2 = new vector<double>;
  };
  if(digit_beta3 == NULL) {
    digit_beta3 = new vector<double>;
  };
  if(digit_beta4 == NULL) {
    digit_beta4 = new vector<double>;
  };
  if(digit_beta5 == NULL) {
    digit_beta5 = new vector<double>;
  };
  if(digit_3hitangle == NULL) {
    digit_3hitangle = new vector<double>;
  };
  if(digit_solar_anisotropy == NULL) {
    digit_solar_anisotropy = new vector<double>;
  };
}

void tmva_tools::ClearVectors()
{
  digit_angle->clear();
  digit_cosangle->clear();
  digit_beta2->clear();
  digit_beta3->clear();
  digit_beta4->clear();
  digit_beta5->clear();
  digit_3hitangle->clear();
  digit_solar_anisotropy->clear();
}

void tmva_tools::FillTree(double duration)
{
  if(verbosity > 1)
    cout << "Starting tmva_tools::FillTree()" << endl;
  ClearVectors();

  if(!digit_times->size()) {
    cerr << "Not filling event; no digits found!" << endl;
    return;
  }

  if(verbosity > 2)
    cout << "tmva_tools::FillTree() Size of all, physics & mix time vectors are "
	 << digit_times->size() << "\t" << digit_times_physics->size() << "\t" << digit_times_mix->size() << endl;

  //Find the start of the window to make distributions from
  double initial_time = digit_times->at(0) + 10; //add a little, noise "smeared out" of window, but not "smeared in"
  //if it's a physics event, initial_time is the time of first physics
  if(digit_times_physics->size() && digit_times_mix->size()) {
    initial_time = TMath::Min(digit_times_physics->at(0), digit_times_mix->at(0));
  }
  else if( digit_times_physics->size() && !digit_times_mix->size()) {
    initial_time = digit_times_physics->at(0);
  }
  else if(!digit_times_physics->size() &&  digit_times_mix->size()) {
    initial_time = digit_times_mix->at(0);
  }

  if(verbosity > 2)
    cout << "tmva_tools::FillTree() Initial time is " << initial_time << endl;

  if(verbosity > 1)
    cout << "tmva_tools::FillTree() Going to erase information outside the time window" << endl;
  //first, remove entries not in the time window
  size_t counter = 0;
  for(int idigit = digit_times_nosort->size() - 1; idigit >= 0; idigit--) {
    double thistime = digit_times_nosort->at(idigit);
    if(thistime > (initial_time + duration) || thistime < initial_time) {
      //not in the tmva window
      if(verbosity > 3)
	cout << "erasing at position " << idigit << " with time " << thistime << endl;
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

  if(verbosity > 1)
    cout << "tmva_tools::FillTree() " << var_nhits << " digits in the analysis window" << endl;

  //cout << "NHITS (in window) is " << var_nhits << endl;

  FillTruthInfo(initial_time, duration);
  FillAngleVectors();
  FillAnisotropyVector();

  if(verbosity > 1)
    cout << "tmva_tools::FillTree() going to calculate distribution quantities for digit charge."
	 << " Size = " << digit_charges->size() << endl;
  GetDistributionQuantities(digit_charges,      var_distributions[0], false, 101);
  if(verbosity > 1)
    cout << "tmva_tools::FillTree() going to calculate distribution quantities for digit time."
      " Size = " << digit_times_nosort->size() << endl;
  GetDistributionQuantities(digit_times_nosort, var_distributions[1], true,  101);
  if(verbosity > 1)
    cout << "tmva_tools::FillTree() going to calculate distribution quantities for digit position theta."
      " Size = " << digit_position_q->size() << endl;
  GetDistributionQuantities(digit_position_q,   var_distributions[2], false, 101);
  if(verbosity > 1)
    cout << "tmva_tools::FillTree() going to calculate distribution quantities for digit position R."
      " Size = " << digit_position_r->size() << endl;
  GetDistributionQuantities(digit_position_r,   var_distributions[3], false, 101);
  if(verbosity > 1)
    cout << "tmva_tools::FillTree() going to calculate distribution quantities for digit position Z."
      " Size = " << digit_position_z->size() << endl;
  GetDistributionQuantities(digit_position_z,   var_distributions[4], false, 101);
  if(verbosity > 1)
    cout << "tmva_tools::FillTree() going to calculate distribution quantities for digit angle."
      " Size = " << digit_angle->size() << endl;
  GetDistributionQuantities(digit_angle,        var_distributions[5], false, 101);
  if(verbosity > 1)
    cout << "tmva_tools::FillTree() going to calculate distribution quantities for digit cos(angle)."
      " Size = " << digit_cosangle->size() << endl;
  GetDistributionQuantities(digit_cosangle,     var_distributions[6], false, 101);
  if(verbosity > 1)
    cout << "tmva_tools::FillTree() going to calculate distribution quantities for digit betas (2-5)."
      " Size = " << digit_beta2->size() << endl;
  GetDistributionQuantities(digit_beta2,       var_distributions[7], false, 101);
  GetDistributionQuantities(digit_beta3,       var_distributions[8], false, 101);
  GetDistributionQuantities(digit_beta4,       var_distributions[9], false, 101);
  GetDistributionQuantities(digit_beta5,       var_distributions[10], false, 101);
  //https://github.com/mastbaum/gamma_discrimination/blob/master/beta14/python/beta14.py
  var_beta14 = var_distributions[7][0] + 4 * var_distributions[7][1];
  if(verbosity > 2)
    cout << "tmva_tools::FillTree() beta14 is " << var_beta14 << endl;
  if(verbosity > 1)
    cout << "tmva_tools::FillTree() going to calculate distribution quantities for digit 3hitangle."
	 << " Size = " << digit_3hitangle->size() << endl;
  GetDistributionQuantities(digit_3hitangle,    var_distributions[11], false, 101);
  if(verbosity > 1)
    cout << "tmva_tools::FillTree() going to calculate distribution quantities for digit solar anisotropy."
	 << " Size = " << digit_solar_anisotropy->size() << endl;
  GetDistributionQuantities(digit_solar_anisotropy, var_distributions[12], false, 101);

  if(!filledfirstevent) {
    if(verbosity > 1)
      cout << "tmva_tools::FillTree() calling FillTreeVectors()" << endl;
    FillTreeVectors();
    filledfirstevent = true;
  }

  if(verbosity > 1)
    cout << "tmva_tools::FillTree() calling TTree::Fill()" << endl;
  tree->Fill();

  ClearTreeVectors();
}

void tmva_tools::FillTruthInfo(double initial_time, double duration)
{
  if(verbosity > 1)
    cout << "Starting FillTruthInfo()" << endl;
  double last = -9999999;
  int nphysics_window = 0;
  const int nphysics_total = digit_times_physics->size() + digit_times_mix->size();
  if(!nphysics_total) {
    var_TRUE_fraction_physics_hit_in_window = 0;
    var_TRUE_last_physics_hit = initial_time;
    if(verbosity > 1)
      cout << "tmva_tools::FillTruthInfo() No physics hits found." 
	   << " Using defaults for physics fraction & last physics hit time" << endl;
    return;
  }
  for(vector<double>::iterator it = digit_times_physics->begin(); it != digit_times_physics->end(); ++it) {
    double thistime = *it;
    if(thistime < (initial_time + duration) && thistime > initial_time) {
      nphysics_window++;
    }
    if(thistime > last)
      last = thistime;
  }
  for(vector<double>::iterator it = digit_times_mix->begin(); it != digit_times_mix->end(); ++it) {
    double thistime = *it;
    if(thistime < (initial_time + duration) && thistime > initial_time) {
      nphysics_window++;
    }
    if(thistime > last)
      last = thistime;
  }
  var_TRUE_fraction_physics_hit_in_window = (double)nphysics_window / (double)nphysics_total;
  var_TRUE_last_physics_hit = last;
  if(verbosity > 1)
    cout << "tmva_tools::FillTruthInfo() has fraction of physics hits used "
	 << var_TRUE_fraction_physics_hit_in_window << " = " << nphysics_window << " / " << nphysics_total << endl
	 << "tmva_tools::FillTruthInfo() has last physics hit at "
	 << var_TRUE_last_physics_hit << endl;
}

void tmva_tools::FillAngleVectors()
{
  if(verbosity > 1)
    cout << "Starting FillAngleVectors()" << endl;
  for(int idigit = 0; idigit < var_nhits; idigit++) {
    for(int jdigit = 0; jdigit < var_nhits; jdigit++) {
      if(jdigit > idigit) {
	double angle = digit_pmtvec->at(idigit).Angle(digit_pmtvec->at(jdigit));
	digit_angle->push_back(angle); //2pmt angle
	//http://mathworld.wolfram.com/LegendrePolynomial.html
	digit_cosangle->push_back(TMath::Cos(angle)); //cosine of 2pmt angle (1st order Legendre polynomial of cos(angle) 2pmt angle)
	digit_beta2->push_back(0.5   * (3 *pow(angle,2) - 1)); //2nd order Legendre polynomial of cos(angle) 2pmt angle
	digit_beta3->push_back(0.5   * (5 *pow(angle,3) - 3*angle)); //3rd order Legendre polynomial of cos(angle) 2pmt angle
	digit_beta4->push_back(0.125 * (35*pow(angle,4) - 30*pow(angle,2) + 3)); //4th order Legendre polynomial of cos(angle) 2pmt angle
	digit_beta5->push_back(0.125 * (63*pow(angle,5) - 70*pow(angle,3) + 15*angle)); //5th order Legendre polynomial of cos(angle) 2pmt angle
      }
      else if(idigit == jdigit)
	continue;
      for(int kdigit = 0; kdigit < var_nhits; kdigit++) {
	break;
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
  if(verbosity > 1)
    cout << "Starting FillAnisotropyVector()" << endl;
  int nplane_hi = 0, nplane_lo = 0;
  double a  = TRUE_e_direction.X();
  double b  = TRUE_e_direction.Y();
  double c  = TRUE_e_direction.Z();
  double denominator = TMath::Sqrt(a*a + b*b + c*c);
  if(verbosity > 2)
    cout << "tmva_tools::FillAnisotropyVector() has a,b,c "
	 << a << "," << b << "," << c << "  sqrt(a^2 + b^2 + c^2) = "
	 << denominator << endl;
  for(int idigit = 0; idigit < var_nhits; idigit++) {
    //http://mathworld.wolfram.com/Point-PlaneDistance.html
    double x0 = digit_pmtvec->at(idigit).X();
    double y0 = digit_pmtvec->at(idigit).Y();
    double z0 = digit_pmtvec->at(idigit).Z();
    double plane_distance = (a*x0 + b*y0 + c*z0) / denominator;
    if(verbosity > 3)
      cout << "tmva_tools::FillAnisotropyVector() for idigit = " << idigit
	   << " has x0,y0,z0 " << x0 << "," << y0 << "," << z0
	   << " and plane_distance " << plane_distance << endl;
    (plane_distance < 0) ? nplane_lo++ : nplane_hi++;
    digit_solar_anisotropy->push_back(plane_distance);
  }//idigit
  var_solar_anisotropy_ratio = (double)nplane_hi / (double)nplane_lo;
  if(verbosity > 2)
    cout << "tmva_tools::FillAnisotropyVector() Event has ratio "
	 << var_solar_anisotropy_ratio << " = " << nplane_hi << " / " << nplane_lo << endl;
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
