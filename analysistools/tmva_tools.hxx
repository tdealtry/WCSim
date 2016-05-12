#ifndef tmva_tools_h
#define tmva_tools_h

#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

#include "trigger_tools.hxx"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::string;

class tmva_tools : public trigger_tools
{
public:
  tmva_tools(TFile * f, bool onetimeslice);
  ~tmva_tools();

  void FillTree(double duration);
  void GetDistributionQuantities(vector<double> * v, double * quantities, const bool mean_subtraction, const size_t nbins);


  void Write();

private:
  void CleanupVectors();
  void CreateVectors();

  void FillAngleVectors();
  void FillAnisotropyVector();
  void FillTruthInfo(double initial_time, double duration);

  static const int ndistributions = 8;
  static const int quantities_per_distribution = 5;
  string distnames[ndistributions];
  string quantnames[quantities_per_distribution];

  int    var_nhits;
  double var_distributions[ndistributions][quantities_per_distribution];
  double var_solar_anisotropy_ratio; // calculate NHITS ratio between close & far (tank split in two)
  //truth
  double var_TRUE_last_physics_hit; //time of last physics hit, relative to start of window (defined as first TRUE physics hit)
  double var_TRUE_fraction_physics_hit_in_window; //fraction of the physics hits that are used in the histograms

  vector<double> * digit_angle;
  vector<double> * digit_3hitangle;
  vector<double> * digit_solar_anisotropy;

  TTree * tree;
  TFile * f;

  bool one_time_slice;
};

#endif //tmva_tools_h
