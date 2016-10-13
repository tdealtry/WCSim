#ifndef tmva_tools_h
#define tmva_tools_h

#include <string>

#include "trigger_tools.hxx"

#include "TTree.h"

using std::string;

class tmva_tools : public trigger_tools
{
public:
  tmva_tools(TFile * f, bool onetimeslice, int verbosity);
  ~tmva_tools();

  void FillTree(double duration);
  void GetDistributionQuantities(vector<double> * v, double * quantities, const bool mean_subtraction, const size_t nbins);


  void Write();

private:
  void CleanupVectors();
  void CreateVectors();
  void ClearVectors();

  void CreateTreeVectors();
  void CleanupTreeVectors();
  void ClearTreeVectors();
  void FillTreeVectors();

  void FillAngleVectors();
  void FillAnisotropyVector();
  void FillTruthInfo(double initial_time, double duration);

  static const int ndistributions = 13;
  static const int quantities_per_distribution = 5;
  string distnames[ndistributions];
  string quantnames[quantities_per_distribution];

  int    var_nhits;
  double var_distributions[ndistributions][quantities_per_distribution];
  double var_beta14;
  double var_solar_anisotropy_ratio; // calculate NHITS ratio between close & far (tank split in two)
  //truth
  double var_TRUE_last_physics_hit; //time of last physics hit, relative to start of window (defined as first TRUE physics hit)
  double var_TRUE_fraction_physics_hit_in_window; //fraction of the physics hits that are used in the histograms

  vector<double> * digit_angle;
  vector<double> * digit_cosangle;
  vector<double> * digit_beta2;
  vector<double> * digit_beta3;
  vector<double> * digit_beta4;
  vector<double> * digit_beta5;
  vector<double> * digit_3hitangle;
  vector<double> * digit_solar_anisotropy;

  vector<double> * tree_digit_angle;
  vector<double> * tree_digit_cosangle;
  vector<double> * tree_digit_beta2;
  vector<double> * tree_digit_beta3;
  vector<double> * tree_digit_beta4;
  vector<double> * tree_digit_beta5;
  vector<double> * tree_digit_3hitangle;
  vector<double> * tree_digit_solar_anisotropy;
  vector<double> * tree_digit_charges;
  vector<double> * tree_digit_times_nosort;
  vector<double> * tree_digit_position_q;
  vector<double> * tree_digit_position_r;
  vector<double> * tree_digit_position_z;

  TTree * tree;
  bool filledfirstevent;
};

#endif //tmva_tools_h
