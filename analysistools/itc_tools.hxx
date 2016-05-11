#ifndef itc_tools_h
#define itc_tools_h

#include <vector>
#include <algorithm>
#include <iostream>

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"

using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::sort;

class itc_tools
{
public:
  itc_tools(int smallwindow, int largewindow, int offset, bool calculatendigits, bool onetimeslice);
  ~itc_tools();

  void CalcMaxITC(vector<double> * digit_times, vector<double> * digit_times_physics, vector<double> * digit_times_mix);
  void TestCalcITC();

  void WriteToFile(TFile * f);

private:
  TH1D * h1_max_ndigits;
  TH1D * h1_max_ndigits_time;
  TH1D * h1_max_ndigits_1slice;
  TH1D * h1_max_itc;
  TH1D * h1_max_itc_time;
  TH1D * h1_max_itc_1slice;

  double var_max_ndigits;
  double var_max_ndigits_time;
  double var_max_ndigits_1slice;
  double var_max_itc;
  double var_max_itc_time;
  double var_max_itc_1slice;

  TTree * tree;
  
  int smallwindow;
  int largewindow;
  int offset;
  bool calcndigits;
  bool one_time_slice;
};

#endif //itc_tools_h
