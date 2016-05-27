#ifndef itc_tools_h
#define itc_tools_h

#include "TTree.h"

#include "trigger_tools.hxx"

class itc_tools : public trigger_tools
{
public:
  itc_tools(TFile * f, int smallwindow, int largewindow, int offset, bool calculatendigits, bool onetimeslice, int verbosity);
  ~itc_tools();

  void CalcMaxITC();
  //void TestCalcITC();

  void Write();

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
};

#endif //itc_tools_h
