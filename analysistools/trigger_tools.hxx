#ifndef trigger_tools_h
#define trigger_tools_h

#include <vector>

#include "TH1D.h"

#include "itc_tools.hxx"

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimEnumerations.hh"

using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::sort;

class trigger_tools
{
public:

  enum DigiType_t {
    kDigiTypeUndefined = -1,
    kDigiTypePhysics,
    kDigiTypeNoise,
    kDigiTypeMix,
    kDigiTypeError
  };
  
  trigger_tools();
  ~trigger_tools();

  void PopulateDigitTimes(WCSimRootTrigger * trigger, bool append, WCSimRootEvent * event = NULL);
  void CalcMaxITC(itc_tools * itc);
  void WriteToFile(TFile * f);
  void PrintDigitTimes(DigiType_t digitype = kDigiTypeUndefined);

  //vector<double> * GetDigitTimes(DigiType_t digitype);

private:
  void CleanupDigitTimes();
  void CreateDigitTimes();
  void SortDigitTimes();
  void FillDigitTimes(double digitime, DigiType_t digitype);
  DigiType_t GetDigitType(WCSimRootCherenkovDigiHit * wcsimrootcherenkovdigihit, WCSimRootEvent * event);
  vector<double> * digit_times;
  vector<double> * digit_times_physics;
  vector<double> * digit_times_noise;
  vector<double> * digit_times_mix;
};

#endif //trigger_tools_h
