#ifndef trigger_tools_h
#define trigger_tools_h

#include <vector>

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimEnumerations.hh"

using std::vector;

class trigger_tools
{
public:
  trigger_tools();
  ~trigger_tools();
  void PopulateDigitTimes(WCSimRootTrigger * trigger, bool append, WCSimRootEvent * event = NULL);

  enum DigiType_t {
    kDigiTypeUndefined = -1,
    kDigiTypePhysics,
    kDigiTypeNoise,
    kDigiTypeMix,
    kDigiTypeError
  };
  
private:
  void CleanupDigitTimes();
  void CreateDigitTimes();
  void FillDigitTimes(double digitime, DigiType_t digitype);
  DigiType_t GetDigitType(WCSimRootCherenkovDigiHit * wcsimrootcherenkovdigihit, WCSimRootEvent * event);
  vector<double> * digit_times;
  vector<double> * digit_times_physics;
  vector<double> * digit_times_noise;
  vector<double> * digit_times_mix;
};

#endif //trigger_tools_h
