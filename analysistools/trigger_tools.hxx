#ifndef trigger_tools_h
#define trigger_tools_h

#include <vector>

class trigger_tools
{
public:
  trigger_tools();
  ~trigger_tools();
  void PopulateDigitTimes(WCSimRootEvent * event, bool append, WCSimRootEvent * event = null);

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
  void FillDigitTimes(double digitime, int filltype);
  DigiType_t GetDigitType(WCSimRootCherenkovDigiHit * wcsimrootcherenkovdigihit, WCSimRootEvent * event);
  vector<double> * digit_times;
  vector<double> * digit_times_physics;
  vector<double> * digit_times_noise;
  vector<double> * digit_times_mix;
};

#endif //trigger_tools_h
