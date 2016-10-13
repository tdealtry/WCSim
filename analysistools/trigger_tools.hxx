#ifndef trigger_tools_h
#define trigger_tools_h

#include <vector>

#include "TH1D.h"
#include "TVector3.h"
#include "TFile.h"

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimEnumerations.hh"

using std::vector;

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
  
  trigger_tools(TFile * f, bool onetimeslice, int verbosity);
  virtual ~trigger_tools();

  void PopulateDigitTimes(WCSimRootTrigger * trigger, bool append, WCSimRootEvent * event = NULL);
  void PopulateVectors(WCSimRootTrigger * trigger, bool append, WCSimRootGeom * geo);
  void PopulateTruthGun(WCSimRootTrigger * trigger);

  void PrintDigitTimes(DigiType_t digitype = kDigiTypeUndefined);

  vector<double> * GetDigitTimes(DigiType_t digitype);
  void CopyDigitTimes(trigger_tools * t);

private:
  void CleanupVectors();
  void CreateVectors();
  void ClearVectors();

  void CleanupDigitTimes();
  void CreateDigitTimes();
  void ClearDigitTimes();
  void SortDigitTimes();
  void FillDigitTimes(double digitime, DigiType_t digitype);
  DigiType_t GetDigitType(WCSimRootCherenkovDigiHit * wcsimrootcherenkovdigihit, WCSimRootEvent * event);

protected:
  TFile * f;
  bool one_time_slice;
  int verbosity;
  
  vector<double> * digit_times;
  vector<double> * digit_times_physics;
  vector<double> * digit_times_noise;
  vector<double> * digit_times_mix;

  vector<double>   * digit_charges;
  vector<double>   * digit_times_nosort;
  vector<double>   * digit_ipmt;
  vector<TVector3> * digit_pmtvec;
  vector<double>   * digit_position_q;
  vector<double>   * digit_position_r;
  vector<double>   * digit_position_z;

  double TRUE_e_vertex_x;
  double TRUE_e_vertex_y;
  double TRUE_e_vertex_z;
  double TRUE_e_direction_x;
  double TRUE_e_direction_y;
  double TRUE_e_direction_z;
  double TRUE_e_energy;
};

#endif //trigger_tools_h
