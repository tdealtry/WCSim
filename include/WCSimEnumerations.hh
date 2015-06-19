#ifndef WCSimEnumerations_h
#define WCSimEnumerations_h 1

#include <string>

typedef enum ETriggerType {
  kTriggerUndefined = -1,
  kTriggerNHits,
  kTriggerNHitsSKDETSIM,
  kTriggerNHitsTest,
  kTriggerITCRatio
} TriggerType_t;

class WCSimEnumerations
{
public:

  static std::string EnumAsString(TriggerType_t t)
  {
    switch(t) {
    case (kTriggerNHits) :
      return "NHits";
      break;
    case (kTriggerNHitsSKDETSIM) :
      return "NHits_SKDETSIM";
      break;
    case (kTriggerNHitsTest) :
      return "NHits_TEST";
      break;
    case (kTriggerITCRatio) :
      return "ITCRatio";
      break;
    default:
      return "";
      break;
    }
    return "";
  }

};

#endif
