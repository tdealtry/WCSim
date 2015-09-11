#ifndef WCSimEnumerations_h
#define WCSimEnumerations_h 1

#include <string>
#include <iostream>

typedef enum ETriggerType {
  kTriggerUndefined = -1,
  kTriggerNHits,
  kTriggerLocalNHits,
  kTriggerNHitsSKDETSIM,
  kTriggerNHitsTest,
  kTriggerFailure // this should always be the last entry (for looping)
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
    case (kTriggerLocalNHits) :
      return "Local_NHits";
      break;
    case (kTriggerNHitsSKDETSIM) :
      return "NHits_SKDETSIM";
      break;
    case (kTriggerNHitsTest) :
      return "NHits_TEST";
      break;
    case (kTriggerFailure) :
      return "No_trigger_passed";
      break;
    default:
      return "";
      break;
    }
    return "";
  }

  static TriggerType_t TriggerTypeFromString(std::string s)
  {
    for(int i = int(kTriggerUndefined)+1; i <= kTriggerFailure; i++) {
      if(s.compare(WCSimEnumerations::EnumAsString((TriggerType_t)i)) == 0) {
	return (TriggerType_t)i;
      }
    }
    std::cerr << "WCSimEnumerations::TriggerTypeFromString() Unknown string value " << s << std::endl;
    return kTriggerUndefined;
  }

};

#endif
