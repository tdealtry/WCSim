#ifndef WCSimWCTriggerNHitsThenITC_h
#define WCSimWCTriggerNHitsThenITC_h 1

#include "WCSimWCTriggerBase.hh"

class WCSimWCTriggerNHitsThenITC : public WCSimWCTriggerBase
{
public:

  //not recommended to override these methods
  WCSimWCTriggerNHitsThenITC(G4String name, WCSimDetectorConstruction*, WCSimWCDAQMessenger*);
  ~WCSimWCTriggerNHitsThenITC();
  
private:
  void DoTheWork(WCSimWCDigitsCollection* WCDCPMT);

};

#endif








