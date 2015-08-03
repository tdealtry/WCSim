/**
 * \class WCSimWCTriggerNHitsThenLocalNHits
 *
 * \brief A trigger class which first looks for a global NHits trigger, then searches for a local NHits trigger
 *
 */

#ifndef WCSimWCTriggerNHitsThenLocalNHits_h
#define WCSimWCTriggerNHitsThenLocalNHits_h 1

#include "WCSimWCTriggerBase.hh"

class WCSimWCTriggerNHitsThenLocalNHits : public WCSimWCTriggerBase
{
public:

  ///Create WCSimWCTriggerNHitsThenLocalNHits instance with knowledge of the detector and DAQ options
  WCSimWCTriggerNHitsThenLocalNHits(G4String name, WCSimDetectorConstruction*, WCSimWCDAQMessenger*);

  ~WCSimWCTriggerNHitsThenLocalNHits();
  
private:
  ///Calls the workhorse of this class: AlgNHitsThenLocalNHits
  void DoTheWork(WCSimWCDigitsCollection* WCDCPMT);

};

#endif








