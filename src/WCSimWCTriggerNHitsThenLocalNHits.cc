#include "WCSimWCTriggerNHitsThenLocalNHits.hh"
#include "WCSimWCPMT.hh"
#include "WCSimWCDigi.hh"
#include "WCSimWCHit.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4ios.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "WCSimDetectorConstruction.hh"
#include "WCSimPmtInfo.hh"
#include "WCSimDarkRateMessenger.hh"

#include <vector>
// for memset
#include <cstring>
#include <iostream>

WCSimWCTriggerNHitsThenLocalNHits::WCSimWCTriggerNHitsThenLocalNHits(G4String name,
					 WCSimDetectorConstruction* myDetector,
					 WCSimWCDAQMessenger* myMessenger)
  :WCSimWCTriggerBase(name, myDetector, myMessenger)
{
  triggerClassName = "NHitsThenLocalNHits";
}

WCSimWCTriggerNHitsThenLocalNHits::~WCSimWCTriggerNHitsThenLocalNHits(){
}

void WCSimWCTriggerNHitsThenLocalNHits::DoTheWork(WCSimWCDigitsCollection* WCDCPMT) {
  //Apply an NHitsThenLocalNHits trigger
  bool remove_hits = false;
  AlgNHitsThenLocalNHits(WCDCPMT, remove_hits);
}
