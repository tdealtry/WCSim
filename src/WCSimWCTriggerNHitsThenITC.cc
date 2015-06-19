#include "WCSimWCTriggerNHitsThenITC.hh"
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

WCSimWCTriggerNHitsThenITC::WCSimWCTriggerNHitsThenITC(G4String name,
					 WCSimDetectorConstruction* myDetector,
					 WCSimWCDAQMessenger* myMessenger)
  :WCSimWCTriggerBase(name, myDetector, myMessenger)
{
}

WCSimWCTriggerNHitsThenITC::~WCSimWCTriggerNHitsThenITC(){
}

void WCSimWCTriggerNHitsThenITC::DoTheWork(WCSimWCDigitsCollection* WCDCPMT) {
  //Apply an NHitsThenITC trigger
  bool remove_hits = false;
  AlgNHitsThenITC(WCDCPMT, remove_hits);
}
