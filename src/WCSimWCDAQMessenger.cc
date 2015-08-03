#include "WCSimWCDAQMessenger.hh"
#include "WCSimEventAction.hh"
#include "WCSimWCDigitizerBase.hh"
#include "WCSimWCTriggerBase.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"

#include <string>

WCSimWCDAQMessenger::WCSimWCDAQMessenger()
{
  WCSimDAQDir = new G4UIdirectory("/DAQ/");
  WCSimDAQDir->SetGuidance("Commands to select DAQ options");

  DigitizerChoice = new G4UIcmdWithAString("/DAQ/Digitizer", this);
  DigitizerChoice->SetGuidance("Set the Digitizer type");
  DigitizerChoice->SetGuidance("Available choices are:\n"
			       "SKI\n"
			       "SKIV\n"
			       "SKI_SKDETSIM (combined trigger & digitization (therefore ignores /DAQ/Trigger); buggy) \n"
			       );
  DigitizerChoice->SetParameterName("Digitizer", false);
  DigitizerChoice->SetCandidates(
				 "SKI "
				 "SKIV "
				 "SKI_SKDETSIM "
				 );
  DigitizerChoice->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  TriggerChoice = new G4UIcmdWithAString("/DAQ/Trigger", this);
  TriggerChoice->SetGuidance("Set the Trigger type");
  TriggerChoice->SetGuidance("Available choices are:\n"
			     "NHits\n"
			     "NHits2\n"
			     "NHitsThenITC\n"
			     "SKI_SKDETSIM (combined trigger & digitization (therefore ignores /DAQ/Digitization); buggy) \n"
			     );
  TriggerChoice->SetParameterName("Trigger", false);
  TriggerChoice->SetCandidates(
			       "NHits "
			       "NHits2 "
			       "NHitsThenITC "
			       "SKI_SKDETSIM "
			       );
  TriggerChoice->AvailableForStates(G4State_PreInit, G4State_Idle);

  //Save failure trigger specific options
  SaveFailuresTriggerDir = new G4UIdirectory("/DAQ/TriggerSaveFailures/");
  SaveFailuresTriggerDir->SetGuidance("Commands specific to the Save Failures trigger");

  SaveFailuresTriggerMode = new G4UIcmdWithAnInteger("/DAQ/TriggerSaveFailures/Mode", this);
  SaveFailuresTriggerMode->SetGuidance("0: save only triggered events; 1: save both triggered and failed events; 2: save only failed events");
  SaveFailuresTriggerMode->SetParameterName("SaveFailuresMode",true);
  SaveFailuresTriggerMode->SetDefaultValue(0);
  StoreSaveFailuresMode = 0;

  SaveFailuresTriggerTime = new G4UIcmdWithADouble("/DAQ/TriggerSaveFailures/TriggerTime", this);
  SaveFailuresTriggerTime->SetGuidance("The trigger time for the events which failed other triggers");
  SaveFailuresTriggerTime->SetParameterName("SaveFailuresTime",true);
  SaveFailuresTriggerTime->SetDefaultValue(100);
  StoreSaveFailuresTime = 100;

  //NHits trigger specifc options
  NHitsTriggerDir = new G4UIdirectory("/DAQ/TriggerNHits/");
  NHitsTriggerDir->SetGuidance("Commands specific to the NHits trigger");
  
  NHitsTriggerThreshold = new G4UIcmdWithAnInteger("/DAQ/TriggerNHits/Threshold", this);
  NHitsTriggerThreshold->SetGuidance("Set the NHits trigger threshold");
  NHitsTriggerThreshold->SetParameterName("NHitsThreshold",false);
  NHitsTriggerThreshold->SetDefaultValue(25);
  StoreNHitsThreshold = 25;

  NHitsTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerNHits/Window", this);
  NHitsTriggerWindow->SetGuidance("Set the NHits trigger window (in ns)");
  NHitsTriggerWindow->SetParameterName("NHitsWindow",false);
  NHitsTriggerWindow->SetDefaultValue(200);
  StoreNHitsWindow = 200;

  NHitsTriggerAdjustForNoise = new G4UIcmdWithABool("/DAQ/TriggerNHits/AdjustForNoise", this);
  NHitsTriggerAdjustForNoise->SetGuidance("Adjust the NHits trigger threshold automatically dependent on the average noise rate");
  NHitsTriggerAdjustForNoise->SetParameterName("NHitsAdjustForNoise",true);
  NHitsTriggerAdjustForNoise->SetDefaultValue(false);
  StoreNHitsAdjustForNoise = false;

  //ITC Ratio trigger specifc options
  ITCRatioTriggerDir = new G4UIdirectory("/DAQ/TriggerITCRatio/");
  ITCRatioTriggerDir->SetGuidance("Commands specific to the ITCRatio trigger");
  
  ITCRatioTriggerThreshold = new G4UIcmdWithADouble("/DAQ/TriggerITCRatio/Threshold", this);
  ITCRatioTriggerThreshold->SetGuidance("Set the ITCRatio trigger threshold");
  ITCRatioTriggerThreshold->SetParameterName("ITCRatioThreshold",true);
  ITCRatioTriggerThreshold->SetDefaultValue(0.3);
  StoreITCRatioTriggerThreshold = 0.3;

  ITCRatioTriggerSmallWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerITCRatio/SmallWindow", this);
  ITCRatioTriggerSmallWindow->SetGuidance("Set the ITCRatio small trigger window (in ns) i.e. for the numerator");
  ITCRatioTriggerSmallWindow->SetParameterName("ITCRatioSmallWindow",true);
  ITCRatioTriggerSmallWindow->SetDefaultValue(200);
  StoreITCRatioTriggerSmallWindow = 200;

  ITCRatioTriggerLargeWindowLow = new G4UIcmdWithAnInteger("/DAQ/TriggerITCRatio/LargeWindowLow", this);
  ITCRatioTriggerLargeWindowLow->SetGuidance("Set the ITCRatio large trigger window low edge (in ns) i.e. for the denominator");
  ITCRatioTriggerLargeWindowLow->SetParameterName("ITCRatioLargeWindowLow",true);
  ITCRatioTriggerLargeWindowLow->SetDefaultValue(200);
  StoreITCRatioTriggerLargeWindowLow = 200;

  ITCRatioTriggerLargeWindowHigh = new G4UIcmdWithAnInteger("/DAQ/TriggerITCRatio/LargeWindowHigh", this);
  ITCRatioTriggerLargeWindowHigh->SetGuidance("Set the ITCRatio large trigger window high edge (in ns) i.e. for the denominator");
  ITCRatioTriggerLargeWindowHigh->SetParameterName("ITCRatioLargeWindowHigh",true);
  ITCRatioTriggerLargeWindowHigh->SetDefaultValue(1000);
  StoreITCRatioTriggerLargeWindowHigh = 1000;
}

WCSimWCDAQMessenger::~WCSimWCDAQMessenger()
{
  delete SaveFailuresTriggerDir;
  delete SaveFailuresTriggerMode;
  delete SaveFailuresTriggerTime;

  delete NHitsTriggerDir;
  delete NHitsTriggerThreshold;
  delete NHitsTriggerWindow;
  delete NHitsTriggerAdjustForNoise;

  delete ITCRatioTriggerThreshold;
  delete ITCRatioTriggerSmallWindow;
  delete ITCRatioTriggerLargeWindowLow;
  delete ITCRatioTriggerLargeWindowHigh;

  delete DigitizerChoice;
  delete TriggerChoice;
  delete WCSimDAQDir;
}

void WCSimWCDAQMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  //Because this Messenger class contains options for classes that don't exist when options are
  // read in (Trigger and Digitizer class options) we need to store each options' value
  // for use in the Tell*() methods later

  if (command == DigitizerChoice) {
    G4cout << "Digitizer choice set to " << newValue << G4endl;
    StoreDigitizerChoice = newValue;
  }
  else if (command == TriggerChoice) {
    G4cout << "Trigger choice set to " << newValue << G4endl;
    StoreTriggerChoice = newValue;
  }

  //Save failures "trigger"
  else if (command == SaveFailuresTriggerMode) {
    StoreSaveFailuresMode = SaveFailuresTriggerMode->GetNewIntValue(newValue);
    std::string failuremode;
    if(StoreSaveFailuresMode == 0)
      failuremode = "Saving only triggered events";
    else if(StoreSaveFailuresMode == 1)
      failuremode = "Saving both triggered and failed events";
    else if(StoreSaveFailuresMode == 2)
      failuremode = "Saving only failed events";
    else {
      G4cerr << "Unknown value of /DAQ/TriggerSaveFailures/Mode " << StoreSaveFailuresMode << " Exiting..." << G4endl;
      exit(-1);
    }      
    G4cout << failuremode << G4endl;
  }
  else if (command == SaveFailuresTriggerTime) {
    G4cout << "Trigger time for events which fail all triggers will be set to " << newValue << G4endl;
    StoreSaveFailuresTime = SaveFailuresTriggerTime->GetNewDoubleValue(newValue);
  }

  //NHits trigger
  else if (command == NHitsTriggerThreshold) {
    G4cout << "NHits trigger threshold set to " << newValue << G4endl;
    StoreNHitsThreshold = NHitsTriggerThreshold->GetNewIntValue(newValue);
  }
  else if (command == NHitsTriggerAdjustForNoise) {
    StoreNHitsAdjustForNoise = NHitsTriggerAdjustForNoise->GetNewBoolValue(newValue);
    if(StoreNHitsAdjustForNoise)
      G4cout << "Will adjust NHits trigger threshold using average dark noise rate" << G4endl;
  }
  else if (command == NHitsTriggerWindow) {
    G4cout << "NHits trigger window set to " << newValue << G4endl;
    StoreNHitsWindow = NHitsTriggerWindow->GetNewIntValue(newValue);
  }

  //ITC ratio trigger
  else if(command == ITCRatioTriggerThreshold) {
    G4cout << "ITC ratio threshold set to " << newValue << G4endl;
    StoreITCRatioTriggerThreshold = ITCRatioTriggerThreshold->GetNewDoubleValue(newValue);
  }
  else if(command == ITCRatioTriggerSmallWindow) {
    G4cout << "ITC ratio small trigger window set to " << newValue << G4endl;
    StoreITCRatioTriggerSmallWindow = ITCRatioTriggerSmallWindow->GetNewIntValue(newValue);
  }
  else if(command == ITCRatioTriggerLargeWindowLow) {
    G4cout << "ITC ratio large trigger window low edge set to " << newValue << G4endl;
    StoreITCRatioTriggerLargeWindowLow = ITCRatioTriggerLargeWindowLow->GetNewIntValue(newValue);
  }
  else if(command == ITCRatioTriggerLargeWindowHigh) {
    G4cout << "ITC ratio large trigger window high edge set to " << newValue << G4endl;
    StoreITCRatioTriggerLargeWindowHigh = ITCRatioTriggerLargeWindowHigh->GetNewIntValue(newValue);
  }
}

void WCSimWCDAQMessenger::SetEventActionOptions()
{
  G4cout << "Passing DAQ options to the event action class instance" << G4endl;
  WCSimEvent->SetDigitizerChoice(StoreDigitizerChoice);
  G4cout << "\tDigitizer choice set to " << StoreDigitizerChoice << G4endl;
  WCSimEvent->SetTriggerChoice(StoreTriggerChoice);
  G4cout << "\tTrigger choice set to " << StoreTriggerChoice << G4endl;
}

void WCSimWCDAQMessenger::SetTriggerOptions()
{
  G4cout << "Passing Trigger options to the trigger class instance" << G4endl;

  WCSimTrigger->SetSaveFailuresMode(StoreSaveFailuresMode);
  std::string failuremode;
  if(StoreSaveFailuresMode == 0)
    failuremode = "Saving only triggered events";
  else if(StoreSaveFailuresMode == 1)
    failuremode = "Saving both triggered and failed events";
  else if(StoreSaveFailuresMode == 2)
    failuremode = "Saving only failed events";
  G4cout << "\t" << failuremode << G4endl;
  WCSimTrigger->SetSaveFailuresTime(StoreSaveFailuresTime);
  G4cout << "\tTrigger time for events which fail all triggers will be set to " << StoreSaveFailuresTime << G4endl;

  WCSimTrigger->SetNHitsThreshold(StoreNHitsThreshold);
  G4cout << "\tNHits trigger threshold set to " << StoreNHitsThreshold << G4endl;
  WCSimTrigger->SetNHitsAdjustForNoise(StoreNHitsAdjustForNoise);
  if(StoreNHitsAdjustForNoise)
    G4cout << "\tWill adjust NHits trigger threshold using average dark noise rate" << G4endl;
  WCSimTrigger->SetNHitsWindow(StoreNHitsWindow);
  G4cout << "\tNHits trigger window set to " << StoreNHitsWindow << G4endl;

  WCSimTrigger->SetITCRatioThreshold(StoreITCRatioTriggerThreshold);
  G4cout << "\tITC ratio threshold set to " << StoreITCRatioTriggerThreshold << G4endl;
  WCSimTrigger->SetITCRatioSmallWindow(StoreITCRatioTriggerSmallWindow);
  G4cout << "\tITC ratio small window set to " << StoreITCRatioTriggerSmallWindow << G4endl;
  WCSimTrigger->SetITCRatioLargeWindowLow(StoreITCRatioTriggerLargeWindowLow);
  G4cout << "\tITC ratio large window low edge set to " << StoreITCRatioTriggerLargeWindowLow << G4endl;
  WCSimTrigger->SetITCRatioLargeWindowHigh(StoreITCRatioTriggerLargeWindowHigh);
  G4cout << "\tITC ratio large window high edge set to " << StoreITCRatioTriggerLargeWindowHigh << G4endl;
}

void WCSimWCDAQMessenger::SetDigitizerOptions()
{
  G4cout << "Passing Digitizer options to the digitizer class instance" << G4endl;
  WCSimDigitize->SKDigitizerType(StoreDigitizerChoice);
  G4cout << "\tDigitizer choice set to " << StoreDigitizerChoice << G4endl;
}
