#include "WCSimWCDAQMessenger.hh"
#include "WCSimEventAction.hh"
#include "WCSimWCDigitizer.hh"
#include "WCSimWCTrigger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"

#include <string>

WCSimWCDAQMessenger::WCSimWCDAQMessenger(WCSimEventAction* eventaction) :
  WCSimEvent(eventaction)
{
  initialiseString = " (this is a default set; it may be overwritten by user commands)";
  initialised = false;

  WCSimDAQDir = new G4UIdirectory("/DAQ/");
  WCSimDAQDir->SetGuidance("Commands to select DAQ options");

  G4String defaultDigitizer = "SKI";
  DigitizerChoice = new G4UIcmdWithAString("/DAQ/Digitizer", this);
  DigitizerChoice->SetGuidance("Set the Digitizer type");
  DigitizerChoice->SetGuidance("Available choices are:\n"
			       "SKI\n"
			       );
  DigitizerChoice->SetParameterName("Digitizer", false);
  DigitizerChoice->SetCandidates(
				 "SKI "
				 );
  DigitizerChoice->AvailableForStates(G4State_PreInit, G4State_Idle);
  DigitizerChoice->SetDefaultValue(defaultDigitizer);
  StoreDigitizerChoice = defaultDigitizer;
  SetNewValue(DigitizerChoice, defaultDigitizer);

  G4String defaultTrigger = "NDigits";
  TriggerChoice = new G4UIcmdWithAString("/DAQ/Trigger", this);
  TriggerChoice->SetGuidance("Set the Trigger type");
  TriggerChoice->SetGuidance("Available choices are:\n"
			     "NDigits\n"
			     "NDigits2\n"
			     "NHitsThenLocalNHits\n"
			     "NHitsThenRegions\n"
			     "NHitsThenITC\n"
			     "NHitsThenAnisotropy\n"
			     "NoTrigger \n"
			     );
  TriggerChoice->SetParameterName("Trigger", false);
  TriggerChoice->SetCandidates(
			       "NDigits "
			       "NDigits2 "
			       "NHitsThenLocalNHits "
			       "NHitsThenRegions "
			       "NHitsThenITC "
			       "NHitsThenAnisotropy "
			       "NoTrigger "
			       );
  TriggerChoice->AvailableForStates(G4State_PreInit, G4State_Idle);
  TriggerChoice->SetDefaultValue(defaultTrigger);
  StoreTriggerChoice = defaultTrigger;
  SetNewValue(TriggerChoice, defaultTrigger);

  bool defaultMultiDigitsPerTrigger = false;
  MultiDigitsPerTrigger = new G4UIcmdWithABool("/DAQ/MultiDigitsPerTrigger", this);
  MultiDigitsPerTrigger->SetGuidance("Allow the number of digits per PMT per trigger to be > 1?");
  MultiDigitsPerTrigger->SetParameterName("MultiDigitsPerTrigger",true);
  MultiDigitsPerTrigger->SetDefaultValue(defaultMultiDigitsPerTrigger);
  StoreMultiDigitsPerTrigger = defaultMultiDigitsPerTrigger;
  MultiDigitsPerTriggerSet = false; //this variable is bool & defaults are class specfic; use this to know if the default is overidden
  //don't SetNewValue -> defaults class-specific and taken from GetDefault*()

  bool defaultWriteGeomInfo = false;
  WriteGeomInfo = new G4UIcmdWithABool("/DAQ/WriteGeomInfo", this);
  WriteGeomInfo->SetGuidance("Allow trigger class to write out the geometry info to a file, then exit");
  WriteGeomInfo->SetParameterName("WriteGeomInfo",true);
  WriteGeomInfo->SetDefaultValue(defaultWriteGeomInfo);
  StoreWriteGeomInfo = defaultWriteGeomInfo;
  SetNewValue(WriteGeomInfo, G4UIcommand::ConvertToString(defaultWriteGeomInfo));

  //Generic digitizer specific options
  DigitizerDir = new G4UIdirectory("/DAQ/DigitizerOpt/");
  DigitizerDir->SetGuidance("Generic commands for digitizers");

  int defaultDigitizerDeadTime = -99;
  DigitizerDeadTime = new G4UIcmdWithAnInteger("/DAQ/DigitizerOpt/DeadTime", this);
  DigitizerDeadTime->SetGuidance("The deadtime for the digitizer (in ns)");
  DigitizerDeadTime->SetParameterName("DigitizerDeadTime",true);
  DigitizerDeadTime->SetDefaultValue(defaultDigitizerDeadTime);
  StoreDigitizerDeadTime = defaultDigitizerDeadTime;
  //don't SetNewValue -> defaults class-specific and taken from GetDefault*()

  int defaultDigitizerIntegrationWindow = -99;
  DigitizerIntegrationWindow = new G4UIcmdWithAnInteger("/DAQ/DigitizerOpt/IntegrationWindow", this);
  DigitizerIntegrationWindow->SetGuidance("The integration window for the digitizer (in ns)");
  DigitizerIntegrationWindow->SetParameterName("DigitizerIntegrationWindow",true);
  DigitizerIntegrationWindow->SetDefaultValue(defaultDigitizerIntegrationWindow);
  StoreDigitizerIntegrationWindow = defaultDigitizerIntegrationWindow;
  //don't SetNewValue -> defaults class-specific and taken from GetDefault*()


  //Save failure trigger specific options
  SaveFailuresTriggerDir = new G4UIdirectory("/DAQ/TriggerSaveFailures/");
  SaveFailuresTriggerDir->SetGuidance("Commands specific to the Save Failures trigger");

  int defaultSaveFailuresTriggerMode = 0;
  SaveFailuresTriggerMode = new G4UIcmdWithAnInteger("/DAQ/TriggerSaveFailures/Mode", this);
  SaveFailuresTriggerMode->SetGuidance("0: save only triggered events; 1: save both triggered and failed events; 2: save only failed events");
  SaveFailuresTriggerMode->SetParameterName("SaveFailuresMode",true);
  SaveFailuresTriggerMode->SetDefaultValue(defaultSaveFailuresTriggerMode);
  StoreSaveFailuresMode = defaultSaveFailuresTriggerMode;
  SetNewValue(SaveFailuresTriggerMode, G4UIcommand::ConvertToString(defaultSaveFailuresTriggerMode));

  double defaultSaveFailuresTriggerTime = 100;
  SaveFailuresTriggerTime = new G4UIcmdWithADouble("/DAQ/TriggerSaveFailures/TriggerTime", this);
  SaveFailuresTriggerTime->SetGuidance("The trigger time for the events which failed other triggers");
  SaveFailuresTriggerTime->SetParameterName("SaveFailuresTime",true);
  SaveFailuresTriggerTime->SetDefaultValue(defaultSaveFailuresTriggerTime);
  StoreSaveFailuresTime = defaultSaveFailuresTriggerTime;
  SetNewValue(SaveFailuresTriggerTime, G4UIcommand::ConvertToString(defaultSaveFailuresTriggerTime));

  int defaultSaveFailuresPreTriggerWindow = -400;
  SaveFailuresPreTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerSaveFailures/PreTriggerWindow", this);
  SaveFailuresPreTriggerWindow->SetGuidance("Set the SaveFailures pretrigger window (in ns)");
  SaveFailuresPreTriggerWindow->SetParameterName("SaveFailuresPreTriggerWindow",false);
  SaveFailuresPreTriggerWindow->SetDefaultValue(defaultSaveFailuresPreTriggerWindow);
  StoreSaveFailuresPreWindow = defaultSaveFailuresPreTriggerWindow;
  SetNewValue(SaveFailuresPreTriggerWindow, G4UIcommand::ConvertToString(defaultSaveFailuresPreTriggerWindow));

  int defaultSaveFailuresPostTriggerWindow = 950;
  SaveFailuresPostTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerSaveFailures/PostTriggerWindow", this);
  SaveFailuresPostTriggerWindow->SetGuidance("Set the SaveFailures posttrigger window (in ns)");
  SaveFailuresPostTriggerWindow->SetParameterName("SaveFailuresPostTriggerWindow",false);
  SaveFailuresPostTriggerWindow->SetDefaultValue(defaultSaveFailuresPostTriggerWindow);
  StoreSaveFailuresPostWindow = defaultSaveFailuresPostTriggerWindow;
  SetNewValue(SaveFailuresPostTriggerWindow, G4UIcommand::ConvertToString(defaultSaveFailuresPostTriggerWindow));


  //NDigits trigger specifc options
  NDigitsTriggerDir = new G4UIdirectory("/DAQ/TriggerNDigits/");
  NDigitsTriggerDir->SetGuidance("Commands specific to the NDigits trigger");  

  int defaultNDigitsTriggerThreshold = -99;
  NDigitsTriggerThreshold = new G4UIcmdWithAnInteger("/DAQ/TriggerNDigits/Threshold", this);
  NDigitsTriggerThreshold->SetGuidance("Set the NDigits trigger threshold");
  NDigitsTriggerThreshold->SetParameterName("NDigitsThreshold",false);
  NDigitsTriggerThreshold->SetDefaultValue(defaultNDigitsTriggerThreshold);
  StoreNDigitsThreshold = defaultNDigitsTriggerThreshold;
  //don't SetNewValue -> defaults class-specific and taken from GetDefault*()

  int defaultNDigitsTriggerWindow = -99;
  NDigitsTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerNDigits/Window", this);
  NDigitsTriggerWindow->SetGuidance("Set the NDigits trigger window (in ns)");
  NDigitsTriggerWindow->SetParameterName("NDigitsWindow",false);
  NDigitsTriggerWindow->SetDefaultValue(defaultNDigitsTriggerWindow);
  StoreNDigitsWindow = defaultNDigitsTriggerWindow;
  //don't SetNewValue -> defaults class-specific and taken from GetDefault*()

  bool defaultNDigitsTriggerAdjustForNoise = true;
  NDigitsTriggerAdjustForNoise = new G4UIcmdWithABool("/DAQ/TriggerNDigits/AdjustForNoise", this);
  NDigitsTriggerAdjustForNoise->SetGuidance("Adjust the NDigits trigger threshold automatically dependent on the average noise rate");
  NDigitsTriggerAdjustForNoise->SetParameterName("NDigitsAdjustForNoise",true);
  NDigitsTriggerAdjustForNoise->SetDefaultValue(defaultNDigitsTriggerAdjustForNoise);
  StoreNDigitsAdjustForNoise = defaultNDigitsTriggerAdjustForNoise;
  SetNewValue(NDigitsTriggerAdjustForNoise, G4UIcommand::ConvertToString(defaultNDigitsTriggerAdjustForNoise));

  int defaultNDigitsPreTriggerWindow = -99;
  NDigitsPreTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerNDigits/PreTriggerWindow", this);
  NDigitsPreTriggerWindow->SetGuidance("Set the NDigits pretrigger window (in ns)");
  NDigitsPreTriggerWindow->SetParameterName("NDigitsPreTriggerWindow",false);
  NDigitsPreTriggerWindow->SetDefaultValue(defaultNDigitsPreTriggerWindow);
  StoreNDigitsPreWindow = defaultNDigitsPreTriggerWindow;
  //don't SetNewValue -> defaults class-specific and taken from GetDefault*()

  int defaultNDigitsPostTriggerWindow = -99;
  NDigitsPostTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerNDigits/PostTriggerWindow", this);
  NDigitsPostTriggerWindow->SetGuidance("Set the NDigits posttrigger window (in ns)");
  NDigitsPostTriggerWindow->SetParameterName("NDigitsPostTriggerWindow",false);
  NDigitsPostTriggerWindow->SetDefaultValue(defaultNDigitsPostTriggerWindow);
  StoreNDigitsPostWindow = defaultNDigitsPostTriggerWindow;
  //don't SetNewValue -> defaults class-specific and taken from GetDefault*()


  //Local NHits trigger specifc options
  LocalNHitsTriggerDir = new G4UIdirectory("/DAQ/TriggerLocalNHits/");
  LocalNHitsTriggerDir->SetGuidance("Commands specific to the Local NHits trigger");

  int defaultLocalNHitsTriggerNeighbours = 50;
  LocalNHitsTriggerNeighbours = new G4UIcmdWithAnInteger("/DAQ/TriggerLocalNHits/Neighbours", this);
  LocalNHitsTriggerNeighbours->SetGuidance("Set the Local NHits trigger locality definition - the number of nearest neighbours");
  LocalNHitsTriggerNeighbours->SetParameterName("LocalNHitsNeighbours",false);
  LocalNHitsTriggerNeighbours->SetDefaultValue(defaultLocalNHitsTriggerNeighbours);
  SetNewValue(LocalNHitsTriggerNeighbours, G4UIcommand::ConvertToString(defaultLocalNHitsTriggerNeighbours));
  StoreLocalNHitsNeighbours = defaultLocalNHitsTriggerNeighbours;

  int defaultLocalNHitsTriggerThreshold = 10;
  LocalNHitsTriggerThreshold = new G4UIcmdWithAnInteger("/DAQ/TriggerLocalNHits/Threshold", this);
  LocalNHitsTriggerThreshold->SetGuidance("Set the Local NHits trigger threshold");
  LocalNHitsTriggerThreshold->SetParameterName("LocalNHitsThreshold",false);
  LocalNHitsTriggerThreshold->SetDefaultValue(defaultLocalNHitsTriggerThreshold);
  SetNewValue(LocalNHitsTriggerThreshold, G4UIcommand::ConvertToString(defaultLocalNHitsTriggerThreshold));
  StoreLocalNHitsThreshold = defaultLocalNHitsTriggerThreshold;

  int defaultLocalNHitsTriggerWindow = 50;
  LocalNHitsTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerLocalNHits/Window", this);
  LocalNHitsTriggerWindow->SetGuidance("Set the Local NHits trigger window (in ns)");
  LocalNHitsTriggerWindow->SetParameterName("LocalNHitsWindow",false);
  LocalNHitsTriggerWindow->SetDefaultValue(defaultLocalNHitsTriggerWindow);
  SetNewValue(LocalNHitsTriggerWindow, G4UIcommand::ConvertToString(defaultLocalNHitsTriggerWindow));
  StoreLocalNHitsWindow = defaultLocalNHitsTriggerWindow;

  bool defaultLocalNHitsTriggerAdjustForNoise = false;
  LocalNHitsTriggerAdjustForNoise = new G4UIcmdWithABool("/DAQ/TriggerLocalNHits/AdjustForNoise", this);
  LocalNHitsTriggerAdjustForNoise->SetGuidance("Adjust the Local NHits trigger threshold automatically dependent on the average noise rate");
  LocalNHitsTriggerAdjustForNoise->SetParameterName("LocalNHitsAdjustForNoise",true);
  LocalNHitsTriggerAdjustForNoise->SetDefaultValue(defaultLocalNHitsTriggerAdjustForNoise);
  SetNewValue(LocalNHitsTriggerAdjustForNoise, G4UIcommand::ConvertToString(defaultLocalNHitsTriggerAdjustForNoise));
  StoreLocalNHitsAdjustForNoise = defaultLocalNHitsTriggerAdjustForNoise;


  // NHits regions trigger specifc options
  TriggerRegionsDir = new G4UIdirectory("/DAQ/TriggerRegions/");
  TriggerRegionsDir->SetGuidance("Commands specific to the NHits regions triggers");

  int defaultTriggerRegionsNBinsPhi = 12;
  TriggerRegionsNBinsPhi = new G4UIcmdWithAnInteger("/DAQ/TriggerRegions/NBinsPhi", this);
  TriggerRegionsNBinsPhi->SetGuidance("Set the Trigger Regions side PMT definition - the number of phi bins");
  TriggerRegionsNBinsPhi->SetParameterName("TriggerRegionsNBinsPhi",false);
  TriggerRegionsNBinsPhi->SetDefaultValue(defaultTriggerRegionsNBinsPhi);
  SetNewValue(TriggerRegionsNBinsPhi, G4UIcommand::ConvertToString(defaultTriggerRegionsNBinsPhi));
  StoreLocalNHitsNBinsPhi = defaultTriggerRegionsNBinsPhi;

  int defaultTriggerRegionsNBinsZ = 11;
  TriggerRegionsNBinsZ = new G4UIcmdWithAnInteger("/DAQ/TriggerRegions/NBinsZ", this);
  TriggerRegionsNBinsZ->SetGuidance("Set the Trigger Regions side PMT definition - the number of Z bins");
  TriggerRegionsNBinsZ->SetParameterName("TriggerRegionsNBinsZ",false);
  TriggerRegionsNBinsZ->SetDefaultValue(defaultTriggerRegionsNBinsZ);
  SetNewValue(TriggerRegionsNBinsZ, G4UIcommand::ConvertToString(defaultTriggerRegionsNBinsZ));
  StoreLocalNHitsNBinsZ = defaultTriggerRegionsNBinsZ;

  int defaultTriggerRegionsNRings = 2;
  TriggerRegionsNRings = new G4UIcmdWithAnInteger("/DAQ/TriggerRegions/NRings", this);
  TriggerRegionsNRings->SetGuidance("Set the Trigger Regions top/bottom PMT definition - the number of rings (note the central circle is not a ring, therefore NRings=0 is valid)");
  TriggerRegionsNRings->SetParameterName("TriggerRegionsNRings",false);
  TriggerRegionsNRings->SetDefaultValue(defaultTriggerRegionsNRings);
  SetNewValue(TriggerRegionsNRings, G4UIcommand::ConvertToString(defaultTriggerRegionsNRings));
  StoreLocalNHitsNRings = defaultTriggerRegionsNRings;

  int defaultTriggerRegionsNCentralSectors = 3;
  TriggerRegionsNCentralSectors = new G4UIcmdWithAnInteger("/DAQ/TriggerRegions/NCentralSectors", this);
  TriggerRegionsNCentralSectors->SetGuidance("Set the Trigger Regions top/bottom PMT definition - the number of sectors in the the central circle");
  TriggerRegionsNCentralSectors->SetParameterName("TriggerRegionsNCentralSectors",false);
  TriggerRegionsNCentralSectors->SetDefaultValue(defaultTriggerRegionsNCentralSectors);
  SetNewValue(TriggerRegionsNCentralSectors, G4UIcommand::ConvertToString(defaultTriggerRegionsNCentralSectors));
  StoreLocalNHitsNCentralSectors = defaultTriggerRegionsNCentralSectors;

  int defaultTriggerRegionsNRingSectors = 3;
  TriggerRegionsNRingSectors = new G4UIcmdWithAnInteger("/DAQ/TriggerRegions/NRingSectors", this);
  TriggerRegionsNRingSectors->SetGuidance("Set the Trigger Regions top/bottom PMT definition - the multiplicative number of sectors more in each subsequent ring");
  TriggerRegionsNRingSectors->SetParameterName("TriggerRegionsNRingSectors",false);
  TriggerRegionsNRingSectors->SetDefaultValue(defaultTriggerRegionsNRingSectors);
  SetNewValue(TriggerRegionsNRingSectors, G4UIcommand::ConvertToString(defaultTriggerRegionsNRingSectors));
  StoreLocalNHitsNRingSectors = defaultTriggerRegionsNRingSectors;


  //ITC Ratio trigger specifc options
  ITCRatioTriggerDir = new G4UIdirectory("/DAQ/TriggerITCRatio/");
  ITCRatioTriggerDir->SetGuidance("Commands specific to the ITCRatio trigger");

  double defaultITCRatioTriggerThreshold = 0.3;
  ITCRatioTriggerThreshold = new G4UIcmdWithADouble("/DAQ/TriggerITCRatio/Threshold", this);
  ITCRatioTriggerThreshold->SetGuidance("Set the ITCRatio trigger threshold");
  ITCRatioTriggerThreshold->SetParameterName("ITCRatioThreshold",true);
  ITCRatioTriggerThreshold->SetDefaultValue(defaultITCRatioTriggerThreshold);
  StoreITCRatioTriggerThreshold = defaultITCRatioTriggerThreshold;
  SetNewValue(ITCRatioTriggerThreshold, G4UIcommand::ConvertToString(defaultITCRatioTriggerThreshold));

  int defaultITCRatioTriggerSmallWindow = 200;
  ITCRatioTriggerSmallWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerITCRatio/SmallWindow", this);
  ITCRatioTriggerSmallWindow->SetGuidance("Set the ITCRatio small trigger window (in ns) i.e. for the numerator");
  ITCRatioTriggerSmallWindow->SetParameterName("ITCRatioSmallWindow",true);
  ITCRatioTriggerSmallWindow->SetDefaultValue(defaultITCRatioTriggerSmallWindow);
  StoreITCRatioTriggerSmallWindow = defaultITCRatioTriggerSmallWindow;
  SetNewValue(ITCRatioTriggerSmallWindow, G4UIcommand::ConvertToString(defaultITCRatioTriggerSmallWindow));

  int defaultITCRatioTriggerLargeWindowLow = 200;
  ITCRatioTriggerLargeWindowLow = new G4UIcmdWithAnInteger("/DAQ/TriggerITCRatio/LargeWindowLow", this);
  ITCRatioTriggerLargeWindowLow->SetGuidance("Set the ITCRatio large trigger window low edge (in ns) i.e. for the denominator");
  ITCRatioTriggerLargeWindowLow->SetParameterName("ITCRatioLargeWindowLow",true);
  ITCRatioTriggerLargeWindowLow->SetDefaultValue(defaultITCRatioTriggerLargeWindowLow);
  StoreITCRatioTriggerLargeWindowLow = defaultITCRatioTriggerLargeWindowLow;
  SetNewValue(ITCRatioTriggerLargeWindowLow, G4UIcommand::ConvertToString(defaultITCRatioTriggerLargeWindowLow));

  int defaultITCRatioTriggerLargeWindowHigh = 1000;
  ITCRatioTriggerLargeWindowHigh = new G4UIcmdWithAnInteger("/DAQ/TriggerITCRatio/LargeWindowHigh", this);
  ITCRatioTriggerLargeWindowHigh->SetGuidance("Set the ITCRatio large trigger window high edge (in ns) i.e. for the denominator");
  ITCRatioTriggerLargeWindowHigh->SetParameterName("ITCRatioLargeWindowHigh",true);
  ITCRatioTriggerLargeWindowHigh->SetDefaultValue(defaultITCRatioTriggerLargeWindowHigh);
  StoreITCRatioTriggerLargeWindowHigh = defaultITCRatioTriggerLargeWindowHigh;
  SetNewValue(ITCRatioTriggerLargeWindowHigh, G4UIcommand::ConvertToString(defaultITCRatioTriggerLargeWindowHigh));


  initialiseString = "";
  initialised = true;
}

WCSimWCDAQMessenger::~WCSimWCDAQMessenger()
{
  delete SaveFailuresTriggerDir;
  delete SaveFailuresTriggerMode;
  delete SaveFailuresTriggerTime;
  delete SaveFailuresPreTriggerWindow;
  delete SaveFailuresPostTriggerWindow;

  delete NDigitsTriggerDir;
  delete NDigitsTriggerThreshold;
  delete NDigitsTriggerWindow;
  delete NDigitsTriggerAdjustForNoise;
  delete NDigitsPreTriggerWindow;
  delete NDigitsPostTriggerWindow;

  delete DigitizerDir;
  delete DigitizerDeadTime;
  delete DigitizerIntegrationWindow;

  delete LocalNHitsTriggerDir;
  delete LocalNHitsTriggerNeighbours;
  delete LocalNHitsTriggerThreshold;
  delete LocalNHitsTriggerWindow;
  delete LocalNHitsTriggerAdjustForNoise;

  delete TriggerRegionsNBinsPhi;
  delete TriggerRegionsNBinsZ;
  delete TriggerRegionsNRings;
  delete TriggerRegionsNCentralSectors;
  delete TriggerRegionsNRingSectors;

  delete ITCRatioTriggerThreshold;
  delete ITCRatioTriggerSmallWindow;
  delete ITCRatioTriggerLargeWindowLow;
  delete ITCRatioTriggerLargeWindowHigh;

  delete DigitizerChoice;
  delete TriggerChoice;
  delete MultiDigitsPerTrigger;
  delete WriteGeomInfo;
  delete WCSimDAQDir;
}

void WCSimWCDAQMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  //Because this Messenger class contains options for classes that don't exist when options are
  // read in (Trigger and Digitizer class options) we need to store each options' value
  // for use in the Tell*() methods later

  if (command == DigitizerChoice) {
    G4cout << "Digitizer choice set to " << newValue << initialiseString.c_str() << G4endl;
    WCSimEvent->SetDigitizerChoice(newValue);
    StoreDigitizerChoice = newValue;
  }
  else if (command == TriggerChoice) {
    G4cout << "Trigger choice set to " << newValue << initialiseString.c_str() << G4endl;
    WCSimEvent->SetTriggerChoice(newValue);
    StoreTriggerChoice = newValue;
  }

  else if (command == MultiDigitsPerTrigger) {
    StoreMultiDigitsPerTrigger = MultiDigitsPerTrigger->GetNewBoolValue(newValue);
    if(!StoreMultiDigitsPerTrigger)
      G4cout << "Will restrict number of digits per PMT per trigger to <= 1" << initialiseString.c_str() << G4endl;
    else
      G4cout << "Will allow number of digits per PMT per trigger to go > 1" << initialiseString.c_str() << G4endl;
    if(initialised)
      MultiDigitsPerTriggerSet = true;
  }
  else if (command == WriteGeomInfo) {
    StoreWriteGeomInfo = WriteGeomInfo->GetNewBoolValue(newValue);
    if(StoreWriteGeomInfo)
      G4cout << "Will write out geometry information for triggering to a file, then exit" << initialiseString.c_str() << G4endl;
  }

  //Generic digitizer options
  else if (command == DigitizerDeadTime) {
    G4cout << "Digitizer deadtime set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoreDigitizerDeadTime = DigitizerDeadTime->GetNewIntValue(newValue);
  }
  else if (command == DigitizerIntegrationWindow) {
    G4cout << "Digitizer integration window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoreDigitizerIntegrationWindow = DigitizerIntegrationWindow->GetNewIntValue(newValue);
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
    G4cout << failuremode << initialiseString.c_str() << G4endl;
  }
  else if (command == SaveFailuresTriggerTime) {
    G4cout << "Trigger time for events which fail all triggers will be set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoreSaveFailuresTime = SaveFailuresTriggerTime->GetNewDoubleValue(newValue);
  }
  else if (command == SaveFailuresPreTriggerWindow) {
    G4cout << "SaveFailures pretrigger window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoreSaveFailuresPreWindow = SaveFailuresPreTriggerWindow->GetNewIntValue(newValue);
  }
  else if (command == SaveFailuresPostTriggerWindow) {
    G4cout << "SaveFailures posttrigger window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoreSaveFailuresPostWindow = SaveFailuresPostTriggerWindow->GetNewIntValue(newValue);
  }

  //NDigits trigger
  else if (command == NDigitsTriggerThreshold) {
    G4cout << "NDigits trigger threshold set to " << newValue << initialiseString.c_str() << G4endl;
    StoreNDigitsThreshold = NDigitsTriggerThreshold->GetNewIntValue(newValue);
  }
  else if (command == NDigitsTriggerAdjustForNoise) {
    StoreNDigitsAdjustForNoise = NDigitsTriggerAdjustForNoise->GetNewBoolValue(newValue);
    if(StoreNDigitsAdjustForNoise)
      G4cout << "Will adjust NDigits trigger threshold using average dark noise rate" << initialiseString.c_str() << G4endl;
  }
  else if (command == NDigitsTriggerWindow) {
    G4cout << "NDigits trigger window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoreNDigitsWindow = NDigitsTriggerWindow->GetNewIntValue(newValue);
  }
  else if (command == NDigitsPreTriggerWindow) {
    G4cout << "NDigits pretrigger window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoreNDigitsPreWindow = NDigitsPreTriggerWindow->GetNewIntValue(newValue);
  }
  else if (command == NDigitsPostTriggerWindow) {
    G4cout << "NDigits posttrigger window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoreNDigitsPostWindow = NDigitsPostTriggerWindow->GetNewIntValue(newValue);
  }

  //Local NHits trigger
  else if (command == LocalNHitsTriggerNeighbours) {
    G4cout << "Local NHits trigger neighbours set to " << newValue << initialiseString.c_str() << G4endl;
    StoreLocalNHitsNeighbours = LocalNHitsTriggerNeighbours->GetNewIntValue(newValue);
  }
  else if (command == LocalNHitsTriggerThreshold) {
    G4cout << "Local NHits trigger threshold set to " << newValue << initialiseString.c_str() << G4endl;
    StoreLocalNHitsThreshold = LocalNHitsTriggerThreshold->GetNewIntValue(newValue);
  }
  else if (command == LocalNHitsTriggerAdjustForNoise) {
    StoreLocalNHitsAdjustForNoise = LocalNHitsTriggerAdjustForNoise->GetNewBoolValue(newValue);
    if(StoreLocalNHitsAdjustForNoise)
      G4cout << "Will adjust Local NHits trigger threshold using average dark noise rate" << initialiseString.c_str() << G4endl;
  }
  else if (command == LocalNHitsTriggerWindow) {
    G4cout << "Local NHits trigger window set to " << newValue << initialiseString.c_str() << G4endl;
    StoreLocalNHitsWindow = LocalNHitsTriggerWindow->GetNewIntValue(newValue);
  }

  //Trigger regions
  else if (command == TriggerRegionsNBinsPhi) {
    G4cout << "Number of phi bins in side PMT region set to " << newValue << initialiseString.c_str() << G4endl;
    StoreLocalNHitsNBinsPhi = TriggerRegionsNBinsPhi->GetNewIntValue(newValue);
  }
  else if (command == TriggerRegionsNBinsZ) {
    G4cout << "Number of Z bins in side PMT region set to " << newValue << initialiseString.c_str() << G4endl;
    StoreLocalNHitsNBinsZ = TriggerRegionsNBinsZ->GetNewIntValue(newValue);
  }
  else if (command == TriggerRegionsNRings) {
    G4cout << "Number of rings in the top/bottom region set to " << newValue << initialiseString.c_str() << G4endl;
    StoreLocalNHitsNRings = TriggerRegionsNRings->GetNewIntValue(newValue);
  }
  else if (command == TriggerRegionsNCentralSectors) {
    G4cout << "Number of sectors in the central circle sector set to " << newValue << initialiseString.c_str() << G4endl;
    StoreLocalNHitsNCentralSectors = TriggerRegionsNCentralSectors->GetNewIntValue(newValue);
  }
  else if (command == TriggerRegionsNRingSectors) {
    G4cout << "Multiplicative number more sectors per ring set to " << newValue << initialiseString.c_str() << G4endl;
    StoreLocalNHitsNRingSectors = TriggerRegionsNRingSectors->GetNewIntValue(newValue);
  }

  //ITC ratio trigger
  else if(command == ITCRatioTriggerThreshold) {
    G4cout << "ITC ratio threshold set to " << newValue << initialiseString.c_str() << G4endl;
    StoreITCRatioTriggerThreshold = ITCRatioTriggerThreshold->GetNewDoubleValue(newValue);
  }
  else if(command == ITCRatioTriggerSmallWindow) {
    G4cout << "ITC ratio small trigger window set to " << newValue << initialiseString.c_str() << G4endl;
    StoreITCRatioTriggerSmallWindow = ITCRatioTriggerSmallWindow->GetNewIntValue(newValue);
  }
  else if(command == ITCRatioTriggerLargeWindowLow) {
    G4cout << "ITC ratio large trigger window low edge set to " << newValue << initialiseString.c_str() << G4endl;
    StoreITCRatioTriggerLargeWindowLow = ITCRatioTriggerLargeWindowLow->GetNewIntValue(newValue);
  }
  else if(command == ITCRatioTriggerLargeWindowHigh) {
    G4cout << "ITC ratio large trigger window high edge set to " << newValue << initialiseString.c_str() << G4endl;
    StoreITCRatioTriggerLargeWindowHigh = ITCRatioTriggerLargeWindowHigh->GetNewIntValue(newValue);
  }

}

void WCSimWCDAQMessenger::SetTriggerOptions()
{
  G4cout << "Passing Trigger options to the trigger class instance" << G4endl;

  if(MultiDigitsPerTriggerSet) {
    WCSimTrigger->SetMultiDigitsPerTrigger(StoreMultiDigitsPerTrigger);
    if(!StoreMultiDigitsPerTrigger)
      G4cout << "\tWill restrict number of digits per PMT per trigger to <= 1" << G4endl;
    else
      G4cout << "\tWill allow number of digits per PMT per trigger to go > 1" << G4endl;
  }

  WCSimTrigger->SetWriteGeomInfo(StoreWriteGeomInfo);
  if(StoreWriteGeomInfo)
    G4cout << "\tWill write out geometry information for triggering to a file, then exit" << G4endl;

  //SaveFailures
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
  G4cout << "\tTrigger time for events which fail all triggers will be set to " << StoreSaveFailuresTime << " ns" << G4endl;
  if(StoreSaveFailuresPreWindow >= -1E6) {
    WCSimTrigger->SetSaveFailuresPreTriggerWindow(StoreSaveFailuresPreWindow);
    G4cout << "\tSaveFailures pretrigger window set to " << StoreSaveFailuresPreWindow << " ns" << G4endl;
  }
  if(StoreSaveFailuresPostWindow >= 0) {
    WCSimTrigger->SetSaveFailuresPostTriggerWindow(StoreSaveFailuresPostWindow);
    G4cout << "\tSaveFailures posttrigger window set to " << StoreSaveFailuresPostWindow << " ns" << G4endl;
  }

  //NDigits
  if(StoreNDigitsThreshold >= 0) {
    WCSimTrigger->SetNDigitsThreshold(StoreNDigitsThreshold);
    G4cout << "\tNDigits trigger threshold set to " << StoreNDigitsThreshold << G4endl;
  }
  WCSimTrigger->SetNDigitsAdjustForNoise(StoreNDigitsAdjustForNoise);
  if(StoreNDigitsAdjustForNoise)
    G4cout << "\tWill adjust NDigits trigger threshold using average dark noise rate" << G4endl;
  if(StoreNDigitsWindow >= 0) {
    WCSimTrigger->SetNDigitsWindow(StoreNDigitsWindow);
    G4cout << "\tNDigits trigger window set to " << StoreNDigitsWindow << " ns" << G4endl;
  }
  if(StoreNDigitsPreWindow >= 0) {
    WCSimTrigger->SetNDigitsPreTriggerWindow(StoreNDigitsPreWindow);
    G4cout << "\tNDigits pretrigger window set to " << StoreNDigitsPreWindow << " ns" << G4endl;
  }
  if(StoreNDigitsPostWindow >= 0) {
    WCSimTrigger->SetNDigitsPostTriggerWindow(StoreNDigitsPostWindow);
    G4cout << "\tNDigits posttrigger window set to " << StoreNDigitsPostWindow << " ns" << G4endl;
  }

  //Local NHits
  WCSimTrigger->SetLocalNHitsNeighbours(StoreLocalNHitsNeighbours);
  G4cout << "\tLocalNHits trigger neighbours set to " << StoreLocalNHitsNeighbours << G4endl;
  WCSimTrigger->SetLocalNHitsThreshold(StoreLocalNHitsThreshold);
  G4cout << "\tLocalNHits trigger threshold set to " << StoreLocalNHitsThreshold << G4endl;
  WCSimTrigger->SetLocalNHitsAdjustForNoise(StoreLocalNHitsAdjustForNoise);
  if(StoreLocalNHitsAdjustForNoise)
    G4cout << "\tWill adjust Local NHits trigger threshold using average dark noise rate" << G4endl;
  WCSimTrigger->SetLocalNHitsWindow(StoreLocalNHitsWindow);
  G4cout << "\tLocalNHits trigger window set to " << StoreLocalNHitsWindow << G4endl;

  //Regions
  WCSimTrigger->SetRegionsNBinsP(StoreLocalNHitsNBinsPhi);
  G4cout << "\tNumber of phi bins in side PMT region set to " << StoreLocalNHitsNBinsPhi << G4endl;
  WCSimTrigger->SetRegionsNBinsZ(StoreLocalNHitsNBinsZ);
  G4cout << "\tNumber of Z bins in side PMT region set to " << StoreLocalNHitsNBinsZ << G4endl;
  WCSimTrigger->SetRegionsNRings(StoreLocalNHitsNRings);
  G4cout << "\tNumber of rings in the top/bottom region set to " << StoreLocalNHitsNRings << G4endl;
  WCSimTrigger->SetRegionsNCentralSectors(StoreLocalNHitsNCentralSectors);
  G4cout << "\tNumber of sectors in the central circle sector set to " << StoreLocalNHitsNCentralSectors << G4endl;
  WCSimTrigger->SetRegionsNRingSectors(StoreLocalNHitsNRingSectors);
  G4cout << "\tMultiplicative number more sectors per ring set to " << StoreLocalNHitsNRingSectors << G4endl;

  //ITC ratio
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
  if(StoreDigitizerDeadTime >= 0) {
    WCSimDigitize->SetDigitizerDeadTime(StoreDigitizerDeadTime);
    G4cout << "\tDigitizer deadtime set to " << StoreDigitizerDeadTime << " ns"  << G4endl;
  }
  if(StoreDigitizerIntegrationWindow >= 0) {
    WCSimDigitize->SetDigitizerIntegrationWindow(StoreDigitizerIntegrationWindow);
    G4cout << "\tDigitizer integration window set to " << StoreDigitizerIntegrationWindow << " ns" << G4endl;
  }
}
