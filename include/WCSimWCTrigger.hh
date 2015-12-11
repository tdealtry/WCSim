#ifndef WCSimWCTrigger_h
#define WCSimWCTrigger_h 1

#include "WCSimEnumerations.hh"
#include "WCSimWCDAQMessenger.hh"
#include "WCSimDetectorConstruction.hh"
#include "G4VDigitizerModule.hh"
#include "WCSimWCDigi.hh"
#include "WCSimWCHit.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <map>
#include <vector>
#include "TVector3.h"

class WCSimWCDigiTrigger;
typedef G4TDigiCollection<WCSimWCDigiTrigger> WCSimWCTriggeredDigitsCollection;

// *******************************************
// BASE CLASS
// *******************************************

/**
 * \class WCSimWCTriggerBase
 *
 * \brief The base class for WCSim triggering algorithms
 *
 * Concrete implementations of a trigger class should inherit from this class.
 * Minimally, only DoTheWork() needs to be implemented in the implementation class.
 *
 */

class WCSimWCTriggerBase : public G4VDigitizerModule
{

public:

  ///Print out the number of occurances of each unique element in a vector
  template<typename T> void PrintVectorCount(std::vector<T> & v);
  void PrintVector3Count(std::vector<std::pair<int, TVector3> > & v, bool r, bool p, bool z);

  ///Create WCSimWCTriggerBase instance with knowledge of the detector and DAQ options
  WCSimWCTriggerBase(G4String name, WCSimDetectorConstruction*, WCSimWCDAQMessenger*);

  virtual ~WCSimWCTriggerBase();

  /**
   * \brief The main user-callable routine of the class. Gets the input & creates the output WCSimWCTriggeredDigitsCollection's, then calls DoTheWork()
   *
   * The virtual keyword for this method is DEPRECATED
   * It is only defined virtual now because it is overridden in the old class (WCSimWCDigitizer)
   */
  virtual void Digitize();

  ///Returns the number of trigger gates in the event (i.e. the number of triggers passed)
  int NumberOfGatesInThisEvent() { return TriggerTimes.size(); }
  ///Get the time of the ith trigger
  Float_t              GetTriggerTime(int i) { return TriggerTimes[i];}
  ///Get the trigger type of the ith trigger
  TriggerType_t        GetTriggerType(int i) { return TriggerTypes[i];}
  ///Get the additional trigger information associated with the ith trigger
  std::vector<Float_t> GetTriggerInfo(int i) { return TriggerInfos[i];}

  //
  // Trigger algorithm option set methods
  //

  ///Set whether to allow the number of digits per PMT per trigger to go > 1
  void SetMultiDigitsPerTrigger(G4bool allow_multi) { multiDigitsPerTrigger = allow_multi; }
  ///Set whether to write the geometry information into a file, then exit
  void SetWriteGeomInfo(G4bool write_geom) { writeGeom = write_geom; }

  // NDigits options
  ///Set the threshold for the NDigits trigger
  void SetNDigitsThreshold(G4int threshold) { ndigitsThreshold = threshold; }
  ///Set the time window for the NDigits trigger
  void SetNDigitsWindow(G4int window) { ndigitsWindow = window; }
  ///Automatically adjust the NDigits threshold based on the average noise occupancy?
  void SetNDigitsAdjustForNoise    (G4bool adjust)      { ndigitsAdjustForNoise = adjust; }
  ///Set the pretrigger window for the NDigits trigger (value will be forced negative)
  void SetNDigitsPreTriggerWindow(G4int window)  { ndigitsPreTriggerWindow  = - abs(window); }
  ///Set the posttrigger window for the NDigits trigger (value will be forced positive)
  void SetNDigitsPostTriggerWindow(G4int window) { ndigitsPostTriggerWindow = + abs(window); }

  // Local NHits options
  ///Set the number of nearest neighbours to use for 'local' in the Local NHits trigger
  void SetLocalNHitsNeighbours(G4int neighbours) { localNHitsNeighbours = neighbours; }
  ///Set the threshold for the Local NHits trigger
  void SetLocalNHitsThreshold(G4int threshold) { localNHitsThreshold = threshold; }
  ///Set the time window for the Local NHits trigger
  void SetLocalNHitsWindow(G4int window) { localNHitsWindow = window; }
  ///Automatically adjust the NHits threshold based on the average noise occupancy?
  void SetLocalNHitsAdjustForNoise(G4bool adjust) { localNHitsAdjustForNoise = adjust; }

  // Regions options
  /// Set the number of phi bins to use when splitting the side PMTs into regions
  void SetRegionsNBinsP(G4int nbinsp) { regionsNBinsP = nbinsp; }
  /// Set the number of Z bins to use when splitting the side PMTs into regions
  void SetRegionsNBinsZ(G4int nbinsz) { regionsNBinsZ = nbinsz; }
  /// Set the number of rings to use when splitting the top/bottom PMTs into regions. Note the central circle is not a ring (i.e. regionsNRings=0 is valid)
  void SetRegionsNRings(G4int NRings) { regionsNRings = NRings; }
  /// Set the number of sectors to split the central circle into when splitting the top/bottom PMTs into regions
  void SetRegionsNCentralSectors(G4int NCentralSectors) { regionsNCentralSectors = NCentralSectors; }
  /// Set the number of additional segments (relative to each inner segment) to split the ring into when splitting the top/bottom PMTs into regions
  void SetRegionsNRingSectors(G4int NRingSectors) { regionsNRingSectors = NRingSectors; }

  // Save trigger failures options
  ///Set the mode for saving failed triggers (0:save only triggered events, 1:save both triggered events & failed events, 2:save only failed events)
  void SetSaveFailuresMode       (G4int mode )        { saveFailuresMode = mode; }
  ///Set the dummy trigger time for the failed triggers
  void SetSaveFailuresTime       (G4double time )     { saveFailuresTime = time; }
  ///Set the pretrigger window for the SaveFailures trigger (value will be forced negative)
  void SetSaveFailuresPreTriggerWindow(G4int window)  { saveFailuresPreTriggerWindow  = - abs(window); }
  ///Set the posttrigger window for the SaveFailures trigger (value will be forced positive)
  void SetSaveFailuresPostTriggerWindow(G4int window) { saveFailuresPostTriggerWindow = + abs(window); }

  ///Knowledge of the dark rate (use for automatically adjusting for noise)
  void SetDarkRate(double idarkrate){ PMTDarkRate = idarkrate; }

  ///DEPRECATED function used in old class (WCSimWCDigitizer), and called in WCSimEventAction
  virtual void SetPMTSize(G4float /*inputSize*/) {};



protected:

  ///This should call the trigger algorithms, and handle any temporary DigitsCollection's
  virtual void DoTheWork(WCSimWCDigitsCollection* WCDCPMT) = 0;

  /// Get the default threshold, etc. from the derived class, and override with read from the .mac file
  void GetVariables();

  ///Set the default trigger class specific decision of whether to save multiple digits per PMT per trigger (overridden by .mac)
  virtual bool GetDefaultMultiDigitsPerTrigger()   { return true; }
  ///Set the default trigger class specific NDigits window (in ns) (overridden by .mac)
  virtual int GetDefaultNDigitsWindow()            { return 200; }
  ///Set the default trigger class specific NDigits threshold (in ns) (overridden by .mac)
  virtual int GetDefaultNDigitsThreshold()         { return 25; }
  ///Set the default trigger class specific NDigits pretrigger window (in ns) (overridden by .mac)
  virtual int GetDefaultNDigitsPreTriggerWindow()  { return -400; }
  ///Set the default trigger class specific NDigits posttrigger window (in ns) (overridden by .mac)
  virtual int GetDefaultNDigitsPostTriggerWindow() { return 950; }

  ///Get the pretrigger window for a given trigger algorithm
  int GetPreTriggerWindow(TriggerType_t t);
  ///Get the posttrigger window for a given trigger algorithm
  int GetPostTriggerWindow(TriggerType_t t);

  //these are the algorithms that perform triggering
  //they are stored here so that different trigger classes can use the same algorithms without copying code
  /**
   * \brief An NDigits trigger algorithm
   *
   * Looks through the input WCSimWCDigitsCollection and integrates the number of hits in a (specified) time window
   * If the integral passes above a (specified) threshold, a trigger is issued
   *
   * The trigger type is kTriggerNDigits
   *
   * The trigger time is the time of the first digit above threshold
   *
   * The trigger information is the number of hits in the time window (i.e. the number of hits that caused the trigger to fire)
   *
   * Currently setup with the optional 'test' argument which runs the algorithm with half the hit threshold
   * for testing purposes. Triggers issued in this mode have type kTriggerNDigitsTest
   */
  void AlgNDigits(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits, bool test=false);

  /**
   * \brief An NHits then local NHits (using nearest neighbours) trigger algorithm
   *
   * Looks through the input WCSimWCDigitsCollection and integrates the number of hits in a (specified) time window
   * If the integral passes above a (specified) threshold, a trigger is issued
   * Else looks at cuts that are tighter in both (specified) time and (specified) space
   * If the integral passes above a (specified) threshold, a trigger is issued
   *
   * The trigger type is kTriggerNHits or kTriggerLocalNHits
   *
   * The trigger time is the time of the first digit above threshold
   *
   * The trigger information is the number of hits in the time window (i.e. the number of hits that caused the trigger to fire) (NHits & local NHits)
   * and the seed PMT ID (local NHits)
   */
  void AlgNHitsThenLocalNHits(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits);

  /**
   * \brief An NHits then local NHits (using regions) trigger algorithm
   *
   * Looks through the input WCSimWCDigitsCollection and integrates the number of hits in a (specified) time window
   * If the integral passes above a (specified) threshold, a trigger is issued
   * Else ???
   *
   * The trigger type is kTriggerNHits or kTriggerLocalNHits
   *
   * The trigger time is the time of the first digit above threshold
   *
   * The trigger information is the number of hits in the time window (i.e. the number of hits that caused the trigger to fire) (NHits & local NHits)
   * and ???
   */
  void AlgNHitsThenRegions(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits);

  // Geometry methods
  ///Write out the geometry information required for this trigger
  virtual void WriteGeomInfo() = 0;
  ///Read in the geometry information required for this trigger
  virtual void ReadGeomInfo() = 0;
  ///Split the tank into areas & get the PMTs IDs for each area
  void PopulatePMTAreas();
  ///Get the region(s) the PMT belongs to
  std::vector<int> FindRegion(const int tubeid);
  ///Find the nearest neighbours of a PMT
  std::vector<int> FindPMTNearestNeighbours(int i);
  ///Find the nearest neighbours of all PMTs
  void FindAllPMTNearestNeighbours();
  std::vector<WCSimPmtInfo*> *    myPMTs;        ///< Vector with the position/orientation of every PMT in the geometry. myPMTs[i] is the PMT with tubeID=i+1
  std::vector< std::vector<int> > pmtNeighbours; ///< Vector of vectors that give the nearest neighbour for each PMT. pmtNeighbours[i] contains the tubeIDs of the neighbours of the PMT with tubeID=i+1
  std::vector<std::pair<std::vector<int>, TVector3> > pmtBlocks; ///< Vector of pairs with first entry: vector of PMT ids in the block, second entry: TVector3 of the average position

  WCSimWCTriggeredDigitsCollection*   DigitsCollection; ///< The main output of the class - collection of digits in the trigger window
  std::map<int,int>          DigiHitMap; ///< Keeps track of the PMTs that have been added to the output WCSimWCTriggeredDigitsCollection

  std::vector<Float_t>                TriggerTimes; ///< The times of the triggers
  std::vector<TriggerType_t>          TriggerTypes; ///< The type of the triggers
  std::vector< std::vector<Float_t> > TriggerInfos; ///< Additional information associated with each trigger

  WCSimWCDAQMessenger*       DAQMessenger; ///< Get the options from the .mac file
  WCSimDetectorConstruction* myDetector;   ///< Know about the detector, so can add appropriate PMT time smearing

  /// Clear the Trigger* vectors and DigiHitMap
  void ReInitialize() {
    TriggerTimes.clear();
    TriggerTypes.clear();
    TriggerInfos.clear();
    DigiHitMap.clear();
  }

  double PMTDarkRate;    ///< Dark noise rate of the PMTs

  // Trigger algorithm options
  G4bool writeGeom;                ///< Write the geometry information for triggering, then exit?
  G4bool multiDigitsPerTrigger;    ///< Allow the number of digits per PMT saved in each trigger window to go > 1?
  //NDigits
  G4int  ndigitsThreshold;         ///< The threshold for the NDigits trigger
  G4int  ndigitsWindow;            ///< The time window for the NDigits trigger
  G4bool ndigitsAdjustForNoise;    ///< Automatically adjust the NDigits trigger threshold based on the average dark noise rate?
  G4int  ndigitsPreTriggerWindow;  ///< The pretrigger window to save before an NDigits trigger
  G4int  ndigitsPostTriggerWindow; ///< The posttrigger window to save after an NDigits trigger
  //local NHits
  G4int localNHitsNeighbours;      ///< The number of nearest neighbours that defines 'local' in the Local NHits trigger
  G4int localNHitsThreshold;       ///< The threshold for the Local NHits trigger
  G4int localNHitsWindow;          ///< The time window for the Local NHits trigger
  G4bool localNHitsAdjustForNoise; ///< Automatically adjust the Local NHits trigger threshold based on the average dark noise rate?
  //regions
  G4int  regionsNBinsP;          ///< Number of phi bins to use when splitting the side PMTs into regions
  G4int  regionsNBinsZ;          ///< Number of Z bins to use when splitting the side PMTs into regions
  G4int  regionsNRings;          ///< Number of rings to use when splitting the top/bottom PMTs into regions. Note the central circle is not a ring (i.e. regionsNRings=0 is valid)
  G4int  regionsNCentralSectors; ///< Number of sectors to split the central circle into when splitting the top/bottom PMTs into regions
  G4int  regionsNRingSectors;    ///< Number of additional segments (relative to each inner segment) to split the ring into when splitting the top/bottom PMTs into regions
  //Save failures
  G4int    saveFailuresMode;              ///< The mode for saving events which don't pass triggers
  G4double saveFailuresTime;              ///< The dummy trigger time for failed events
  G4int    saveFailuresPreTriggerWindow;  ///< The pretrigger window to save before an SaveFailures trigger
  G4int    saveFailuresPostTriggerWindow; ///< The posttrigger window to save after an SaveFailures trigger

  G4String triggerClassName; ///< Save the name of the trigger class

  unsigned int nPMTs;    ///< Store the number of PMTs in the geometry

  int * localNHitsHits; ///< Array to store the number of hits in a time window for each PMT tube id
  int * regionsHits;    ///< Array to store the number of hits in a time window for each PMT region

private:

  ///calculate the average dark noise occupancy (used to modify the NHits threshold)
  int CalculateAverageDarkNoiseOccupancy(int npmts, int window);
  ///modify the NDigits threshold based on the average dark noise rate
  void AdjustNDigitsThresholdForNoise();
  ///modify the LocalNDigits threshold based on the average dark noise rate
  void AdjustLocalNDigitsThresholdForNoise();

  ///takes all trigger times, then loops over all Digits & fills the output DigitsCollection
  void FillDigitsCollection(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits, TriggerType_t save_triggerType);

  static const double offset;        ///< Hit time offset (ns)
  static const double LongTime;      ///< An arbitrary long time to use in loops (ns)

  bool   digitizeCalled; ///< Has Digitize() been called yet?
};

// *******************************************
// CONTAINER CLASS
// *******************************************

class WCSimWCDigiTrigger : public G4VDigi
{
public:

  WCSimWCDigiTrigger();
  ~WCSimWCDigiTrigger();
  WCSimWCDigiTrigger(const WCSimWCDigiTrigger&);
  const WCSimWCDigiTrigger& operator=(const WCSimWCDigiTrigger&);
  int operator==(const WCSimWCDigiTrigger&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void Draw() {}
  void Print();

  inline void SetTubeID(G4int tube) { tubeID = tube; }
  inline void AddGate  (G4int gate) { Gates.insert(gate); }
  inline void AddPe    ()           { totalPe++; }
  inline void SetPe    (G4int gate, G4float Q) {   pe.insert(std::pair<int,float>(gate,Q)); }
  inline void SetTime  (G4int gate, G4float T) { time.insert(std::pair<int,float>(gate,T)); }

  /// Add a whole vector for one digit to fDigiComp. Clear input vector once added.
  void AddDigiCompositionInfo(G4int gate, std::vector<int> &digi_comp){
    fDigiComp.insert(std::pair<int, std::vector<int> >(gate, digi_comp));
    digi_comp.clear();
  }

  inline G4int   GetTubeID() {return tubeID;}
  inline std::vector<G4float> GetPe      (int gate) { return FindInMultimap(gate, pe); }
  inline std::vector<G4float> GetTime    (int gate) { return FindInMultimap(gate, time); }
  std::vector<std::vector<int> > GetDigiCompositionInfo(int gate)
  {
    std::vector<std::vector<int> > v;
    std::multimap<int, std::vector<int> >::iterator it = fDigiComp.begin();
    for (; it != fDigiComp.end(); ++it) {
      if((it->first) == gate)
	v.push_back(it->second);
    }
    return v;
  }

  inline int NumberOfGates()     { return Gates.size();      }
  inline int NumberOfSubEvents() { return Gates.size() - 1;  }
  inline bool HasHitsInGate(int gate) { return Gates.count(gate) > 0; }

private:
  G4int tubeID; ///< PMT id of the digit

  std::set<int> Gates;  ///< 'Gates' specifies subevent

  //lists (meaning multimap) of information for each digit created on the PMT
  std::multimap<int,float> pe;   ///< Digit charge
  std::multimap<int,float> time; ///< Digit time
  std::multimap<int, std::vector<int> > fDigiComp;   ///< Stores the unique IDs of each photon making up a digit

  //integrated hit/digit parameters
  G4int                 totalPe; ///< Total charge on digit

  template <typename T> std::vector<T> FindInMultimap(const int compare, typename std::multimap<int,T> &map)
  {
    typename std::vector<T> v;
    typename std::multimap<int,T>::iterator it = map.begin();
    for (; it != map.end(); ++it) {
      if((it->first) == compare)
	v.push_back(it->second);
    }
    return v;
  }

};

extern G4Allocator<WCSimWCDigiTrigger> WCSimWCDigiTriggerAllocator;

inline void* WCSimWCDigiTrigger::operator new(size_t)
{
  void* aDigi;
  aDigi = (void*) WCSimWCDigiTriggerAllocator.MallocSingle();
  return aDigi;
}

inline void WCSimWCDigiTrigger::operator delete(void* aDigi)
{
  WCSimWCDigiTriggerAllocator.FreeSingle((WCSimWCDigiTrigger*) aDigi);
}



// *******************************************
// DERIVED CLASSES
// *******************************************


/**
 * \class WCSimWCTriggerNDigits
 *
 * \brief A simple NDigits trigger class
 *
 */

class WCSimWCTriggerNDigits : public WCSimWCTriggerBase
{
public:

  ///Create WCSimWCTriggerNDigits instance with knowledge of the detector and DAQ options
  WCSimWCTriggerNDigits(G4String name, WCSimDetectorConstruction*, WCSimWCDAQMessenger*);

  ~WCSimWCTriggerNDigits();

private:
  ///Calls the workhorse of this class: AlgNDigits
  void DoTheWork(WCSimWCDigitsCollection* WCDCPMT);

  bool GetDefaultMultiDigitsPerTrigger()    { return false; } ///< SKI saves only earliest digit on a PMT in the trigger window
  int  GetDefaultNDigitsWindow()            { return 200;   } ///< SK max light travel time ~200 ns
  int  GetDefaultNDigitsThreshold()         { return 25;    } ///< SK NDigits threshold ~25
  int  GetDefaultNDigitsPreTriggerWindow()  { return -400;  } ///< SK SLE trigger window ~-400
  int  GetDefaultNDigitsPostTriggerWindow() { return 950;   } ///< SK SLE trigger window ~+950

  void WriteGeomInfo() {}
  void ReadGeomInfo()  {}
};


/**
 * \class WCSimWCTriggerNDigits2
 *
 * \brief An (incomplete) example of running two trigger algorithms, one after the other
 *
 */

class WCSimWCTriggerNDigits2 : public WCSimWCTriggerBase
{
public:

  //not recommended to override these methods
  WCSimWCTriggerNDigits2(G4String name, WCSimDetectorConstruction*, WCSimWCDAQMessenger*);
  ~WCSimWCTriggerNDigits2();

private:
  void DoTheWork(WCSimWCDigitsCollection* WCDCPMT);

  bool GetDefaultMultiDigitsPerTrigger()    { return false; } ///< SKI saves only earliest digit on a PMT in the trigger window
  int  GetDefaultNDigitsWindow()            { return 200;   } ///< SK max light travel time ~200 ns
  int  GetDefaultNDigitsThreshold()         { return 50;    } ///< 2 * SK NDigits threshold ~25
  int  GetDefaultNDigitsPreTriggerWindow()  { return -400;  } ///< SK SLE trigger window ~-400
  int  GetDefaultNDigitsPostTriggerWindow() { return 950;   } ///< SK SLE trigger window ~+950

  void WriteGeomInfo() {}
  void ReadGeomInfo()  {}
};

/**
 * \class WCSimWCTriggerNHitsThenLocalNHits
 *
 * \brief A trigger class which first looks for a global NHits trigger, then searches for a local NHits trigger using fraction of nearest neighbours
 *
 */

class WCSimWCTriggerNHitsThenLocalNHits : public WCSimWCTriggerBase
{
public:

  ///Create WCSimWCTriggerNHitsThenLocalNHits instance with knowledge of the detector and DAQ options
  WCSimWCTriggerNHitsThenLocalNHits(G4String name, WCSimDetectorConstruction*, WCSimWCDAQMessenger*);

  ~WCSimWCTriggerNHitsThenLocalNHits();

private:
  ///Calls the workhorse of this class: AlgNHitsThenLocalNHits
  void DoTheWork(WCSimWCDigitsCollection* WCDCPMT);

  bool GetDefaultMultiDigitsPerTrigger()    { return false; } ///< SKI saves only earliest digit on a PMT in the trigger window
  //int  GetDefaultNDigitsWindow()            { return 200;   } ///< SK max light travel time ~200 ns
  //int  GetDefaultNDigitsThreshold()         { return 50;    } ///< 2 * SK NDigits threshold ~25
  int  GetDefaultNDigitsPreTriggerWindow()  { return -400;  } ///< SK SLE trigger window ~-400
  int  GetDefaultNDigitsPostTriggerWindow() { return 950;   } ///< SK SLE trigger window ~+950

  void WriteGeomInfo();
  void ReadGeomInfo();
};

/**
 * \class WCSimWCTriggerNHitsThenRegions
 *
 * \brief A trigger class which first looks for a global NHits trigger, then searches for a local NHits trigger using regions
 *
 */

class WCSimWCTriggerNHitsThenRegions : public WCSimWCTriggerBase
{
public:

  ///Create WCSimWCTriggerNHitsThenRegions instance with knowledge of the detector and DAQ options
  WCSimWCTriggerNHitsThenRegions(G4String name, WCSimDetectorConstruction*, WCSimWCDAQMessenger*);

  ~WCSimWCTriggerNHitsThenRegions();

private:
  ///Calls the workhorse of this class: AlgNHitsThenRegions
  void DoTheWork(WCSimWCDigitsCollection* WCDCPMT);

  bool GetDefaultMultiDigitsPerTrigger()    { return false; } ///< SKI saves only earliest digit on a PMT in the trigger window
  //int  GetDefaultNDigitsWindow()            { return 200;   } ///< SK max light travel time ~200 ns
  //int  GetDefaultNDigitsThreshold()         { return 50;    } ///< 2 * SK NDigits threshold ~25
  int  GetDefaultNDigitsPreTriggerWindow()  { return -400;  } ///< SK SLE trigger window ~-400
  int  GetDefaultNDigitsPostTriggerWindow() { return 950;   } ///< SK SLE trigger window ~+950

  void WriteGeomInfo() {};
  void ReadGeomInfo()  {};
};



#endif //WCSimWCTrigger_h
