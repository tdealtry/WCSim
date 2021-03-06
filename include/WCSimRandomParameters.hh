#ifndef WCSimRandomParameters_h
#define WCSimRandomParameters_h 1
#include "WCSimRandomMessenger.hh"
#include "WCSimRootOptions.hh"
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RanluxEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RanecuEngine.h"
#include "TRandom3.h"

class WCSimRandomParameters
{

public:
  WCSimRandomParameters() 
    {
      seed=31415; 
      generator=RANDOM_E_HEPJAMES;
      RandomMessenger = new WCSimRandomMessenger(this);
    }
  ~WCSimRandomParameters() {delete RandomMessenger;}

  WCSimRandomGenerator_t GetGenerator() {return generator; }
  void SetGenerator(WCSimRandomGenerator_t rng)
  {
    switch (rng)
      {
        case RANDOM_E_RANLUX:
	  {
	    G4cout << "Setting the random number generator to RANLUX" << G4endl;
	    CLHEP::RanluxEngine *newluxengine = new CLHEP::RanluxEngine(31415,4);  // highest luxury level
	    CLHEP::HepRandom::setTheEngine(newluxengine);
	  }
	  break;
      case RANDOM_E_RANECU:
	  {
	    G4cout << "Setting the random number generator to RANECU" << G4endl;
	    CLHEP::RanecuEngine *newecuengine = new CLHEP::RanecuEngine();
	    CLHEP::HepRandom::setTheEngine(newecuengine);
	  }
	  break;
	  
        case RANDOM_E_HEPJAMES:
	  {
	    G4cout << "Setting the random number generator to HEPJAMES" << G4endl;
	    CLHEP::HepJamesRandom *newjamesengine = new CLHEP::HepJamesRandom();
	    CLHEP::HepRandom::setTheEngine(newjamesengine);
	  }  
	  break;
      default:
	{
	  G4cout << "Random number generator type not understood: " << rng << G4endl;
	  exit(0);
	}
      }
  };
  int GetSeed() {return CLHEP::HepRandom::getTheSeed();}
  void SetSeed(int iseed) 
    { 
      CLHEP::HepRandom::setTheSeed(iseed);
      gRandom->SetSeed(iseed);
      G4cout << "Setting the Random Seed to: " << iseed << G4endl;
      seed = iseed;
    }

  void SaveOptionsToOutput(WCSimRootOptions * wcopt)
  {
    wcopt->SetRandomSeed(seed);
    wcopt->SetRandomGenerator(generator);
  }

private:
  WCSimRandomGenerator_t generator;
  int seed;
  WCSimRandomMessenger *RandomMessenger;
};

#endif
