#ifndef gen_Pythia6Hadronizer_h
#define gen_Pythia6Hadronizer_h

// -*- C++ -*-

// class Pythia6Hadronizer is an example of a class that models the
// Hadronizer concept.

#include <memory>

#include <boost/shared_ptr.hpp>

#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"

#include "GeneratorInterface/Core/interface/BaseHadronizer.h"

namespace lhef
{
class LHERunInfo;
class LHEEvent;
}

class LHEEventProduct;

namespace HepMC
{
class GenEvent;
}

namespace CLHEP {
class HepRandomEngine;
class RandFlat;
}


namespace gen
{

class Pythia6Service;
class JetMatching;

  class Pythia6Hadronizer : public BaseHadronizer
  {
  
  public:
     Pythia6Hadronizer(edm::ParameterSet const& ps);
     ~Pythia6Hadronizer();

     // bool generatePartons();
     bool generatePartonsAndHadronize();
     bool hadronize();
     bool decay();
     bool residualDecay();
     bool initializeForExternalPartons();
     bool initializeForInternalPartons();
     bool declareStableParticles( const std::vector<int> );
     
     static JetMatching* getJetMatching() { return fJetMatching; }
          
     void finalizeEvent();

     void statistics();

     const char* classname() const;
     
  private:
     
     Pythia6Service* fPy6Service;

     // the following 3 params are common for all generators(interfaces)
     // probably better to wrap them up in a class and reuse ?
     // (the event/run pointers are already moved to BaseHadronizer)
     //
     double fCOMEnergy ;  // this one is irrelevant for setting py6 as hadronizer
                          // or if anything, it should be picked up from LHERunInfoProduct !

     static JetMatching* fJetMatching; 

     bool            fHepMCVerbosity;
     unsigned int    fMaxEventsToPrint ;
           
     // this is the only one specific to Pythia6
     //
     unsigned int    fPythiaListVerbosity ;
     
     // these two params control stop- and r-handron features,
     // that are "custom" add-ons to Py6; 
     // I doubt they should drag along Py6Int main library...
     //
     bool fStopHadronsEnabled;
     bool fGluinoHadronsEnabled;
     
     // and final touch - conversion of Py6 PID's into PDG convension 
     bool fConvertToPDG;
               
  };
}

#endif
