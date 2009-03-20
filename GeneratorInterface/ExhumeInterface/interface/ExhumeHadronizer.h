#ifndef gen_ExhumeHadronizer_h
#define gen_ExhumeHadronizer_h

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

namespace Exhume{
class Event;
class CrossSection;
}

namespace gen
{
  //class Pythia6Hadronizer;

  class ExhumeHadronizer : public BaseHadronizer
  {
  
  public:
     ExhumeHadronizer(edm::ParameterSet const& ps);
     ~ExhumeHadronizer();

     // bool generatePartons();
     bool generatePartonsAndHadronize();
     bool hadronize();
     bool decay();
     bool residualDecay();
     bool initializeForExternalPartons();
     bool initializeForInternalPartons();
     bool declareStableParticles( const std::vector<int> );
     
     void finalizeEvent();

     void statistics();

     const char* classname() const;
     
  private:
     double comEnergy_;
                    
     //edm::ParameterSet processPSet_;
     //edm::ParameterSet paramsPSet_;
     edm::ParameterSet myPSet_;
    
     bool            hepMCVerbosity_;
     unsigned int    maxEventsToPrint_;
     unsigned int    pythiaListVerbosity_;
     
     bool convertToPDG_;

     //Pythia6Hadronizer* pythia6Hadronizer_;
     Exhume::Event* exhumeEvent_;	
     Exhume::CrossSection* exhumeProcess_;
  };
}

#endif
