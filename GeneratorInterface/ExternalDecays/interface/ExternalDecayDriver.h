#ifndef gen_ExternalDecayDriver_h
#define gen_ExternalDecayDriver_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"

namespace HepMC
{
class GenEvent;
}

namespace gen {

class TauolaInterface;
class EvtGenInterface;
class PhotosInterface;

   class ExternalDecayDriver 
   {
      public:
         
	 // ctor & dtor
	 ExternalDecayDriver( const edm::ParameterSet& );
	 ~ExternalDecayDriver();
	 
	 void init( const edm::EventSetup& );

	 const std::vector<int>& operatesOnParticles() { return fPDGs; }
	 
	 HepMC::GenEvent* decay( HepMC::GenEvent* );
	 
	 void statistics() const;
      
      private:
      	 
	 bool             fIsInitialized;
	 TauolaInterface* fTauolaInterface;
	 EvtGenInterface* fEvtGenInterface;
	 PhotosInterface* fPhotosInterface;
	 std::vector<int> fPDGs;
         
   };

}

#endif
