
#include "GeneratorInterface/ExternalDecays/interface/ExternalDecayDriver.h"

#include "GeneratorInterface/ExternalDecays/interface/EvtGenInterface.h"
#include "GeneratorInterface/ExternalDecays/interface/TauolaInterface.h"
#include "GeneratorInterface/ExternalDecays/interface/PhotosInterface.h"

#include "HepMC/GenEvent.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"

using namespace gen;
using namespace edm;

ExternalDecayDriver::ExternalDecayDriver( const ParameterSet& pset )
   : fIsInitialized(false),
     fTauolaInterface(0),
     fEvtGenInterface(0),
     fPhotosInterface(0)
{
    
    std::vector<std::string> extGenNames =
       pset.getParameter< std::vector<std::string> >("parameterSets");

    Service<RandomNumberGenerator> rng;
    if(!rng.isAvailable()) {
       throw cms::Exception("Configuration")
       << "The RandomNumberProducer module requires the RandomNumberGeneratorService\n"
          "which appears to be absent.  Please add that service to your configuration\n"
          "or remove the modules that require it." << std::endl;
    }      
    
    for (unsigned int ip=0; ip<extGenNames.size(); ++ip )
    {
      std::string curSet = extGenNames[ip];
      if ( curSet == "EvtGen" )
      {
         fEvtGenInterface = new gen::EvtGenInterface(pset.getUntrackedParameter< ParameterSet >(curSet));
      }
      else if ( curSet == "Tauola" )
      {
         fTauolaInterface = new gen::TauolaInterface(pset.getUntrackedParameter< ParameterSet >(curSet));
	 tauolaRandomEngine = &rng->getEngine();
      }

      else if ( curSet == "Photos" )
      {
         fPhotosInterface = new gen::PhotosInterface();
      }

    }

}

ExternalDecayDriver::~ExternalDecayDriver()
{
   if ( fEvtGenInterface ) delete fEvtGenInterface;
   if ( fTauolaInterface ) delete fTauolaInterface;
   if ( fPhotosInterface ) delete fPhotosInterface;
}

HepMC::GenEvent* ExternalDecayDriver::decay( HepMC::GenEvent* evt )
{
   
   if ( !fIsInitialized ) return evt;
   
   if ( fEvtGenInterface )
   {  
      evt = fEvtGenInterface->decay( evt ); 
      if ( !evt ) return 0;
   }

   if ( fTauolaInterface ) 
   {
      evt = fTauolaInterface->decay( evt ); 
      if ( !evt ) return 0;
   }
   

   if ( fPhotosInterface )
   {
      evt = fPhotosInterface->apply( evt );
      if ( !evt ) return 0;
   }

         
   return evt;
}

void ExternalDecayDriver::init( const edm::EventSetup& es )
{

   if ( fIsInitialized ) return;
   
   if ( fTauolaInterface ) 
   {
      fTauolaInterface->init( es );
      for ( std::vector<int>::const_iterator i=fTauolaInterface->operatesOnParticles().begin();
            i!=fTauolaInterface->operatesOnParticles().end(); i++ ) 
               fPDGs.push_back( *i );
   }
   
   if ( fEvtGenInterface ) 
   {
      fEvtGenInterface->init();
      for ( std::vector<int>::const_iterator i=fEvtGenInterface->operatesOnParticles().begin();
            i!=fEvtGenInterface->operatesOnParticles().end(); i++ )
               fPDGs.push_back( *i );
   }
   

   if ( fPhotosInterface )
   {
      fPhotosInterface->init();
   }

   
   fIsInitialized = true;
   
   return;
}

void ExternalDecayDriver::statistics() const
{
   if ( fTauolaInterface ) fTauolaInterface->statistics();
   // similar for EvtGen and/or Photos, if needs be
   return;
}
