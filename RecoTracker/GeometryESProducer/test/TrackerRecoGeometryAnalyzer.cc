// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"

using namespace std;

//
//
// class decleration
//

class TrackerRecoGeometryAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TrackerRecoGeometryAnalyzer( const edm::ParameterSet& );
      ~TrackerRecoGeometryAnalyzer();


      virtual void analyze( const edm::Event&, const edm::EventSetup& );
   private:
      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrackerRecoGeometryAnalyzer::TrackerRecoGeometryAnalyzer( const edm::ParameterSet& iConfig )
{
   //now do what ever initialization is needed

}


TrackerRecoGeometryAnalyzer::~TrackerRecoGeometryAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TrackerRecoGeometryAnalyzer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
   using namespace edm;

   std::cout << "Here I am " << std::endl;
   //
   // get the GeometricSearchDet
   //
   edm::ESHandle<GeometricSearchTracker> track;
   iSetup.get<TrackerRecoGeometryRecord>().get( track );     
   
   //---- testing access to barrelLayers ----
   vector<BarrelDetLayer*> theBarrelLayers = track->barrelLayers();
   cout << "number of BarrelLayers: " << theBarrelLayers.size() << endl;

   for(unsigned int i=0; i<3; i++){
     BarrelDetLayer* theLayer = theBarrelLayers[i];   
     cout << "theLayer[" << i << "]->position().perp(): " 
	  << theLayer->components().front()->surface().position().perp() << endl;     
     //cout << "theLayer[" << i << "]->module(): " 
     //<< theLayer->module() << endl;     
   }   
   //--------------------------------------

}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackerRecoGeometryAnalyzer)
 
 
