//////////////////////////
//  Producer by Anders  //
//     and Emmanuele    //
//    july 2012 @ CU    //
//////////////////////////


#ifndef L1TTRACK_PRDC_H
#define L1TTRACK_PRDC_H

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//
// #include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Ref.h"
//
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 
#include "DataFormats/Common/interface/DetSetVector.h"
//
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
#include "SimDataFormats/SLHC/interface/slhcevent.hh"
//#include "SimDataFormats/SLHC/interface/L1TRod.hh"
//#include "SimDataFormats/SLHC/interface/L1TSector.hh"
#include "SimDataFormats/SLHC/interface/L1TBarrel.hh"
#include "SimDataFormats/SLHC/interface/L1TDisk.hh"
#include "SimDataFormats/SLHC/interface/L1TStub.hh"
//#include "SimDataFormats/SLHC/interface/L1TWord.hh"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
//
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
//
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/GeomPropagators/interface/HelixArbitraryPlaneCrossing.h"
////////////////////////
// FAST SIMULATION STUFF
#include "FastSimulation/Particle/interface/RawParticle.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
//
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerDetUnit.h"
#include "DataFormats/SiPixelDetId/interface/StackedTrackerDetId.h"
//#include "SLHCUpgradeSimulations/Utilities/interface/classInfo.h" REMOVE

////////////////
// PHYSICS TOOLS
#include "RecoTracker/TkSeedGenerator/interface/FastHelix.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
//
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementVector.h"
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
//#include "SLHCUpgradeSimulations/Utilities/interface/constants.h" REMOVE

//////////////
// STD HEADERS
#include <memory>
#include <string>
#include <iostream>
#include <fstream>

//////////////
// NAMESPACES
// using namespace std;
// using namespace reco;
using namespace edm;
//using namespace cmsUpgrades;


//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

/////////////////////////////////////
// this class is needed to make a map
// between different types of stubs
class L1TStubCompare 
{
public:
  bool operator()(const L1TStub& x, const L1TStub& y) {
    if (x.layer() != y.layer()) return (y.layer()-x.layer())>0;
    else {
      if (x.ladder() != y.ladder()) return (y.ladder()-x.ladder())>0;
      else {
	if (x.module() != y.module()) return (y.module()-x.module())>0;
	else {
	  if (x.iz() != y.iz()) return (y.iz()-x.iz())>0;
	  else return (x.iphi()-y.iphi())>0;
	}
      }
    }
  }
};


class L1TrackProducer : public edm::EDProducer
{
public:

  typedef L1TkStub_PixelDigi_                           L1TkStubType;
  typedef std::vector< L1TkStubType >                                L1TkStubCollectionType;
  typedef edm::Ptr< L1TkStubType >                                   L1TkStubPtrType;
  typedef std::vector< L1TkStubPtrType >                             L1TkStubPtrCollection;
  typedef std::vector< L1TkStubPtrCollection >                       L1TkStubPtrCollVectorType;

  //typedef L1TkTracklet_PixelDigi_                       L1TkTrackletType;
  //typedef std::vector< L1TkTrackletType >                            L1TkTrackletCollectionType;
  //typedef edm::Ptr< L1TkTrackletType >                               L1TkTrackletPtrType;

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;
  typedef std::vector< L1TTrack >                                    L1TrackCollectionType;

  /// Constructor/destructor
  explicit L1TrackProducer(const edm::ParameterSet& iConfig);
  virtual ~L1TrackProducer();

protected:
                     
private:
  GeometryMap geom;
  int eventnum;

  /// Containers of parameters passed by python configuration file
  edm::ParameterSet config;

  /// ///////////////// ///
  /// MANDATORY METHODS ///
  virtual void beginRun( edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );
};


//////////////
// CONSTRUCTOR
L1TrackProducer::L1TrackProducer(edm::ParameterSet const& iConfig) // :   config(iConfig)
{
  //produces<L1TrackCollectionType>( "Level1Tracks" ).setBranchAlias("Level1Tracks");
  produces< L1TkStubPtrCollVectorType >( "L1TkStubs" ).setBranchAlias("L1TkStubs");
  produces< L1TkTrackCollectionType >( "Level1TkTracks" ).setBranchAlias("Level1TkTracks");
  // produces<L1TkStubMapType>( "L1TkStubMap" ).setBranchAlias("L1TkStubMap");
  // produces< L1TkTrackletCollectionType >( "L1TkTracklets" ).setBranchAlias("L1TkTracklets");
}

/////////////
// DESTRUCTOR
L1TrackProducer::~L1TrackProducer()
{
  /// Insert here what you need to delete
  /// when you close the class instance
}  

//////////
// END JOB
void L1TrackProducer::endRun(edm::Run& run, const edm::EventSetup& iSetup)
{
  /// Things to be done at the exit of the event Loop 

}

////////////
// BEGIN JOB
void L1TrackProducer::beginRun(edm::Run& run, const edm::EventSetup& iSetup )
{
  eventnum=0;
  std::cout << "L1TrackProducer" << std::endl;
}

//////////
// PRODUCE
void L1TrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  typedef std::map< L1TStub, L1TkStubPtrType, L1TStubCompare > stubMapType;

  /// Prepare output
  //std::auto_ptr< L1TrackCollectionType > L1TracksForOutput( new L1TrackCollectionType );
  std::auto_ptr< L1TkStubPtrCollVectorType > L1TkStubsForOutput( new L1TkStubPtrCollVectorType );
  //std::auto_ptr< L1TkTrackletCollectionType > L1TkTrackletsForOutput( new L1TkTrackletCollectionType );
  std::auto_ptr< L1TkTrackCollectionType > L1TkTracksForOutput( new L1TkTrackCollectionType );

  /// Geometry handles etc
  edm::ESHandle<TrackerGeometry>                               geometryHandle;
  const TrackerGeometry*                                       theGeometry;
  edm::ESHandle<StackedTrackerGeometry>           stackedGeometryHandle;
  const StackedTrackerGeometry*                   theStackedGeometry;
  StackedTrackerGeometry::StackContainerIterator  StackedTrackerIterator;

  /// Geometry setup
  /// Set pointers to Geometry
  iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);
  theGeometry = &(*geometryHandle);
  /// Set pointers to Stacked Modules
  iSetup.get<StackedTrackerGeometryRecord>().get(stackedGeometryHandle);
  theStackedGeometry = stackedGeometryHandle.product(); /// Note this is different 
                                                        /// from the "global" geometry
  ////////////////////////
  // GET MAGNETIC FIELD //
  edm::ESHandle<MagneticField> magneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
  const MagneticField* theMagneticField = magneticFieldHandle.product();
  double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();

  ////////////
  // GET BS //
  ////////////
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel("BeamSpotFromSim","BeamSpot",recoBeamSpotHandle);
  math::XYZPoint bsPosition=recoBeamSpotHandle->position();

  cout << "L1TrackProducer: B="<<mMagneticFieldStrength
       <<" vx reco="<<bsPosition.x()
       <<" vy reco="<<bsPosition.y()
       <<" vz reco="<<bsPosition.z()
       <<endl;

  SLHCEvent ev;
  ev.setIPx(bsPosition.x());
  ev.setIPy(bsPosition.y());
  eventnum++;


  ///////////////////
  // GET SIMTRACKS //
  edm::Handle<edm::SimTrackContainer>   simTrackHandle;
  edm::Handle<edm::SimVertexContainer>  simVtxHandle;
  //iEvent.getByLabel( "famosSimHits", simTrackHandle );
  //iEvent.getByLabel( "famosSimHits", simVtxHandle );
  iEvent.getByLabel( "g4SimHits", simTrackHandle );
  iEvent.getByLabel( "g4SimHits", simVtxHandle );

  //////////////////////
  // GET MC PARTICLES //
  edm::Handle<reco::GenParticleCollection> genpHandle;
  iEvent.getByLabel( "genParticles", genpHandle );

  /////////////////////
  // GET PIXEL DIGIS //
  edm::Handle<edm::DetSetVector<PixelDigi> >         pixelDigiHandle;
  edm::Handle<edm::DetSetVector<PixelDigiSimLink> >  pixelDigiSimLinkHandle;
  iEvent.getByLabel("simSiPixelDigis", pixelDigiHandle);
  iEvent.getByLabel("simSiPixelDigis", pixelDigiSimLinkHandle);

  ////////////////////////
  // GET THE PRIMITIVES //
  edm::Handle<L1TkCluster_PixelDigi_Collection>  pixelDigiL1TkClusterHandle;
  edm::Handle<L1TkStub_PixelDigi_Collection>     pixelDigiL1TkStubHandle;
  iEvent.getByLabel("L1TkClustersFromPixelDigis", pixelDigiL1TkClusterHandle);
  iEvent.getByLabel("L1TkStubsFromPixelDigis", "StubsPass", pixelDigiL1TkStubHandle);

  ////////////////////////
  /// LOOP OVER SimTracks
  SimTrackContainer::const_iterator iterSimTracks;
  for ( iterSimTracks = simTrackHandle->begin();
	iterSimTracks != simTrackHandle->end();
	++iterSimTracks ) {

    /// Get the corresponding vertex
    int vertexIndex = iterSimTracks->vertIndex();
    const SimVertex& theSimVertex = (*simVtxHandle)[vertexIndex];
    math::XYZTLorentzVectorD trkVtxPos = theSimVertex.position();
    GlobalPoint trkVtxCorr = GlobalPoint( trkVtxPos.x() - bsPosition.x(), 
					  trkVtxPos.y() - bsPosition.y(), 
					  trkVtxPos.z() - bsPosition.z() );
    
    double pt=iterSimTracks->momentum().pt();
    if (pt!=pt) pt=9999.999;
    ev.addL1SimTrack(iterSimTracks->trackId(),iterSimTracks->type(),pt,
		     iterSimTracks->momentum().eta(), 
		     iterSimTracks->momentum().phi(), 
		     trkVtxCorr.x(),
		     trkVtxCorr.y(),
		     trkVtxCorr.z());
   
    
  } /// End of Loop over SimTracks


  std::cout << "Will loop over digis:"<<std::endl;

  DetSetVector<PixelDigi>::const_iterator iterDet;
  for ( iterDet = pixelDigiHandle->begin();
        iterDet != pixelDigiHandle->end();
        iterDet++ ) {

    /// Build Detector Id
    DetId tkId( iterDet->id );
    StackedTrackerDetId stdetid(tkId);
    /// Check if it is Pixel
    if ( tkId.subdetId() == 2 ) {

      cout << "Will create pxfId"<<endl;
      PXFDetId pxfId(tkId);
      DetSetVector<PixelDigiSimLink>::const_iterator itDigiSimLink1=pixelDigiSimLinkHandle->find(pxfId.rawId());
      if (itDigiSimLink1!=pixelDigiSimLinkHandle->end()){
	cout << "Found forward digisim link"<<endl;
	DetSet<PixelDigiSimLink> digiSimLink = *itDigiSimLink1;
	//DetSet<PixelDigiSimLink> digiSimLink = (*pixelDigiSimLinkHandle)[ pxfId.rawId() ];
	DetSet<PixelDigiSimLink>::const_iterator iterSimLink;
	/// Renormalize layer number from 5-14 to 0-9 and skip if inner pixels

	int disk = pxfId.disk();
	
	if (disk<4) {
	  continue;
	}

	disk-=3;
	
	// Layer 0-20
	//DetId digiDetId = iterDet->id;
	//int sensorLayer = 0.5*(2*PXFDetId(digiDetId).layer() + (PXFDetId(digiDetId).ladder() + 1)%2 - 8);
	
	/// Loop over PixelDigis within Module and select those above threshold
	DetSet<PixelDigi>::const_iterator iterDigi;
	for ( iterDigi = iterDet->data.begin();
	      iterDigi != iterDet->data.end();
	      iterDigi++ ) {
      
	  /// Threshold (here it is NOT redundant)
	  if ( iterDigi->adc() <= 30 ) continue;
	    
	  /// Try to learn something from PixelDigi position
	  cout << "Here 010" <<endl;
	  const GeomDetUnit* gDetUnit = theGeometry->idToDetUnit( tkId );
	  cout << "Here 011" <<endl;
	  MeasurementPoint mp( iterDigi->row() + 0.5, iterDigi->column() + 0.5 );
	  GlobalPoint pdPos = gDetUnit->surface().toGlobal( gDetUnit->topology().localPosition( mp ) ) ;
	    
	  cout << "pxfId.side():"<<pxfId.side()<<endl;

	  int offset=1000;

	  if (pxfId.side()==1) {
	    offset=2000;
	  }

	  assert(pxfId.panel()==1);

	  vector<int> simtrackids;
	  /// Loop over PixelDigiSimLink to find the
	  /// correct link to the SimTrack collection
	  for ( iterSimLink = digiSimLink.data.begin();
		iterSimLink != digiSimLink.data.end();
		iterSimLink++) {
	        
	    /// When the channel is the same, the link is found
	    if ( (int)iterSimLink->channel() == iterDigi->channel() ) {
	            
	      /// Map wrt SimTrack Id
	      unsigned int simTrackId = iterSimLink->SimTrackId();
	      simtrackids.push_back(simTrackId); 
	    }
	  }
	ev.addDigi(offset+disk,iterDigi->row(),iterDigi->column(),
		   pxfId.blade(),pxfId.panel(),pxfId.module(),
		   pdPos.x(),pdPos.y(),pdPos.z(),simtrackids);
	}
      }
      cout << "Done with pxfId"<<endl;
    }

    if ( tkId.subdetId() == 1 ) {
      /// Get the PixelDigiSimLink corresponding to this one
      PXBDetId pxbId(tkId);
      DetSetVector<PixelDigiSimLink>::const_iterator itDigiSimLink=pixelDigiSimLinkHandle->find(pxbId.rawId());
      if (itDigiSimLink==pixelDigiSimLinkHandle->end()){
	continue;
      }
      DetSet<PixelDigiSimLink> digiSimLink = *itDigiSimLink;
      //DetSet<PixelDigiSimLink> digiSimLink = (*pixelDigiSimLinkHandle)[ pxbId.rawId() ];
      DetSet<PixelDigiSimLink>::const_iterator iterSimLink;
      /// Renormalize layer number from 5-14 to 0-9 and skip if inner pixels
      if ( pxbId.layer() < 5 ) {
	continue;
	
      }

      // Layer 0-20
      DetId digiDetId = iterDet->id;
      int sensorLayer = 0.5*(2*PXBDetId(digiDetId).layer() + (PXBDetId(digiDetId).ladder() + 1)%2 - 8);
      
      /// Loop over PixelDigis within Module and select those above threshold
      DetSet<PixelDigi>::const_iterator iterDigi;
      for ( iterDigi = iterDet->data.begin();
	    iterDigi != iterDet->data.end();
	    iterDigi++ ) {
	
	/// Threshold (here it is NOT redundant)
	if ( iterDigi->adc() <= 30 ) continue;
	
	/// Try to learn something from PixelDigi position
	cout << "Here 010" <<endl;
	const GeomDetUnit* gDetUnit = theGeometry->idToDetUnit( tkId );
	cout << "Here 011" <<endl;
	MeasurementPoint mp( iterDigi->row() + 0.5, iterDigi->column() + 0.5 );
	GlobalPoint pdPos = gDetUnit->surface().toGlobal( gDetUnit->topology().localPosition( mp ) ) ;
	
	/// Loop over PixelDigiSimLink to find the
	/// correct link to the SimTrack collection
	vector<int > simtrackids;
	for ( iterSimLink = digiSimLink.data.begin();
	      iterSimLink != digiSimLink.data.end();
	      iterSimLink++) {
	    
	  /// When the channel is the same, the link is found
	  if ( (int)iterSimLink->channel() == iterDigi->channel() ) {
	        
	    /// Map wrt SimTrack Id
	    unsigned int simTrackId = iterSimLink->SimTrackId();
	    simtrackids.push_back(simTrackId);
	  }
	}
	ev.addDigi(sensorLayer,iterDigi->row(),iterDigi->column(),
		   pxbId.layer(),pxbId.ladder(),pxbId.module(),
		   pdPos.x(),pdPos.y(),pdPos.z(),simtrackids);
      }
    }
  }    




  /// Loop over L1TkStubs
  L1TkStub_PixelDigi_Collection::const_iterator iterL1TkStub;
  for ( iterL1TkStub = pixelDigiL1TkStubHandle->begin();
	iterL1TkStub != pixelDigiL1TkStubHandle->end();
	++iterL1TkStub ) {

    double stubPt = theStackedGeometry->findRoughPt(mMagneticFieldStrength,&(*iterL1TkStub));
        
    if (stubPt>10000.0) stubPt=9999.99;
    GlobalPoint stubPosition = theStackedGeometry->findGlobalPosition(&(*iterL1TkStub));

    StackedTrackerDetId stubDetId = iterL1TkStub->getDetId();
    unsigned int iStack = stubDetId.iLayer();
    unsigned int iRing = stubDetId.iRing();
    unsigned int iPhi = stubDetId.iPhi();
    unsigned int iZ = stubDetId.iZ();

    std::vector<bool> innerStack;
    std::vector<int> irphi;
    std::vector<int> iz;
    std::vector<int> iladder;
    std::vector<int> imodule;


    if (iStack==999999) {
      iStack=1000+iRing;
    }

    /// Get the Inner and Outer L1TkCluster
    edm::Ptr<L1TkCluster_PixelDigi_> innerCluster = iterL1TkStub->getClusterPtr(0);
    edm::Ptr<L1TkCluster_PixelDigi_> outerCluster = iterL1TkStub->getClusterPtr(1);
      
          
    /// Loop over Detector Modules
    DetSetVector<PixelDigi>::const_iterator iterDet;
    for ( iterDet = pixelDigiHandle->begin();
	  iterDet != pixelDigiHandle->end();
	  iterDet++ ) {
      
      /// Build Detector Id
      DetId tkId( iterDet->id );
      
      const DetId innerDetId = theStackedGeometry->idToDet( innerCluster->getDetId(), 0 )->geographicalId();
      const DetId outerDetId = theStackedGeometry->idToDet( outerCluster->getDetId(), 1 )->geographicalId();
      
      if (innerDetId.rawId()==outerDetId.rawId()) {
	std::cerr<<"STUB DEBUGGING INNER LAYER == OUTER LAYER RAW ID"<<std::endl;
      }
    
      if (innerDetId.rawId()==tkId.rawId()) {
	// Layer 1-10.5
	//int sensorLayer = (2*PXBDetId(tkId).layer() + (PXBDetId(tkId).ladder() + 1)%2 - 8);

	/// Loop over PixelDigis within Module and select those above threshold
	DetSet<PixelDigi>::const_iterator iterDigi;
	for ( iterDigi = iterDet->data.begin();
	      iterDigi != iterDet->data.end();
	      iterDigi++ ) {

	  /// Threshold (here it is NOT redundant)
	  if ( iterDigi->adc() <= 30 ) continue;
          for (unsigned int ihit=0;ihit<innerCluster->getHits().size();ihit++){
	    if (iterDigi->channel() == innerCluster->getHits().at(ihit)->channel()) {
	      if (iStack<1000) {
		innerStack.push_back(true);
		irphi.push_back(iterDigi->row());
		iz.push_back(iterDigi->column());
		iladder.push_back(PXBDetId(tkId).ladder());
		imodule.push_back(PXBDetId(tkId).module());
	      }
	      else {
		innerStack.push_back(true);
		irphi.push_back(iterDigi->row());
		iz.push_back(iterDigi->column());
		iladder.push_back(PXFDetId(tkId).disk());
		imodule.push_back(PXFDetId(tkId).module());
	      }
	    }
	  }
	}
      }
    
      if (outerDetId.rawId()==tkId.rawId()) {
	// Layer 0-20
	//int sensorLayer = (2*PXBDetId(tkId).layer() + (PXBDetId(tkId).ladder() + 1)%2 - 8);

	/// Loop over PixelDigis within Module and select those above threshold
	DetSet<PixelDigi>::const_iterator iterDigi;
	for ( iterDigi = iterDet->data.begin();
	      iterDigi != iterDet->data.end();
	      iterDigi++ ) {

	  /// Threshold (here it is NOT redundant)
	  if ( iterDigi->adc() <= 30 ) continue;
          for (unsigned int ihit=0;ihit<outerCluster->getHits().size();ihit++){
	    if (iterDigi->channel() == outerCluster->getHits().at(ihit)->channel()) {
	      if (iStack<1000) {
		innerStack.push_back(false);
		irphi.push_back(iterDigi->row());
		iz.push_back(iterDigi->column());
		iladder.push_back(PXBDetId(tkId).ladder());
		imodule.push_back(PXBDetId(tkId).module());
	      }
	      else{
		innerStack.push_back(false);
		irphi.push_back(iterDigi->row());
		iz.push_back(iterDigi->column());
		iladder.push_back(PXFDetId(tkId).disk());
		imodule.push_back(PXFDetId(tkId).module());
	      }
	    }
	  }     
	}
      }
    }      


    ev.addStub(iStack-1,iPhi,iZ,stubPt,
	       stubPosition.x(),stubPosition.y(),stubPosition.z(),
	       innerStack,irphi,iz,iladder,imodule);
  }



  std::cout << "Will actually do L1 tracking:"<<std::endl;


  //////////////////////////
  // NOW RUN THE L1 tracking


  const int NSector=24;

  int NBarrel=6;

  L1TBarrel* L[6];

  if (1) {
    L[0]=new L1TBarrel(20.0,30.0,125.0,NSector);
    L[1]=new L1TBarrel(30.0,40.0,125.0,NSector);
    L[2]=new L1TBarrel(40.0,60.0,125.0,NSector);
    L[3]=new L1TBarrel(60.0,80.0,125.0,NSector);
    L[4]=new L1TBarrel(80.0,100.0,125.0,NSector);
    L[5]=new L1TBarrel(100.0,120.0,125.0,NSector);
  }
  else{
    L[0]=new L1TBarrel(20.0,30.0,125.0,NSector);
    L[1]=new L1TBarrel(30.0,40.0,125.0,NSector);
    L[2]=new L1TBarrel(40.0,65.0,125.0,NSector);
    L[3]=new L1TBarrel(65.0,80.0,125.0,NSector);
    L[4]=new L1TBarrel(80.0,100.0,125.0,NSector);
    L[5]=new L1TBarrel(100.0,120.0,125.0,NSector);
  }

  int NDisk=7;

  L1TDisk* F[7];

  F[0]=new L1TDisk(125.0,140.0,NSector);
  F[1]=new L1TDisk(140.0,157.0,NSector);
  F[2]=new L1TDisk(157.0,175.0,NSector);
  F[3]=new L1TDisk(175.0,195.0,NSector);
  F[4]=new L1TDisk(195.0,220.0,NSector);
  F[5]=new L1TDisk(220.0,245.0,NSector);
  F[6]=new L1TDisk(245.0,279.0,NSector);

  L1TDisk* B[7];

  B[0]=new L1TDisk(-125.0,-140.0,NSector);
  B[1]=new L1TDisk(-140.0,-157.0,NSector);
  B[2]=new L1TDisk(-157.0,-175.0,NSector);
  B[3]=new L1TDisk(-175.0,-195.0,NSector);
  B[4]=new L1TDisk(-195.0,-220.0,NSector);
  B[5]=new L1TDisk(-220.0,-245.0,NSector);
  B[6]=new L1TDisk(-245.0,-279.0,NSector);
      
  for (int j=0;j<ev.nstubs();j++){

    Stub stub=ev.stub(j);

    int simtrackid=ev.simtrackid(stub);
    //int simtrackid=-1;

    //cout << "Simtrack id:"<<simtrackid<<endl;
      

    //int simtrackid=ev.simtrackid(aStub);
    //int simtrackid=-1;

    //cout << "Simtrack id:"<<simtrackid<<endl;

    int layer=stub.layer()+1;
    int ladder=stub.ladder()+1;
    int module=stub.module();

    double sigmax=0.0029;
    double sigmaz=0.036;

    if (hypot(stub.x(),stub.y())>52) sigmaz=1.44;

    L1TStub aStub(simtrackid,-1,-1,layer,ladder,module,stub.x(),stub.y(),stub.z(),sigmax,sigmaz);

    //cout << "layer ladder module:"<<layer<<" "<<ladder<<" "<<module
    //   << " "<<aStub.r()<<" "<<aStub.phi()<<" "<<aStub.z()<<endl;
      

    int count_match=0;
    for (int ilayer=0;ilayer<NBarrel;ilayer++){
      if (L[ilayer]->addStub(aStub)) count_match++;
    }
    for (int idisk=0;idisk<NDisk;idisk++){
      if (F[idisk]->addStub(aStub)) count_match++;
      if (B[idisk]->addStub(aStub)) count_match++;
    }
      
    if (count_match!=1) {
      cout << "layer ladder module:"<<layer<<" "<<ladder<<" "<<module
	   << " "<<aStub.r()<<" "<<aStub.phi()<<" "<<aStub.z()
	   << " "<<count_match<<endl;
    }
    assert(count_match==1);

  }


  L1TTracks allTracks;


  L[0]->findTracklets(L[1]);
  L[0]->findMatches(L[2]);
  L[0]->findMatches(L[3]);
  L[0]->findMatches(L[4]);
  L[0]->findMatches(L[5]);
  L[0]->findMatches(F[0]);
  L[0]->findMatches(F[1]);
  L[0]->findMatches(F[2]);
  L[0]->findMatches(F[3]);
  L[0]->findMatches(F[4]);
  L[0]->findMatches(F[5]);
  L[0]->findMatches(F[6]);
  L[0]->findMatches(B[0]);
  L[0]->findMatches(B[1]);
  L[0]->findMatches(B[2]);
  L[0]->findMatches(B[3]);
  L[0]->findMatches(B[4]);
  L[0]->findMatches(B[5]);
  L[0]->findMatches(B[6]);
  L[0]->fitTracks();
    
  allTracks.addTracks(L[0]->allTracks());

  L[1]->findTracklets(L[2]);
  L[1]->findMatches(L[0]);
  L[1]->findMatches(L[3]);
  L[1]->findMatches(L[4]);
  L[1]->findMatches(L[5]);
  L[1]->findMatches(F[0]);
  L[1]->findMatches(F[1]);
  L[1]->findMatches(F[2]);
  L[1]->findMatches(F[3]);
  L[1]->findMatches(F[4]);
  L[1]->findMatches(F[5]);
  L[1]->findMatches(F[6]);
  L[1]->findMatches(B[0]);
  L[1]->findMatches(B[1]);
  L[1]->findMatches(B[2]);
  L[1]->findMatches(B[3]);
  L[1]->findMatches(B[4]);
  L[1]->findMatches(B[5]);
  L[1]->findMatches(B[6]);
  L[1]->fitTracks();

  allTracks.addTracks(L[1]->allTracks());


  F[0]->findTracklets(F[1]);
  F[0]->findMatches(F[2]);
  F[0]->findMatches(F[3]);
  F[0]->findMatches(F[4]);
  F[0]->findMatches(F[5]);
  F[0]->findMatches(F[6]);
  F[0]->findBarrelMatches(L[0]);
  F[0]->findBarrelMatches(L[1]);
  F[0]->findBarrelMatches(L[2]);
  F[0]->findBarrelMatches(L[3]);
  F[0]->findBarrelMatches(L[4]);
  F[0]->findBarrelMatches(L[5]);
  F[0]->fitTracks();

  allTracks.addTracks(F[0]->allTracks());

  F[1]->findTracklets(F[2]);
  F[1]->findMatches(F[0]);
  F[1]->findMatches(F[3]);
  F[1]->findMatches(F[4]);
  F[1]->findMatches(F[5]);
  F[1]->findMatches(F[6]);
  F[1]->findBarrelMatches(L[0]);
  F[1]->findBarrelMatches(L[1]);
  F[1]->findBarrelMatches(L[2]);
  F[1]->findBarrelMatches(L[3]);
  F[1]->findBarrelMatches(L[4]);
  F[1]->findBarrelMatches(L[5]);
  F[1]->fitTracks();

  allTracks.addTracks(F[1]->allTracks());

  F[2]->findTracklets(F[3]);
  F[2]->findMatches(F[0]);
  F[2]->findMatches(F[1]);
  F[2]->findMatches(F[4]);
  F[2]->findMatches(F[5]);
  F[2]->findMatches(F[6]);
  F[2]->findBarrelMatches(L[0]);
  F[2]->findBarrelMatches(L[1]);
  F[2]->findBarrelMatches(L[2]);
  F[2]->findBarrelMatches(L[3]);
  F[2]->findBarrelMatches(L[4]);
  F[2]->findBarrelMatches(L[5]);
  F[2]->fitTracks();

  allTracks.addTracks(F[2]->allTracks());

  F[3]->findTracklets(F[4]);
  F[3]->findMatches(F[0]);
  F[3]->findMatches(F[1]);
  F[3]->findMatches(F[2]);
  F[3]->findMatches(F[5]);
  F[3]->findMatches(F[6]);
  F[3]->findBarrelMatches(L[0]);
  F[3]->findBarrelMatches(L[1]);
  F[3]->findBarrelMatches(L[2]);
  F[3]->findBarrelMatches(L[3]);
  F[3]->findBarrelMatches(L[4]);
  F[3]->findBarrelMatches(L[5]);
  F[3]->fitTracks();

  allTracks.addTracks(F[3]->allTracks());

  F[4]->findTracklets(F[5]);
  F[4]->findMatches(F[0]);
  F[4]->findMatches(F[1]);
  F[4]->findMatches(F[2]);
  F[4]->findMatches(F[3]);
  F[4]->findMatches(F[6]);
  F[4]->findBarrelMatches(L[0]);
  F[4]->findBarrelMatches(L[1]);
  F[4]->findBarrelMatches(L[2]);
  F[4]->findBarrelMatches(L[3]);
  F[4]->findBarrelMatches(L[4]);
  F[4]->findBarrelMatches(L[5]);
  F[4]->fitTracks();

  allTracks.addTracks(F[4]->allTracks());

  F[5]->findTracklets(F[6]);
  F[5]->findMatches(F[0]);
  F[5]->findMatches(F[1]);
  F[5]->findMatches(F[2]);
  F[5]->findMatches(F[3]);
  F[5]->findMatches(F[4]);
  F[5]->findBarrelMatches(L[0]);
  F[5]->findBarrelMatches(L[1]);
  F[5]->findBarrelMatches(L[2]);
  F[5]->findBarrelMatches(L[3]);
  F[5]->findBarrelMatches(L[4]);
  F[5]->findBarrelMatches(L[5]);
  F[5]->fitTracks();

  allTracks.addTracks(F[5]->allTracks());

  B[0]->findTracklets(B[1]);
  B[0]->findMatches(B[2]);
  B[0]->findMatches(B[3]);
  B[0]->findMatches(B[4]);
  B[0]->findMatches(B[5]);
  B[0]->findMatches(B[6]);
  B[0]->findBarrelMatches(L[0]);
  B[0]->findBarrelMatches(L[1]);
  B[0]->findBarrelMatches(L[2]);
  B[0]->findBarrelMatches(L[3]);
  B[0]->findBarrelMatches(L[4]);
  B[0]->findBarrelMatches(L[5]);
  B[0]->fitTracks();

  allTracks.addTracks(B[0]->allTracks());

  B[1]->findTracklets(B[2]);
  B[1]->findMatches(B[0]);
  B[1]->findMatches(B[3]);
  B[1]->findMatches(B[4]);
  B[1]->findMatches(B[5]);
  B[1]->findMatches(B[6]);
  B[1]->findBarrelMatches(L[0]);
  B[1]->findBarrelMatches(L[1]);
  B[1]->findBarrelMatches(L[2]);
  B[1]->findBarrelMatches(L[3]);
  B[1]->findBarrelMatches(L[4]);
  B[1]->findBarrelMatches(L[5]);
  B[1]->fitTracks();

  allTracks.addTracks(B[1]->allTracks());

  B[2]->findTracklets(B[3]);
  B[2]->findMatches(B[0]);
  B[2]->findMatches(B[1]);
  B[2]->findMatches(B[4]);
  B[2]->findMatches(B[5]);
  B[2]->findMatches(B[6]);
  B[2]->findBarrelMatches(L[0]);
  B[2]->findBarrelMatches(L[1]);
  B[2]->findBarrelMatches(L[2]);
  B[2]->findBarrelMatches(L[3]);
  B[2]->findBarrelMatches(L[4]);
  B[2]->findBarrelMatches(L[5]);
  B[2]->fitTracks();

  allTracks.addTracks(B[2]->allTracks());

  B[3]->findTracklets(B[4]);
  B[3]->findMatches(B[0]);
  B[3]->findMatches(B[1]);
  B[3]->findMatches(B[2]);
  B[3]->findMatches(B[5]);
  B[3]->findMatches(B[6]);
  B[3]->findBarrelMatches(L[0]);
  B[3]->findBarrelMatches(L[1]);
  B[3]->findBarrelMatches(L[2]);
  B[3]->findBarrelMatches(L[3]);
  B[3]->findBarrelMatches(L[4]);
  B[3]->findBarrelMatches(L[5]);
  B[3]->fitTracks();

  allTracks.addTracks(B[3]->allTracks());

  B[4]->findTracklets(B[5]);
  B[4]->findMatches(B[0]);
  B[4]->findMatches(B[1]);
  B[4]->findMatches(B[2]);
  B[4]->findMatches(B[3]);
  B[4]->findMatches(B[6]);
  B[4]->findBarrelMatches(L[0]);
  B[4]->findBarrelMatches(L[1]);
  B[4]->findBarrelMatches(L[2]);
  B[4]->findBarrelMatches(L[3]);
  B[4]->findBarrelMatches(L[4]);
  B[4]->findBarrelMatches(L[5]);
  B[4]->fitTracks();

  allTracks.addTracks(B[4]->allTracks());

  B[5]->findTracklets(B[6]);
  B[5]->findMatches(B[0]);
  B[5]->findMatches(B[1]);
  B[5]->findMatches(B[2]);
  B[5]->findMatches(B[3]);
  B[5]->findMatches(B[4]);
  B[5]->findBarrelMatches(L[0]);
  B[5]->findBarrelMatches(L[1]);
  B[5]->findBarrelMatches(L[2]);
  B[5]->findBarrelMatches(L[3]);
  B[5]->findBarrelMatches(L[4]);
  B[5]->findBarrelMatches(L[5]);
  B[5]->fitTracks();

  allTracks.addTracks(B[5]->allTracks());


 
 
    
  L1TTracks purgedTracks=allTracks.purged();
    
  cout << "allTracks purgedTracks :"<<allTracks.size()<<" "
       <<purgedTracks.size()<<endl;

  for (unsigned int i=0;i<purgedTracks.size();i++){
    const L1TTrack& aTrack=purgedTracks.get(i);
    double frac=0.0;
    int simtrackid=aTrack.simtrackid(frac);
    int isimtrack=ev.getSimtrackFromSimtrackid(simtrackid);
    if (isimtrack!=-1) {
      L1SimTrack simtrack=ev.simtrack(isimtrack);
      //cout << "Track "<<i<<" with simtrackid="<<simtrackid<<" pt="<<simtrack.pt()<<" phi="<<simtrack.phi()<<" eta="<<simtrack.eta()<<endl;
      //cout << "Track pt="<<aTrack.pt(mMagneticFieldStrength)<<" phi="<<aTrack.phi0()<<" eta="<<aTrack.eta()<<endl;
      if (1) {
	static ofstream out("tracks.txt");
	out <<simtrack.pt()<<" "<<simtrack.phi()<<" "
	    <<simtrack.eta()<<" "<<simtrack.vz()<<" "
	    <<aTrack.pt(mMagneticFieldStrength)<<" "<<aTrack.phi0()<<" "
	    <<aTrack.eta()<<" "<<aTrack.z0()<<" "
	    <<aTrack.ptseed(mMagneticFieldStrength)<<" "<<aTrack.phi0seed()<<" "
	    <<aTrack.etaseed()<<" "<<aTrack.z0seed()<<" "
	    <<aTrack.chisqdof()<<endl;
      }
    }
    else {
      cout << "Track "<<i<<" had no matched simtrack"<<endl;
    }
  }


  if (1) {

    static ofstream out("trackeff.txt");
      
    int nsim=ev.nsimtracks();

    //cout << "Number of sim tracks="<<nsim<<endl;

    for (int jj=0;jj<nsim;jj++){
      
      L1SimTrack aSimTrack=ev.simtrack(jj);
      
      int simtrackid=aSimTrack.id();

      double r=sqrt(aSimTrack.vx()*aSimTrack.vx()+
		    aSimTrack.vy()*aSimTrack.vy());

      if (simtrackid>10000) continue; //hack
      
      //cout << "SimTrackID="<<aSimTrack.id()<<" "<<aSimTrack.type()<<endl;

      bool match=false;

      for (unsigned itrack=0;itrack<purgedTracks.size();itrack++) {
	L1TTrack track=purgedTracks.get(itrack);
	  
	double frac;
	int simtrackidmatch=track.simtrackid(frac);

	//cout << "simtrackidmatch frac:"<<simtrackidmatch<<" "<<frac<<endl;

	if (simtrackidmatch==simtrackid&&frac>0.7) {
	  //cout << "Found Track"<<endl;
	  match=true;
	}

      }
      
      out << aSimTrack.pt() << " " << aSimTrack.eta() 
	  << " " << aSimTrack.phi() << " " 
	  << aSimTrack.type() << " " << r << " ";
      
      if (match) {
	out <<"1"<<endl;
      }
      else {
	out <<"0"<<endl;  
      }

    }

  }

  
  for (unsigned itrack=0; itrack<purgedTracks.size(); itrack++) {
    L1TTrack track=purgedTracks.get(itrack);

    //L1TkTrackType TkTrack(TkStubs, aSeedTracklet);
    L1TkTrackType TkTrack;
    //double frac;
    //TkTrack.setSimTrackId(track.simtrackid(frac));  FIXME
    //TkTrack.setRadius(1./track.rinv());  FIXME
    //GlobalPoint bsPosition(recoBeamSpotHandle->position().x(),
    //			   recoBeamSpotHandle->position().y(),
    //			   track.z0()
    //			   ); //store the L1 track vertex position 
    GlobalPoint bsPosition(0.0,
			   0.0,
			   track.z0()
			   ); //store the L1 track vertex position 
    //TkTrack.setVertex(bsPosition);  FIXME
    //TkTrack.setChi2RPhi(track.chisq1()); FIXME
    //TkTrack.setChi2ZPhi(track.chisq2()); FIXME
    cout << "L1TrackProducer::analyze Track with pt="<<track.pt(mMagneticFieldStrength)<<endl;
    TkTrack.setMomentum( GlobalVector ( GlobalVector::Cylindrical(fabs(track.pt(mMagneticFieldStrength)), 
								  track.phi0(), 
								  fabs(track.pt(mMagneticFieldStrength))*sinh(track.eta())) ) );

    L1TkTracksForOutput->push_back(TkTrack);

    vector<L1TkStubPtrType> TkStubs;
    L1TTracklet tracklet = track.getSeed();
    vector<L1TStub> stubComponents;// = tracklet.getStubComponents();
    vector<L1TStub> stubs = track.getStubs();
    //L1TkTrackletType TkTracklet;

    stubMapType::iterator it;
    //for (it = stubMap.begin(); it != stubMap.end(); it++) {
      //if (it->first == stubComponents[0] || it->first == stubComponents[1]) {
      //L1TkStubPtrType TkStub = it->second;
	//if (TkStub->getStack()%2 == 0)
	//  TkTracklet.addStub(0, TkStub);
	//else
	//  TkTracklet.addStub(1, TkStub);
      //}
      
      //for (int j=0; j<(int)stubs.size(); j++) {
    //	if (it->first == stubs[j])
    //  TkStubs.push_back(it->second);
    //}
    //}

    L1TkStubsForOutput->push_back( TkStubs );
    //TkTracklet.checkSimTrack();
    //TkTracklet.fitTracklet(mMagneticFieldStrength, GlobalPoint(bsPosition.x(), bsPosition.y(), 0.0), true);
    //L1TkTrackletsForOutput->push_back( TkTracklet );
  }



  // }

  iEvent.put( L1TkStubsForOutput, "L1TkStubs");
  //iEvent.put( L1TkTrackletsForOutput, "L1TkTracklets" );
  iEvent.put( L1TkTracksForOutput, "Level1TkTracks");

  //for (unsigned int j=0;j<NSECTORS;j++){
  //  delete Sectors[j];
  //}


} /// End of produce()


// ///////////////////////////
// // DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(L1TrackProducer);

#endif
