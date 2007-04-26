#ifndef CtfSpecialSeedGenerator_H
#define CtfSpecialSeedGenerator_H

/** \class CombinatorialSeedGeneratorForCOsmics
 *  A concrete seed generator providing seeds constructed 
 *  from combinations of hits in pairs of strip layers 
 */
//FWK
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
//DataFormats
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"    
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
//RecoTracker
#include "RecoTracker/TkSeedingLayers/interface/SeedingHit.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingLayerSetsBuilder.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h"
#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGenerator.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducer.h"
#include "RecoTracker/SpecialSeedGenerators/interface/SeedFromGenericPairOrTriplet.h"
//#include "RecoTracker/SpecialSeedGenerators/interface/GenericPairOrTripletGenerator.h"
//#include "RecoTracker/SpecialSeedGenerators/interface/SeedCleaner.h"
//MagneticField
#include "MagneticField/Engine/interface/MagneticField.h"
//TrackingTools
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/DetLayers/interface/NavigationDirection.h"
//Geometry
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include <map>

class CtfSpecialSeedGenerator : public edm::EDProducer
{
 public:
  typedef TrajectoryStateOnSurface TSOS;
  

  CtfSpecialSeedGenerator(const edm::ParameterSet& conf);

  virtual ~CtfSpecialSeedGenerator();//{};

  virtual void beginJob(const edm::EventSetup& c);	

  virtual void produce(edm::Event& e, const edm::EventSetup& c);

 private:
  
  void  run(const edm::EventSetup& c, 
	    const edm::Event& e, 
            TrajectorySeedCollection& output);

  void buildSeeds(const edm::EventSetup& iSetup,
                  const edm::Event& e,
		  const OrderedSeedingHits& osh,
		  const NavigationDirection& navdir,
                  const PropagationDirection& dir,
                  TrajectorySeedCollection& output);
  //checks that the hits used are at positive y and are on different layers
  bool preliminaryCheck(const SeedingHitSet& shs);
  //We can check if the seed  points in a region covered by scintillators. To be used only in noB case
  //because it uses StraightLinePropagation
  bool postCheck(const TrajectorySeed& seed);
 
 private:
  edm::ParameterSet conf_;
  edm::ESHandle<MagneticField> magfield;
  edm::ESHandle<TrackerGeometry> tracker;
  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
  //edm::ESHandle<SeedCleaner> theCleaner;
  OrderedHitsGenerator*  hitsGeneratorOutIn;
  OrderedHitsGenerator*  hitsGeneratorInOut;
  PropagationDirection inOutPropagationDirection;
  PropagationDirection outInPropagationDirection;
  //GenericPairOrTripletGenerator* hitsGeneratorOutIn;
  //GenericPairOrTripletGenerator* hitsGeneratorInOut;	
  TrackingRegionProducer* theRegionProducer;	
  TrajectoryStateTransform transformer;
  SeedFromGenericPairOrTriplet* theSeedBuilder; 
  bool useScintillatorsConstraint;
  BoundPlane* upperScintillator;
  BoundPlane* lowerScintillator;
};
#endif


