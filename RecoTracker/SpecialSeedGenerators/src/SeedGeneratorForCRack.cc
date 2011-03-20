#include "RecoTracker/SpecialSeedGenerators/interface/SeedGeneratorForCRack.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoTracker/TkHitPairs/interface/CosmicLayerPairs.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "RecoTracker/TkSeedGenerator/interface/SeedFromConsecutiveHits.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 
void 
SeedGeneratorForCRack::init(const SiStripRecHit2DCollection &collstereo,
			      const SiStripRecHit2DCollection &collrphi ,
			      const SiStripMatchedRecHit2DCollection &collmatched,
			      const edm::EventSetup& iSetup)
{
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
  thePropagatorAl=    new PropagatorWithMaterial(alongMomentum,0.1057,&(*magfield) );
  thePropagatorOp=    new PropagatorWithMaterial(oppositeToMomentum,0.1057,&(*magfield) );
  theUpdator=         new KFUpdator();
  
  // get the transient builder

  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;

  iSetup.get<TransientRecHitRecord>().get(builderName,theBuilder);
  TTTRHBuilder = theBuilder.product();
  CosmicLayerPairs cosmiclayers(geometry);
  cosmiclayers.init(collstereo,collrphi,collmatched,iSetup);
  thePairGenerator=new CosmicHitPairGenerator(cosmiclayers,iSetup);
  HitPairs.clear();
  thePairGenerator->hitPairs(region,HitPairs,iSetup);
  LogDebug("CosmicSeedFinder") <<"Initialized with " << HitPairs.size() << " hit pairs" << std::endl;
}

SeedGeneratorForCRack::SeedGeneratorForCRack(edm::ParameterSet const& conf):
  conf_(conf)
{  
  float ptmin=conf_.getParameter<double>("ptMin");
  float originradius=conf_.getParameter<double>("originRadius");
  float halflength=conf_.getParameter<double>("originHalfLength");
  float originz=conf_.getParameter<double>("originZPosition");
  seedpt = conf_.getParameter<double>("SeedPt");

  builderName = conf_.getParameter<std::string>("TTRHBuilder");   
  geometry=conf_.getUntrackedParameter<std::string>("GeometricStructure","STANDARD");
  region=GlobalTrackingRegion(ptmin,originradius,
 			      halflength,originz);
  multipleScatteringFactor=conf_.getUntrackedParameter<double>("multipleScatteringFactor", 1.0);
  seedMomentum =conf_.getUntrackedParameter<double>("SeedMomentum",1);
  edm::LogInfo("SeedGeneratorForCRack")<<" PtMin of track is "<<ptmin<< 
    " The Radius of the cylinder for seeds is "<<originradius <<"cm"  << " The set Seed Momentum" <<  seedpt;
}

void SeedGeneratorForCRack::run(TrajectorySeedCollection &output,const edm::EventSetup& iSetup){
  seeds(output,iSetup,region);
  delete thePairGenerator;
  delete thePropagatorAl;
  delete thePropagatorOp;
  delete theUpdator; 
}

void SeedGeneratorForCRack::seeds(TrajectorySeedCollection &output,
				    const edm::EventSetup& iSetup,
				    const TrackingRegion& region){
  for(unsigned int is=0;is<HitPairs.size();is++){

    GlobalPoint inner = tracker->idToDet((*(HitPairs[is].inner())).geographicalId())->surface().toGlobal((*(HitPairs[is].inner())).localPosition());
    GlobalPoint outer = tracker->idToDet((*(HitPairs[is].outer())).geographicalId())->surface().toGlobal((*(HitPairs[is].outer())).localPosition());
    
    LogDebug("CosmicSeedFinder") <<"inner point of the seed "<<inner <<" outer point of the seed "<<outer; 
    TransientTrackingRecHit::ConstRecHitPointer inrhit=TTTRHBuilder->build(HitPairs[is].inner()->hit());
    TransientTrackingRecHit::ConstRecHitPointer outrhit = TTTRHBuilder->build(HitPairs[is].outer()->hit());

    edm::OwnVector<TrackingRecHit> hits;
    hits.push_back(HitPairs[is].outer()->hit()->clone());

    for (int i=0;i<2;i++){
      //FIRST STATE IS CALCULATED CONSIDERING THAT THE CHARGE CAN BE POSITIVE OR NEGATIVE
      int predsign=(2*i)-1;
      if((outer.y()-inner.y())>0){
	GlobalVector momentum = GlobalVector(inner-outer);
	momentum = momentum.unit()*seedMomentum;
	GlobalTrajectoryParameters Gtp(inner,
				       momentum,
				       predsign, 
				       &(*magfield));
	AlgebraicSymMatrix55 errMatrix = ROOT::Math::SMatrixIdentity();
	TSOS innerState = TSOS(Gtp, CurvilinearTrajectoryError(errMatrix), tracker->idToDet((HitPairs[is].inner()->hit())->geographicalId())->surface());
	const TSOS innerUpdated = theUpdator->update(innerState, *inrhit);
	//Cosmic Seed update inner...
	LogDebug("CosmicSeedFinder") << " FirstTSOS " << innerUpdated;
	
	//First propagation
	const TSOS outerState =
	  thePropagatorOp->propagate(innerUpdated,
				     tracker->idToDet((*(HitPairs[is].outer())).geographicalId())->surface());
	if ( outerState.isValid()) {
	  LogDebug("CosmicSeedFinder") <<"outerState "<<outerState;
	  TSOS outerUpdated= theUpdator->update( outerState,*outrhit);
	  //fudge factor for multiple scattering
	  outerUpdated.rescaleError( multipleScatteringFactor);
	  if ( outerUpdated.isValid()) {
	    LogDebug("CosmicSeedFinder") <<"outerUpdated "<<outerUpdated;
	    
	    PTrajectoryStateOnDet *PTraj=  
	      transformer.persistentState(outerUpdated, (*(HitPairs[is].outer())).geographicalId().rawId());
	    output.push_back(TrajectorySeed(*PTraj,hits,alongMomentum));
            delete PTraj;
	    
	  }else      edm::LogWarning("CosmicSeedFinder") << " SeedForCosmics first update failed ";
	}else      edm::LogWarning("CosmicSeedFinder") << " SeedForCosmics first propagation failed ";
      
      
      }
      else{
        GlobalVector momentum= GlobalVector(outer-inner);
	momentum=momentum.unit()*seedMomentum;
	GlobalTrajectoryParameters Gtp(inner,
				       momentum,
				       predsign, 
				       &(*magfield));
	AlgebraicSymMatrix55 errMatrix = ROOT::Math::SMatrixIdentity();
	TSOS innerState = TSOS(Gtp, CurvilinearTrajectoryError(errMatrix), tracker->idToDet((HitPairs[is].inner()->hit())->geographicalId())->surface());
	const TSOS innerUpdated = theUpdator->update(innerState, *inrhit);
	LogDebug("CosmicSeedFinder") << " FirstTSOS "<< innerState;
	
	//First propagation
	const TSOS outerState =
	  thePropagatorOp->propagate(innerUpdated,
				     tracker->idToDet((*(HitPairs[is].outer())).geographicalId())->surface());
	if ( outerState.isValid()) {
	  
	  LogDebug("CosmicSeedFinder") <<"outerState "<<outerState;
	  TSOS outerUpdated = theUpdator->update( outerState,*outrhit);
	  //fudge factor for multiple scattering
	  outerUpdated.rescaleError( multipleScatteringFactor);
	  if ( outerUpdated.isValid()) {
	  LogDebug("CosmicSeedFinder") <<"outerUpdated "<<outerUpdated;
	  PTrajectoryStateOnDet *PTraj=  
	    transformer.persistentState(outerUpdated,(*(HitPairs[is].outer())).geographicalId().rawId());
          output.push_back(TrajectorySeed(*PTraj,hits,oppositeToMomentum));
          delete PTraj;
	
	  }else      edm::LogWarning("CosmicSeedFinder") << " SeedForCosmics first update failed ";
	}else      edm::LogWarning("CosmicSeedFinder") << " SeedForCosmics first propagation failed ";
      }
    }
  }
}

