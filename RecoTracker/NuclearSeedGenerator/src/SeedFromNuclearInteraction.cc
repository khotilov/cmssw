#include "RecoTracker/NuclearSeedGenerator/interface/SeedFromNuclearInteraction.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

SeedFromNuclearInteraction::SeedFromNuclearInteraction(const edm::EventSetup& es, const edm::ParameterSet& iConfig) : 
rescaleDirectionFactor(iConfig.getParameter<double>("rescaleDirectionFactor")),
rescalePositionFactor(iConfig.getParameter<double>("rescalePositionFactor")),
rescaleCurvatureFactor(iConfig.getParameter<double>("rescaleCurvatureFactor")) {

  edm::ESHandle<Propagator>  thePropagatorHandle;
  es.get<TrackingComponentsRecord>().get("PropagatorWithMaterial",thePropagatorHandle);
  thePropagator = &(*thePropagatorHandle);
  isValid_=true;

  es.get<TrackerDigiGeometryRecord> ().get (pDD);
}

//----------------------------------------------------------------------
void SeedFromNuclearInteraction::setMeasurements(const TM& inner_TM, const TM& outer_TM) {

       // delete pointer to TrackingRecHits
       theHits.clear();

       // get the inner and outer transient TrackingRecHits
       outerHit = outer_TM.recHit().get();

       theHits.push_back(  inner_TM.recHit() ); // put temporarily - TODO: remove this line
       theHits.push_back(  outer_TM.recHit() );

       innerTM = &inner_TM;
       outerTM = &outer_TM;

       // calculate the initial FreeTrajectoryState from the inner and outer TM.
       freeTS  = stateWithError();

       // convert freeTS into a persistent TSOS on the outer surface
       isValid_ = construct();
}
//----------------------------------------------------------------------
void SeedFromNuclearInteraction::setMeasurements(TangentHelix& helix, const TM& inner_TM, const TM& outer_TM) {

       // delete pointer to TrackingRecHits
       theHits.clear();

       // get the inner and outer transient TrackingRecHits
       outerHit = outer_TM.recHit().get();

       theHits.push_back( inner_TM.recHit() );
       theHits.push_back( outer_TM.recHit() );

       innerTM = &inner_TM;
       outerTM = &outer_TM;

       // calculate the initial FreeTrajectoryState from the inner and outer TM assuming that the helix equation is already known.
       freeTS = stateWithError(helix); 

       // convert freeTS into a persistent TSOS on the outer surface
       isValid_ = construct();
}

//----------------------------------------------------------------------
FreeTrajectoryState SeedFromNuclearInteraction::stateWithError() const {

   // Calculation of the helix assuming that the secondary track has the same direction
   // than the primary track and pass through the inner and outer hits.
   GlobalVector direction = innerTM->updatedState().globalDirection();
   GlobalPoint inner = innerTM->updatedState().globalPosition();
   GlobalPoint outer = pDD->idToDet(outerHit->geographicalId())->surface().toGlobal(outerHit->localPosition());
   TangentHelix helix(direction, inner, outer);
   LogDebug("NuclearSeedGenerator") << "First vtx position : " << helix.vertexPoint() << "\n"
                                    << "Rho = " << helix.rho() << "\n";

   return stateWithError(helix);
}
//----------------------------------------------------------------------
FreeTrajectoryState SeedFromNuclearInteraction::stateWithError(TangentHelix& helix) const {

//   typedef TkRotation<float> Rotation;

   GlobalVector dirAtVtx = helix.directionAtVertex();
   const MagneticField& mag = innerTM->updatedState().globalParameters().magneticField();
   // Get the global parameters of the trajectory
   // we assume that the magnetic field at the vertex is equal to the magnetic field at the inner TM.
   GlobalTrajectoryParameters gtp(helix.vertexPoint(), dirAtVtx , helix.charge(mag.inTesla(helix.vertexPoint()).z())/helix.rho(), 0, &mag);
   LogDebug("NuclearSeedGenerator") << "Momentum = " << gtp.momentum() << "\n"
                                    << "Charge = " << helix.charge(mag.inTesla(helix.vertexPoint()).z()) << "\n"
                                    << "Inner = " << helix.innerPoint() << "\n"
                                    << "Outer = " << helix.outerPoint() << "\n";

   // Error matrix in a frame where z is in the direction of the track at the vertex
   AlgebraicSymMatrix55 m = ROOT::Math::SMatrixIdentity();
   double vtxerror = helix.circle().vertexError();
   double curvatureError = helix.curvatureError();
   m(0,0)=curvatureError*curvatureError;
   m(1,1)=1E-5;
   m(2,2)=1E-5;
   m(3,3)=1E-5;
   m(4,4)=1E-4;

   LogDebug("NuclearSeedGenerator") << "vtxError : " << vtxerror << "\n";

/*
   //rotation around the z-axis by  -phi
   Rotation tmpRotz ( cos(dirAtVtx.phi()), -sin(dirAtVtx.phi()), 0., 
                        sin(dirAtVtx.phi()), cos(dirAtVtx.phi()), 0.,
                         0.,              0.,              1. );

   //rotation around y-axis by -theta
   Rotation tmpRoty ( cos(dirAtVtx.theta()), 0.,sin(dirAtVtx.theta()),
                               0.,              1.,              0.,
                              -sin(dirAtVtx.theta()), 0., cos(dirAtVtx.theta()) );

   Rotation position(m(0,0), 0, 0, 0, m(1,1), 0, 0, 0, m(2,2) );
   Rotation momentum(m(3,3), 0, 0, 0, m(4,4), 0, 0, 0, m(5,5) ); 

   // position = position * tmpRoty * tmpRotz
   // momentum = momentum * tmpRoty * tmpRotz
   position *= tmpRoty;   momentum *= tmpRoty; 
   position *= tmpRotz;   momentum *= tmpRotz; 

   m(0,0) = position.xx();
   m(1,0) = position.yx();
   m(2,0) = position.zx();
   m(0,1) = position.xy();
   m(1,1) = position.yy();
   m(2,1) = position.zy();
   m(0,2) = position.xz();
   m(1,2) = position.yz();
   m(2,2) = position.zz();
   m(3,3) = momentum.xx();
   m(4,3) = momentum.yx();
   m(5,3) = momentum.zx();
   m(3,4) = momentum.xy();
   m(4,4) = momentum.yy();
   m(5,4) = momentum.zy();
   m(3,5) = momentum.xz();
   m(4,5) = momentum.yz();
   m(5,5) = momentum.zz();
 */  

   FreeTrajectoryState result( gtp, CurvilinearTrajectoryError(m) );

   LogDebug("NuclearSeedGenerator") << "FreeTrajectoryState used for seeds : " << result << "\n";

   return result;
}

//----------------------------------------------------------------------
bool SeedFromNuclearInteraction::construct() {

   // loop on all hits in theHits
   KFUpdator                 theUpdator;
   TrajectoryStateOnSurface  updatedTSOS;

   const TrackingRecHit* hit = 0;

   for ( unsigned int iHit = 0; iHit < theHits.size(); iHit++) {
     hit = theHits[iHit]->hit();
     TrajectoryStateOnSurface state = (iHit==0) ? 
        thePropagator->propagate(freeTS,pDD->idToDet(hit->geographicalId())->surface())
       : thePropagator->propagate(updatedTSOS, pDD->idToDet(hit->geographicalId())->surface());

     if (!state.isValid()) return false; 
 
     const TransientTrackingRecHit::ConstRecHitPointer& tth = theHits[iHit]; 
     updatedTSOS =  theUpdator.update(state, *tth);

     LogDebug("NuclearSeedGenerator") << "TrajectorySeed updated with hit at position : " << pDD->idToDet(hit->geographicalId())->surface().toGlobal(hit->localPosition()) << "\n"
                                      << "state = " << updatedTSOS << "\n";
       
   } 

   TrajectoryStateTransform transformer;

   pTraj = boost::shared_ptr<PTrajectoryStateOnDet>( transformer.persistentState(updatedTSOS, outerHit->geographicalId().rawId()) );
   return true;
}

//----------------------------------------------------------------------
edm::OwnVector<TrackingRecHit>  SeedFromNuclearInteraction::hits() const { 
    recHitContainer      _hits;
    for( ConstRecHitContainer::const_iterator it = theHits.begin(); it!=theHits.end(); it++ ){
           _hits.push_back( it->get()->hit()->clone() );
           LogDebug("NuclearSeedGenerator") << "Hit put in TrajectorySeed at position : " << pDD->idToDet(it->get()->hit()->geographicalId())->surface().toGlobal(it->get()->hit()->localPosition()) << "\n";
    }
    return _hits; 
}
