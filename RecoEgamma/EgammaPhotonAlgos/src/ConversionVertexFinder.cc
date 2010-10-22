#include <iostream>
#include <vector>
#include <memory>
#include "RecoEgamma/EgammaPhotonAlgos/interface/ConversionVertexFinder.h"
#include "RecoEgamma/EgammaPhotonAlgos/interface/TangentApproachInRPhi.h"
#include "TrackingTools/GeomPropagators/interface/TrackerBounds.h"
// Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
//Kinematic constraint vertex fitter
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>

#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/ColinearityKinematicConstraint.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

// new templated one
#ifdef TemplateKineFit
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitterT.h"
#include "RecoVertex/KinematicFit/interface/ColinearityKinematicConstraintT.h"
#endif


//
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include <TMath.h>





ConversionVertexFinder::ConversionVertexFinder(const edm::ParameterSet& config ):
  conf_(config)
{ 
  LogDebug("ConversionVertexFinder") << "ConversionVertexFinder CTOR  " <<  "\n";  
  maxDelta_ = conf_.getParameter<double>("maxDelta");
  maxReducedChiSq_ = conf_.getParameter<double>("maxReducedChiSq");
  minChiSqImprovement_  = conf_.getParameter<double>("minChiSqImprovement");
  maxNbrOfIterations_ = conf_.getParameter<int>("maxNbrOfIterations");
  kcvFitter_ = new KinematicConstrainedVertexFitter();
  kcvFitter_->setParameters(conf_);

}

ConversionVertexFinder::~ConversionVertexFinder() {

  LogDebug("ConversionVertexFinder") << "ConversionVertexFinder DTOR " <<  "\n";  
  delete   kcvFitter_;
 
}

bool  ConversionVertexFinder::run(std::vector<reco::TransientTrack>  pair, reco::Vertex& the_vertex) {
  bool found= false;

  if ( pair.size() < 2) return found;
  
  TangentApproachInRPhi tangent;
  tangent.calculate(pair[0].innermostMeasurementState(),pair[1].innermostMeasurementState());
  if (!tangent.status()) {
    return false;
  }
  
  GlobalPoint tangentPoint = tangent.crossingPoint();
  if (!TrackerBounds::isInside(tangentPoint)) {
    return false;
  }
  
  float sigma = 0.00000001;
  float chi = 0.;
  float ndf = 0.;
  float mass = 0.000511;
  

  KinematicParticleFactoryFromTransientTrack pFactory;
  
  std::vector<RefCountedKinematicParticle> particles;
  
  particles.push_back(pFactory.particle (pair[0],mass,chi,ndf,sigma,*pair[0].innermostMeasurementState().freeState()));
  particles.push_back(pFactory.particle (pair[1],mass,chi,ndf,sigma,*pair[1].innermostMeasurementState().freeState()));
  

#ifndef TemplateKineFit
  ColinearityKinematicConstraint constr(ColinearityKinematicConstraint::PhiTheta);
  
  RefCountedKinematicTree myTree = kcvFitter_->fit(particles, &constr, &tangentPoint);
#else

  // bizzare way to the get the field...
  const MagneticField* mf = pair[0].field();

  ColinearityKinematicConstraintT<colinearityKinematic::PhiTheta> constr;
  KinematicConstrainedVertexFitterT<2,2> kcvFitter(mf);
  kcvFitter.setParameters(conf_);
  RefCountedKinematicTree myTree =  kcvFitter.fit(particles, &constr, &tangentPoint);

#ifdef TemplateKineFitDebug

  ColinearityKinematicConstraint oldconstr(ColinearityKinematicConstraint::PhiTheta);
  
  RefCountedKinematicTree oldTree = kcvFitter_->fit(particles, &oldconstr, &tangentPoint);


  if( oldTree->isValid() ) {
    std::cout << "old" << std::endl;
    std::cout <<  oldTree->currentParticle()->currentState().globalMomentum() <<
      oldTree->currentParticle()->currentState().globalPosition()<< std::endl;
    std::vector<RefCountedKinematicParticle> fStates=oldTree->finalStateParticles();
    for (unsigned int kk=0; kk<fStates.size(); kk++) {
      std::cout <<  fStates[kk]->currentState().globalMomentum() << 
	fStates[kk]->currentState().globalPosition() << std::endl;
    }
  } else       std::cout << "old invalid" << std::endl;
  
  if( myTree->isValid() ) {
    std::cout << "new" << std::endl;
    std::cout <<  myTree->currentParticle()->currentState().globalMomentum() <<
      myTree->currentParticle()->currentState().globalPosition()<< std::endl;
    std::vector<RefCountedKinematicParticle> fStates=myTree->finalStateParticles();
    for (unsigned int kk=0; kk<fStates.size(); kk++) {
      std::cout <<  fStates[kk]->currentState().globalMomentum() << 
	fStates[kk]->currentState().globalPosition() << std::endl;
    }
  } else       std::cout << "new invalid" << std::endl;

#endif // TemplateKineFitDebug

#endif

  if( myTree->isValid() ) {
    myTree->movePointerToTheTop();                                                                                
    RefCountedKinematicParticle the_photon = myTree->currentParticle();                                           
    if (the_photon->currentState().isValid()){                                                                    
      //const ParticleMass photon_mass = the_photon->currentState().mass();                                       
      RefCountedKinematicVertex gamma_dec_vertex;                                                               
      gamma_dec_vertex = myTree->currentDecayVertex();                                                          
      if( gamma_dec_vertex->vertexIsValid() ){                                                                  
	const float chi2Prob = ChiSquaredProbability(gamma_dec_vertex->chiSquared(), gamma_dec_vertex->degreesOfFreedom());
	if (chi2Prob>0.){// no longer cut here, only ask positive probability here 
	  //const math::XYZPoint vtxPos(gamma_dec_vertex->position());                                           
	  the_vertex = *gamma_dec_vertex;
	  found = true;
	}
      }
    }
  }
 
  return found;
}


TransientVertex  ConversionVertexFinder::run(std::vector<reco::TransientTrack>  pair) {
  LogDebug("ConversionVertexFinder") << "ConversionVertexFinder run pair size " << pair.size() <<  "\n";  
  
  //for ( std::vector<reco::TransientTrack>::const_iterator iTk=pair.begin(); iTk!=pair.end(); ++iTk) {
  // LogDebug("ConversionVertexFinder") << "  ConversionVertexFinder  Tracks in the pair  charge " << iTk->charge() << " Num of RecHits " << iTk->recHitsSize() << " inner momentum " << iTk->track().innerMomentum() << "\n";  
  //}


  reco::Vertex theVertex;  
  KalmanVertexFitter fitter(true);
  TransientVertex transientVtx;

  const std::string metname =  "ConversionVertexFinder| ConversionVertexFinder";
  try{

    transientVtx = fitter.vertex(pair); 

  }  catch ( cms::Exception& e ) {


    edm::LogWarning(metname) << "cms::Exception caught in ConversionVertexFinder::run\n"
			     << e.explainSelf();
    
  }
  

  return transientVtx;

    
    
}

