#ifndef CkfDebugger_H
#define CkfDebugger_H

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "TrackingTools/PatternTools/interface/MeasurementExtractor.h"

#include <vector>
#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

class Trajectory;
class TrajectoryMeasurement;
class PSimHit;
class TransientTrackingRecHit;
class MeasurementTracker;
class TrajectoryStateOnSurface;
class MagneticField;
class Chi2MeasurementEstimator;
class Propagator;

typedef TransientTrackingRecHit::ConstRecHitPointer CTTRHp;

using namespace std;

class CkfDebugger {
 public:
  CkfDebugger( edm::EventSetup const & es );

  ~CkfDebugger();
  
  void printSimHits( const edm::Event& iEvent);

  void countSeed(){totSeeds++;}

  void fillSeedHist(CTTRHp h1,CTTRHp h2, TrajectoryStateOnSurface t) {
    hchi2seedAll->Fill( testSeed(h1,h2,t) );
  }

  bool analyseCompatibleMeasurements( const Trajectory&,
				      const std::vector<TrajectoryMeasurement>&,
				      const MeasurementTracker*,
				      const Propagator*,
				      const Chi2MeasurementEstimatorBase*,
				      const TransientTrackingRecHitBuilder*);

  void deleteHitAssociator(){
    edm::LogVerbatim("CkfPattern") << "deleting hitAssociator " << hitAssociator ;
    delete hitAssociator;
  }

 private:
  typedef TrajectoryMeasurement        TM;
  typedef TrajectoryStateOnSurface     TSOS;

  class SimHit {
  public:

    SimHit( const PSimHit* phit, const GeomDetUnit* gdu) : thePHit( phit), theDet(gdu) {}
    LocalPoint localPosition() const {return thePHit->localPosition();}
    GlobalPoint globalPosition() const {return theDet->toGlobal( thePHit->localPosition());}
    const GeomDetUnit* det() const {return theDet;}
    unsigned int trackId()      const {return thePHit->trackId();}
    LocalVector  localDirection()  const {return thePHit->localDirection();}
    Geom::Theta<float> thetaAtEntry() const {return thePHit->thetaAtEntry();}
    Geom::Phi<float>   phiAtEntry()   const {return thePHit->phiAtEntry();}
    float        pabs()         const {return thePHit->pabs();}
    float        timeOfFlight() const {return thePHit->timeOfFlight();}
    float        energyLoss()   const {return thePHit->energyLoss();}
    int          particleType() const {return thePHit->particleType();}
    unsigned int detUnitId()    const {return thePHit->detUnitId();}
    unsigned short processType() const {return thePHit->processType();}
    const PSimHit& psimHit() const { return *thePHit;}

  private:

    const PSimHit*     thePHit;
    const GeomDetUnit* theDet;

  };

  const TrackerGeometry*           theTrackerGeom;
  const MagneticField*             theMagField;
  const GeometricSearchTracker*    theGeomSearchTracker;
  const MeasurementEstimator*  theChi2;
  const Propagator*                theForwardPropagator;
  TrackerHitAssociator*      hitAssociator;
  const MeasurementTracker*        theMeasurementTracker;
  const TransientTrackingRecHitBuilder* theTTRHBuilder;

  std::map<unsigned int, std::vector<PSimHit*> > idHitsMap;

  void dumpSimHit( const SimHit& hit) const;

  bool correctTrajectory( const Trajectory&, unsigned int&) const;

  int assocTrackId(CTTRHp rechit) const;

  //const PSimHit* nextCorrectHit( const Trajectory&, unsigned int&) ;
  vector<const PSimHit*> nextCorrectHits( const Trajectory&, unsigned int&) ;

  bool associated(CTTRHp rechit, const PSimHit& sh) const;

  bool goodSimHit(const PSimHit& sh) const;

  bool correctMeas( const TM& tm, const PSimHit* correctHit) const;

  std::pair<CTTRHp, double> analyseRecHitExistance( const PSimHit& sh, const TSOS& startingState);

  int analyseRecHitNotFound(const Trajectory&,CTTRHp);

  double testSeed(CTTRHp,CTTRHp, TrajectoryStateOnSurface);

  const PSimHit* pSimHit(unsigned int tkId, DetId detId);

  bool hasDelta(const PSimHit* correctHit){
    bool delta = false;
    for (vector<PSimHit>::iterator isim = hitAssociator->SimHitMap[correctHit->detUnitId()].begin();
	 isim != hitAssociator->SimHitMap[correctHit->detUnitId()].end(); ++isim){ 
/*       edm::LogVerbatim("CkfDebugger") << "SimHit on this det at pos="<< position(&*isim)  */
/* 	     << " det=" << isim->detUnitId() << " process=" << isim->processType() ; */
      if (isim->processType() == 9) delta = true;
    }
    return delta;
  }

  Global3DPoint position(const PSimHit* sh) const {
    return theTrackerGeom->idToDetUnit(DetId(sh->detUnitId()))->toGlobal(sh->localPosition());
  };
  const GeomDetUnit* det(const PSimHit* sh) const {return theTrackerGeom->idToDetUnit(DetId(sh->detUnitId()));};

  int layer(const GeomDetUnit* det){return ((int)(((det->geographicalId().rawId() >>16) & 0xF)));}
  int layer(const GeomDet* det){return ((int)(((det->geographicalId().rawId() >>16) & 0xF)));}

  template<unsigned int D>  
  pair<double,double> computePulls(CTTRHp recHit, TSOS startingState){
    typedef typename AlgebraicROOTObject<D>::Vector VecD;
    typedef typename AlgebraicROOTObject<D,D>::SymMatrix SMatDD;
    TSOS detState = theForwardPropagator->propagate(startingState,recHit->det()->surface());
    edm::LogVerbatim("CkfDebugger") << "parameters=" << recHit->parameters() ;
    edm::LogVerbatim("CkfDebugger") << "parametersError=" << recHit->parametersError() ;
    MeasurementExtractor me(detState);
    VecD r = asSVector<D>(recHit->parameters()) - me.measuredParameters<D>(*recHit);
    edm::LogVerbatim("CkfDebugger") << "me.measuredParameters=" << me.measuredParameters<D>(*recHit) ;
    edm::LogVerbatim("CkfDebugger") << "me.measuredError=" << me.measuredError<D>(*recHit) ;
    SMatDD R = asSMatrix<D>(recHit->parametersError()) + me.measuredError<D>(*recHit);
    edm::LogVerbatim("CkfDebugger") << "r=" << r ;
    edm::LogVerbatim("CkfDebugger") << "R=" << R ;
    R.Invert();
    edm::LogVerbatim("CkfDebugger") << "R(-1)=" << R ;
    edm::LogVerbatim("CkfDebugger") << "chi2=" << ROOT::Math::Similarity(r,R) ;
    double pullX=(-r[0])*sqrt(R(0,0));
    double pullY=(-r[1])*sqrt(R(1,1));
    edm::LogVerbatim("CkfDebugger") << "pullX=" << pullX ;
    edm::LogVerbatim("CkfDebugger") << "pullY=" << pullY ;
    return  pair<double,double>(pullX,pullY);
  }
  pair<double,double> computePulls(CTTRHp recHit, TSOS startingState) {
        switch (recHit->dimension()) {
                case 1: return computePulls<1>(recHit,startingState);
                case 2: return computePulls<2>(recHit,startingState);
                case 3: return computePulls<3>(recHit,startingState);
                case 4: return computePulls<4>(recHit,startingState);
                case 5: return computePulls<5>(recHit,startingState);
        }
        throw cms::Exception("CkfDebugger error: rechit of dimension not 1,2,3,4,5");
  }

  vector<int> dump;
  map<pair<int,int>, int> dump2;
  map<pair<int,int>, int> dump3;
  map<pair<int,int>, int> dump4;
  map<pair<int,int>, int> dump5;

  TFile file;
  TH1F* hchi2seedAll, *hchi2seedProb;

  map<string,TH1F*> hPullX_shrh;
  map<string,TH1F*> hPullY_shrh;
  map<string,TH1F*> hPullX_shst;
  map<string,TH1F*> hPullY_shst;
  map<string,TH1F*> hPullX_strh;
  map<string,TH1F*> hPullY_strh;

  map<string,TH1F*> hPullM_shrh;
  map<string,TH1F*> hPullS_shrh;
  map<string,TH1F*> hPullM_shst;
  map<string,TH1F*> hPullS_shst;
  map<string,TH1F*> hPullM_strh;
  map<string,TH1F*> hPullS_strh;

  map<string,TH1F*> hPullGP_X_shst;
  map<string,TH1F*> hPullGP_Y_shst;
  map<string,TH1F*> hPullGP_Z_shst;

  TH2F* hPullGPXvsGPX_shst;
  TH2F* hPullGPXvsGPY_shst;
  TH2F* hPullGPXvsGPZ_shst;
  TH2F* hPullGPXvsGPr_shst;
  TH2F* hPullGPXvsGPeta_shst;
  TH2F* hPullGPXvsGPphi_shst;

  int seedWithDelta;
  int problems;
  int no_sim_hit;
  int no_layer;
  int layer_not_found;
  int det_not_found;
  int chi2gt30;
  int chi2gt30delta;
  int chi2gt30deltaSeed;
  int chi2ls30;
  int simple_hit_not_found;
  int no_component;
  int only_one_component;
  int matched_not_found;
  int matched_not_associated;
  int partner_det_not_fuond;
  int glued_det_not_fuond;
  int propagation;
  int other;
  int totchi2gt30;

  int totSeeds;
};

class less_mag : public std::binary_function<PSimHit*, PSimHit*, bool> {
 public:
  less_mag(){ }
  bool operator()(const PSimHit* a,const PSimHit* b) { 
    return 
      ( a->timeOfFlight()< b->timeOfFlight() );
  }
};

#endif
