/** \file
 *
 *  $Date: 2006/04/25 17:03:24 $
 *  $Revision: 1.1 $
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"


#include "RecoMuon/DetLayers/interface/MuRodBarrelLayer.h"
#include "RecoMuon/DetLayers/interface/MuDetRod.h"

//#include "Geometry/Vector/interface/CoordinateSets.h"


#include "CLHEP/Random/RandFlat.h"

using namespace std;
using namespace edm;

class MuonRecoGeometryAnalyzer : public EDAnalyzer {
 public:

  MuonRecoGeometryAnalyzer( const ParameterSet& pset);

  virtual void analyze( const Event& ev, const EventSetup& es);

  void testDTLayers(const MuonDetLayerGeometry*, const MagneticField* field);

  MeasurementEstimator *theEstimator;
};


  
MuonRecoGeometryAnalyzer::MuonRecoGeometryAnalyzer(const ParameterSet& iConfig) 
{
  float theMaxChi2=25.;
  float theNSigma=3.;
  theEstimator = new Chi2MeasurementEstimator(theMaxChi2,theNSigma);
  
}


void MuonRecoGeometryAnalyzer::analyze( const Event& ev,
				       const EventSetup& es ) {

  ESHandle<MuonDetLayerGeometry> geo;
  es.get<MuonRecoGeometryRecord>().get(geo);

  ESHandle<MagneticField> magfield;
  es.get<IdealMagneticFieldRecord>().get(magfield);
  
  testDTLayers(geo.product(),magfield.product());
}


void MuonRecoGeometryAnalyzer::testDTLayers(const MuonDetLayerGeometry* geo,const MagneticField* field) {

  const vector<DetLayer*>& layers = geo->allDTLayers();
  
  // Get a layer
  const MuRodBarrelLayer* layer = (const MuRodBarrelLayer*) layers.front();

  
  const BoundCylinder& cyl = layer->specificSurface();  

  double halfZ = cyl.bounds().length()/2.;

  // Generate a random point on the cylinder
  double aPhi = RandFlat::shoot(-Geom::pi(),Geom::pi());
  double aZ = RandFlat::shoot(-halfZ, halfZ);
  GlobalPoint gp(GlobalPoint::Cylindrical(cyl.radius(), aPhi, aZ));  

  // Momentum: 10 GeV, straight from the origin
  GlobalVector gv(GlobalVector::Spherical(gp.theta(), aPhi, 10.));

  //FIXME: only negative charge
  int charge = -1;
  cout << "A" <<endl;

  GlobalTrajectoryParameters gtp(gp,gv,charge,field);
  cout << "B" <<endl;
  TrajectoryStateOnSurface tsos(gtp, cyl);
  cout << "testDTLayers: at " << tsos.globalPosition()
       << " R=" << tsos.globalPosition().perp()
       << " Z=" << tsos.globalPosition().z()
       << " p = " << tsos.globalMomentum()
       << endl;


  SteppingHelixPropagator prop(field,alongMomentum);

  pair<bool, TrajectoryStateOnSurface> comp = layer->compatible(tsos,prop,*theEstimator);
  cout << "is compatible: " << comp.first
       << " at: R=" << comp.second.globalPosition().perp()
       << " Z=" <<  comp.second.globalPosition().z()
       << endl;
  
  vector<DetLayer::DetWithState> compDets = layer->compatibleDets(tsos,prop,*theEstimator);
  cout << "compatibleDets: " << compDets.size()
       << endl;    

}



//define this as a plug-in
#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(MuonRecoGeometryAnalyzer)
