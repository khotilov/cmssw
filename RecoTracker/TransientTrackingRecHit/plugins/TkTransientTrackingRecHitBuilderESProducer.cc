#include "RecoTracker/TransientTrackingRecHit/plugins/TkTransientTrackingRecHitBuilderESProducer.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"


#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include <string>
#include <memory>

using namespace edm;

TkTransientTrackingRecHitBuilderESProducer::TkTransientTrackingRecHitBuilderESProducer(const edm::ParameterSet & p) 
{
  std::string myname = p.getParameter<std::string>("ComponentName");
  pset_ = p;
  setWhatProduced(this,myname);
}

TkTransientTrackingRecHitBuilderESProducer::~TkTransientTrackingRecHitBuilderESProducer() {}

boost::shared_ptr<TransientTrackingRecHitBuilder> 
TkTransientTrackingRecHitBuilderESProducer::produce(const TransientRecHitRecord & iRecord){ 
//   if (_propagator){
//     delete _propagator;
//     _propagator = 0;
//   }

  std::string sname = pset_.getParameter<std::string>("StripCPE");
  std::string pname = pset_.getParameter<std::string>("PixelCPE");
  std::string mname = pset_.getParameter<std::string>("Matcher");
  
  edm::ESHandle<StripClusterParameterEstimator> se; 
  edm::ESHandle<PixelClusterParameterEstimator> pe; 
  edm::ESHandle<SiStripRecHitMatcher>           me; 
  const StripClusterParameterEstimator * sp ;
  const PixelClusterParameterEstimator * pp ;
  const SiStripRecHitMatcher           * mp ;
    
  if (sname == "Fake") {
    sp = 0;
  }else{
    iRecord.getRecord<TkStripCPERecord>().get( sname, se );     
    sp = se.product();
  }
  
  if (pname == "Fake") {
    pp = 0;
  }else{
    iRecord.getRecord<TkPixelCPERecord>().get( pname, pe );     
    pp = pe.product();
  }
  
  if (mname == "Fake") {
    mp = 0;
  }else{
    iRecord.getRecord<TkStripCPERecord>().get( mname, me );     
    mp = me.product();
  }
  

  edm::ESHandle<TrackerGeometry> pDD;
  iRecord.getRecord<TrackerDigiGeometryRecord>().get( pDD );     
  
  _builder  = boost::shared_ptr<TransientTrackingRecHitBuilder>(new TkTransientTrackingRecHitBuilder(pDD.product(), pp, sp, mp));
  return _builder;
}


