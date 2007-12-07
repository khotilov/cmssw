#ifndef RecoEgamma_EgammaPhotonProducers_ConvertedPhotonProducer_h
#define RecoEgamma_EgammaPhotonProducers_ConvertedPhotonProducer_h
/** \class ConvertedPhotonProducer
 **  
 **
 **  $Id: ConvertedPhotonProducer.h,v 1.10 2007/09/25 10:46:56 nancy Exp $ 
 **  $Date: 2007/09/25 10:46:56 $ 
 **  $Revision: 1.10 $
 **  \author Nancy Marinelli, U. of Notre Dame, US
 **
 ***/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"

#include "TrackingTools/MeasurementDet/interface/LayerMeasurements.h"
#include "TrackingTools/DetLayers/interface/NavigationSetter.h"
#include "TrackingTools/DetLayers/interface/NavigationSchool.h"
#include "RecoTracker/TkNavigation/interface/SimpleNavigationSchool.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"

class ConversionTrackEcalImpactPoint;
class ConversionTrackPairFinder;
class ConversionVertexFinder;
class ConvertedPhotonProducer : public edm::EDProducer {

 public:

  ConvertedPhotonProducer (const edm::ParameterSet& ps);
  ~ConvertedPhotonProducer();


  virtual void beginJob (edm::EventSetup const & es);
  virtual void endJob ();
  virtual void produce(edm::Event& evt, const edm::EventSetup& es);

 private:

  
  std::string conversionOITrackProducerBarrel_;
  std::string conversionIOTrackProducerBarrel_;

  std::string conversionOITrackProducerEndcap_;
  std::string conversionIOTrackProducerEndcap_;

  std::string outInTrackSCBarrelAssociationCollection_;
  std::string inOutTrackSCBarrelAssociationCollection_;

  std::string outInTrackSCEndcapAssociationCollection_;
  std::string inOutTrackSCEndcapAssociationCollection_;


  std::string ConvertedPhotonCollection_;
  std::string PhotonExtraCollection_;


  std::string photonProducer_   ;
  std::string photonCollection_ ;
  std::string photonCorrCollection_ ;

  
  std::string bcProducer_;
  std::string bcBarrelCollection_;
  std::string bcEndcapCollection_;
  std::string scHybridBarrelProducer_;
  std::string scIslandEndcapProducer_;
  std::string scHybridBarrelCollection_;
  std::string scIslandEndcapCollection_;
  edm::ParameterSet conf_;

  std::string barrelClusterShapeMapProducer_;
  std::string barrelClusterShapeMapCollection_;
  std::string endcapClusterShapeMapProducer_;
  std::string endcapClusterShapeMapCollection_;

 

  edm::ESHandle<MagneticField> theMF_;
  edm::ESHandle<GeometricSearchTracker>       theGeomSearchTracker_;
 
  const MeasurementTracker*     theMeasurementTracker_;
  


  ConversionTrackPairFinder*      theTrackPairFinder_;
  ConversionVertexFinder*         theVertexFinder_;
  const LayerMeasurements*      theLayerMeasurements_;
  const NavigationSchool*       theNavigationSchool_;

  ConversionTrackEcalImpactPoint* theEcalImpactPositionFinder_;


 
  
  bool isInitialized;
  int nEvt_;


};
#endif
