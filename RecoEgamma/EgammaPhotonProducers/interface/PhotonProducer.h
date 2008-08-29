#ifndef RecoEgamma_EgammaPhotonProducers_PhotonProducer_h
#define RecoEgamma_EgammaPhotonProducers_PhotonProducer_h
/** \class PhotonProducer
 **  
 **
 **  $Id: PhotonProducer.h,v 1.24 2008/08/17 12:47:38 nancy Exp $ 
 **  $Date: 2008/08/17 12:47:38 $ 
 **  $Revision: 1.24 $
 **  \author Nancy Marinelli, U. of Notre Dame, US
 **
 ***/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeedFwd.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "RecoEgamma/EgammaTools/interface/HoECalculator.h"
#include "RecoEgamma/EgammaTools/interface/ConversionLikelihoodCalculator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"


// PhotonProducer inherits from EDProducer, so it can be a module:
class PhotonProducer : public edm::EDProducer {

 public:

  PhotonProducer (const edm::ParameterSet& ps);
  ~PhotonProducer();

  virtual void beginJob (edm::EventSetup const & es);
  virtual void produce(edm::Event& evt, const edm::EventSetup& es);

 private:

  void fillPhotonCollection(const edm::Handle<reco::SuperClusterCollection> & scHandle,
			    const CaloSubdetectorGeometry *geometry,
			    const CaloSubdetectorGeometry *geometryES,
                            const CaloTopology *topology,
			    const EcalRecHitCollection* hits,
			    HBHERecHitMetaCollection *mhbhe,
			    const edm::Handle<reco::ConversionCollection> & conversionHandle,
			    const reco::ElectronPixelSeedCollection& pixelSeeds,
			    math::XYZPoint & vtx,
			    reco::PhotonCollection & outputCollection,
			    int& iSC);

  reco::ConversionRef solveAmbiguity( const edm::Handle<reco::ConversionCollection> & conversionHandle, reco::SuperClusterRef& sc);

  double hOverE(const reco::SuperClusterRef & scRef, HBHERecHitMetaCollection *mhbhe);

  std::string PhotonCollection_;
  edm::InputTag scHybridBarrelProducer_;
  edm::InputTag scIslandEndcapProducer_;
  edm::InputTag scHybridBarrelCollection_;
  edm::InputTag scIslandEndcapCollection_;

  edm::InputTag barrelEcalHits_;
  edm::InputTag endcapEcalHits_;


  std::string conversionProducer_;
  std::string conversionCollection_;


  std::string hbheLabel_;
  std::string hbheInstanceName_;
  double hOverEConeSize_;
  double maxHOverE_;
  double minSCEt_;
  double minR9_;
  bool validConversions_;
  std::string pixelSeedProducer_;
  std::string vertexProducer_;
  bool usePrimaryVertex_;
  bool risolveAmbiguity_;
  edm::ParameterSet conf_;

  PositionCalc posCalculator_;
  std::string likelihoodWeights_;

  edm::ESHandle<CaloGeometry> theCaloGeom_;
  edm::ESHandle<CaloTopology> theCaloTopo_;
  HoECalculator  theHoverEcalc_;
  ConversionLikelihoodCalculator* theLikelihoodCalc_;

  bool validPixelSeeds_;

  //int nEvt_;

};
#endif
