// -*- C++ -*-
//
// Package:    EgammaElectronAlgos
// Class:      GsfElectronAlgo
//
/**\class GsfElectronAlgo EgammaElectronAlgos/GsfElectronAlgo

 Description: top algorithm producing TrackCandidate and Electron objects from supercluster
              driven pixel seeded Ckf tracking

*/
//
// Original Author:  Ursula Berthon, Claude Charlot
//         Created:  Thu july 6 13:22:06 CEST 2006
// $Id: GsfElectronAlgo.cc,v 1.84 2009/12/14 23:40:21 chamont Exp $
//
//


#include "RecoEgamma/EgammaElectronAlgos/interface/GsfElectronAlgo.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/EgAmbiguityTools.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronClassification.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronMomentumCorrector.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronEnergyCorrector.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionFactory.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Provenance/interface/Provenance.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
#include "TrackingTools/TransientTrackingRecHit/interface/RecHitComparatorByPosition.h"

#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/Registry.h"


#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include <TMath.h>
#include <Math/VectorUtil.h>
#include <Math/Point3D.h>
#include <sstream>
#include <algorithm>


using namespace edm ;
using namespace std ;
using namespace reco ;

GsfElectronAlgo::GsfElectronAlgo
 ( const edm::ParameterSet & conf,
   double minSCEtBarrel, double minSCEtEndcaps,
   double maxEOverPBarrel, double maxEOverPEndcaps,
   double minEOverPBarrel, double minEOverPEndcaps,
   double maxDeltaEtaBarrel, double maxDeltaEtaEndcaps,
   double maxDeltaPhiBarrel,double maxDeltaPhiEndcaps,
   double maxSigmaIetaIetaBarrel, double maxSigmaIetaIetaEndcaps,
   double maxFbremBarrel, double maxFbremEndcaps,
   bool isBarrel, bool isEndcaps, bool isFiducial,
   bool seedFromTEC,
   double minMVA, double maxTIP,
   double minSCEtBarrelPflow, double minSCEtEndcapsPflow,
   double maxEOverPBarrelPflow, double maxEOverPEndcapsPflow,
   double minEOverPBarrelPflow, double minEOverPEndcapsPflow,
   double maxDeltaEtaBarrelPflow, double maxDeltaEtaEndcapsPflow,
   double maxDeltaPhiBarrelPflow,double maxDeltaPhiEndcapsPflow,
//   double hOverEConeSizePflow, double hOverEPtMinPflow,
//   double maxHOverEDepth1BarrelPflow, double maxHOverEDepth1EndcapsPflow,
//   double maxHOverEDepth2Pflow,
   double maxSigmaIetaIetaBarrelPflow, double maxSigmaIetaIetaEndcapsPflow,
   double maxFbremBarrelPflow, double maxFbremEndcapsPflow,
   bool isBarrelPflow, bool isEndcapsPflow, bool isFiducialPflow,
   double minMVAPflow, double maxTIPPflow,
   bool applyPreselection, bool applyEtaCorrection,
   bool applyAmbResolution, unsigned ambSortingStrategy, unsigned ambClustersOverlapStrategy,
   bool addPflowElectrons,
   double intRadiusTk, double ptMinTk, double maxVtxDistTk, double maxDrbTk,
   double intRadiusHcal, double etMinHcal,
   double intRadiusEcalBarrel, double intRadiusEcalEndcaps, double jurassicWidth,
   double etMinBarrel, double eMinBarrel, double etMinEndcaps, double eMinEndcaps,
   bool vetoClustered, bool useNumCrystals
 )
 : minSCEtBarrel_(minSCEtBarrel), minSCEtEndcaps_(minSCEtEndcaps), maxEOverPBarrel_(maxEOverPBarrel), maxEOverPEndcaps_(maxEOverPEndcaps),
   minEOverPBarrel_(minEOverPBarrel), minEOverPEndcaps_(minEOverPEndcaps),
   maxDeltaEtaBarrel_(maxDeltaEtaBarrel), maxDeltaEtaEndcaps_(maxDeltaEtaEndcaps),
   maxDeltaPhiBarrel_(maxDeltaPhiBarrel),maxDeltaPhiEndcaps_(maxDeltaPhiEndcaps),
   hcalHelper_(0), hcalHelperPflow_(0),
   maxSigmaIetaIetaBarrel_(maxSigmaIetaIetaBarrel), maxSigmaIetaIetaEndcaps_(maxSigmaIetaIetaEndcaps),
   maxFbremBarrel_(maxFbremBarrel), maxFbremEndcaps_(maxFbremEndcaps),
   isBarrel_(isBarrel), isEndcaps_(isEndcaps), isFiducial_(isFiducial),
   seedFromTEC_(seedFromTEC),
   minMVA_(minMVA), maxTIP_(maxTIP),
   minSCEtBarrelPflow_(minSCEtBarrelPflow), minSCEtEndcapsPflow_(minSCEtEndcapsPflow), maxEOverPBarrelPflow_(maxEOverPBarrelPflow), maxEOverPEndcapsPflow_(maxEOverPEndcapsPflow),
   minEOverPBarrelPflow_(minEOverPBarrelPflow), minEOverPEndcapsPflow_(minEOverPEndcapsPflow),
   maxDeltaEtaBarrelPflow_(maxDeltaEtaBarrelPflow), maxDeltaEtaEndcapsPflow_(maxDeltaEtaEndcapsPflow),
   maxDeltaPhiBarrelPflow_(maxDeltaPhiBarrelPflow),maxDeltaPhiEndcapsPflow_(maxDeltaPhiEndcapsPflow),
//   hOverEConeSizePflow_(hOverEConeSizePflow), hOverEPtMinPflow_(hOverEPtMinPflow),
//   maxHOverEDepth1BarrelPflow_(maxHOverEDepth1BarrelPflow), maxHOverEDepth1EndcapsPflow_(maxHOverEDepth1EndcapsPflow),
//   maxHOverEDepth2Pflow_(maxHOverEDepth2Pflow),
   maxSigmaIetaIetaBarrelPflow_(maxSigmaIetaIetaBarrelPflow), maxSigmaIetaIetaEndcapsPflow_(maxSigmaIetaIetaEndcapsPflow),
   maxFbremBarrelPflow_(maxFbremBarrelPflow), maxFbremEndcapsPflow_(maxFbremEndcapsPflow),
   isBarrelPflow_(isBarrelPflow), isEndcapsPflow_(isEndcapsPflow), isFiducialPflow_(isFiducialPflow),
   minMVAPflow_(minMVAPflow), maxTIPPflow_(maxTIPPflow),
   applyPreselection_(applyPreselection), applyEtaCorrection_(applyEtaCorrection),
   applyAmbResolution_(applyAmbResolution), ambSortingStrategy_(ambSortingStrategy), ambClustersOverlapStrategy_(ambClustersOverlapStrategy),
   addPflowElectrons_(addPflowElectrons),
   intRadiusTk_(intRadiusTk), ptMinTk_(ptMinTk),  maxVtxDistTk_(maxVtxDistTk),  maxDrbTk_(maxDrbTk),
   intRadiusHcal_(intRadiusHcal), etMinHcal_(etMinHcal), intRadiusEcalBarrel_(intRadiusEcalBarrel),  intRadiusEcalEndcaps_(intRadiusEcalEndcaps),  jurassicWidth_(jurassicWidth),
   etMinBarrel_(etMinBarrel),  eMinBarrel_(eMinBarrel),  etMinEndcaps_(etMinEndcaps),  eMinEndcaps_(eMinEndcaps),
   vetoClustered_(vetoClustered), useNumCrystals_(useNumCrystals), ctfTracksCheck_(false),
   cacheIDGeom_(0),cacheIDTopo_(0),cacheIDTDGeom_(0),cacheIDMagField_(0),
   superClusterErrorFunction_(0),
   pfTranslatorParametersChecked_(false), ecalSeedingParametersChecked_(false)
 {
   // this is the new version allowing to configurate the algo
   // interfaces still need improvement!!
   mtsTransform_ = 0 ;
   constraintAtVtx_ = 0;

   // get nested parameter set for the TransientInitialStateEstimator
   ParameterSet tise_params = conf.getParameter<ParameterSet>("TransientInitialStateEstimatorParameters") ;

   // hcal strategy
   hcalHelper_ = new ElectronHcalHelper(conf) ;
   hcalHelperPflow_ = new ElectronHcalHelper(conf,true,true) ;
   //hOverEConeSize_ = conf.getParameter<double>("hOverEConeSize") ;
   hcalTowers_ = conf.getParameter<edm::InputTag>("hcalTowers") ;
   //hOverEPtMin_ = conf.getParameter<double>("hOverEPtMin") ;
//     maxHOverEDepth1Barrel_ = conf.getParameter<double>("maxHOverEDepth1Barrel") ;
//     maxHOverEDepth1Endcaps_ = conf.getParameter<double>("maxHOverEDepth1Endcaps") ;
//     maxHOverEDepth2_ = conf.getParameter<double>("maxHOverEDepth2") ;
//     maxHOverEDepth1BarrelPflow_ = conf.getParameter<double>("maxHOverEDepth1BarrelPflow") ;
//     maxHOverEDepth1EndcapsPflow_ = conf.getParameter<double>("maxHOverEDepth1EndcapsPflow") ;
//     maxHOverEDepth2Pflow_ = conf.getParameter<double>("maxHOverEDepth2Pflow") ;
    maxHOverEBarrel_ = conf.getParameter<double>("maxHOverEBarrel") ;
    maxHOverEEndcaps_ = conf.getParameter<double>("maxHOverEEndcaps") ;
    maxHBarrel_ = conf.getParameter<double>("maxHBarrel") ;
    maxHEndcaps_ = conf.getParameter<double>("maxHEndcaps") ;
    maxHOverEBarrelPflow_ = conf.getParameter<double>("maxHOverEBarrelPflow") ;
    maxHOverEEndcapsPflow_ = conf.getParameter<double>("maxHOverEEndcapsPflow") ;
    maxHBarrelPflow_ = conf.getParameter<double>("maxHBarrelPflow") ;
    maxHEndcapsPflow_ = conf.getParameter<double>("maxHEndcapsPflow") ;

  // get input collections
  //tracks_ = conf.getParameter<edm::InputTag>("tracks");
  gsfElectronCores_ = conf.getParameter<edm::InputTag>("gsfElectronCores");
  reducedBarrelRecHitCollection_ = conf.getParameter<edm::InputTag>("reducedBarrelRecHitCollection") ;
  reducedEndcapRecHitCollection_ = conf.getParameter<edm::InputTag>("reducedEndcapRecHitCollection") ;
  pfMVA_ = conf.getParameter<edm::InputTag>("pfMVA") ;
  ctfTracks_ = conf.getParameter<edm::InputTag>("ctfTracks");
  seedsTag_ = conf.getParameter<edm::InputTag>("seedsTag");

  // for backward compatibility
  ctfTracksCheck_ = conf.getParameter<bool>("ctfTracksCheck");

  // function for corrector
  std::string superClusterErrorFunctionName
   = conf.getParameter<std::string>("superClusterErrorFunction") ;
  if (superClusterErrorFunctionName!="")
   {
    superClusterErrorFunction_
     = EcalClusterFunctionFactory::get()->create(superClusterErrorFunctionName,conf) ;
   }
}

GsfElectronAlgo::~GsfElectronAlgo() {
	delete constraintAtVtx_;
	delete mtsTransform_;
  delete hcalHelper_ ;
  delete hcalHelperPflow_ ;
}

void GsfElectronAlgo::setupES(const edm::EventSetup& es) {

  // get EventSetupRecords if needed
  bool updateField(false);
  if (cacheIDMagField_!=es.get<IdealMagneticFieldRecord>().cacheIdentifier()){
    updateField = true;
    cacheIDMagField_=es.get<IdealMagneticFieldRecord>().cacheIdentifier();
    es.get<IdealMagneticFieldRecord>().get(theMagField);
  }

  bool updateGeometry(false);
  if (cacheIDTDGeom_!=es.get<TrackerDigiGeometryRecord>().cacheIdentifier()){
    updateGeometry = true;
    cacheIDTDGeom_=es.get<TrackerDigiGeometryRecord>().cacheIdentifier();
    es.get<TrackerDigiGeometryRecord>().get(trackerHandle_);
  }

  if ( updateField || updateGeometry ) {
    delete mtsTransform_;
    mtsTransform_ = new MultiTrajectoryStateTransform(trackerHandle_.product(),theMagField.product());
  }

  if ( updateField || updateGeometry ) {
    delete constraintAtVtx_;
    constraintAtVtx_ = new GsfConstraintAtVertex(es);
  }

  if (cacheIDGeom_!=es.get<CaloGeometryRecord>().cacheIdentifier()){
    cacheIDGeom_=es.get<CaloGeometryRecord>().cacheIdentifier();
    es.get<CaloGeometryRecord>().get(theCaloGeom);
  }

  if (cacheIDTopo_!=es.get<CaloTopologyRecord>().cacheIdentifier()){
    cacheIDTopo_=es.get<CaloTopologyRecord>().cacheIdentifier();
    es.get<CaloTopologyRecord>().get(theCaloTopo);
  }

//  if (useHcalRecHits_)
//   {
  hcalHelper_->checkSetup(es) ;
  hcalHelperPflow_->checkSetup(es) ;
//   }
  if (superClusterErrorFunction_)
   { superClusterErrorFunction_->init(es) ; }
 }

void  GsfElectronAlgo::run(Event& e, GsfElectronCollection & outEle) {

  // get the input
  edm::Handle<GsfElectronCoreCollection> coresH;
  e.getByLabel(gsfElectronCores_,coresH);
  edm::Handle<TrackCollection> ctfTracksH;
  e.getByLabel(ctfTracks_, ctfTracksH);
  edm::Handle< EcalRecHitCollection > pEBRecHits;
  e.getByLabel( reducedBarrelRecHitCollection_, pEBRecHits );
  edm::Handle< EcalRecHitCollection > pEERecHits;
  e.getByLabel( reducedEndcapRecHitCollection_, pEERecHits ) ;
  edm::Handle<CaloTowerCollection> towersH;
  e.getByLabel(hcalTowers_, towersH);
  edm::Handle<edm::ValueMap<float> > pfMVAH;
  e.getByLabel(pfMVA_,pfMVAH);
  edm::Handle<ElectronSeedCollection> seedsCollection;
  e.getByLabel(seedsTag_,seedsCollection);

  if (!pfTranslatorParametersChecked_)
   {
    pfTranslatorParametersChecked_ = true ;
    checkPfTranslatorParameters(pfMVAH.provenance()->psetID()) ;
   }
  if (!seedsCollection.isValid())
   {
    edm::LogWarning("GsfElectronAlgo")
      <<"Cannot check consistency of parameters with ecal seeding ones,"
      <<" because the original collection of seeds is not any more available." ;
   }
  else if (!ecalSeedingParametersChecked_)
   {
    ecalSeedingParametersChecked_ = true ;
    checkEcalSeedingParameters(seedsCollection.provenance()->psetID()) ;
   }

  // get the beamspot from the Event:
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  e.getByType(recoBeamSpotHandle);
  const BeamSpot bs = *recoBeamSpotHandle;

  // prepare access to hcal data
  hcalHelper_->readEvent(e) ;
  hcalHelperPflow_->readEvent(e) ;

  // temporay array for electrons before preselection and before amb. solving
  GsfElectronPtrCollection tempEle, tempEle1;

  // create electrons
  process(coresH,ctfTracksH,pfMVAH,towersH,pEBRecHits,pEERecHits,bs,tempEle);

  std::ostringstream str;

  str << "\n========== GsfElectronAlgo Info (before preselection) ==========";
  str << "\nEvent " << e.id();
  str << "\nNumber of electron cores: " << coresH.product()->size();
  str << "\nNumber of electrons: " << tempEle.size();
  for (GsfElectronPtrCollection::const_iterator it = tempEle.begin(); it != tempEle.end(); it++) {
    str << "\nNew electron with charge, pt, eta, phi : "  << (*it)->charge() << " , "
        << (*it)->pt() << " , " << (*it)->eta() << " , " << (*it)->phi();
  }

  if (applyPreselection_)
   {

    preselectElectrons(tempEle, tempEle1, bs);

    std::ostringstream str1 ;

    str1 << "\n=================================================";
    LogDebug("GsfElectronAlgo") << str.str();

    str1 << "\n========== GsfElectronAlgo Info (before amb. solving) ==========";
    str1 << "\nEvent " << e.id();
    str1 << "\nNumber of preselected electrons: " << tempEle1.size();
    GsfElectronPtrCollection::const_iterator it ;
    for ( it = tempEle1.begin(); it != tempEle1.end(); it++)
     {
      str1 << "\nNew electron with charge, pt, eta, phi : "  << (*it)->charge() << " , "
          << (*it)->pt() << " , " << (*it)->eta() << " , " << (*it)->phi();
     }
   }
  else
   {
    for ( GsfElectronPtrCollection::const_iterator it = tempEle.begin() ; it != tempEle.end() ; it++ )
     { tempEle1.push_back(*it) ; }
   }

  str << "\n=================================================";
  LogDebug("GsfElectronAlgo") << str.str();

  if (applyAmbResolution_)
   {
    //resolveElectrons(tempEle1, outEle);
    resolveElectrons(tempEle1, outEle, pEBRecHits, pEERecHits);

    std::ostringstream str2 ;

    str2 << "\n========== GsfElectronAlgo Info (after amb. solving) ==========";
    str2 << "\nEvent " << e.id();
    str2 << "\nNumber of preselected and resolved electrons: " << outEle.size();
    for ( GsfElectronCollection::const_iterator it = outEle.begin(); it != outEle.end(); it++) {
      str2 << "\nNew electron with charge, pt, eta, phi : "  << it->charge() << " , "
          << it->pt() << " , " << it->eta() << " , " << it->phi();
    }
    str2 << "\n=================================================";
    LogDebug("GsfElectronAlgo") << str2.str() ;

   }
  else
   {
    for ( GsfElectronPtrCollection::const_iterator it = tempEle1.begin() ; it != tempEle1.end() ; it++ )
     { outEle.push_back(**it) ; }
   }

  // delete temporary electrons
  GsfElectronPtrCollection::const_iterator it ;
  for ( it = tempEle.begin() ; it != tempEle.end() ; it++ )
   { delete (*it) ; }

  return ;
}

void GsfElectronAlgo::process(
  edm::Handle<GsfElectronCoreCollection> coresH,
  edm::Handle<TrackCollection> ctfTracksH,
  edm::Handle<edm::ValueMap<float> > pfMVAH,
  edm::Handle<CaloTowerCollection> towersH,
  edm::Handle<EcalRecHitCollection> reducedEBRecHits,
  edm::Handle<EcalRecHitCollection> reducedEERecHits,
  const BeamSpot & bs,
  GsfElectronPtrCollection & outEle )
 {

  // Isolation algos

  float extRadiusSmall=0.3, extRadiusLarge=0.4, intRadius=intRadiusTk_;
  float ptMin=ptMinTk_, maxVtxDist=maxVtxDistTk_, drb=maxDrbTk_;
  ElectronTkIsolation tkIsolation03(extRadiusSmall,intRadius,ptMin,maxVtxDist,drb,ctfTracksH.product(),bs.position()) ;
  ElectronTkIsolation tkIsolation04(extRadiusLarge,intRadius,ptMin,maxVtxDist,drb,ctfTracksH.product(),bs.position()) ;

  float egHcalIsoConeSizeOutSmall=0.3, egHcalIsoConeSizeOutLarge=0.4;
  float egHcalIsoConeSizeIn=intRadiusHcal_,egHcalIsoPtMin=etMinHcal_;
  int egHcalDepth1=1, egHcalDepth2=2;
  EgammaTowerIsolation hadDepth1Isolation03(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth1,towersH.product()) ;
  EgammaTowerIsolation hadDepth2Isolation03(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth2,towersH.product()) ;
  EgammaTowerIsolation hadDepth1Isolation04(egHcalIsoConeSizeOutLarge,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth1,towersH.product()) ;
  EgammaTowerIsolation hadDepth2Isolation04(egHcalIsoConeSizeOutLarge,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth2,towersH.product()) ;

  float egIsoConeSizeOutSmall=0.3, egIsoConeSizeOutLarge=0.4, egIsoJurassicWidth=jurassicWidth_;
  float egIsoPtMinBarrel=etMinBarrel_,egIsoEMinBarrel=eMinBarrel_, egIsoConeSizeInBarrel=intRadiusEcalBarrel_;
  float egIsoPtMinEndcap=etMinEndcaps_,egIsoEMinEndcap=eMinEndcaps_, egIsoConeSizeInEndcap=intRadiusEcalEndcaps_;
  EcalRecHitMetaCollection ecalBarrelHits(*reducedEBRecHits);
  EcalRecHitMetaCollection ecalEndcapHits(*reducedEERecHits);
  EgammaRecHitIsolation ecalBarrelIsol03(egIsoConeSizeOutSmall,egIsoConeSizeInBarrel,egIsoJurassicWidth,egIsoPtMinBarrel,egIsoEMinBarrel,theCaloGeom,&ecalBarrelHits,DetId::Ecal);
  EgammaRecHitIsolation ecalBarrelIsol04(egIsoConeSizeOutLarge,egIsoConeSizeInBarrel,egIsoJurassicWidth,egIsoPtMinBarrel,egIsoEMinBarrel,theCaloGeom,&ecalBarrelHits,DetId::Ecal);
  EgammaRecHitIsolation ecalEndcapIsol03(egIsoConeSizeOutSmall,egIsoConeSizeInEndcap,egIsoJurassicWidth,egIsoPtMinEndcap,egIsoEMinEndcap,theCaloGeom,&ecalEndcapHits,DetId::Ecal);
  EgammaRecHitIsolation ecalEndcapIsol04(egIsoConeSizeOutLarge,egIsoConeSizeInEndcap,egIsoJurassicWidth,egIsoPtMinEndcap,egIsoEMinEndcap,theCaloGeom,&ecalEndcapHits,DetId::Ecal);
  ecalBarrelIsol03.setUseNumCrystals(useNumCrystals_);
  ecalBarrelIsol03.setVetoClustered(vetoClustered_);
  ecalBarrelIsol04.setUseNumCrystals(useNumCrystals_);
  ecalBarrelIsol04.setVetoClustered(vetoClustered_);
  ecalEndcapIsol03.setUseNumCrystals(useNumCrystals_);
  ecalEndcapIsol03.setVetoClustered(vetoClustered_);
  ecalEndcapIsol04.setUseNumCrystals(useNumCrystals_);
  ecalEndcapIsol04.setVetoClustered(vetoClustered_);

  //const GsfTrackCollection * gsfTrackCollection = gsfTracksH.product() ;
  const GsfElectronCoreCollection * coreCollection = coresH.product() ;
  for (unsigned int i=0;i<coreCollection->size();++i) {

    // retreive core, track and scl
    const GsfElectronCoreRef coreRef = edm::Ref<GsfElectronCoreCollection>(coresH,i);
    const GsfTrackRef gsfTrackRef = coreRef->gsfTrack() ; //edm::Ref<GsfTrackCollection>(gsfTracksH,i);

    // don't add pflow only electrons if one so wish
    if (!coreRef->ecalDrivenSeed() && !addPflowElectrons_) continue ;

    // Get the super cluster
    SuperClusterRef scRef = coreRef->superCluster() ;
    SuperCluster theClus = *scRef ;

    // mva
    //const edm::ValueMap<float> & pfmvas = *pfMVAH.product() ;
    //float mva=std::numeric_limits<float>::infinity();
    //if (coreRef->trackerDrivenSeed()) mva = pfmvas[gsfTrackRef];

    // MVA
    // we check that the value is never inferior to the no-cut value
    // we generally use in the configuration file for minMVA.
    float mva = (*pfMVAH.product())[gsfTrackRef] ;
    double noCutMin = -999999999 ;
    if (mva<noCutMin) { edm::LogError("GsfElectronAlgo")<<"unexpected MVA value: "<<mva ; }

    // electron basic cluster
    CaloClusterPtr elbcRef = getEleBasicCluster(gsfTrackRef,&theClus) ;

    // calculate Trajectory StatesOnSurface....
    if (!calculateTSOS(*gsfTrackRef,theClus, bs)) continue ;
    mtsMode_->momentumFromModeCartesian(vtxTSOS_,vtxMom_) ;
    sclPos_=sclTSOS_.globalPosition() ;

    // H/E
    double HoE1 = 0. ;
    double HoE2 = 0. ;
    if (coreRef->ecalDrivenSeed())
     {
      HoE1 = hcalHelper_->hcalESumDepth1(theClus)/theClus.energy() ;
      HoE2 = hcalHelper_->hcalESumDepth2(theClus)/theClus.energy() ;
     }
    else
     {
      HoE1 = hcalHelperPflow_->hcalESumDepth1(theClus)/theClus.energy() ;
      HoE2 = hcalHelperPflow_->hcalESumDepth2(theClus)/theClus.energy() ;
     }

    // ctf track and charge ID
    TrackRef ctfTrackRef = coreRef->ctfTrack() ;
    float fracShHits = coreRef->ctfGsfOverlap() ;
    if (ctfTracksCheck_&&ctfTrackRef.isNull())
     {
      pair<TrackRef,float> ctfpair = getCtfTrackRef(gsfTrackRef,ctfTracksH) ;
      ctfTrackRef = ctfpair.first ;
      fracShHits = ctfpair.second ;
     }
    int eleCharge ;
    GsfElectron::ChargeInfo eleChargeInfo ;
    computeCharge(gsfTrackRef,ctfTrackRef,scRef,bs,eleCharge,eleChargeInfo) ;

    createElectron(coreRef,eleCharge,eleChargeInfo,
      elbcRef,ctfTrackRef,fracShHits,HoE1,HoE2,tkIsolation03,tkIsolation04,
      hadDepth1Isolation03,hadDepth2Isolation03,hadDepth1Isolation04,hadDepth2Isolation04,
      ecalBarrelIsol03,ecalEndcapIsol03,ecalBarrelIsol04,ecalEndcapIsol04,
      reducedEBRecHits,reducedEERecHits,mva,outEle) ;

     LogInfo("")<<"Constructed new electron with energy  "<< theClus.energy();

  } // loop over tracks
}

void GsfElectronAlgo::preselectElectrons( GsfElectronPtrCollection & inEle, GsfElectronPtrCollection & outEle, const reco::BeamSpot& bs )
 {
  GsfElectronPtrCollection::size_type ei, emax = inEle.size() ;
  GsfElectronPtrCollection::iterator e1 ;
  for( ei=1, e1=inEle.begin() ;  e1!=inEle.end() ; ++ei, ++e1 )
   {
    LogDebug("")<<"========== preSelection "<<ei<<"/"<<emax<<"==========" ;

    // kind of construction algorithm
    bool eg = (*e1)->core()->ecalDrivenSeed();
    bool pf = (*e1)->core()->trackerDrivenSeed() && !(*e1)->core()->ecalDrivenSeed() ;
    if (eg&&pf) { edm::LogError("GsfElectronAlgo")<<"An electron cannot be both egamma and purely pflow" ; }
    if ((!eg)&&(!pf)) { edm::LogError("GsfElectronAlgo")<<"An electron cannot be neither egamma nor purely pflow" ; }

    // Et cut
    LogDebug("") << "Et : " << (*e1)->superCluster()->energy()/cosh((*e1)->superCluster()->eta());
    if (eg && (*e1)->isEB() && ((*e1)->superCluster()->energy()/cosh((*e1)->superCluster()->eta()) < minSCEtBarrel_)) continue;
    if (eg && (*e1)->isEE() && ((*e1)->superCluster()->energy()/cosh((*e1)->superCluster()->eta()) < minSCEtEndcaps_)) continue;
    if (pf && (*e1)->isEB() && ((*e1)->superCluster()->energy()/cosh((*e1)->superCluster()->eta()) < minSCEtBarrelPflow_)) continue;
    if (pf && (*e1)->isEE() && ((*e1)->superCluster()->energy()/cosh((*e1)->superCluster()->eta()) < minSCEtEndcapsPflow_)) continue;

    // E/p cut
    LogDebug("") << "E/p : " << (*e1)->eSuperClusterOverP();
    if (eg && (*e1)->isEB() && ((*e1)->eSuperClusterOverP() > maxEOverPBarrel_)) continue;
    if (eg && (*e1)->isEE() && ((*e1)->eSuperClusterOverP() > maxEOverPEndcaps_)) continue;
    if (eg && (*e1)->isEB() && ((*e1)->eSuperClusterOverP() < minEOverPBarrel_)) continue;
    if (eg && (*e1)->isEE() && ((*e1)->eSuperClusterOverP() < minEOverPEndcaps_)) continue;
    if (pf && (*e1)->isEB() && ((*e1)->eSuperClusterOverP() > maxEOverPBarrelPflow_)) continue;
    if (pf && (*e1)->isEE() && ((*e1)->eSuperClusterOverP() > maxEOverPEndcapsPflow_)) continue;
    if (pf && (*e1)->isEB() && ((*e1)->eSuperClusterOverP() < minEOverPBarrelPflow_)) continue;
    if (pf && (*e1)->isEE() && ((*e1)->eSuperClusterOverP() < minEOverPEndcapsPflow_)) continue;
    LogDebug("") << "E/p criteria is satisfied ";

    // HoE cuts
    LogDebug("") << "HoE1 : " << (*e1)->hcalDepth1OverEcal() << "HoE2 : " << (*e1)->hcalDepth2OverEcal();
//     if ( eg && (*e1)->isEB() && ((*e1)->hcalDepth1OverEcal() > maxHOverEDepth1Barrel_) ) continue;
//     if ( eg && (*e1)->isEE() && ((*e1)->hcalDepth1OverEcal() > maxHOverEDepth1Endcaps_) ) continue;
//     if ( eg && ((*e1)->hcalDepth2OverEcal() > maxHOverEDepth2_) ) continue;
//     if ( pf && (*e1)->isEB() && ((*e1)->hcalDepth1OverEcal() > maxHOverEDepth1BarrelPflow_) ) continue;
//     if ( pf && (*e1)->isEE() && ((*e1)->hcalDepth1OverEcal() > maxHOverEDepth1EndcapsPflow_) ) continue;
//     if ( pf && ((*e1)->hcalDepth2OverEcal() > maxHOverEDepth2Pflow_) ) continue;
    double had = (*e1)->hcalOverEcal()*(*e1)->superCluster()->energy();
    bool HoEveto = false;
    const reco::CaloCluster & seedCluster = *((*e1)->superCluster()->seed()) ;
    int detector = seedCluster.hitsAndFractions()[0].first.subdetId() ;
    if (detector==EcalBarrel && (had<maxHBarrel_ || (had/(*e1)->superCluster()->energy())<maxHOverEBarrel_)) HoEveto=true;
    else if (detector==EcalEndcap && (had<maxHEndcaps_ || (had/(*e1)->superCluster()->energy())<maxHOverEEndcaps_)) HoEveto=true;
    if ( eg && !HoEveto ) continue;
    bool HoEvetoPflow = false;
    if (detector==EcalBarrel && (had<maxHBarrelPflow_ || (had/(*e1)->superCluster()->energy())<maxHOverEBarrelPflow_)) HoEvetoPflow=true;
    else if (detector==EcalEndcap && (had<maxHEndcapsPflow_ || (had/(*e1)->superCluster()->energy())<maxHOverEEndcapsPflow_)) HoEvetoPflow=true;
    if ( pf && !HoEvetoPflow ) continue;
    LogDebug("") << "H/E criteria is satisfied ";

    // delta eta criteria
    double deta = (*e1)->deltaEtaSuperClusterTrackAtVtx();
    LogDebug("") << "delta eta : " << deta;
    if (eg && (*e1)->isEB() && (fabs(deta) > maxDeltaEtaBarrel_)) continue;
    if (eg && (*e1)->isEE() && (fabs(deta) > maxDeltaEtaEndcaps_)) continue;
    if (pf && (*e1)->isEB() && (fabs(deta) > maxDeltaEtaBarrelPflow_)) continue;
    if (pf && (*e1)->isEE() && (fabs(deta) > maxDeltaEtaEndcapsPflow_)) continue;
    LogDebug("") << "Delta eta criteria is satisfied ";

    // delta phi criteria
    double dphi = (*e1)->deltaPhiSuperClusterTrackAtVtx();
    LogDebug("") << "delta phi : " << dphi;
    if (eg && (*e1)->isEB() && (fabs(dphi) > maxDeltaPhiBarrel_)) continue;
    if (eg && (*e1)->isEE() && (fabs(dphi) > maxDeltaPhiEndcaps_)) continue;
    if (pf && (*e1)->isEB() && (fabs(dphi) > maxDeltaPhiBarrelPflow_)) continue;
    if (pf && (*e1)->isEE() && (fabs(dphi) > maxDeltaPhiEndcapsPflow_)) continue;
    LogDebug("") << "Delta phi criteria is satisfied ";

    //
    if (eg && (*e1)->isEB() && ((*e1)->sigmaIetaIeta() > maxSigmaIetaIetaBarrel_)) continue;
    if (eg && (*e1)->isEE() && ((*e1)->sigmaIetaIeta() > maxSigmaIetaIetaEndcaps_)) continue;
    if (pf && (*e1)->isEB() && ((*e1)->sigmaIetaIeta() > maxSigmaIetaIetaBarrelPflow_)) continue;
    if (pf && (*e1)->isEE() && ((*e1)->sigmaIetaIeta() > maxSigmaIetaIetaEndcapsPflow_)) continue;

    // fiducial
    if (eg && !(*e1)->isEB() && isBarrel_) continue;
    if (eg && !(*e1)->isEE() && isEndcaps_) continue;
    if (eg && ((*e1)->isEBEEGap() || (*e1)->isEBEtaGap() || (*e1)->isEBPhiGap() || (*e1)->isEERingGap() || (*e1)->isEEDeeGap())
     && isFiducial_) continue;
    if (pf && !(*e1)->isEB() && isBarrelPflow_) continue;
    if (pf && !(*e1)->isEE() && isEndcapsPflow_) continue;
    if (pf && ((*e1)->isEBEEGap() || (*e1)->isEBEtaGap() || (*e1)->isEBPhiGap() || (*e1)->isEERingGap() || (*e1)->isEEDeeGap())
     && isFiducialPflow_) continue;

    // seed in TEC
    edm::RefToBase<TrajectorySeed> seed = (*e1)->gsfTrack()->extra()->seedRef() ;
    ElectronSeedRef elseed = seed.castTo<ElectronSeedRef>() ;
    if (eg && !seedFromTEC_)
     {
      if (elseed.isNull())
	   { edm::LogError("GsfElectronAlgo")<<"The GsfTrack seed is not an ElectronSeed ?!" ; }
	  else
	   { if (elseed->subDet2()==6) continue ; }
     }

    // BDT output
    if (eg && ((*e1)->mva()<minMVA_)) continue ;
    if (pf && ((*e1)->mva()<minMVAPflow_)) continue ;

    // transverse impact parameter
    if (eg && fabs((*e1)->gsfTrack()->dxy(bs.position()))>maxTIP_) continue;
    if (pf && fabs((*e1)->gsfTrack()->dxy(bs.position()))>maxTIPPflow_) continue;

    LogDebug("") << "electron has passed preselection criteria ";
    LogDebug("") << "=================================================";

    outEle.push_back(*e1) ;

   }
 }

// utilities for constructor
float normalized_dphi( float dphi )
 {
  if (fabs(dphi)>CLHEP::pi) return (dphi<0?CLHEP::twopi+dphi:dphi-CLHEP::twopi) ;
  else return dphi ;
 }

math::XYZPoint convert( const GlobalPoint & gp )
 { return math::XYZPoint(gp.x(),gp.y(),gp.z()) ; }

math::XYZVector convert( const GlobalVector & gv )
 { return math::XYZVector(gv.x(),gv.y(),gv.z()) ; }

// interface to be improved...
void GsfElectronAlgo::createElectron
 ( const GsfElectronCoreRef & coreRef,
   int charge, const reco::GsfElectron::ChargeInfo & chargeInfo,
   const CaloClusterPtr & elbcRef,
   const TrackRef & ctfTrackRef, const float shFracInnerHits,
   double HoE1, double HoE2,
   ElectronTkIsolation & tkIso03, ElectronTkIsolation & tkIso04,
   EgammaTowerIsolation & had1Iso03, EgammaTowerIsolation & had2Iso03,
   EgammaTowerIsolation & had1Iso04, EgammaTowerIsolation & had2Iso04,
   EgammaRecHitIsolation & ecalBarrelIso03,EgammaRecHitIsolation & ecalEndcapsIso03,
   EgammaRecHitIsolation & ecalBarrelIso04,EgammaRecHitIsolation & ecalEndcapsIso04,
   edm::Handle<EcalRecHitCollection> reducedEBRecHits,edm::Handle<EcalRecHitCollection> reducedEERecHits,
   float mva, GsfElectronPtrCollection & outEle )

 {
  GsfTrackRef trackRef = coreRef->gsfTrack() ;
  SuperClusterRef scRef = coreRef->superCluster() ;
  if (scRef.isNull()) return ;

  // Seed info
  const reco::CaloCluster & seedCluster = *(scRef->seed()) ;
  DetId seedXtalId = seedCluster.hitsAndFractions()[0].first ;

  // various useful positions and momemtums
  GlobalVector innMom, seedMom, eleMom, sclMom, outMom ;
  mtsMode_->momentumFromModeCartesian(innTSOS_,innMom) ;
  GlobalPoint innPos=innTSOS_.globalPosition() ;
  mtsMode_->momentumFromModeCartesian(seedTSOS_,seedMom) ;
  GlobalPoint  seedPos=seedTSOS_.globalPosition() ;
  mtsMode_->momentumFromModeCartesian(eleTSOS_,eleMom) ;
  GlobalPoint  elePos=eleTSOS_.globalPosition() ;
  mtsMode_->momentumFromModeCartesian(sclTSOS_,sclMom) ;
  // sclPos_ already here
  // vtxMom_ already here
  GlobalPoint vtxPos=vtxTSOS_.globalPosition() ;
  mtsMode_->momentumFromModeCartesian(outTSOS_,outMom);
  GlobalPoint outPos=outTSOS_.globalPosition() ;
  GlobalVector vtxMomWithConstraint;
  mtsMode_->momentumFromModeCartesian(constrainedVtxTSOS_,vtxMomWithConstraint);


  //====================================================
  // Candidate attributes
  //====================================================

  double scale = (*scRef).energy()/vtxMom_.mag() ;
  Candidate::LorentzVector momentum = Candidate::LorentzVector
   ( vtxMom_.x()*scale,vtxMom_.y()*scale,vtxMom_.z()*scale,
	 (*scRef).energy() ) ;


  //====================================================
  // Track-Cluster Matching
  //====================================================

  reco::GsfElectron::TrackClusterMatching tcMatching ;
  tcMatching.electronCluster = elbcRef ;
  tcMatching.eSuperClusterOverP = (vtxMom_.mag()>0)?(scRef->energy()/vtxMom_.mag()):(-1.) ;
  tcMatching.eSeedClusterOverP = (vtxMom_.mag()>0.)?(seedCluster.energy()/vtxMom_.mag()):(-1) ;
  tcMatching.eSeedClusterOverPout = (seedMom.mag()>0.)?(seedCluster.energy()/seedMom.mag()):(-1.) ;
  tcMatching.eEleClusterOverPout = (eleMom.mag()>0.)?(elbcRef->energy()/eleMom.mag()):(-1.) ;
  tcMatching.deltaEtaSuperClusterAtVtx = scRef->eta()-sclPos_.eta() ;
  tcMatching.deltaEtaSeedClusterAtCalo = seedCluster.eta() - seedPos.eta() ;
  tcMatching.deltaEtaEleClusterAtCalo = elbcRef->eta() - seedPos.eta() ;
  tcMatching.deltaPhiEleClusterAtCalo = normalized_dphi(elbcRef->phi()-elePos.phi()) ;
  tcMatching.deltaPhiSuperClusterAtVtx = normalized_dphi(scRef->phi()-sclPos_.phi()) ;
  tcMatching.deltaPhiSeedClusterAtCalo = normalized_dphi(seedCluster.phi()-seedPos.phi()) ;


  //=======================================================
  // Track extrapolations
  //=======================================================

  reco::GsfElectron::TrackExtrapolations tkExtra ;
  tkExtra.positionAtVtx = convert(vtxPos) ;
  tkExtra.positionAtCalo = convert(sclPos_) ;
  tkExtra.momentumAtVtx = convert(vtxMom_) ;
  tkExtra.momentumAtCalo = convert(sclMom) ;
  tkExtra.momentumOut = convert(seedMom) ;
  tkExtra.momentumAtEleClus = convert(eleMom) ;
  tkExtra.momentumAtVtxWithConstraint = convert(vtxMomWithConstraint);


  //=======================================================
  // Closest Ctf Track
  //=======================================================

  reco::GsfElectron::ClosestCtfTrack ctfInfo ;
  ctfInfo.ctfTrack = ctfTrackRef  ;
  ctfInfo.shFracInnerHits = shFracInnerHits ;


  //====================================================
  // FiducialFlags, using nextToBoundary definition of gaps
  //====================================================

  reco::GsfElectron::FiducialFlags fiducialFlags ;
  int detector = seedXtalId.subdetId() ;
  double feta=fabs(scRef->position().eta()) ;
  if (detector==EcalBarrel)
   {
	fiducialFlags.isEB = true ;
	if (EBDetId::isNextToEtaBoundary(EBDetId(seedXtalId)))
	 {
	  if (fabs(feta-1.479)<.1)
	   { fiducialFlags.isEBEEGap = true ; }
	  else
	   { fiducialFlags.isEBEtaGap = true ; }
	 }
	if (EBDetId::isNextToPhiBoundary(EBDetId(seedXtalId)))
	 { fiducialFlags.isEBPhiGap = true ; }
   }
  else if (detector==EcalEndcap)
   {
	fiducialFlags.isEE = true ;
	if (EEDetId::isNextToRingBoundary(EEDetId(seedXtalId)))
	 {
	  if (fabs(feta-1.479)<.1)
	   { fiducialFlags.isEBEEGap = true ; }
	  else
	   { fiducialFlags.isEERingGap = true ; }
	 }
	if (EEDetId::isNextToDBoundary(EEDetId(seedXtalId)))
	 { fiducialFlags.isEEDeeGap = true ; }
   }
  else
   { edm::LogWarning("GsfElectronAlgo")<<"createElectron(): do not know if it is a barrel or endcap seed cluster !!!!" ; }


  //====================================================
  // ShowerShape
  //====================================================

  reco::GsfElectron::ShowerShape showerShape ;
  const CaloTopology * topology = theCaloTopo.product() ;
  const CaloGeometry * geometry = theCaloGeom.product() ;
  const EcalRecHitCollection * reducedRecHits = 0 ;
  if (fiducialFlags.isEB)
   { reducedRecHits = reducedEBRecHits.product() ; }
  else
   { reducedRecHits = reducedEERecHits.product() ; }
  std::vector<float> covariances = EcalClusterTools::covariances(seedCluster,reducedRecHits,topology,geometry) ;
  std::vector<float> localCovariances = EcalClusterTools::localCovariances(seedCluster,reducedRecHits,topology) ;
  showerShape.sigmaEtaEta = sqrt(covariances[0]) ;
  showerShape.sigmaIetaIeta = sqrt(localCovariances[0]) ;
  showerShape.e1x5 = EcalClusterTools::e1x5(seedCluster,reducedRecHits,topology)  ;
  showerShape.e2x5Max = EcalClusterTools::e2x5Max(seedCluster,reducedRecHits,topology)  ;
  showerShape.e5x5 = EcalClusterTools::e5x5(seedCluster,reducedRecHits,topology) ;
  showerShape.hcalDepth1OverEcal = HoE1 ;
  showerShape.hcalDepth2OverEcal = HoE2 ;


  //====================================================
  // brems fraction
  //====================================================

  float fbrem = (outMom.mag()>0.)?((innMom.mag()-outMom.mag())/innMom.mag()):1.e30 ;


  //====================================================
  // Go !
  //====================================================

  GsfElectron * ele = new
	GsfElectron
	 ( momentum,charge,chargeInfo,coreRef,
	   tcMatching, tkExtra, ctfInfo,
	   fiducialFlags,showerShape,
	   fbrem,mva) ;

  // set corrections + classification
  ElectronClassification theClassifier;
  theClassifier.correct(*ele);
  // energy corrections only for ecalDriven electrons
  if (ele->core()->ecalDrivenSeed()) {
    ElectronEnergyCorrector theEnCorrector(superClusterErrorFunction_);
    theEnCorrector.correct(*ele, applyEtaCorrection_);
    ElectronMomentumCorrector theMomCorrector;
    theMomCorrector.correct(*ele,vtxTSOS_);
  }

  // now isolation variables
  reco::GsfElectron::IsolationVariables dr03, dr04 ;
  dr03.tkSumPt = tkIso03.getPtTracks(ele);
  dr03.hcalDepth1TowerSumEt = had1Iso03.getTowerEtSum(ele);
  dr03.hcalDepth2TowerSumEt = had2Iso03.getTowerEtSum(ele);
  dr03.ecalRecHitSumEt = ecalBarrelIso03.getEtSum(ele)+ecalEndcapsIso03.getEtSum(ele);
  dr04.tkSumPt = tkIso04.getPtTracks(ele);
  dr04.hcalDepth1TowerSumEt = had1Iso04.getTowerEtSum(ele);
  dr04.hcalDepth2TowerSumEt = had2Iso04.getTowerEtSum(ele);
  dr04.ecalRecHitSumEt = ecalBarrelIso04.getEtSum(ele)+ecalEndcapsIso04.getEtSum(ele);
  ele->setIsolation03(dr03);
  ele->setIsolation04(dr04);

  outEle.push_back(ele) ;
 }


const CaloClusterPtr GsfElectronAlgo::getEleBasicCluster(const GsfTrackRef &t, const SuperCluster *scRef) {

    CaloClusterPtr eleRef;
    TrajectoryStateOnSurface tempTSOS;
    TrajectoryStateOnSurface outTSOS = mtsTransform_->outerStateOnSurface(*t);
    float dphimin = 1.e30;
    for (CaloCluster_iterator bc=scRef->clustersBegin(); bc!=scRef->clustersEnd(); bc++) {
      GlobalPoint posclu((*bc)->position().x(),(*bc)->position().y(),(*bc)->position().z());
      tempTSOS = mtsTransform_->extrapolatedState(outTSOS,posclu) ;
      if (!tempTSOS.isValid()) tempTSOS=outTSOS;
      GlobalPoint extrap = tempTSOS.globalPosition();
      float dphi = posclu.phi() - extrap.phi();
      if (fabs(dphi)>CLHEP::pi) dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
      if (fabs(dphi)<dphimin) {
        dphimin = fabs(dphi);
	eleRef = (*bc);
	eleTSOS_ = tempTSOS;
      }
    }
    return eleRef;

}

bool  GsfElectronAlgo::calculateTSOS(const GsfTrack &t,const SuperCluster & theClus, const BeamSpot & bs){

    //at innermost point
    innTSOS_ = mtsTransform_->innerStateOnSurface(t);
    if (!innTSOS_.isValid()) return false;

    //at vertex
    // innermost state propagation to the beam spot position
    vtxTSOS_ = mtsTransform_->extrapolatedState(innTSOS_,
	GlobalPoint(bs.position().x(),bs.position().y(),bs.position().z()));
    if (!vtxTSOS_.isValid()) vtxTSOS_=innTSOS_;

    //at seed
    outTSOS_ = mtsTransform_->outerStateOnSurface(t);
    if (!outTSOS_.isValid()) return false;

    //    TrajectoryStateOnSurface seedTSOS
    seedTSOS_ = mtsTransform_->extrapolatedState(outTSOS_,
						 GlobalPoint(theClus.seed()->position().x(),
							     theClus.seed()->position().y(),
							     theClus.seed()->position().z()));
    if (!seedTSOS_.isValid()) seedTSOS_=outTSOS_;

    //at scl
    sclTSOS_ = mtsTransform_->extrapolatedState(innTSOS_,GlobalPoint(theClus.x(),theClus.y(),theClus.z()));
    if (!sclTSOS_.isValid()) sclTSOS_=outTSOS_;

    // constrained momentum
    constrainedVtxTSOS_ = constraintAtVtx_->constrainAtBeamSpot(t,bs);

    return true;
}


//=======================================================================================
// Ambiguity solving
//=======================================================================================

bool better_electron( const reco::GsfElectron * e1, const reco::GsfElectron * e2 )
 { return (fabs(e1->eSuperClusterOverP()-1)<fabs(e2->eSuperClusterOverP()-1)) ; }

//void GsfElectronAlgo::resolveElectrons( GsfElectronPtrCollection & inEle, reco::GsfElectronCollection & outEle)
void GsfElectronAlgo::resolveElectrons( GsfElectronPtrCollection & inEle, reco::GsfElectronCollection & outEle,
       edm::Handle<EcalRecHitCollection> & reducedEBRecHits,
       edm::Handle<EcalRecHitCollection> & reducedEERecHits )
 {
  GsfElectronPtrCollection::iterator e1, e2 ;
  if (ambSortingStrategy_==0)
   { inEle.sort(EgAmbiguityTools::isBetter) ; }
  else if (ambSortingStrategy_==1)
   { inEle.sort(EgAmbiguityTools::isInnerMost(trackerHandle_)) ; }
  else
   { edm::LogError("GsfElectronAlgo")<<"unknown ambSortingStrategy "<<ambSortingStrategy_ ; }

  // resolve when e/g SC is found
  for( e1 = inEle.begin() ;  e1 != inEle.end() ; ++e1 )
   {
    SuperClusterRef scRef1 = (*e1)->superCluster();
    CaloClusterPtr eleClu1 = (*e1)->electronCluster();
    LogDebug("GsfElectronAlgo")
      << "Blessing electron with E/P " << (*e1)->eSuperClusterOverP()
      << ", cluster " << scRef1.get()
      << " & track " << (*e1)->gsfTrack().get() ;

    for( e2 = e1, ++e2 ;  e2 != inEle.end() ; )
     {
      SuperClusterRef scRef2 = (*e2)->superCluster();
      CaloClusterPtr eleClu2 = (*e2)->electronCluster();

      // search if same cluster
      bool sameCluster = false ;
      if (ambClustersOverlapStrategy_==0)
       { sameCluster = (scRef1==scRef2) ; }
      else if (ambClustersOverlapStrategy_==1)
       {
    	sameCluster =
         ( (EgAmbiguityTools::sharedEnergy(&(*eleClu1),&(*eleClu2),reducedEBRecHits,reducedEERecHits)>=1.*cosh(scRef1->eta())) ||
    	   (EgAmbiguityTools::sharedEnergy(&(*scRef1->seed()),&(*eleClu2),reducedEBRecHits,reducedEERecHits)>=1.*cosh(scRef1->eta())) ||
    	   (EgAmbiguityTools::sharedEnergy(&(*eleClu1),&(*scRef2->seed()),reducedEBRecHits,reducedEERecHits)>=1.*cosh(scRef1->eta())) ||
    	   (EgAmbiguityTools::sharedEnergy(&(*scRef1->seed()),&(*scRef2->seed()),reducedEBRecHits,reducedEERecHits)>=1.*cosh(scRef1->eta())) ) ;
       }
      else
       { edm::LogError("GsfElectronAlgo")<<"unknown ambClustersOverlapStrategy_ "<<ambClustersOverlapStrategy_ ; }

      // main instructions
      if (sameCluster)
       {
        LogDebug("GsfElectronAlgo")
          << "Discarding electron with E/P " << (*e2)->eSuperClusterOverP()
          << ", cluster " << scRef2.get()
          << " and track " << (*e2)->gsfTrack().get() ;
        (*e1)->addAmbiguousGsfTrack((*e2)->gsfTrack()) ;
        e2 = inEle.erase(e2) ;
       }
      else if ((*e1)->gsfTrack()==(*e2)->gsfTrack())
       {
        LogDebug("GsfElectronAlgo")
          << "Forgetting electron with E/P " << (*e2)->eSuperClusterOverP()
          << ", cluster " << scRef2.get()
          << " and track " << (*e2)->gsfTrack().get() ;
        e2 = inEle.erase(e2) ;
       }
      else
       { ++e2 ; }
     }
    outEle.push_back(**e1) ;
   }

 }


//=======================================================================================
// Code from Puneeth Kalavase
//=======================================================================================

pair<TrackRef,float> GsfElectronAlgo::getCtfTrackRef(const GsfTrackRef& gsfTrackRef, edm::Handle<reco::TrackCollection> ctfTracksH ) {

  float maxFracShared = 0;
  TrackRef ctfTrackRef = TrackRef() ;
  const TrackCollection * ctfTrackCollection = ctfTracksH.product() ;

  // get the Hit Pattern for the gsfTrack
  const HitPattern& gsfHitPattern = gsfTrackRef->hitPattern();

  unsigned int counter ;
  TrackCollection::const_iterator ctfTkIter ;
  for ( ctfTkIter = ctfTrackCollection->begin() , counter = 0 ;
        ctfTkIter != ctfTrackCollection->end() ; ctfTkIter++, counter++ )
   {

    double dEta = gsfTrackRef->eta() - ctfTkIter->eta();
    double dPhi = gsfTrackRef->phi() - ctfTkIter->phi();
    double pi = acos(-1.);
    if(fabs(dPhi) > pi) dPhi = 2*pi - fabs(dPhi);

    // dont want to look at every single track in the event!
    if(sqrt(dEta*dEta + dPhi*dPhi) > 0.3) continue;

    unsigned int shared = 0 ;
    int gsfHitCounter = 0 ;
    int numGsfInnerHits = 0 ;
    int numCtfInnerHits = 0 ;
    // get the CTF Track Hit Pattern
    const HitPattern& ctfHitPattern = ctfTkIter->hitPattern() ;

    trackingRecHit_iterator elHitsIt ;
    for ( elHitsIt = gsfTrackRef->recHitsBegin() ;
          elHitsIt != gsfTrackRef->recHitsEnd() ;
          elHitsIt++, gsfHitCounter++ )
     {
      if(!((**elHitsIt).isValid()))  //count only valid Hits
       { continue ; }

      // look only in the pixels/TIB/TID
      uint32_t gsfHit = gsfHitPattern.getHitPattern(gsfHitCounter) ;
      if (!(gsfHitPattern.pixelHitFilter(gsfHit) ||
	        gsfHitPattern.stripTIBHitFilter(gsfHit) ||
	        gsfHitPattern.stripTIDHitFilter(gsfHit) ) )
       { continue ; }

      numGsfInnerHits++ ;

      int ctfHitsCounter = 0 ;
      numCtfInnerHits = 0 ;
      trackingRecHit_iterator ctfHitsIt ;
      for ( ctfHitsIt = ctfTkIter->recHitsBegin() ;
            ctfHitsIt != ctfTkIter->recHitsEnd() ;
            ctfHitsIt++, ctfHitsCounter++ )
       {
        if(!((**ctfHitsIt).isValid())) //count only valid Hits!
         { continue ; }

	    uint32_t ctfHit = ctfHitPattern.getHitPattern(ctfHitsCounter);
	    if( !(ctfHitPattern.pixelHitFilter(ctfHit) ||
	          ctfHitPattern.stripTIBHitFilter(ctfHit) ||
	          ctfHitPattern.stripTIDHitFilter(ctfHit) ) )
	     { continue ; }

	    numCtfInnerHits++ ;

        if( (**elHitsIt).sharesInput(&(**ctfHitsIt),TrackingRecHit::all) )
         {
          shared++ ;
          break ;
         }

       } //ctfHits iterator

     } //gsfHits iterator

    if ((numGsfInnerHits==0)||(numCtfInnerHits==0))
     { continue ; }

    if ( static_cast<float>(shared)/min(numGsfInnerHits,numCtfInnerHits) > maxFracShared )
     {
      maxFracShared = static_cast<float>(shared)/min(numGsfInnerHits, numCtfInnerHits);
      ctfTrackRef = TrackRef(ctfTracksH,counter);
     }

   } //ctfTrack iterator

  return make_pair(ctfTrackRef,maxFracShared) ;
 }

void GsfElectronAlgo::computeCharge
 ( const GsfTrackRef & gsf, const TrackRef & ctf, const SuperClusterRef & sc, const BeamSpot & bs,
   int & charge, GsfElectron::ChargeInfo & info )
 {
  // determine charge from SC
  GlobalPoint orig(bs.position().x(), bs.position().y(), bs.position().z()) ;
  GlobalPoint scpos(sc->position().x(), sc->position().y(), sc->position().z()) ;
  GlobalVector scvect(scpos-orig) ;
  GlobalPoint inntkpos = innTSOS_.globalPosition() ;
  GlobalVector inntkvect = GlobalVector(inntkpos-orig) ;
  float dPhiInnEle=normalized_dphi(scvect.phi()-inntkvect.phi()) ;
  if(dPhiInnEle>0) info.scPixCharge = -1 ;
  else info.scPixCharge = 1 ;

  // flags
  int chargeGsf = gsf->charge() ;
  info.isGsfScPixConsistent = ((chargeGsf*info.scPixCharge)>0) ;
  info.isGsfCtfConsistent = (ctf.isNonnull()&&((chargeGsf*ctf->charge())>0)) ;
  info.isGsfCtfScPixConsistent = (info.isGsfScPixConsistent&&info.isGsfCtfConsistent) ;

  // default charge
  if (info.isGsfScPixConsistent||ctf.isNull())
   { charge = info.scPixCharge ; }
  else
   { charge = ctf->charge() ; }
 }

void GsfElectronAlgo::checkPfTranslatorParameters( edm::ParameterSetID const & psetid )
 {
  edm::ParameterSet pset ;
  edm::pset::Registry::instance()->getMapped(psetid,pset) ;
  edm::ParameterSet mvaBlock = pset.getParameter<edm::ParameterSet>("MVACutBlock") ;
  double pfTranslatorMinMva = mvaBlock.getParameter<double>("MVACut") ;
  double pfTranslatorUndefined = -99. ;
  if (minMVAPflow_<pfTranslatorMinMva)
   {
    // For pure tracker seeded electrons, if MVA is under translatorMinMva, there is no supercluster
    // of any kind available, so GsfElectronCoreProducer has already discarded the electron.
    edm::LogWarning("GsfElectronAlgo")
      <<"Parameter minMVAPflow will have no effect on purely tracker seeded electrons."
      <<" It is inferior to the cut already applied by PFlow translator." ;
   }
  if (minMVA_<pfTranslatorMinMva)
   {
    // For ecal seeded electrons, there is a cluster and GsfElectronCoreProducer has kept all electrons,
    // but when MVA is under translatorMinMva, the translator has not stored the supercluster and
    // forced the MVA value to translatorUndefined
    if (minMVA_<pfTranslatorUndefined)
     {
      edm::LogWarning("GsfElectronAlgo")
        <<"Parameter minMVA is inferior to the cut applied by PFlow translator."
        <<" Some ecal (and eventually tracker) seeded electrons may lack their MVA value and PFlow supercluster." ;
     }
    else
     {
      // the MVA value has been forced to translatorUndefined, inferior minMVAPflow
      // so the cut actually applied is the PFlow one
      edm::LogWarning("GsfElectronAlgo")
        <<"Parameter minMVA will have no effect on ecal (and eventually tracker) seeded electrons."
        <<" It is inferior to the cut already defined by PFlow translator." ;
     }
   }
 }

void GsfElectronAlgo::checkEcalSeedingParameters( edm::ParameterSetID const & psetid )
 {
  edm::ParameterSet pset ;
  edm::pset::Registry::instance()->getMapped(psetid,pset) ;
  edm::ParameterSet seedConfiguration = pset.getParameter<edm::ParameterSet>("SeedConfiguration") ;
//  edm::ParameterSet seedParameters = seedConfiguration.getParameter<edm::ParameterSet>("ecalDrivenElectronSeedsParameters") ;
//
//  if (seedParameters.getParameter<bool>("applyHOverECut"))
//   {
//    if (hOverEConeSize_!=seedParameters.getParameter<double>("hOverEConeSize"))
//     { edm::LogWarning("GsfElectronAlgo") <<"The H/E cone size is different from ecal seeding." ; }
//    if (maxHOverEBarrel_<seedParameters.getParameter<double>("maxHOverEBarrel"))
//     { edm::LogWarning("GsfElectronAlgo") <<"The max barrel H/E is lower than during ecal seeding." ; }
//    if (maxHOverEEndcaps_<seedParameters.getParameter<double>("maxHOverEEndcaps"))
//     { edm::LogWarning("GsfElectronAlgo") <<"The max endcaps H/E is lower than during ecal seeding." ; }
//   }
//
//  if (minSCEtBarrel_<seedParameters.getParameter<double>("SCEtCut"))
//   { edm::LogWarning("GsfElectronAlgo") <<"The minimum super-cluster Et in barrel is lower than during ecal seeding." ; }
//  if (minSCEtEndcaps_<seedParameters.getParameter<double>("SCEtCut"))
//   { edm::LogWarning("GsfElectronAlgo") <<"The minimum super-cluster Et in endcaps is lower than during ecal seeding." ; }
 }
