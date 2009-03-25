#ifndef PhotonIsolationCalculator_H
#define PhotonIsolationCalculator_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"


#include <string>

class PhotonIsolationCalculator {

public:

  PhotonIsolationCalculator(){};

  virtual ~PhotonIsolationCalculator(){};

  void setup(const edm::ParameterSet& conf);

  void calculate(const reco::Photon*, 
		 const edm::Event&, const edm::EventSetup& es,
		 reco::Photon::FiducialFlags& phofid, 
		 reco::Photon::IsolationVariables& phoisolR03, 
		 reco::Photon::IsolationVariables& phoisolR04 );

  void classify(const reco::Photon* photon, 
		bool &isEBPho,
		bool &isEEPho,
		bool &isEBGap,
		bool &isEEGap,
		bool &isEBEEGap);
  void calculateTrackIso(const reco::Photon* photon,
			 const edm::Event &e,
			 double &trkCone,
			 int &ntrkCone,
			 double pTThresh=0,
			 double RCone=.4,
			 double RinnerCone=.1);


  double calculateEcalRecHitIso(const reco::Photon* photon,
				const edm::Event& iEvent,
				const edm::EventSetup& iSetup,
				double RCone,
				double RConeInner,
                                double etaSlice,
				double eMin,
				double etMin);
  double calculateHcalTowerIso(const reco::Photon* photon,
			       const edm::Event& iEvent,
			       const edm::EventSetup& iSetup,
			       double RCone,
			       double RConeInner,
			       double eMin,
                               signed int depth);


  
 protected:


  std::string barrelecalCollection_;
  std::string barrelecalProducer_;
  std::string endcapecalCollection_;
  std::string endcapecalProducer_;
  std::string hcalCollection_;
  std::string hcalProducer_;

  edm::InputTag trackInputTag_;
  //  edm::InputTag gsfRecoInputTag_;

  double modulePhiBoundary_;
  std::vector<double> moduleEtaBoundary_;

  std::vector<double>  trkIsoBarrelRadiusA_;
  std::vector<double>  ecalIsoBarrelRadiusA_;
  std::vector<double>  hcalIsoBarrelRadiusA_;
  std::vector<double>  trkIsoBarrelRadiusB_;
  std::vector<double>  ecalIsoBarrelRadiusB_;
  std::vector<double>  hcalIsoBarrelRadiusB_;

  std::vector<double>  trkIsoEndcapRadiusA_;
  std::vector<double>  ecalIsoEndcapRadiusA_;
  std::vector<double>  hcalIsoEndcapRadiusA_;
  std::vector<double>  trkIsoEndcapRadiusB_;
  std::vector<double>  ecalIsoEndcapRadiusB_;
  std::vector<double>  hcalIsoEndcapRadiusB_;

  

  //Isolation parameters variables
  double photonEcalRecHitConeInnerRadiusA_;
  double photonEcalRecHitConeOuterRadiusA_;
  double photonEcalRecHitEtaSliceA_;
  double photonEcalRecHitThreshEA_;
  double photonEcalRecHitThreshEtA_;
  double photonHcalTowerConeInnerRadiusA_;
  double photonHcalTowerConeOuterRadiusA_;
  double photonHcalTowerThreshEA_;
  double photonHcalDepth1TowerConeInnerRadiusA_;
  double photonHcalDepth1TowerConeOuterRadiusA_;
  double photonHcalDepth1TowerThreshEA_;
  double photonHcalDepth2TowerConeInnerRadiusA_;
  double photonHcalDepth2TowerConeOuterRadiusA_;
  double photonHcalDepth2TowerThreshEA_;
  double trackConeOuterRadiusA_;
  double trackConeInnerRadiusA_;
  double isolationtrackThresholdA_;

  double photonEcalRecHitConeInnerRadiusB_;
  double photonEcalRecHitConeOuterRadiusB_;
  double photonEcalRecHitEtaSliceB_;
  double photonEcalRecHitThreshEB_;
  double photonEcalRecHitThreshEtB_;
  double photonHcalTowerConeInnerRadiusB_;
  double photonHcalTowerConeOuterRadiusB_;
  double photonHcalTowerThreshEB_;
  double photonHcalDepth1TowerConeInnerRadiusB_;
  double photonHcalDepth1TowerConeOuterRadiusB_;
  double photonHcalDepth1TowerThreshEB_;
  double photonHcalDepth2TowerConeInnerRadiusB_;
  double photonHcalDepth2TowerConeOuterRadiusB_;
  double photonHcalDepth2TowerThreshEB_;
  double trackConeOuterRadiusB_;
  double trackConeInnerRadiusB_;
  double isolationtrackThresholdB_;



  };

#endif // PhotonIsolationCalculator_H
