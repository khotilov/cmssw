#ifndef TauAnalysis_Core_ElectronHistManager_h  
#define TauAnalysis_Core_ElectronHistManager_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "TauAnalysis/Core/interface/HistManagerBase.h"

#include "TauAnalysis/RecoTools/interface/ParticlePFIsolationExtractor.h"

#include <vector>
#include <string>

class ElectronHistManager : public HistManagerBase 
{
 public:
  explicit ElectronHistManager(const edm::ParameterSet&);
  ~ElectronHistManager();
  
 private:
//--- histogram booking and filling functions 
//    inherited from HistManagerBase class
  void bookHistogramsImp();
  void fillHistogramsImp(const edm::Event&, const edm::EventSetup&, double);

//--- auxiliary functions
  void bookElectronHistograms(MonitorElement*&, MonitorElement*&, MonitorElement*&, const char*);
  void bookElectronIsoConeSizeDepHistograms();

  double getElectronWeight(const pat::Electron&);

  void fillElectronHistograms(const pat::Electron&, MonitorElement*, MonitorElement*, MonitorElement*, double);
  void fillElectronIsoHistograms(const pat::Electron&, const reco::PFCandidateCollection&, double);
  void fillElectronIsoConeSizeDepHistograms(const pat::Electron&, double);
  void fillConversionHistograms(const pat::Electron&, edm::Handle<reco::TrackCollection>&, double, double );

//--- configuration parameters
  edm::InputTag electronSrc_;
  edm::InputTag vertexSrc_;
  edm::InputTag jetSrc_;
  edm::InputTag pfCandidateSrc_;
  edm::InputTag genParticleSrc_;
  edm::InputTag conversionPartnerTrackSrc_;
  edm::InputTag dcsTag_;

  bool requireGenElectronMatch_;

  typedef std::vector<int> vint;
  vint skipPdgIdsGenParticleMatch_;

  PATElectronPFIsolationExtractor* pfCombIsoExtractor_;
  PATElectronPFIsolationExtractor* pfChargedHadronIsoExtractor_;
  PATElectronPFIsolationExtractor* pfNeutralHadronIsoExtractor_;
  PATElectronPFIsolationExtractor* pfPhotonIsoExtractor_;

  bool makeIsoPtCtrlHistograms_;
  bool makeIsoPtConeSizeDepHistograms_;
  bool makeConversionHistograms_;
  bool isData_;

  unsigned numElectronIsoConeSizes_;
  double electronIsoConeSizeIncr_;
  unsigned numElectronIsoPtThresholds_;
  double electronIsoPtThresholdIncr_;
  double electronEtaMaxBarrel_;
  double electronEtaMinEndcap_;

//--- histograms
  MonitorElement* hNumElectrons_; 

  MonitorElement* hElectronPt_; 
  MonitorElement* hElectronEta_;
  MonitorElement* hElectronPtVsEta_;
  MonitorElement* hElectronPhi_;
  MonitorElement* hElectronCharge_;

  MonitorElement* hElectronWeightPosLog_;
  MonitorElement* hElectronWeightNegLog_;
  MonitorElement* hElectronWeightZero_;
  MonitorElement* hElectronWeightLinear_;

  MonitorElement* hGenElectronPt_;
  MonitorElement* hGenElectronEta_;
  MonitorElement* hGenElectronPhi_;

  MonitorElement* hElectronEnCompToGen_;
  MonitorElement* hElectronThetaCompToGen_;
  MonitorElement* hElectronPhiCompToGen_;

  MonitorElement* hElectronMatchingGenParticlePdgId_;
  MonitorElement* hElectronMatchingFinalStateGenParticlePdgId_;

  MonitorElement* hElectronTrackPt_;
  MonitorElement* hElectronTrackIPxy_;
  MonitorElement* hElectronTrackIPz_;

  MonitorElement* hElectronSuperclEnOverTrackMomBarrel_;
  MonitorElement* hElectronSuperclEnOverTrackMomEndcap_;
  MonitorElement* hElectronSeedClEnOverTrackMomBarrel_;
  MonitorElement* hElectronSeedClEnOverTrackMomEndcap_;
  MonitorElement* hElectronSuperclShapeSigmaEtaEta_;
  MonitorElement* hElectronBremsFraction_;
  MonitorElement* hElectronHadEnOverEmEnBarrel_;
  MonitorElement* hElectronSuperclShapeSigmaIetaIetaBarrel_;
  MonitorElement* hElectronDeltaPhiSuperclToExtrapolTrackBarrel_;
  MonitorElement* hElectronDeltaEtaSuperclToExtrapolTrackBarrel_;
  MonitorElement* hElectronHadEnOverEmEnEndcap_;
  MonitorElement* hElectronSuperclShapeSigmaIetaIetaEndcap_;
  MonitorElement* hElectronDeltaPhiSuperclToExtrapolTrackEndcap_;
  MonitorElement* hElectronDeltaEtaSuperclToExtrapolTrackEndcap_;

//--- IsoDeposits reconstructed from ECAL and HCAL recHits/CaloTowers and reco::Tracks
  MonitorElement* hElectronTrkIsoPt_;
  MonitorElement* hElectronTrkIsoPtRel_;
  MonitorElement* hElectronTrkIsoPtVsElectronPt_;
  MonitorElement* hElectronEcalIsoPt_;
  MonitorElement* hElectronEcalIsoPtRel_;
  MonitorElement* hElectronEcalIsoPtBarrel_;
  MonitorElement* hElectronEcalIsoPtEndcap_;
  MonitorElement* hElectronEcalIsoPtBarrelRel_;
  MonitorElement* hElectronEcalIsoPtEndcapRel_;
  MonitorElement* hElectronEcalIsoPtVsElectronPt_;
  MonitorElement* hElectronHcalIsoPt_;
  MonitorElement* hElectronHcalIsoPtRel_;
  MonitorElement* hElectronIsoSumPt_;
  MonitorElement* hElectronIsoSumPtRel_;
  MonitorElement* hElectronIsoSumPtVsElectronPt_;
  
  MonitorElement* hElectronDeltaRnearestJet_;

  MonitorElement* hElectronTrkIsoValProfile_;
  MonitorElement* hElectronTrkIsoEtaDistProfile_;
  MonitorElement* hElectronTrkIsoPhiDistProfile_;

  MonitorElement* hElectronEcalIsoValProfile_;
  MonitorElement* hElectronEcalIsoEtaDistProfile_;
  MonitorElement* hElectronEcalIsoPhiDistProfile_;

  MonitorElement* hElectronHcalIsoValProfile_;
  MonitorElement* hElectronHcalIsoEtaDistProfile_;
  MonitorElement* hElectronHcalIsoPhiDistProfile_;

  reco::isodeposit::AbsVetos electronTrkIsoParam_;
  reco::isodeposit::AbsVetos electronEcalIsoParam_;
  reco::isodeposit::AbsVetos electronHcalIsoParam_;

  std::vector<MonitorElement*> hElectronTrkIsoPtConeSizeDep_;
  std::vector<MonitorElement*> hElectronEcalIsoPtConeSizeDep_;
  std::vector<MonitorElement*> hElectronHcalIsoPtConeSizeDep_;

//--- IsoDeposits reconstructed from Partcile Flow
  MonitorElement* hElectronParticleFlowIsoPt_;
  MonitorElement* hElectronParticleFlowIsoPtRel_;
  MonitorElement* hElectronPFChargedHadronIsoPt_;
  MonitorElement* hElectronPFChargedHadronIsoPtRel_;
  MonitorElement* hElectronPFNeutralHadronIsoPt_;
  MonitorElement* hElectronPFNeutralHadronIsoPtRel_;
  MonitorElement* hElectronPFGammaIsoPt_;
  MonitorElement* hElectronPFGammaIsoPtRel_;
  
  MonitorElement* hElectronPFCombIsoPt_;
  MonitorElement* hElectronPFCombIsoPtBarrel_;
  MonitorElement* hElectronPFCombIsoPtEndcap_;
  MonitorElement* hElectronPFCombIsoPtRel_;
  MonitorElement* hElectronPFCombIsoPtRelBarrel_;
  MonitorElement* hElectronPFCombIsoPtRelEndcap_;

  MonitorElement* hElectronPFChargedHadronIsoPtCtrl_;
  MonitorElement* hElectronPFGammaIsoPtCtrl_;

  MonitorElement* hElectronConvDoca_;
  MonitorElement* hElectronConvDeltaCotTheta_;
  MonitorElement* hElectronConvRadius_;
  MonitorElement* hElectronMissExpInnerHits_;

  std::vector<MonitorElement*> hElectronParticleFlowIsoPtConeSizeDep_;
  std::vector<MonitorElement*> hElectronPFChargedHadronIsoPtConeSizeDep_;
  std::vector<MonitorElement*> hElectronPFNeutralHadronIsoPtConeSizeDep_;
  std::vector<MonitorElement*> hElectronPFGammaIsoPtConeSizeDep_;

  reco::isodeposit::AbsVetos electronParticleFlowIsoParam_;
};

#endif  


