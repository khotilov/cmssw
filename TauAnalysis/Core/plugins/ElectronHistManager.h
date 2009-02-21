#ifndef TauAnalysis_Core_ElectronHistManager_h  
#define TauAnalysis_Core_ElectronHistManager_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TauAnalysis/Core/interface/HistManagerBase.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include <vector>
#include <string>

class ElectronHistManager : public HistManagerBase 
{
 public:
  
  explicit ElectronHistManager(const edm::ParameterSet&);
  ~ElectronHistManager();
  
  void bookHistograms(const edm::EventSetup&);
  void fillHistograms(const edm::Event&, const edm::EventSetup&);

 private:

//--- private functions
  void bookElectronHistograms(DQMStore&, MonitorElement*&, MonitorElement*&, MonitorElement*&, const char*);
  
  void fillElectronHistograms(const pat::Electron&, MonitorElement*, MonitorElement*, MonitorElement*);
  void fillElectronIsoHistograms(const pat::Electron&);
  void fillElectronIsoConeSizeDepHistograms(const pat::Electron&);

//--- configuration parameters
  edm::InputTag electronSrc_;
  edm::InputTag vertexSrc_;

  std::string dqmDirectory_store_;

  bool requireGenElectronMatch_;

  unsigned numElectronIsoConeSizes_;
  double electronIsoConeSizeIncr_;
  unsigned numElectronIsoPtThresholds_;
  double electronIsoPtThresholdIncr_;
  double electronEtaMaxBarrel_;
  double electronEtaMinEndcap_;

//--- histograms
  MonitorElement* hElectronPt_; 
  MonitorElement* hElectronEta_;
  MonitorElement* hElectronPtVsEta_;
  MonitorElement* hElectronPhi_;

  MonitorElement* hElectronEnCompToGen_;
  MonitorElement* hElectronThetaCompToGen_;
  MonitorElement* hElectronPhiCompToGen_;

  MonitorElement* hElectronTrackPt_;
  MonitorElement* hElectronTrackIPxy_;
  MonitorElement* hElectronTrackIPz_;

  MonitorElement* hElectronSuperclEnOverTrackMomBarrel_;
  MonitorElement* hElectronSuperclEnOverTrackMomEndcap_;

  MonitorElement* hElectronIdRobust_;
  
  MonitorElement* hElectronTrkIsoPt_;
  MonitorElement* hElectronEcalIsoPt_;
  MonitorElement* hElectronEcalIsoPtBarrel_;
  MonitorElement* hElectronEcalIsoPtEndcap_;
  MonitorElement* hElectronHcalIsoPt_;
  MonitorElement* hElectronIsoSumPt_;
  
  MonitorElement* hElectronTrkIsoValProfile_;
  MonitorElement* hElectronTrkIsoEtaDistProfile_;
  MonitorElement* hElectronTrkIsoPhiDistProfile_;

  MonitorElement* hElectronEcalIsoValProfile_;
  MonitorElement* hElectronEcalIsoEtaDistProfile_;
  MonitorElement* hElectronEcalIsoPhiDistProfile_;

  MonitorElement* hElectronHcalIsoValProfile_;
  MonitorElement* hElectronHcalIsoEtaDistProfile_;
  MonitorElement* hElectronHcalIsoPhiDistProfile_;

  std::vector<MonitorElement*> hElectronTrkIsoPtConeSizeDep_;
  std::vector<MonitorElement*> hElectronEcalIsoPtConeSizeDep_;
  std::vector<MonitorElement*> hElectronHcalIsoPtConeSizeDep_;
};

#endif  


