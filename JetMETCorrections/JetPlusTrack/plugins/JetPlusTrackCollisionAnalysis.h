#ifndef JetPlusTrackCollisionAnalysis_h
#define JetPlusTrackCollisionAnalysis_h

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <string>
#include <memory>
#include <map>

namespace edm {
  class Event;
  class EventSetup;
  class ParameterSet;
}
#include "TFile.h"
#include "TTree.h"

///
/// jet energy corrections from MCjet calibration
///
namespace cms
{

class JetPlusTrackCollisionAnalysis : public edm::EDAnalyzer
{
public:  

  JetPlusTrackCollisionAnalysis(const edm::ParameterSet& fParameters);

  virtual ~JetPlusTrackCollisionAnalysis();

      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void beginJob(const edm::EventSetup& ) ;
      virtual void endJob() ;
   
private:
// Histograms/Tree
     std::string fOutputFileName ;
     bool allowMissingInputs_;
     TFile*      hOutputFile ;
     TTree*      myTree;
     int NumRecoJetsCaloTower, NumRecoJetsZSPCorrected, NumRecoJetsRecHit, NumRecoJetsJPTCorrected;
     float JetRecoEtCaloTower[10],JetRecoEtaCaloTower[10],JetRecoPhiCaloTower[10];
     float JetRecoEtZSPCorrected[10],JetRecoEtaZSPCorrected[10],JetRecoPhiZSPCorrected[10];
     float JetRecoEtRecHit[10],EcalEnergyCone[10],HcalEnergyConeZSP[10],HcalEnergyConeNZSP[10];
     float JetRecoEtJPTCorrected[10],JetRecoEtaJPTCorrected[10],JetRecoPhiJPTCorrected[10];
     int NumRecoTrack,NumRecoCone;
     float TrackRecoEt[5000],TrackRecoEta[5000],TrackRecoPhi[5000]; 
     int run, event;

  double mCone;
  edm::InputTag mInputCaloTower;
  edm::InputTag mInputJetsCaloTower;
  edm::InputTag mInputJetsZSPCorrected; 
  edm::InputTag mInputJetsJPTCorrected; 

  std::string m_inputTrackLabel;
  std::vector<edm::InputTag> ecalLabels_;
  edm::InputTag ebrechit;
  edm::InputTag eerechit;
  edm::InputTag hbhelabel_;
  edm::InputTag hbhelabelNZS_;
};
}
#endif
