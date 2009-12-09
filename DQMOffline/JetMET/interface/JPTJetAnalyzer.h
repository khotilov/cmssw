#ifndef JPTJetAnalyzer_H
#define JPTJetAnalyzer_H

/** \class JPTJetAnalyzer
 *
 *  DQM monitoring source for JPT Jets
 *
 *  $Date: 2009/12/04 11:00:25 $
 *  $Revision: 1.4 $
 *  \author N. Cripps - Imperial
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMOffline/JetMET/interface/JetAnalyzerBase.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DQMServices/Core/interface/MonitorElement.h"
// forward declare classes which do not need to be defined for interface
class DQMStore;
namespace reco {
  class CaloJet;
}
namespace jptJetAnalysis {
  class TrackPropagatorToCalo;
  class StripSignalOverNoiseCalculator;
}
namespace jpt {
  class MatchedTracks;
}
class JetPlusTrackCorrector;
class JetCorrector;
class TrackingRecHit;
class SiStripRecHit2D;


/// JPT jet analyzer class definition
class JPTJetAnalyzer : public JetAnalyzerBase {
 public:
  /// Constructor
  JPTJetAnalyzer(const edm::ParameterSet& config);
  
  /// Destructor
  virtual ~JPTJetAnalyzer();
  
  /// Inizialize parameters for histo binning
  void beginJob(const edm::EventSetup& eventSetup, DQMStore* dqmStore);
  
  /// Do the analysis
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup, const reco::CaloJet& rawJet, double& pt1, double& pt2, double& pt3);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup, const reco::CaloJetCollection& rawJets);
  
  /// Finish up a job
  virtual void endJob();
  
 private:
   
  // Helpper classes
  /// Helpper class to hold the configuration for a histogram
  struct HistogramConfig {
    bool enabled;
    unsigned int nBins;
    double min;
    double max;
    unsigned int nBinsY;
    double minY;
    double maxY;
    HistogramConfig();
    HistogramConfig(const unsigned int theNBins, const double theMin, const double theMax);
    HistogramConfig(const unsigned int theNBinsX, const double theMinX, const double theMaxX,
                    const unsigned int theNBinsY, const double theMinY, const double theMaxY);
  };
  /// Helpper class for grouping histograms belowing to a set of tracks
  struct TrackHistograms {
    MonitorElement* nTracksHisto;
    MonitorElement* ptHisto;
    MonitorElement* phiHisto;
    MonitorElement* etaHisto;
    MonitorElement* nHitsHisto;
    MonitorElement* nLayersHisto;
    MonitorElement* ptVsEtaHisto;
    MonitorElement* dzHisto;
    MonitorElement* dxyHisto;
    MonitorElement* trackDirectionJetDRHisto;
    MonitorElement* trackImpactPointJetDRHisto;
    TrackHistograms();
    TrackHistograms(MonitorElement* theNTracksHisto, MonitorElement* thePtHisto, MonitorElement* thePhiHisto, MonitorElement* theEtaHisto,
                    MonitorElement* theNHitsHisto, MonitorElement* theNLayersHisto, MonitorElement* thePtVsEtaHisto,
                    MonitorElement* dzHisto, MonitorElement* dxyHisto,
                    MonitorElement* theTrackDirectionJetDRHisto, MonitorElement* theTrackImpactPointJetDRHisto);
  };
  
  // Private methods
  /// Load the config for a hitogram
  void getConfigForHistogram(const std::string& configName, const edm::ParameterSet& psetContainingConfigPSet, std::ostringstream* pDebugStream = NULL);
  /// Load the configs for histograms associated with a set of tracks
  void getConfigForTrackHistograms(const std::string& tag, const edm::ParameterSet& psetContainingConfigPSet,std::ostringstream* pDebugStream = NULL);
  /// Book histograms and profiles
  MonitorElement* bookHistogram(const std::string& name, const std::string& title, const std::string& xAxisTitle, DQMStore* dqm);
  MonitorElement* book2DHistogram(const std::string& name, const std::string& title, const std::string& xAxisTitle, const std::string& yAxisTitle, DQMStore* dqm);
  MonitorElement* bookProfile(const std::string& name, const std::string& title, const std::string& xAxisTitle, const std::string& yAxisTitle, DQMStore* dqm);
  /// Book all histograms
  void bookHistograms(DQMStore* dqm);
  /// Book the histograms for a track
  void bookTrackHistograms(TrackHistograms* histos, const std::string& tag, const std::string& titleTag,
                           MonitorElement* trackDirectionJetDRHisto, MonitorElement* trackImpactPointJetDRHisto, DQMStore* dqm);
  /// Fill histogram or profile if it has been booked
  void fillHistogram(MonitorElement* histogram, const double value);
  void fillHistogram(MonitorElement* histogram, const double valueX, const double valueY);
  /// Fill all track histograms
  void fillTrackHistograms(TrackHistograms& allTracksHistos, TrackHistograms& inCaloInVertexHistos,
                           TrackHistograms& inCaloOutVertexHistos, TrackHistograms& outCaloInVertexHistos,
                           const jpt::MatchedTracks& tracks, const reco::CaloJet& rawJet);
  void fillTrackHistograms(TrackHistograms& histos, const reco::TrackRefVector& tracks, const reco::CaloJet& rawJet);
  /// Fill the SoN hisotgram for hits on tracks
  void fillSiStripSoNForTracks(const reco::TrackRefVector& tracks);
  void fillSiStripHitSoN(const TrackingRecHit& hit);
  void fillSiStripHitSoNForSingleHit(const SiStripRecHit2D& hit);
  /// Utility function to calculate the fraction of track Pt in cone
  static double findPtFractionInCone(const reco::TrackRefVector& inConeTracks, const reco::TrackRefVector& outOfConeTracks);
  
  /// String constant for message logger category
  static const char* messageLoggerCatregory;
  
  // Config
  /// Path of directory used to store histograms in DQMStore
  const std::string histogramPath_;
  /// Create verbose debug messages
  const bool verbose_;
  /// JPT corrector name
  const std::string jptCorrectorName_;
  /// ZSP corrector name
  const std::string zspCorrectorName_;
  /// Histogram configuration (nBins etc)
  std::map<std::string,HistogramConfig> histogramConfig_;
  
  /// JPT Corrector
  const JetPlusTrackCorrector* jptCorrector_;
  /// ZSP Corrector
  const JetCorrector* zspCorrector_;
  
  /// Write DQM store to a file?
  const bool writeDQMStore_;
  /// DQM store file name
  std::string dqmStoreFileName_;
  
  /// Helpper object to propagate tracks to the calo surface
  jptJetAnalysis::TrackPropagatorToCalo* trackPropagator_;
  /// Helpper object to calculate strip SoN for tracks
  jptJetAnalysis::StripSignalOverNoiseCalculator* sOverNCalculator_;
  
  // Histograms
  MonitorElement *JetE_, *JetEt_, *JetP_, *JetMass_, *JetPt_;
  MonitorElement *JetPt1_, *JetPt2_, *JetPt3_;
  MonitorElement *JetPx_, *JetPy_, *JetPz_;
  MonitorElement *JetEta_, *JetPhi_, *JetDeltaEta_, *JetDeltaPhi_, *JetPhiVsEta_;
  MonitorElement *TrackSiStripHitStoNHisto_;
  MonitorElement *InCaloTrackDirectionJetDRHisto_, *OutCaloTrackDirectionJetDRHisto_;
  MonitorElement *InVertexTrackImpactPointJetDRHisto_, *OutVertexTrackImpactPointJetDRHisto_;
  MonitorElement *NTracksPerJetHisto_, *NTracksPerJetVsJetEtHisto_, *NTracksPerJetVsJetEtaHisto_;
  MonitorElement *PtFractionInConeHisto_, *PtFractionInConeVsJetRawEtHisto_, *PtFractionInConeVsJetEtaHisto_;
  MonitorElement *CorrFactorHisto_, *CorrFactorVsJetEtHisto_, *CorrFactorVsJetEtaHisto_;
  MonitorElement *ZSPCorrFactorHisto_, *ZSPCorrFactorVsJetEtHisto_, *ZSPCorrFactorVsJetEtaHisto_;
  MonitorElement *JPTCorrFactorHisto_, *JPTCorrFactorVsJetEtHisto_, *JPTCorrFactorVsJetEtaHisto_;
  TrackHistograms allPionHistograms_, inCaloInVertexPionHistograms_, inCaloOutVertexPionHistograms_, outCaloInVertexPionHistograms_;
  TrackHistograms allMuonHistograms_, inCaloInVertexMuonHistograms_, inCaloOutVertexMuonHistograms_, outCaloInVertexMuonHistograms_;
  TrackHistograms allElectronHistograms_, inCaloInVertexElectronHistograms_, inCaloOutVertexElectronHistograms_, outCaloInVertexElectronHistograms_;
  
  ///DQMStore. Used to write out to file
  DQMStore* dqm_;
};

inline void JPTJetAnalyzer::fillHistogram(MonitorElement* histogram, const double value)
{
  if (histogram) histogram->Fill(value);
}

inline void JPTJetAnalyzer::fillHistogram(MonitorElement* histogram, const double valueX, const double valueY)
{
  if (histogram) histogram->Fill(valueX,valueY);
}

#endif
