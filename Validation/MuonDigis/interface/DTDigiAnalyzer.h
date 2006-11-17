#ifndef SimMuon_DTDigiAnalyzer_h
#define SimMuon_DTDigiAnalyzer_h

/** \class DTDigiAnalyzer
 *  Analyse the the muon-drift-tubes digitizer. 
 *  
 *  $Date: 2006/08/07 15:47:36 $
 *  $Revision: 1.3 $
 *  \authors: R. Bellan
 */

#include<vector>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <DataFormats/MuonDetId/interface/DTWireId.h>
#include<TH1F.h>

class TH1F;
class TH2F;
class TFile;
class PSimHit;
//class DTMCStatistics;
//class DTMuonDigiStatistics;
//class DTHitsAnalysis;
class   hDigis;

namespace edm {
  class ParameterSet; class Event; class EventSetup;}

class DTDigiAnalyzer : public edm::EDAnalyzer{
  
 public:
  explicit DTDigiAnalyzer(const edm::ParameterSet& pset);
  virtual ~DTDigiAnalyzer();
  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  hDigis* WheelHistos(int wheel);
  
 private:
  typedef std::map<DTWireId, std::vector<const PSimHit*> > DTWireIdMap; 

  TH1F *DigiTimeBox;
  TH1F *DigiTimeBox_wheel2m;
  TH1F *DigiTimeBox_wheel1m;
  TH1F *DigiTimeBox_wheel0;
  TH1F *DigiTimeBox_wheel1p;
  TH1F *DigiTimeBox_wheel2p;
  TH1F *DigiEfficiency;
  TH1F *DigiEfficiencyMu;
  TH1F *DoubleDigi;
  TH2F *SimvsDigi;
  TH1F *Wire_DoubleDigi; 

  TH1F *MB1_sim_occup;
  TH1F *MB1_digi_occup;
  TH1F *MB2_sim_occup;
  TH1F *MB2_digi_occup;
  TH1F *MB3_sim_occup;
  TH1F *MB3_digi_occup;
  TH1F *MB4_sim_occup;
  TH1F *MB4_digi_occup;
 
  std::vector<TH1F*> DigiTimeBox_SL;
  TH1F *DigiHisto;
  /*
  int sl_num;
  char stringcham[20];
  //for ( slnum = 1; slnum < 61 ; ++slnum )
  //{     
    sl_num = 21;
     sprintf(stringcham, "DigiTimeBox_sl_%d", sl_num) ; 
     TH1F *stringcham;
  //}
  */

  TFile *file;
  std::string label;
  //  DTMCStatistics        *MCStatistics;
  // DTMuonDigiStatistics  *MuonDigiStatistics;
  // DTHitsAnalysis        *HitsAnalysis;
  
};

#endif    
