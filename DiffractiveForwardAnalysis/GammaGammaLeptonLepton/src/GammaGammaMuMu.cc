// -*- C++ -*-
//
// Package:    GammaGammaMuMu
// Class:      GammaGammaMuMu
// 
/**\class GammaGammaMuMu GammaGammaMuMu.cc GammaGammaLeptonLepton/GammaGammaMuMu/src/GammaGammaMuMu.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jonathan Hollar
//         Created:  Wed Sep 20 10:08:38 BST 2006
// $Id: GammaGammaMuMu.cc,v 1.65 2010/05/07 09:43:38 jjhollar Exp $
//
//


// system include files
#include "DataFormats/PatCandidates/interface/Muon.h" 
#include "DataFormats/PatCandidates/interface/Jet.h" 
#include "DataFormats/PatCandidates/interface/Electron.h" 
#include "DataFormats/PatCandidates/interface/Tau.h" 
#include "DataFormats/PatCandidates/interface/Photon.h" 
#include "DataFormats/PatCandidates/interface/MET.h" 
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h" 
#include "DataFormats/CastorReco/interface/CastorTower.h" 
#include "DataFormats/HcalRecHit/interface/CastorRecHit.h"
#include "DataFormats/Common/interface/Ref.h"   
#include "DataFormats/Common/interface/TriggerResults.h"   
#include "DataFormats/HLTReco/interface/TriggerEvent.h" 
 
#include "FWCore/ParameterSet/interface/ParameterSet.h"  
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterDescriptionNode.h"
#include "FWCore/Framework/interface/TriggerNames.h"  
   
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "DataFormats/RecoCandidate/interface/CaloRecHitCandidate.h"
#include "DataFormats/HcalRecHit/interface/ZDCRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"

#include "FWCore/Framework/interface/ESHandle.h" 
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h" 
#include "DataFormats/JetReco/interface/CaloJetCollection.h" 
#include "DataFormats/JetReco/interface/CaloJet.h"  
#include "DataFormats/EgammaCandidates/interface/Electron.h"  
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"   
#include "DataFormats/TrackReco/interface/Track.h"  
#include "DataFormats/TrackReco/interface/TrackFwd.h" 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"  
#include "DataFormats/MuonReco/interface/Muon.h"  
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"   
#include "DataFormats/METReco/interface/CaloMET.h"  
#include "DataFormats/METReco/interface/CaloMETFwd.h"   
#include "DataFormats/METReco/interface/CaloMETCollection.h"  
#include "DataFormats/EgammaCandidates/interface/Photon.h"  
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"  
#include "DataFormats/CaloTowers/interface/CaloTower.h"  
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"   
#include "DataFormats/Candidate/interface/Candidate.h"  
#include "DataFormats/Candidate/interface/CandidateFwd.h"   
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h" 
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/src/L1GctJetCounts.cc"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h" 
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h" 
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h" 
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapFwd.h" 
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h" 


#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h" // for PU

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaMuMu.h"

// user include files
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 
// Muons:
#include <DataFormats/TrackReco/interface/Track.h>
// Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

// Vertexing 
#include "DataFormats/VertexReco/interface/Vertex.h" 
#include "DataFormats/VertexReco/interface/VertexFwd.h" 
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" 
#include "TrackingTools/Records/interface/TransientTrackRecord.h" 
#include "TrackingTools/TransientTrack/interface/TransientTrack.h" 
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h" 
#include "RecoVertex/VertexPrimitives/interface/ConvertError.h" 
#include "SimTracker/Records/interface/TrackAssociatorRecord.h" 
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h" 

//#include "MuonAnalysis/TagAndProbe/interface/MuonPerformanceReadback.h" 
//#include "MuonAnalysis/TagAndProbe/interface/MuonPerformance.h"  


#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AcceptanceTableHelper.h"  


// C++
#include <memory>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GammaGammaMuMu::GammaGammaMuMu(const edm::ParameterSet& pset)
{
   //now do what ever initialization is needed
  recTrackLabel      = pset.getParameter<edm::InputTag>("RecoTrackLabel");
  recVertexLabel     = pset.getParameter<edm::InputTag>("RecoVertexLabel");
  theGLBMuonLabel    = pset.getParameter<edm::InputTag>("GlobalMuonCollectionLabel");
  thePixelGsfELabel  = pset.getParameter<edm::InputTag>("ElectronCollectionLabel");
  theJetLabel        = pset.getParameter<edm::InputTag>("JetCollectionLabel");
  theMetLabel        = pset.getParameter<edm::InputTag>("MetLabel");
  thePhotonLabel     = pset.getParameter<edm::InputTag>("PhotonCollectionLabel");
  theCaloTowLabel    = pset.getParameter<edm::InputTag>("CaloTowerLabel");
  recCastorTowerLabel = pset.getParameter<edm::InputTag>("CastorTowerLabel"); 
  recZDCRecHitsLabel = pset.getParameter<edm::InputTag>("ZDCRecHitsLabel");
  recCastorRecHitsLabel = pset.getParameter<edm::InputTag>("CastorRecHitsLabel");
  hltMenuLabel       = pset.getParameter<std::string>("HLTMenuLabel");

  mudptmax           = pset.getParameter<double>("DimuonMaxdpt");
  mudphimin          = pset.getParameter<double>("DimuonMindphi");
  drisocalo          = pset.getParameter<double>("CaloTowerdR");
  keepsamesign       = pset.getParameter<bool>("KeepSameSignDimuons");

  algonames          =  pset.getParameter< std::vector<std::string> >("AlgoNames"); 

  rootfilename       = pset.getUntrackedParameter<std::string>("outfilename","test.root");

  //  edm::FileInPath myDataFile("FastSimulation/ProtonTaggers/data/acceptance_420_220.root");  
  edm::FileInPath myDataFile("FastSimulation/ForwardDetectors/data/acceptance_420_220.root");
  std::string fullPath = myDataFile.fullPath();  
  std::cout << "Opening " << fullPath << std::endl;  
  TFile f(fullPath.c_str());  
  if (f.Get("description") != NULL)  
    std::cout << "Description found: " << f.Get("description")->GetTitle() << std::endl;  
    
  std::cout << "Reading acceptance tables " << std::endl;  
 
  helper420beam1.Init(f, "a420");  
  helper420beam2.Init(f, "a420_b2");  
  helper220beam1.Init(f, "a220");  
  helper220beam2.Init(f, "a220_b2");  
  helper420a220beam1.Init(f, "a420a220");  
  helper420a220beam2.Init(f, "a420a220_b2");  

  nEvt=0;
  MUONMAX=10;
  JETMAX=30;
  TRACKMAX=500;
  PHOTONMAX=50;
  GENPHOTONMAX=5;
  GENMUONMAX=10;

  thefile = new TFile(rootfilename.c_str(),"recreate");
  thefile->cd();
  thetree= new TTree("ntp1","ntp1");

  thetree->Branch("nJetCand",&nJetCand,"nJetCand/I");
  thetree->Branch("JetCand_px",JetCand_px,"JetCand_px[nJetCand]/D");
  thetree->Branch("JetCand_py",JetCand_py,"JetCand_py[nJetCand]/D");
  thetree->Branch("JetCand_pz",JetCand_pz,"JetCand_pz[nJetCand]/D");
  thetree->Branch("JetCand_e",JetCand_e,"JetCand_e[nJetCand]/D");
  thetree->Branch("JetCand_eta",JetCand_eta,"JetCand_eta[nJetCand]/D");
  thetree->Branch("JetCand_phi",JetCand_phi,"JetCand_phi[nJetCand]/D");
  thetree->Branch("HighestJet_e",&HighestJet_e,"HighestJet_e/D");
  thetree->Branch("HighestJet_eta",&HighestJet_eta,"HighestJet_eta/D"); 
  thetree->Branch("HighestJet_phi",&HighestJet_phi,"HighestJet_phi/D"); 
  thetree->Branch("SumJet_e",&SumJet_e,"SumJet_e/D");

  thetree->Branch("nMuonCand",&nMuonCand,"nMuonCand/I");
  thetree->Branch("MuonCand_px",MuonCand_px,"MuonCand_px[nMuonCand]/D");
  thetree->Branch("MuonCand_py",MuonCand_py,"MuonCand_py[nMuonCand]/D");
  thetree->Branch("MuonCand_pz",MuonCand_pz,"MuonCand_pz[nMuonCand]/D");
  thetree->Branch("MuonCand_p",MuonCand_p,"MuonCand_p[nMuonCand]/D");
  thetree->Branch("MuonCand_pt",MuonCand_pt,"MuonCand_pt[nMuonCand]/D");
  thetree->Branch("MuonCand_eta",MuonCand_eta,"MuonCand_eta[nMuonCand]/D");
  thetree->Branch("MuonCand_phi",MuonCand_phi,"MuonCand_phi[nMuonCand]/D");
  thetree->Branch("MuonCand_charge",MuonCand_charge,"MuonCand_charge[nMuonCand]/I");
  thetree->Branch("MuonCand_tmlsloosemuonid",MuonCand_tmlsloosemuonid,"MuonCand_tmlsloosemuonid[nMuonCand]/I");
  thetree->Branch("MuonCand_tmlsOptLowPtloosemuonid",MuonCand_tmlsOptLowPtloosemuonid,"MuonCand_tmlsOptLowPtloosemuonid[nMuonCand]/I");
  thetree->Branch("MuonCand_tm2dloosemuid",MuonCand_tm2dloosemuid,"MuonCand_tm2dloosemuid[nMuonCand]/I");
  thetree->Branch("MuonCand_tmlsAngloosemuonid",MuonCand_tmlsAngloosemuonid,"MuonCand_tmlsAngloosemuonid[nMuonCand]/I"); 
  thetree->Branch("MuonCand_tmlsAngtightmuonid",MuonCand_tmlsAngtightmuonid,"MuonCand_tmlsAngtightmuonid[nMuonCand]/I");
  thetree->Branch("MuonCand_tmosAngloosemuonid",MuonCand_tmosAngloosemuonid,"MuonCand_tmosAngloosemuonid[nMuonCand]/I");  
  thetree->Branch("MuonCand_tmosAngtightmuonid",MuonCand_tmosAngtightmuonid,"MuonCand_tmosAngtightmuonid[nMuonCand]/I");
  thetree->Branch("MuonCand_arbmuid",MuonCand_arbmuid,"MuonCand_arbmuid[nMuonCand]/I");
  thetree->Branch("MuonCand_isglobal",MuonCand_isglobal,"MuonCand_isglobal[nMuonCand]/I");
  thetree->Branch("MuonCand_istracker",MuonCand_istracker,"MuonCand_istracker[nMuonCand]/I"); 
  thetree->Branch("MuonCand_isstandalone",MuonCand_isstandalone,"MuonCand_isstandalone[nMuonCand]/I"); 
  thetree->Branch("MuonCand_hcalisor3",MuonCand_hcalisor3,"MuonCand_hcalisor3[nMuonCand]/D"); 
  thetree->Branch("MuonCand_ecalisor3",MuonCand_ecalisor3,"MuonCand_ecalisor3[nMuonCand]/D"); 
  thetree->Branch("MuonCand_hoisor3",MuonCand_hoisor3,"MuonCand_hoisor3[nMuonCand]/D"); 
  thetree->Branch("MuonCand_trkisor3",MuonCand_trkisor3,"MuonCand_trkisor3[nMuonCand]/D");  
  thetree->Branch("MuonCand_hcalisor5",MuonCand_hcalisor5,"MuonCand_hcalisor5[nMuonCand]/D");  
  thetree->Branch("MuonCand_ecalisor5",MuonCand_ecalisor5,"MuonCand_ecalisor5[nMuonCand]/D");  
  thetree->Branch("MuonCand_hoisor5",MuonCand_hoisor5,"MuonCand_hoisor5[nMuonCand]/D");
  thetree->Branch("MuonCand_trkisor5",MuonCand_trkisor5,"MuonCand_trkisor5[nMuonCand]/D"); 
  thetree->Branch("MuonCand_timein", MuonCand_timein, "MuonCand_timein[nMuonCand]/D"); 	 
  thetree->Branch("MuonCand_timeout", MuonCand_timeout, "MuonCand_timeout[nMuonCand]/D"); 
  thetree->Branch("MuonCand_timeinerr", MuonCand_timeinerr, "MuonCand_timeinerr[nMuonCand]/D");  
  thetree->Branch("MuonCand_timeouterr", MuonCand_timeouterr, "MuonCand_timeouterr[nMuonCand]/D"); 	 
  thetree->Branch("MuonCand_efficiency", MuonCand_efficiency, "MuonCand_efficiency[nMuonCand]/D");        
  thetree->Branch("MuonCand_validtrackhits", MuonCand_validtrackhits, "MuonCand_validtrackhits[nMuonCand]/I");        
  thetree->Branch("MuonCand_validhits", MuonCand_validhits, "MuonCand_validhits[nMuonCand]/I");        
  thetree->Branch("MuonCand_normchi2", MuonCand_normchi2, "MuonCand_normchi2[nMuonCand]/D");        
  thetree->Branch("MuonCand_normtrackchi2", MuonCand_normtrackchi2, "MuonCand_normtrackchi2[nMuonCand]/D");         


  thetree->Branch("nHLTMu3MuonCand",&nHLTMu3MuonCand,"nHLTMu3MuonCand/I"); 
  thetree->Branch("HLT_Mu3_MuonCand_pt",&HLT_Mu3_MuonCand_pt,"HLT_Mu3_MuonCand_pt[nHLTMu3MuonCand]/D"); 
  thetree->Branch("HLT_Mu3_MuonCand_eta",&HLT_Mu3_MuonCand_eta,"HLT_Mu3_MuonCand_eta[nHLTMu3MuonCand]/D"); 
  thetree->Branch("HLT_Mu3_MuonCand_phi",&HLT_Mu3_MuonCand_phi,"HLT_Mu3_MuonCand_phi[nHLTMu3MuonCand]/D"); 
  thetree->Branch("HLT_Mu3_MuonCand_charge",&HLT_Mu3_MuonCand_charge,"HLT_Mu3_MuonCand_charge[nHLTMu3MuonCand]/I");  

  thetree->Branch("nHLTDiMu0MuonCand",&nHLTDiMu0MuonCand,"nHLTDiMu0MuonCand/I");  
  thetree->Branch("HLT_DoubleMu0_MuonCand_pt",&HLT_DoubleMu0_MuonCand_pt,"HLT_DoubleMu0_MuonCand_pt[nHLTDiMu0MuonCand]/D");  
  thetree->Branch("HLT_DoubleMu0_MuonCand_eta",&HLT_DoubleMu0_MuonCand_eta,"HLT_DoubleMu0_MuonCand_eta[nHLTDiMu0MuonCand]/D");  
  thetree->Branch("HLT_DoubleMu0_MuonCand_phi",&HLT_DoubleMu0_MuonCand_phi,"HLT_DoubleMu0_MuonCand_phi[nHLTDiMu0MuonCand]/D");  
  thetree->Branch("HLT_DoubleMu0_MuonCand_charge",&HLT_DoubleMu0_MuonCand_charge,"HLT_DoubleMu0_MuonCand_charge[nHLTDiMu0MuonCand]/I");   

  thetree->Branch("nHLTDiMu3MuonCand",&nHLTDiMu3MuonCand,"nHLTDiMu3MuonCand/I");  
  thetree->Branch("HLT_DoubleMu3_MuonCand_pt",&HLT_DoubleMu3_MuonCand_pt,"HLT_DoubleMu3_MuonCand_pt[nHLTDiMu3MuonCand]/D");  
  thetree->Branch("HLT_DoubleMu3_MuonCand_eta",&HLT_DoubleMu3_MuonCand_eta,"HLT_DoubleMu3_MuonCand_eta[nHLTDiMu3MuonCand]/D");  
  thetree->Branch("HLT_DoubleMu3_MuonCand_phi",&HLT_DoubleMu3_MuonCand_phi,"HLT_DoubleMu3_MuonCand_phi[nHLTDiMu3MuonCand]/D");  
  thetree->Branch("HLT_DoubleMu3_MuonCand_charge",&HLT_DoubleMu3_MuonCand_charge,"HLT_DoubleMu3_MuonCand_charge[nHLTDiMu3MuonCand]/I");   

  thetree->Branch("nCaloCand",&nCaloCand,"nCaloCand/I");
  thetree->Branch("CaloTower_e",CaloTower_e,"CaloTower_e[nCaloCand]/D");
  thetree->Branch("CaloTower_et",CaloTower_et,"CaloTower_et[nCaloCand]/D");
  thetree->Branch("CaloTower_eta",CaloTower_eta,"CaloTower_eta[nCaloCand]/D"); 
  thetree->Branch("CaloTower_phi",CaloTower_phi,"CaloTower_phi[nCaloCand]/D"); 
  thetree->Branch("CaloTower_dr",CaloTower_dr,"CaloTower_dr[nCaloCand]/D");
  thetree->Branch("CaloTower_emE",CaloTower_emE,"CaloTower_emE[nCaloCand]/D");
  thetree->Branch("CaloTower_hadE",CaloTower_hadE,"CaloTower_hadE[nCaloCand]/D");
  thetree->Branch("CaloTower_outE",CaloTower_outE,"CaloTower_outE[nCaloCand]/D");
  thetree->Branch("CaloTower_ID",CaloTower_ID,"CaloTower_ID[nCaloCand]/I");
  thetree->Branch("CaloTower_x",CaloTower_x,"CaloTower_x[nCaloCand]/D");
  thetree->Branch("CaloTower_y",CaloTower_y,"CaloTower_y[nCaloCand]/D");
  thetree->Branch("CaloTower_z",CaloTower_z,"CaloTower_z[nCaloCand]/D");
  thetree->Branch("CaloTower_t",CaloTower_t,"CaloTower_t[nCaloCand]/D");
  thetree->Branch("CaloTower_badhcalcells",CaloTower_badhcalcells,"CaloTower_badhcalcells[nCaloCand]/I"); 
  thetree->Branch("CaloTower_problemhcalcells",CaloTower_problemhcalcells,"CaloTower_problemhcalcells[nCaloCand]/I"); 
  thetree->Branch("CaloTower_badecalcells",CaloTower_badecalcells,"CaloTower_badecalcells[nCaloCand]/I");  
  thetree->Branch("CaloTower_problemecalcells",CaloTower_problemecalcells,"CaloTower_problemecalcells[nCaloCand]/I");  
  thetree->Branch("HighestCaloTower_e",&HighestCaloTower_e,"HighestCaloTower_e/D");
  thetree->Branch("HighestCaloTower_eta",&HighestCaloTower_eta,"HighestCaloTower_eta/D");
  thetree->Branch("HighestCaloTower_phi",&HighestCaloTower_phi,"HighestCaloTower_phi/D"); 
  thetree->Branch("HighestCaloTower_dr",&HighestCaloTower_dr,"HighestCaloTower_dr/D");
  thetree->Branch("HighestEtCaloTower_et",&HighestEtCaloTower_et,"HighestEtCaloTower_et/D");
  thetree->Branch("HighestEtCaloTower_eta",&HighestEtCaloTower_eta,"HighestEtCaloTower_eta/D");
  thetree->Branch("HighestEtCaloTower_phi",&HighestEtCaloTower_phi,"HighestEtCaloTower_phi/D"); 
  thetree->Branch("HighestEtCaloTower_dr",&HighestEtCaloTower_dr,"HighestEtCaloTower_dr/D");
  thetree->Branch("SumCalo_e",&SumCalo_e,"SumCalo_e/D");

  thetree->Branch("nCastorTowerCand",&nCastorTowerCand,"nCastorTowerCand/I");  
  thetree->Branch("CastorTower_e",CastorTower_e,"CastorTower_e[nCastorTowerCand]/D");  
  thetree->Branch("CastorTower_eta",CastorTower_eta,"CastorTower_eta[nCastorTowerCand]/D");   
  thetree->Branch("CastorTower_phi",CastorTower_phi,"CastorTower_phi[nCastorTowerCand]/D");  
  thetree->Branch("CastorTower_emratio",CastorTower_emratio,"CastorTower_emratio[nCastorTowerCand]/D");  
  thetree->Branch("HighestCastorTowerFwd_e",&HighestCastorTowerFwd_e,"HighestCastorTowerFwd_e/D"); 
  thetree->Branch("HighestCastorTowerBwd_e",&HighestCastorTowerBwd_e,"HighestCastorTowerBwd_e/D"); 
  thetree->Branch("SumCastorFwd_e",&SumCastorFwd_e,"SumCastorFwd_e/D");
  thetree->Branch("SumCastorBwd_e",&SumCastorBwd_e,"SumCastorBwd_e/D"); 

  thetree->Branch("nZDChitCand", &nZDChitCand, "nZDChitCand/I");
  thetree->Branch("ZDChit_section", ZDChit_section, "ZDChit_section[nZDChitCand]/I");
  thetree->Branch("ZDChit_energy", ZDChit_energy, "ZDChit_energy[nZDChitCand]/D");
  thetree->Branch("ZDChit_time", ZDChit_time, "ZDChit_time[nZDChitCand]/D");
  thetree->Branch("ZDChit_side", ZDChit_side, "ZDChit_side[nZDChitCand]/I");
  thetree->Branch("ZDCsumEMplus", &ZDCsumEMplus, "ZDCsumEMplus/D");
  thetree->Branch("ZDCsumHADplus", &ZDCsumHADplus, "ZDCsumHADplus/D");
  thetree->Branch("ZDCsumEMminus", &ZDCsumEMminus, "ZDCsumEMminus/D");
  thetree->Branch("ZDCsumHADminus", &ZDCsumHADminus, "ZDCsumHADminus/D");
  thetree->Branch("CASTORsumRecHitsE", &CASTORsumRecHitsE, "CASTORsumRecHitsE/D");

  thetree->Branch("nExtraCaloTowersE1",&nExtraCaloTowersE1,"nExtraCaloTowersE1/I"); 
  thetree->Branch("nExtraCaloTowersE2",&nExtraCaloTowersE2,"nExtraCaloTowersE2/I"); 
  thetree->Branch("nExtraCaloTowersE3",&nExtraCaloTowersE3,"nExtraCaloTowersE3/I");  
  thetree->Branch("nExtraCaloTowersE4",&nExtraCaloTowersE4,"nExtraCaloTowersE4/I");  
  thetree->Branch("nExtraCaloTowersE5",&nExtraCaloTowersE5,"nExtraCaloTowersE5/I");  
  thetree->Branch("nExtraCaloTowersE6",&nExtraCaloTowersE6,"nExtraCaloTowersE6/I");   
  thetree->Branch("nExtraCaloTowersE7",&nExtraCaloTowersE7,"nExtraCaloTowersE7/I");   
  thetree->Branch("nExtraCaloTowersE8",&nExtraCaloTowersE8,"nExtraCaloTowersE8/I");   
  thetree->Branch("nExtraCaloTowersE9",&nExtraCaloTowersE9,"nExtraCaloTowersE9/I");   

  thetree->Branch("nExtraCaloTowersE0hf", &nExtraCaloTowersE0hf, "nExtraCaloTowersE0hf/I"); 	 
  thetree->Branch("nExtraCaloTowersE1hf", &nExtraCaloTowersE1hf, "nExtraCaloTowersE1hf/I"); 	 
  thetree->Branch("nExtraCaloTowersE2hf", &nExtraCaloTowersE2hf, "nExtraCaloTowersE2hf/I"); 	 
  thetree->Branch("nExtraCaloTowersE3hf", &nExtraCaloTowersE3hf, "nExtraCaloTowersE3hf/I");  
  thetree->Branch("nExtraCaloTowersE4hf", &nExtraCaloTowersE4hf, "nExtraCaloTowersE4hf/I");  
  thetree->Branch("nExtraCaloTowersE5hf", &nExtraCaloTowersE5hf, "nExtraCaloTowersE5hf/I");  

  thetree->Branch("nExtraCaloTowersE1he", &nExtraCaloTowersE1he, "nExtraCaloTowersE1he/I"); 	 
  thetree->Branch("nExtraCaloTowersE2he", &nExtraCaloTowersE2he, "nExtraCaloTowersE2he/I"); 	 
  thetree->Branch("nExtraCaloTowersE3he", &nExtraCaloTowersE3he, "nExtraCaloTowersE3he/I"); 	 
  thetree->Branch("nExtraCaloTowersE4he", &nExtraCaloTowersE4he, "nExtraCaloTowersE4he/I");   
  thetree->Branch("nExtraCaloTowersE5he", &nExtraCaloTowersE5he, "nExtraCaloTowersE5he/I");   

  thetree->Branch("nExtraCaloTowersE2hb", &nExtraCaloTowersE2hb, "nExtraCaloTowersE2hb/I"); 	 
  thetree->Branch("nExtraCaloTowersE3hb", &nExtraCaloTowersE3hb, "nExtraCaloTowersE3hb/I"); 	 
  thetree->Branch("nExtraCaloTowersE4hb", &nExtraCaloTowersE4hb, "nExtraCaloTowersE4hb/I");
  thetree->Branch("nExtraCaloTowersE5hb", &nExtraCaloTowersE5hb, "nExtraCaloTowersE5hb/I");   

  thetree->Branch("nExtraCaloTowersEt0pt1",&nExtraCaloTowersEt0pt1,"nExtraCaloTowersEt0pt1/I"); 
  thetree->Branch("nExtraCaloTowersEt0pt2",&nExtraCaloTowersEt0pt2,"nExtraCaloTowersEt0pt2/I"); 
  thetree->Branch("nExtraCaloTowersEt0pt5",&nExtraCaloTowersEt0pt5,"nExtraCaloTowersEt0pt5/I"); 
  thetree->Branch("nExtraCaloTowersEt1",&nExtraCaloTowersEt1,"nExtraCaloTowersEt1/I"); 
  thetree->Branch("nExtraCaloTowersEt2",&nExtraCaloTowersEt2,"nExtraCaloTowersEt2/I"); 
  thetree->Branch("nExtraCaloTowersEt3",&nExtraCaloTowersEt3,"nExtraCaloTowersEt3/I");  
  thetree->Branch("nExtraCaloTowersEt4",&nExtraCaloTowersEt4,"nExtraCaloTowersEt4/I");  

  thetree->Branch("nTrackCand",&nTrackCand,"nTrackCand/I");
  thetree->Branch("nQualityTrackCand",&nQualityTrackCand,"nQualityTrackCand/I"); 
  thetree->Branch("TrackCand_px",TrackCand_px,"TrackCand_px[nTrackCand]/D");
  thetree->Branch("TrackCand_py",TrackCand_py,"TrackCand_py[nTrackCand]/D");
  thetree->Branch("TrackCand_pz",TrackCand_pz,"TrackCand_pz[nTrackCand]/D");
  thetree->Branch("TrackCand_p",TrackCand_p,"TrackCand_p[nTrackCand]/D");
  thetree->Branch("TrackCand_pt",TrackCand_pt,"TrackCand_pt[nTrackCand]/D");
  thetree->Branch("TrackCand_eta",TrackCand_eta,"TrackCand_eta[nTrackCand]/D");
  thetree->Branch("TrackCand_phi",TrackCand_phi,"TrackCand_phi[nTrackCand]/D");
  thetree->Branch("TrackCand_vtxdxyz",TrackCand_vtxdxyz,"TrackCand_vtxdxyz[nTrackCand]/D");
  thetree->Branch("TrackCand_purity",TrackCand_purity,"TrackCand_purity[nTrackCand]/D");
  thetree->Branch("TrackCand_nhits",TrackCand_nhits,"TrackCand_nhits[nTrackCand]/I");

  thetree->Branch("TrackCand_vtxZ",TrackCand_vtxZ,"TrackCand_vtxZ[nTrackCand]/D");
  thetree->Branch("TrackCand_vtxT",TrackCand_vtxT,"TrackCand_vtxT[nTrackCand]/D");
  thetree->Branch("ClosestExtraTrack_vtxdxyz",&ClosestExtraTrack_vtxdxyz,"ClosestExtraTrack_vtxdxyz/D");
  
  thetree->Branch("nPFPhotonCand",&nPFPhotonCand,"nPFPhotonCand/I");
  thetree->Branch("PFPhotonCand_pt",PFPhotonCand_pt,"PFPhotonCand_pt[nPFPhotonCand]/D");
  thetree->Branch("PFPhotonCand_eta",PFPhotonCand_eta,"PFPhotonCand_eta[nPFPhotonCand]/D");
  thetree->Branch("PFPhotonCand_phi",PFPhotonCand_phi,"PFPhotonCand_phi[nPFPhotonCand]/D");
  thetree->Branch("PFPhotonCand_drtrue",PFPhotonCand_drtrue,"PFPhotonCand_drtrue[nPFPhotonCand]/D"); 

  thetree->Branch("MuMu_mass",&MuMu_mass,"MuMu_mass/D");
  thetree->Branch("MuMu_dphi",&MuMu_dphi,"MuMu_dphi/D");
  thetree->Branch("MuMu_dpt",&MuMu_dpt,"MuMu_dpt/D");
  thetree->Branch("MuMu_vtxx",&MuMu_vtxx,"MuMu_vtxx/D");
  thetree->Branch("MuMu_vtxy",&MuMu_vtxy,"MuMu_vtxy/D"); 
  thetree->Branch("MuMu_vtxz",&MuMu_vtxz,"MuMu_vtxz/D");
  thetree->Branch("MuMu_vtxT",&MuMu_vtxT,"MuMu_vtxT/D"); 
  thetree->Branch("MuMu_vtxchi2dof",&MuMu_vtxchi2dof,"MuMu_vtxchi2dof/D");
  thetree->Branch("MuMu_vtxisvalid",&MuMu_vtxisvalid,"MuMu_vtxisvalid/I");
  thetree->Branch("MuMu_extratracks1mm",&MuMu_extratracks1mm,"MuMu_extratracks1mm/I");
  thetree->Branch("MuMu_extratracks3mm",&MuMu_extratracks3mm,"MuMu_extratracks3mm/I");
  thetree->Branch("MuMu_extratracks5mm",&MuMu_extratracks5mm,"MuMu_extratracks5mm/I");
  thetree->Branch("MuMu_extratracks1cm",&MuMu_extratracks1cm,"MuMu_extratracks1cm/I");
  thetree->Branch("MuMu_extratracks3cm",&MuMu_extratracks3cm,"MuMu_extratracks3cm/I");
  thetree->Branch("MuMu_extratracks5cm",&MuMu_extratracks5cm,"MuMu_extratracks5cm/I"); 
  thetree->Branch("MuMu_extratracks10cm",&MuMu_extratracks10cm,"MuMu_extratracks10cm/I"); 
  thetree->Branch("MuMuGamma_mass",MuMuGamma_mass,"MuMuGamma_mass[nPFPhotonCand]/D"); 

  thetree->Branch("HitInZDC",&HitInZDC,"HitInZDC/I");
  thetree->Branch("HitInCastor",&HitInCastor,"HitInCastor/I");

  thetree->Branch("nGenPhotCand",&nGenPhotCand,"nGenPhotCand/I"); 
  thetree->Branch("GenPhotCand_pt",GenPhotCand_pt,"GenPhotCand_pt[nGenPhotCand]/D"); 
  thetree->Branch("GenPhotCand_eta",GenPhotCand_eta,"GenPhotCand_eta[nGenPhotCand]/D");  
  thetree->Branch("GenPhotCand_phi",GenPhotCand_phi,"GenPhotCand_phi[nGenPhotCand]/D");  

  thetree->Branch("nGenMuonCand",&nGenMuonCand,"nGenMuonCand/I");  
  thetree->Branch("GenMuonCand_px",GenMuonCand_px,"GenMuonCand_px[nGenMuonCand]/D");  
  thetree->Branch("GenMuonCand_py",GenMuonCand_py,"GenMuonCand_py[nGenMuonCand]/D");   
  thetree->Branch("GenMuonCand_pz",GenMuonCand_pz,"GenMuonCand_pz[nGenMuonCand]/D");   

  thetree->Branch("GenMuMu_eta",&GenMuMu_eta,"GenMuMu_eta/D");    
  thetree->Branch("GenMuMu_pt",&GenMuMu_pt,"GenMuMu_pt/D");     
  
  thetree->Branch("Etmiss",&Etmiss,"Etmiss/D");

  thetree->Branch("HLT_DoubleMu3",&HLT_DoubleMu3,"HLT_DoubleMu3/I");
  thetree->Branch("HLT_DoubleMu0",&HLT_DoubleMu0,"HLT_DoubleMu0/I");
  thetree->Branch("HLT_Mu3",&HLT_Mu3,"HLT_Mu3/I");
  thetree->Branch("HLT_L2Mu0", &HLT_L2Mu0, "HLT_L2Mu0/I"); 
  thetree->Branch("HLT_L1DoubleMuOpen", &HLT_L1DoubleMuOpen, "HLT_L1DoubleMuOpen/I"); 

  thetree->Branch("Run",&Run,"Run/I");
  thetree->Branch("LumiSection",&LumiSection,"LumiSection/I");
  thetree->Branch("BX",&BX,"BX/I");
  thetree->Branch("EventNum",&EventNum,"EventNum/I");
  thetree->Branch("L1TechnicalTriggers",L1TechnicalTriggers,"L1TechnicalTriggers[128]/I"); 

  thetree->Branch("nPrimVertexCand",&nPrimVertexCand,"nPrimVertexCand/I");
  thetree->Branch("PrimVertexCand_x",&PrimVertexCand_x,"PrimVertexCand_x[nPrimVertexCand]/D");
  thetree->Branch("PrimVertexCand_y",&PrimVertexCand_y,"PrimVertexCand_y[nPrimVertexCand]/D");
  thetree->Branch("PrimVertexCand_z",&PrimVertexCand_z,"PrimVertexCand_z[nPrimVertexCand]/D");
  thetree->Branch("PrimVertexCand_tracks",&PrimVertexCand_tracks,"PrimVertexCand_tracks[nPrimVertexCand]/I");
  thetree->Branch("PrimVertexCand_chi2",&PrimVertexCand_chi2,"PrimVertexCand_chi2[nPrimVertexCand]/D");
  thetree->Branch("PrimVertexCand_ndof",&PrimVertexCand_ndof,"PrimVertexCand_ndof[nPrimVertexCand]/D");

  thetree->Branch("LowPt_pt",LowPt_pt,"LowPt_pt[nMuonCand]/D");
  thetree->Branch("LowPt_eta",LowPt_eta,"LowPt_eta[nMuonCand]/D");

  thetree->Branch("nPU",&nPU,"nPU/I");

//  thetree->Branch("evweight",&evweight,"evweight/D");
}


GammaGammaMuMu::~GammaGammaMuMu()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GammaGammaMuMu::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
  nPrimVertexCand=0;
  nMuonCand=0;
  nHLTMu3MuonCand=0;
  nHLTDiMu3MuonCand=0;  
  nHLTDiMu0MuonCand=0;
  nJetCand=0;
  nCaloCand=0;
  nTrackCand=0;
  nQualityTrackCand=0;
  nCastorTowerCand=0;
  nZDChitCand=0;
  ZDCsumHADminus=0;
  ZDCsumEMminus=0;
  ZDCsumHADplus=0;
  ZDCsumEMplus=0;
  nExtraCaloTowersE1=0;
  nExtraCaloTowersE2=0;
  nExtraCaloTowersE3=0; 
  nExtraCaloTowersE4=0; 
  nExtraCaloTowersE5=0;
  nExtraCaloTowersE6=0; 
  nExtraCaloTowersE7=0;  
  nExtraCaloTowersE8=0;  
  nExtraCaloTowersE9=0; 
  nExtraCaloTowersEt0pt1=0; 
  nExtraCaloTowersEt0pt2=0;
  nExtraCaloTowersEt0pt5=0;  
  nExtraCaloTowersEt1=0; 
  nExtraCaloTowersEt2=0;  
  nExtraCaloTowersEt3=0;  
  nExtraCaloTowersEt4=0;   
  nExtraCaloTowersE0hf=0; 	 
  nExtraCaloTowersE1hf=0; 	 
  nExtraCaloTowersE2hf=0; 	 
  nExtraCaloTowersE3hf=0; 
  nExtraCaloTowersE4hf=0; 
  nExtraCaloTowersE5hf=0;   
 
  nExtraCaloTowersE1he=0; 	 
  nExtraCaloTowersE2he=0; 	 
  nExtraCaloTowersE3he=0; 	 
  nExtraCaloTowersE4he=0;  
  nExtraCaloTowersE5he=0;   

  nExtraCaloTowersE2hb=0; 	 
  nExtraCaloTowersE3hb=0; 	 
  nExtraCaloTowersE4hb=0;
  nExtraCaloTowersE5hb=0;   

  HitInZDC=0;
  HitInCastor=0;
  nGenPhotCand=0;
  nGenMuonCand=0;

  nPFPhotonCand=0;

  MuMu_mass = -1;
  MuMu_dphi = -1;
  MuMu_dpt = -1;
  MuMu_extratracks1mm = 0; 
  MuMu_extratracks3mm = 0; 
  MuMu_extratracks5mm = 0; 
  MuMu_extratracks1cm = 0; 
  MuMu_extratracks3cm = 0;
  MuMu_extratracks5cm = 0;
  MuMu_extratracks10cm = 0;
  ClosestExtraTrack_vtxdxyz = 999.;

  bool passed = true;
  int LS = 0;
  int LSopt = 0;

  nEvt++;
  using reco::TrackCollection;

  // Run and BX information
  BX = event.bunchCrossing();
  Run = event.id().run();
  LumiSection = event.luminosityBlock();
  EventNum = event.id().event();

  // L1 technical triggers 
  edm::Handle<L1GlobalTriggerReadoutRecord> L1GTRR; 
  edm::Handle<L1GlobalTriggerObjectMapRecord> L1GTOMRec; 
  event.getByLabel(InputTag("gtDigis::RECO"), L1GTRR); 
  event.getByLabel(InputTag("hltL1GtObjectMap::HLT"), L1GTOMRec); 
  if (L1GTRR.isValid()) { 
    DecisionWord gtDecisionWord = L1GTRR->decisionWord(); 
    const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = L1GTRR->technicalTriggerWord(); 
    const unsigned int numberTechnicalTriggerBits(technicalTriggerWordBeforeMask.size()); 
    for (unsigned int iBit = 0; iBit < numberTechnicalTriggerBits; ++iBit) { 
      int techTrigger = (int) technicalTriggerWordBeforeMask.at(iBit); 
      L1TechnicalTriggers[iBit] = techTrigger; 
    } 
  } 

  // Get the trigger information from the event
  edm::Handle<edm::TriggerResults> hltResults ; 
  event.getByLabel(InputTag("TriggerResults","",hltMenuLabel),hltResults) ; 
  trigNames.init(*hltResults) ;
  for (unsigned int i=0; i<trigNames.size(); i++)  
    { 
      if ( trigNames.triggerNames().at(i) == "HLT_Mu3" )       
        {  
          if ( hltResults->accept(i) )  
            {HLT_Mu3 = 1;}
	  else
	    HLT_Mu3 = 0;
        }  
      if ( trigNames.triggerNames().at(i) == "HLT_DoubleMu0" ) 
        {   
          if ( hltResults->accept(i) )  
	    {HLT_DoubleMu0 = 1;}
	  else
	    HLT_DoubleMu0 = 0;
        }  
      if ( trigNames.triggerNames().at(i) == "HLT_DoubleMu3" ) 
        {   
          if ( hltResults->accept(i) )  
	    {HLT_DoubleMu3 = 1;}
	  else
	    HLT_DoubleMu3 = 0;
        }  
      if ( trigNames.triggerNames().at(i) == "HLT_L1DoubleMuOpen" )
        {    
          if ( hltResults->accept(i) )   
            {HLT_L1DoubleMuOpen = 1;}
          else 
            HLT_L1DoubleMuOpen = 0; 
        }   
      if ( trigNames.triggerNames().at(i) == "HLT_L2Mu0" ) 
        {     
          if ( hltResults->accept(i) )    
            {HLT_L2Mu0 = 1;} 
          else  
            HLT_L2Mu0 = 0;  
        }    
    }

  Handle<TriggerEvent> hltObjects;
  event.getByLabel(InputTag("hltTriggerSummaryAOD","",hltMenuLabel),hltObjects);
  if (hltObjects.isValid()) 
    {
      size_type muindex = hltObjects->filterIndex(InputTag("hltSingleMu3L3Filtered3::"+hltMenuLabel)); 
      size_type dimu0index = hltObjects->filterIndex(InputTag("hltDiMuonL3PreFiltered0::"+hltMenuLabel)); 
      size_type dimuindex = hltObjects->filterIndex(InputTag("hltDiMuonL3PreFiltered::"+hltMenuLabel)); 

      if( dimu0index < hltObjects->sizeFilters() )
	{
	  const trigger::Keys& DIMU0KEYS(hltObjects->filterKeys(dimu0index)); 
	  const size_type nK(DIMU0KEYS.size()); 
	  const TriggerObjectCollection& TOC(hltObjects->getObjects());

	  for(int ipart = 0; ipart != nK; ++ipart)   
	    {
	      const TriggerObject& TO = TOC[DIMU0KEYS[ipart]];  
	      
	      if(fabs(TO.id()) == 13)
		{
		  HLT_DoubleMu0_MuonCand_pt[nHLTDiMu0MuonCand] = TO.pt();
		  HLT_DoubleMu0_MuonCand_eta[nHLTDiMu0MuonCand] = TO.eta(); 
		  HLT_DoubleMu0_MuonCand_phi[nHLTDiMu0MuonCand] = TO.phi(); 
		  if(TO.id() > 0)	  
		    HLT_DoubleMu0_MuonCand_charge[nHLTDiMu0MuonCand] = 1; 
		  else
		    HLT_DoubleMu0_MuonCand_charge[nHLTDiMu0MuonCand] = -1;  
		  
		  nHLTDiMu0MuonCand++;
		}
	    }
	}
      if( dimuindex < hltObjects->sizeFilters() )
	{
	  const trigger::Keys& DIMUKEYS(hltObjects->filterKeys(dimuindex)); 
	  const size_type nK(DIMUKEYS.size()); 
	  const TriggerObjectCollection& TOC(hltObjects->getObjects());

	  for(int ipart = 0; ipart != nK; ++ipart)   
	    {
	      const TriggerObject& TO = TOC[DIMUKEYS[ipart]];  
	      
	      if(fabs(TO.id()) == 13)
		{
		  HLT_DoubleMu3_MuonCand_pt[nHLTDiMu3MuonCand] = TO.pt();
		  HLT_DoubleMu3_MuonCand_eta[nHLTDiMu3MuonCand] = TO.eta(); 
		  HLT_DoubleMu3_MuonCand_phi[nHLTDiMu3MuonCand] = TO.phi(); 
		  if(TO.id() > 0)	  
		    HLT_DoubleMu3_MuonCand_charge[nHLTDiMu3MuonCand] = 1; 
		  else
		    HLT_DoubleMu3_MuonCand_charge[nHLTDiMu3MuonCand] = -1;  
		  
		  nHLTDiMu3MuonCand++;
		}
	    }
	}
      if( muindex < hltObjects->sizeFilters() ) 
	{
          const trigger::Keys& MUKEYS(hltObjects->filterKeys(muindex));  
          const size_type nK(MUKEYS.size());  
          const TriggerObjectCollection& TOC(hltObjects->getObjects()); 
	  
          for(int ipart = 0; ipart != nK; ++ipart)    
            { 
              const TriggerObject& TO = TOC[MUKEYS[ipart]];   
               
              if(fabs(TO.id()) == 13) 
                { 
                  HLT_Mu3_MuonCand_pt[nHLTMu3MuonCand] = TO.pt(); 
                  HLT_Mu3_MuonCand_eta[nHLTMu3MuonCand] = TO.eta();  
                  HLT_Mu3_MuonCand_phi[nHLTMu3MuonCand] = TO.phi();  
                  if(TO.id() > 0)          
                    HLT_Mu3_MuonCand_charge[nHLTMu3MuonCand] = 1;  
                  else 
                    HLT_Mu3_MuonCand_charge[nHLTMu3MuonCand] = -1;   
                   
                  nHLTMu3MuonCand++; 
                } 
            } 
	  
	}
    }


  // Get the #PU information
  nPU=0;
  edm::Handle<CrossingFrame<edm::HepMCProduct> > crossingFrameHepMCH;
  event.getByLabel("mix","source",crossingFrameHepMCH);
  
  if((crossingFrameHepMCH.isValid())) 
    {  
      //    unsigned int nbunches = 0;
      //    for(int ibunch = crossingFrameHepMCH->getBunchRange().first; 
      //            ibunch <= crossingFrameHepMCH->getBunchRange().second; ++ibunch,++nbunches){
      //      int nrPileUpB = crossingFrameHepMCH->getNrPileups(ibunch);
      //      std::cout << "  --> Number of added pile-up's for bunch " << ibunch   << ": " << nrPileUpB << std::endl;
      //    }
      int nrPileUp = crossingFrameHepMCH->getNrPileups(0);
      nPU = nrPileUp;
    }

  // Get the muon collection from the event
  // PAT
  edm::Handle<edm::View<pat::Muon> > muons; 
  event.getByLabel(theGLBMuonLabel,muons); 
  edm::View<pat::Muon>::const_iterator muon;

  // AOD
  /*
    Handle<reco::MuonCollection> muons;
    event.getByLabel(theGLBMuonLabel, muons);
    reco::MuonCollection::const_iterator muon;
  */

//cout << "================" << endl;
  if(muons->size() == 2)
    {
      for (muon = muons->begin(); muon != muons->end() && nMuonCand<MUONMAX; ++muon)
	{
	  MuonCand_p[nMuonCand]=muon->p();
	  MuonCand_px[nMuonCand]=muon->px();
	  MuonCand_py[nMuonCand]=muon->py();
	  MuonCand_pz[nMuonCand]=muon->pz();
	  MuonCand_pt[nMuonCand]=muon->pt();
	  MuonCand_eta[nMuonCand]=muon->eta();
	  MuonCand_phi[nMuonCand]=muon->phi();
	  MuonCand_charge[nMuonCand]=muon->charge();

	  // Muon ID - 31X compatible
	  MuonCand_tmlsloosemuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMLastStationLoose);
          MuonCand_tmlsOptLowPtloosemuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMLastStationOptimizedLowPtLoose);
	  MuonCand_tm2dloosemuid[nMuonCand]=muon::isGoodMuon(*muon, muon::TM2DCompatibilityLoose);
	  MuonCand_tmlsAngloosemuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMLastStationAngLoose); 
	  MuonCand_tmlsAngtightmuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMLastStationAngTight);  
	  MuonCand_tmosAngloosemuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMOneStationAngLoose);  
	  MuonCand_tmosAngtightmuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMOneStationAngTight);  
	  MuonCand_arbmuid[nMuonCand]=muon::isGoodMuon(*muon, muon::AllArbitrated);
	  MuonCand_isglobal[nMuonCand]=muon->isGlobalMuon();
          MuonCand_istracker[nMuonCand]=muon->isTrackerMuon(); 
          MuonCand_isstandalone[nMuonCand]=muon->isStandAloneMuon();

	  // Muon ID - 22X compatible
// 	  MuonCand_tmlsloosemuonid[nMuonCand]=muon->isGood(reco::Muon::TMLastStationLoose);
//           MuonCand_tmlsOptLowPtloosemuonid[nMuonCand]=muon->isGood(reco::Muon::TMLastStationOptimizedLowPtLoose);
// 	  //muonid::isGoodMuon(*muon, *muon,muonid::TMLastStationLoose);
// 	  MuonCand_tm2dloosemuid[nMuonCand]=muon->isGood(reco::Muon::TM2DCompatibilityLoose);
// 	  //Muon(*muon,muonid::TM2DCompatibilityLoose);
// 	  MuonCand_arbmuid[nMuonCand]=muon->isGood(reco::Muon::AllArbitrated);
// 	  MuonCand_isglobal[nMuonCand]=muon->isGlobalMuon();
//           MuonCand_istracker[nMuonCand]=muon->isTrackerMuon(); 
//           MuonCand_isstandalone[nMuonCand]=muon->isStandAloneMuon();
 
	  //if(!(muon->isGood(reco::Muon::TMLastStationLoose)) && (muon->isGood(reco::Muon::TMLastStationOptimizedLowPtLoose))){
	  if(!(muon::isGoodMuon(*muon, muon::TMLastStationLoose)) && (muon::isGoodMuon(*muon, muon::TMLastStationOptimizedLowPtLoose))){	  
	    LowPt_eta[nMuonCand]=muon->eta();
	    LowPt_pt[nMuonCand]=muon->pt();
	  }else{  LowPt_eta[nMuonCand]=10.;
	    LowPt_pt[nMuonCand]=-1.;}
	  
	  //	if(muon->isGood(reco::Muon::TMLastStationLoose)) LS++;
	  //	if(muon->isGood(reco::Muon::TMLastStationOptimizedLowPtLoose)) LSopt++;

	  if(muon::isGoodMuon(*muon, muon::TMLastStationLoose)) LS++;
	  if(muon::isGoodMuon(*muon, muon::TMLastStationOptimizedLowPtLoose)) LSopt++;

	  // Isolation 
	  MuonCand_hcalisor3[nMuonCand]=muon->isolationR03().hadEt;
	  MuonCand_ecalisor3[nMuonCand]=muon->isolationR03().emEt;  
          MuonCand_hoisor3[nMuonCand]=muon->isolationR03().hoEt;
	  MuonCand_trkisor3[nMuonCand]=muon->isolationR03().nTracks;  

	  MuonCand_hcalisor5[nMuonCand]=muon->isolationR05().hadEt;  
	  MuonCand_ecalisor5[nMuonCand]=muon->isolationR05().emEt;   
          MuonCand_hoisor5[nMuonCand]=muon->isolationR05().hoEt;
	  MuonCand_trkisor5[nMuonCand]=muon->isolationR05().nTracks;  

          MuonCand_timeout[nMuonCand]=muon->time().timeAtIpInOut; 	 
          MuonCand_timein[nMuonCand]=muon->time().timeAtIpOutIn; 	 
          MuonCand_timeouterr[nMuonCand]=muon->time().timeAtIpInOutErr; 	 
          MuonCand_timeinerr[nMuonCand]=muon->time().timeAtIpOutInErr; 	 
 
	  if(muon->isTrackerMuon() || muon->isGlobalMuon()) 
	    {
	      MuonCandTrack_p[nMuonCand] = muon->innerTrack()->p();
	      MuonCand_validtrackhits[nMuonCand]=muon->innerTrack()->numberOfValidHits();  
	      MuonCand_validhits[nMuonCand]=muon->numberOfValidHits(); 
              MuonCand_normtrackchi2[nMuonCand]=muon->innerTrack()->normalizedChi2();  

	      if(muon->isGlobalMuon())
		MuonCand_normchi2[nMuonCand]=muon->normChi2(); 
	    }

	  std::string algoname;
	  double totalmuoneff = 1.0;

	  // October exercise - testing efficiencies!
	  /*
          for(unsigned int i = 0; i < algonames.size(); ++i)  
            { 
              algoname = algonames[i]; 

	      const MuonPerformance &muonefficiency = effreader->getPerformanceRecord(algoname, iSetup);   
	      double muoneff = effreader->getEff(MuonCand_pt[nMuonCand], MuonCand_eta[nMuonCand], MuonCand_phi[nMuonCand], MuonCand_charge[nMuonCand], muonefficiency); 
	      effreader->getEffError(MuonCand_pt[nMuonCand], MuonCand_eta[nMuonCand], MuonCand_phi[nMuonCand], MuonCand_charge[nMuonCand], muonefficiency);  
	      
	      if(muoneff > -1)
		totalmuoneff *= muoneff;
	    }
	  */
	  MuonCand_efficiency[nMuonCand] = totalmuoneff;

/*	
          if(muon->isStandAloneMuon() )
		{cout << "--> debug : stdalone  µ" << endl;
		 cout << "            n valid hit = " << muon->standAloneMuon()->numberOfValidHits() << endl;
		 cout << "            chi² = " << muon->standAloneMuon()->normalizedChi2() << endl;
		}
          if(muon->isGlobalMuon() )
                {cout << "--> debug : global   µ" << endl;
                 cout << "            n valid hit = " << muon->combinedMuon()->numberOfValidHits() << endl;
                 cout << "            chi² = " << muon->combinedMuon()->normalizedChi2() << endl;
		}
          if(muon->isTrackerMuon() )
                {cout << "--> debug : tracker    µ" << endl;
                 cout << "            n valid hit = " << muon->track()->numberOfValidHits() << endl;
                 cout << "            chi² = " << muon->track()->normalizedChi2() << endl;
		}
*/
	  nMuonCand++;
	}  

      // Calculate invariant mass, delta-phi and delta-pT
      if((MuonCand_charge[0]*MuonCand_charge[1]<0) || (keepsamesign == true))
	{
	  double mass = pow(MuonCand_p[0]+MuonCand_p[1],2);
	  mass-=pow(MuonCand_px[0]+MuonCand_px[1],2);
	  mass-=pow(MuonCand_py[0]+MuonCand_py[1],2);
	  mass-=pow(MuonCand_pz[0]+MuonCand_pz[1],2);
	  MuMu_mass = sqrt(mass);

          MuMu_dpt = fabs(MuonCand_pt[0]-MuonCand_pt[1]);

	  double dphi = fabs(MuonCand_phi[0]-MuonCand_phi[1]);
	  if(dphi < 3.14159)
	    MuMu_dphi = dphi;
	  else
	    MuMu_dphi = (2.0*3.14159)-dphi;
	}
    }
  // Get the Jet collection from the event
  // PAT
  edm::Handle<edm::View<pat::Jet> > jets; 
  event.getByLabel(theJetLabel,jets); 
  edm::View<pat::Jet>::const_iterator jet;

  // AOD
  /*
  edm::Handle<reco::CaloJetCollection> pJets;
  event.getByLabel(theJetLabel,pJets);
  const reco::CaloJetCollection* jets = pJets.product();
  reco::CaloJetCollection::const_iterator jet;
  */

  // Get the MET collection from the event
  // PAT
  /*
  edm::Handle<edm::View<pat::MET> > mets; 
  event.getByLabel(theMetLabel,mets); 
  edm::View<pat::MET>::const_iterator met;
  */

  // AOD
  edm::Handle<reco::CaloMETCollection> pMET;
  event.getByLabel(theMetLabel,pMET);
  const reco::CaloMETCollection* mets = pMET.product();
  reco::CaloMETCollection::const_iterator met;

  // Get the CaloTower collection from the event
  edm::Handle<CaloTowerCollection> caloTowers; 
  event.getByLabel(theCaloTowLabel,caloTowers); 
  const CaloTowerCollection* towers = caloTowers.product(); 
  CaloTowerCollection::const_iterator calo; 

  // Get the track collection from the event
  edm::Handle<reco::TrackCollection> recoTracks;
  event.getByLabel(recTrackLabel, recoTracks);
  const TrackCollection* tracks = recoTracks.product();
  TrackCollection::const_iterator track;

  // Get the vertex collection from the event
  edm::Handle<reco::VertexCollection> recoVertexs;
  event.getByLabel(recVertexLabel, recoVertexs);
  const VertexCollection* vertexs = recoVertexs.product();
  VertexCollection::const_iterator vertex_i;

  for (vertex_i = vertexs->begin(); vertex_i != vertexs->end(); vertex_i++){
    PrimVertexCand_x[nPrimVertexCand] = vertex_i->x();
    PrimVertexCand_y[nPrimVertexCand] = vertex_i->y();
    PrimVertexCand_z[nPrimVertexCand] = vertex_i->z();
    PrimVertexCand_tracks[nPrimVertexCand] = vertex_i->tracksSize();
    PrimVertexCand_chi2[nPrimVertexCand] = vertex_i->chi2();
    PrimVertexCand_ndof[nPrimVertexCand] = vertex_i->ndof();
    nPrimVertexCand++;
  }


  // Get the CASTOR towers collection from the event 
  edm::Handle<reco::CastorTowerCollection> recoCastorTowers;  
  event.getByLabel(recCastorTowerLabel, recoCastorTowers);  

  // Get the ZDC rechits collection from the event
  edm::Handle<ZDCRecHitCollection> recoZDChits;
  event.getByLabel(recZDCRecHitsLabel, recoZDChits);
  const ZDCRecHitCollection* zdchits = recoZDChits.product();
  ZDCRecHitCollection::const_iterator zdchit;

  // Get the PFlow collection from the event
  edm::Handle<reco::PFCandidateCollection> pflows;
  event.getByLabel("particleFlow",pflows);
  reco::PFCandidateCollection::const_iterator pflow;

  double highestejet = -1.0;
  double highestejeteta = -999.0;
  double highestejetphi = -999.0;
  double totalejet = -1.0;
  double highestetower = -1.0; 
  double highestetowerdr = -999.0;
  double highestetowereta = -999.0;
  double highestetowerphi = -999.0;
  double highestettower = -1.0; 
  double highestettowerdr = -999.0;
  double highestettowereta = -999.0;
  double highestettowerphi = -999.0;
  double highestcastortowerfwd = -999.0; 
  double highestcastortowerbwd = -999.0; 
  double totalecastorfwd = 0.0;
  double totalecastorbwd = 0.0;
  double totalecalo = -1.0; 
  double closesttrkdxyz = 999.0;

  // If this event contains a di-mu/e/gamma candidate, look at Jets & MET & CaloTowers & Tracks
  if(nMuonCand == 2)
    {
      for ( jet = jets->begin(); jet != jets->end() && nJetCand<JETMAX; ++jet )
	{
	  JetCand_e[nJetCand]=jet->energy();
	  JetCand_px[nJetCand]=jet->px();
	  JetCand_py[nJetCand]=jet->py();
	  JetCand_pz[nJetCand]=jet->pz();
	  JetCand_phi[nJetCand]=jet->phi();
	  JetCand_eta[nJetCand]=jet->eta();

	  totalejet = totalejet + JetCand_e[nJetCand];
	  if(JetCand_e[nJetCand] > highestejet)
	    {
	      highestejet = JetCand_e[nJetCand];
	      highestejeteta = JetCand_eta[nJetCand];
	      highestejetphi = JetCand_phi[nJetCand];
	    }
	  nJetCand++;
	}
      HighestJet_e = highestejet;
      HighestJet_eta = highestejeteta;
      HighestJet_phi = highestejetphi;
      SumJet_e = totalejet;
      met = mets->begin();
      float e_met = met->energy();
      Etmiss = e_met;
      for (calo = towers->begin(); calo != towers->end(); ++calo )
	{
	  CaloTower_e[nCaloCand]=calo->energy();
	  CaloTower_et[nCaloCand]=calo->et();
	  CaloTower_phi[nCaloCand]=calo->phi(); 
	  CaloTower_eta[nCaloCand]=calo->eta(); 
	  CaloTower_emE[nCaloCand]=calo->emEnergy();
	  CaloTower_hadE[nCaloCand]=calo->hadEnergy();  //hcal
	  CaloTower_outE[nCaloCand]=calo->outerEnergy(); //ho
	  GlobalPoint emPosition=calo->emPosition();
          GlobalPoint hadPosition=calo->hadPosition();
	  if(CaloTower_emE[nCaloCand]>0){ //if ECAL
		CaloTower_x[nCaloCand]=emPosition.x();
                CaloTower_y[nCaloCand]=emPosition.y();
                CaloTower_z[nCaloCand]=emPosition.z();
                CaloTower_t[nCaloCand]=calo->ecalTime();
		if(fabs(CaloTower_eta[nCaloCand])>=0   && fabs(CaloTower_eta[nCaloCand])<1.4) {CaloTower_ID[nCaloCand]=1;} //cout<<" -> EB"<<endl;}//EB
		if(fabs(CaloTower_eta[nCaloCand])>=1.4 && fabs(CaloTower_eta[nCaloCand])<2.95){CaloTower_ID[nCaloCand]=2;} //cout<<" -> EE"<<endl;}//EE
		if(fabs(CaloTower_eta[nCaloCand])>=2.95 && fabs(CaloTower_eta[nCaloCand])<5.2){CaloTower_ID[nCaloCand]=3;} //cout<<" -> EF"<<endl;}//EF ??
	  }

	  else if(CaloTower_hadE[nCaloCand]>0){
                CaloTower_x[nCaloCand]=hadPosition.x();
                CaloTower_y[nCaloCand]=hadPosition.y();
                CaloTower_z[nCaloCand]=hadPosition.z();
                CaloTower_t[nCaloCand]=calo->hcalTime();
		if(fabs(CaloTower_eta[nCaloCand])>=0   && fabs(CaloTower_eta[nCaloCand])<1.4) {CaloTower_ID[nCaloCand]=4;} // cout<<" -> HB"<<endl;}//HB
		if(fabs(CaloTower_eta[nCaloCand])>=1.4 && fabs(CaloTower_eta[nCaloCand])<2.95) {CaloTower_ID[nCaloCand]=5;} // cout<<" -> HE"<<endl;}//HE
		if(fabs(CaloTower_eta[nCaloCand])>=2.95 && fabs(CaloTower_eta[nCaloCand])<5.2) {CaloTower_ID[nCaloCand]=6;} // cout<<" -> HF"<<endl;}//HF ??
	  }

	  CaloTower_badhcalcells[nCaloCand]=calo->numBadHcalCells();
	  CaloTower_problemhcalcells[nCaloCand]=calo->numProblematicHcalCells(); 
          CaloTower_badecalcells[nCaloCand]=calo->numBadEcalCells(); 
          CaloTower_problemecalcells[nCaloCand]=calo->numProblematicEcalCells();  
	  
	  float calodr1 = sqrt(((CaloTower_eta[nCaloCand]-MuonCand_eta[0])*(CaloTower_eta[nCaloCand]-MuonCand_eta[0])) + 
			       ((CaloTower_phi[nCaloCand]-MuonCand_phi[0])*(CaloTower_phi[nCaloCand]-MuonCand_phi[0])));
	  float calodr2 = sqrt(((CaloTower_eta[nCaloCand]-MuonCand_eta[1])*(CaloTower_eta[nCaloCand]-MuonCand_eta[1])) + 
			       ((CaloTower_phi[nCaloCand]-MuonCand_phi[1])*(CaloTower_phi[nCaloCand]-MuonCand_phi[1])));
	  
	  if(calodr1 < calodr2)
	    CaloTower_dr[nCaloCand] = calodr1;
	  else
	    CaloTower_dr[nCaloCand] = calodr2;

	  totalecalo = totalecalo + CaloTower_e[nCaloCand]; 

	  if(CaloTower_e[nCaloCand] > highestetower) 
	    {
	      highestetower = CaloTower_e[nCaloCand]; 
	      highestetowereta = CaloTower_eta[nCaloCand];
	      highestetowerphi = CaloTower_phi[nCaloCand];
	      highestetowerdr = CaloTower_dr[nCaloCand];
	    }

	  if(CaloTower_et[nCaloCand] > highestettower) 
	    {
	      highestettower = CaloTower_et[nCaloCand]; 
	      highestettowereta = CaloTower_eta[nCaloCand];
	      highestettowerphi = CaloTower_phi[nCaloCand];
	      highestettowerdr = CaloTower_dr[nCaloCand];
	    }


	  if(CaloTower_dr[nCaloCand] > drisocalo)
	    {
              if(CaloTower_e[nCaloCand] > 1.0) 
                nExtraCaloTowersE1++; 
              if(CaloTower_e[nCaloCand] > 2.0)  
                nExtraCaloTowersE2++;  
              if(CaloTower_e[nCaloCand] > 3.0)  
                nExtraCaloTowersE3++;  
              if(CaloTower_e[nCaloCand] > 4.0)  
                nExtraCaloTowersE4++;  
              if(CaloTower_e[nCaloCand] > 5.0)  
                nExtraCaloTowersE5++;  
	      if(CaloTower_e[nCaloCand] > 6.0)   
                nExtraCaloTowersE6++;   
              if(CaloTower_e[nCaloCand] > 7.0)   
                nExtraCaloTowersE7++;   
              if(CaloTower_e[nCaloCand] > 8.0)   
                nExtraCaloTowersE8++;   
              if(CaloTower_e[nCaloCand] > 9.0)   
                nExtraCaloTowersE9++;   

	      
              if(CaloTower_et[nCaloCand] > 0.1) 
                nExtraCaloTowersEt0pt1++; 
              if(CaloTower_et[nCaloCand] > 0.2)  
                nExtraCaloTowersEt0pt2++;  
	      if(CaloTower_et[nCaloCand] > 0.5)  
                nExtraCaloTowersEt0pt5++;  
              if(CaloTower_et[nCaloCand] > 1.0)  
                nExtraCaloTowersEt1++;  
              if(CaloTower_et[nCaloCand] > 2.0)  
                nExtraCaloTowersEt2++;  
              if(CaloTower_et[nCaloCand] > 3.0)   
                nExtraCaloTowersEt3++;   
              if(CaloTower_et[nCaloCand] > 4.0)   
                nExtraCaloTowersEt4++;   

	 
              if(CaloTower_e[nCaloCand] > 0.0 && abs(CaloTower_eta[nCaloCand]) > 3.0) 	 
                nExtraCaloTowersE0hf++; 	 
              if(CaloTower_e[nCaloCand] > 1.0 && abs(CaloTower_eta[nCaloCand]) > 3.0) 	 
                nExtraCaloTowersE1hf++; 	 
              if(CaloTower_e[nCaloCand] > 2.0 && abs(CaloTower_eta[nCaloCand]) > 3.0) 	 
                nExtraCaloTowersE2hf++; 	 
	      if(CaloTower_e[nCaloCand] > 3.0 && abs(CaloTower_eta[nCaloCand]) > 3.0)   
                nExtraCaloTowersE3hf++;  
              if(CaloTower_e[nCaloCand] > 4.0 && abs(CaloTower_eta[nCaloCand]) > 3.0)    
                nExtraCaloTowersE4hf++;  
              if(CaloTower_e[nCaloCand] > 5.0 && abs(CaloTower_eta[nCaloCand]) > 3.0)    
                nExtraCaloTowersE5hf++; 

              if(CaloTower_e[nCaloCand] > 1.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5) 	 
                nExtraCaloTowersE1he++; 	 
              if(CaloTower_e[nCaloCand] > 2.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5) 	 
                nExtraCaloTowersE2he++; 	 
              if(CaloTower_e[nCaloCand] > 3.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5) 	 
                nExtraCaloTowersE3he++; 	 
	      if(CaloTower_e[nCaloCand] > 4.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5)    
                nExtraCaloTowersE4he++;   
              if(CaloTower_e[nCaloCand] > 5.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5)    
                nExtraCaloTowersE5he++;    

              if(CaloTower_e[nCaloCand] > 2.0 && abs(CaloTower_eta[nCaloCand]) < 1.5) 	 
                nExtraCaloTowersE2hb++; 	 
              if(CaloTower_e[nCaloCand] > 3.0 && abs(CaloTower_eta[nCaloCand]) < 1.5) 	 
                nExtraCaloTowersE3hb++; 	 
              if(CaloTower_e[nCaloCand] > 4.0 && abs(CaloTower_eta[nCaloCand]) < 1.5) 	 
                nExtraCaloTowersE4hb++;
              if(CaloTower_e[nCaloCand] > 5.0 && abs(CaloTower_eta[nCaloCand]) < 1.5)    
                nExtraCaloTowersE5hb++;    
	    }

	  nCaloCand++;
	}

      SumCalo_e = totalecalo;
      HighestCaloTower_e = highestetower;
      HighestCaloTower_eta = highestetowereta;
      HighestCaloTower_phi = highestetowerphi;
      HighestCaloTower_dr = highestetowerdr;
      HighestEtCaloTower_et = highestettower;
      HighestEtCaloTower_eta = highestettowereta;
      HighestEtCaloTower_phi = highestettowerphi;
      HighestEtCaloTower_dr = highestettowerdr;

      // Now CASTOR rechits
      edm::Handle<CastorRecHitCollection> recoCASTORhits;
      event.getByLabel(recCastorRecHitsLabel, recoCASTORhits);
      const CastorRecHitCollection* castorhits = recoCASTORhits.product();
      CastorRecHitCollection::const_iterator castorhit;
      for ( castorhit = castorhits->begin(); castorhit != castorhits->end(); ++castorhit )
	{
	  CASTORsumRecHitsE += castorhit->energy();
	}

      // Now CASTOR towers 
      if(recoCastorTowers.isValid()) 
	{ 
	  const CastorTowerCollection* castortowers = recoCastorTowers.product();   
	  CastorTowerCollection::const_iterator castortower;   
	  
	  for ( castortower = castortowers->begin(); castortower != castortowers->end(); ++castortower )  
	    { 
	      CastorTower_e[nCastorTowerCand] = castortower->energy(); 
	      CastorTower_eta[nCastorTowerCand] = castortower->eta();  
	      CastorTower_phi[nCastorTowerCand] = castortower->phi();  
	      CastorTower_emratio[nCastorTowerCand] = castortower->fem(); 
	      
	      if(CastorTower_eta[nCastorTowerCand] > 0)
		{
		  totalecastorfwd+=CastorTower_e[nCastorTowerCand];
		  if(CastorTower_e[nCastorTowerCand] > highestcastortowerfwd) 
		    highestcastortowerfwd = CastorTower_e[nCastorTowerCand];
		} 
	      if(CastorTower_eta[nCastorTowerCand] < 0) 
		{
		  totalecastorbwd+=CastorTower_e[nCastorTowerCand];
		  if(CastorTower_e[nCastorTowerCand] > highestcastortowerbwd)  
		    highestcastortowerbwd = CastorTower_e[nCastorTowerCand];  
		}
	      
	      nCastorTowerCand++;  
	    } 
	}

      HighestCastorTowerFwd_e = highestcastortowerfwd; 
      HighestCastorTowerBwd_e = highestcastortowerbwd; 
      SumCastorFwd_e = totalecastorfwd;
      SumCastorBwd_e = totalecastorbwd; 
    }

  // Now ZDC rechits
  for ( zdchit = zdchits->begin(); zdchit != zdchits->end(); ++zdchit )
    {
      HcalZDCDetId id(zdchit->id());
      int Side      = (zdchit->id()).zside();
      int Section   = (zdchit->id()).section();

      ZDChit_section[nZDChitCand] = Section;
      ZDChit_energy[nZDChitCand] = zdchit->energy();
      ZDChit_time[nZDChitCand] = zdchit->time();
      ZDChit_side[nZDChitCand] = Side;
      if((Section == 1) && (Side == 1))
        {
          ZDCsumEMplus = ZDCsumEMplus + ZDChit_energy[nZDChitCand];
        }
      if((Section == 1) && (Side == -1))
        {
          ZDCsumEMminus = ZDCsumEMminus + ZDChit_energy[nZDChitCand];
        }
      if((Section == 2) && (Side == 1))
        {
          ZDCsumHADplus = ZDCsumHADplus + ZDChit_energy[nZDChitCand];
        }
      if((Section == 2) && (Side == -1))
        {
          ZDCsumHADminus = ZDCsumHADminus + ZDChit_energy[nZDChitCand];
        }
      nZDChitCand++;
    }

  // Check for particles in ZDC/Castor acceptance. 
  // Use MC truth for now, replace with real RECO when available
  double MCPar_px,MCPar_py,MCPar_pz,MCPar_e,MCPar_eta,MCPar_mass;
  int MCPar_pdgid;

  Handle<GenParticleCollection> genParticles;
  event.getByLabel( "genParticles", genParticles );
  if(genParticles.isValid())
    {
      for ( size_t i = 0; i < genParticles->size(); ++ i ) 
	{
	  const Candidate & p = (*genParticles)[ i ];
	  MCPar_pdgid=p.pdgId();
	  MCPar_eta=p.eta();
	  MCPar_px=p.px();
	  MCPar_py=p.py();
	  MCPar_pz=p.pz();
	  MCPar_mass=p.mass();
	  MCPar_e = sqrt(MCPar_mass*MCPar_mass + (MCPar_px*MCPar_px + MCPar_py*MCPar_py + MCPar_pz*MCPar_pz));
	  
	  if(MCPar_pdgid == 22)
	    {
	      if(p.status() == 1 && nGenPhotCand < GENPHOTONMAX)
		{
		  GenPhotCand_pt[nGenPhotCand]=p.pt();
		  GenPhotCand_eta[nGenPhotCand]=p.eta(); 
		  GenPhotCand_phi[nGenPhotCand]=p.phi(); 
		  nGenPhotCand++;
		}
	    }
	  
	  if(MCPar_pdgid == 13 || MCPar_pdgid == -13)
	    {
	      if(p.status() == 1 && nGenMuonCand < GENMUONMAX) 
		{ 
		  GenMuonCand_px[nGenMuonCand]=p.px(); 
		  GenMuonCand_py[nGenMuonCand]=p.py();  
		  GenMuonCand_pz[nGenMuonCand]=p.pz();  
		  nGenMuonCand++; 
		} 
	    }
	  
	  if(MCPar_pdgid == 22 && abs(MCPar_eta) > 8.6 && MCPar_e > 20.0) 
	    HitInZDC++;
	  if(MCPar_pdgid == 2112 && abs(MCPar_eta) > 8.6 && MCPar_e > 50.0)
	    HitInZDC++;
	  if((MCPar_pdgid != 22 && MCPar_pdgid != 2112) && (abs(MCPar_eta) > 5.2 && abs(MCPar_eta) < 6.6))
	    HitInCastor++;
	  
	  if(MCPar_pdgid == 2212 && MCPar_pz > 3000.0) 
	    { 
	      double MCPar_pt = sqrt(MCPar_px*MCPar_px + MCPar_py*MCPar_py); 
	      double phi = p.phi(); 
	      double mp = 0.938272029; 
	      // ... compute kinimatical variable  
	      
	      float xi  = 1.0;    // fractional momentum loss  
	      if (MCPar_pz>0)  
		xi -= MCPar_pz/7000.0;  
	      else  
		xi += MCPar_pz/7000.0;  
	      
	      double t   = (-MCPar_pt*MCPar_pt - mp*mp*xi*xi) / (1-xi); // "t"  
	      
	      float acc420b1, acc220b1, acc420and220b1, acc420or220b1; // beam 1 (clockwise)  
	      float acc420b2, acc220b2, acc420and220b2, acc420or220b2; // beam 2 (anti-clockwise)  
	      
	      acc420b1 = acc220b1 = acc420and220b1 = acc420or220b1 = 0;  
	      acc420b2 = acc220b2 = acc420and220b2 = acc420or220b2 = 0;  
	      
	      
	      if(MCPar_pz > 0) 
		{ 
		  acc420b1       = helper420beam1.GetAcceptance(t, xi, phi);  
		  acc220b1       = helper220beam1.GetAcceptance(t, xi, phi);  
		  acc420and220b1 = helper420a220beam1.GetAcceptance(t, xi, phi);  
		  acc420or220b1  = acc420b1 + acc220b1 - acc420and220b1;  
		} 
	      else 
		{ 
		  acc420b2       = helper420beam2.GetAcceptance(t, xi, phi);  
		  acc220b2       = helper220beam2.GetAcceptance(t, xi, phi);  
		  acc420and220b2 = helper420a220beam2.GetAcceptance(t, xi, phi);  
		  acc420or220b2  = acc420b2 + acc220b2 - acc420and220b2;  
		} 
	      
	    } 
	  
	}
    }

  GenMuMu_eta = 0.0;
  GenMuMu_pt = 0.0;
  if(nGenMuonCand == 2)
    {
      TLorentzVector muvec1;
      TLorentzVector muvec2;
      TLorentzVector mumuvec;

      muvec1.SetXYZM(GenMuonCand_px[0],GenMuonCand_py[0],GenMuonCand_pz[0],0.1057);
      muvec2.SetXYZM(GenMuonCand_px[1],GenMuonCand_py[1],GenMuonCand_pz[1],0.1057); 
      mumuvec = muvec1 + muvec2;
      GenMuMu_eta = mumuvec.Eta();
      GenMuMu_pt = mumuvec.Pt();
    }

  // Now ParticleFlow photons 
  double leadingphotpx, leadingphotpy, leadingphotpz, leadingphotp; 
  for(pflow = pflows->begin(); pflow != pflows->end(); ++pflow) 
    { 
      int parttype = PFCandidate::ParticleType (pflow->particleId()); 
      if(parttype == 4 && nPFPhotonCand < PHOTONMAX) 
        { 
          PFPhotonCand_pt[nPFPhotonCand] = pflow->pt(); 
          PFPhotonCand_eta[nPFPhotonCand] = pflow->eta();  
          PFPhotonCand_phi[nPFPhotonCand] = pflow->phi();  
	  PFPhotonCand_drtrue[nPFPhotonCand] = -999.;

	  for(int ntruephot = 0; ntruephot < nGenPhotCand;ntruephot++)
	    {
	      double photdeta = (PFPhotonCand_eta[nPFPhotonCand]-GenPhotCand_eta[ntruephot]);
	      double photdphi = (PFPhotonCand_phi[nPFPhotonCand]-GenPhotCand_phi[ntruephot]); 
	      PFPhotonCand_drtrue[nPFPhotonCand] = sqrt((photdeta*photdeta) + (photdphi*photdphi));
	    }

          if(nPFPhotonCand == 0) 
            { 
              leadingphotpx = pflow->px(); 
              leadingphotpy = pflow->py();  
              leadingphotpz = pflow->pz();  
              leadingphotp = pflow->p();  
            } 

	  if(nMuonCand == 2)  
	    {  
	      double mmgmass = pow(MuonCand_p[0]+MuonCand_p[1]+pflow->p(),2);   
	      mmgmass-=pow(MuonCand_px[0]+MuonCand_px[1]+pflow->px(),2);   
	      mmgmass-=pow(MuonCand_py[0]+MuonCand_py[1]+pflow->py(),2);   
	      mmgmass-=pow(MuonCand_pz[0]+MuonCand_pz[1]+pflow->pz(),2);   
	      MuMuGamma_mass[nPFPhotonCand] = sqrt(mmgmass);   
	    }  

          nPFPhotonCand++; 
        } 
    } 
   
  // Now do vertexing and track counting
  edm::ESHandle<TransientTrackBuilder> theVtx;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theVtx);
  //  vector < reco::TransientTrack > mutrks;
  vector<TransientTrack> transmutrks; 
  reco::TrackCollection * mutrks = new reco::TrackCollection;

  // First get "muon" tracks
  bool isMuon = false;
  for( track = tracks->begin(); track != tracks->end(); ++ track ) 
    { 
      isMuon = false;
      for(int j = 0;j < nMuonCand; j++)
	{
	  if(MuonCandTrack_p[j] == track->p())
	    {/*
              cout << "µ vertex Z"<<j << " = " << track->vertex().z() << endl;*/
	      isMuon = true;
	      mutrks->push_back( *track );
	      TransientTrack tmptrk = (*theVtx).build( *track );
	      transmutrks.push_back( tmptrk );
	    }
	}
    }

  // If 2 muons, make a vertex
  if(transmutrks.size() == 2) 
   {
      KalmanVertexFitter fitter(true); 
      TransientVertex mumuVertex = fitter.vertex(transmutrks); 
      if(mumuVertex.isValid())
	{
	  MuMu_vtxx = mumuVertex.position().x(); 
	  MuMu_vtxy = mumuVertex.position().y(); 
	  MuMu_vtxz = mumuVertex.position().z(); 
	  MuMu_vtxchi2dof = mumuVertex.normalisedChiSquared();
	  MuMu_vtxT = sqrt(mumuVertex.position().x()*mumuVertex.position().x() + mumuVertex.position().y()*mumuVertex.position().y() ); 
	  MuMu_vtxisvalid = 1;/*
	  cout << "µµ vertex Z = " << MuMu_vtxz << endl;
          cout << "          x²= " << MuMu_vtxchi2dof << endl;
          cout << "          reffited  = " << mumuVertex.hasRefittedTracks() << endl;
          cout << "          weight_tk = " << mumuVertex.hasTrackWeight() << endl;
          cout << "          error ZZ   = " << mumuVertex.positionError().czz() << endl;
          cout << "          error ZX   = " << mumuVertex.positionError().czx() << endl;
          cout << "          error ZY   = " << mumuVertex.positionError().czy() << endl;
          cout << "          error XX   = " << mumuVertex.positionError().cxx() << endl;
          cout << "          error YY   = " << mumuVertex.positionError().cyy() << endl;
	  */
	}
      else
	{
	  MuMu_vtxx = 99;  
	  MuMu_vtxy = 99;  
	  MuMu_vtxz = 99;  
	  MuMu_vtxT = 99;
	  MuMu_vtxchi2dof = 9999;
	  MuMu_vtxisvalid = 0;/*
          cout << "µµ vertex X is not valid " << endl;*/
	}
//=================================================
/*      for (vertex_i = vertexs->begin(); vertex_i != vertexs->end(); vertex_i++){
	cout << "?? vertex X = " << vertex_i->x() 
                     << ", Y = " << vertex_i->y() 
                     << ", Z = " << vertex_i->z() << endl;}*/
//===================================

/*
      vector<Int_t> vec_isPU;
      vector<Int_t>::iterator ivec_isPU=vec_isPU.begin();
      for ( ; ivec_isPU!=vec_isPU.end(); ivec_isPU++){
	vec_isPU.push_back(0);
      }
       //  OK, now go back and count "extra" tracks on the dimuon vertex
      // Loop1 = tag a track as off-time PU-track
      for(track = tracks->begin(); track != tracks->end() && nTrackCand<TRACKMAX; ++ track) 
        { 
	  if(track->p() == MuonCandTrack_p[0] || track->p() == MuonCandTrack_p[1])  continue;

      //Reject off-time PU events
          bool is_PU = false;
          bool is_assovtx = false;
          Double_t extratrkZ = track->vertex().z();
//   cout << "extra track @ z = " << extratrkZ << endl;
          for (vertex_i = vertexs->begin(); vertex_i != vertexs->end(); vertex_i++){
             Double_t extravtxZ = vertex_i->z();
//           if the extra track is associated to a reco vertex
             if(sqrt( (extravtxZ - extratrkZ)*(extravtxZ - extratrkZ) ) < 0.1){
		is_assovtx = true;
//   cout << "   > is associated " << endl;
//              if this vertex is far from the µµ vertex and then could be considered as PU vtx
                if(sqrt( (extravtxZ - MuMu_vtxz)*(extravtxZ - MuMu_vtxz)) > 5.) is_PU = true;
             }
          }
// in this sample, there are still PU tracks, but not associated to a RecoVertex


          if(!(is_PU)) vec_isPU.push_back(0);
	  else         vec_isPU.push_back(1);
      }
*/

      //      vector<Int_t>::iterator ivec_isPU_read=vec_isPU.begin();
      // OK, now go back and count "extra" tracks on the dimuon vertex
      // Loop2 = compute "track" quantities
      for(track = tracks->begin(); track != tracks->end() && nTrackCand<TRACKMAX; ++ track)
        {
          if(track->p() == MuonCandTrack_p[0] || track->p() == MuonCandTrack_p[1]) continue;

	  //          if(*ivec_isPU_read == 1) {ivec_isPU_read++; continue;}

	  //          ivec_isPU_read++;
          TrackCand_purity[nTrackCand]=track->quality(TrackBase::highPurity); 
          TrackCand_p[nTrackCand]=track->p(); 
          TrackCand_px[nTrackCand]=track->px(); 
          TrackCand_py[nTrackCand]=track->py(); 
          TrackCand_pz[nTrackCand]=track->pz(); 
          TrackCand_pt[nTrackCand]=track->pt(); 
          TrackCand_eta[nTrackCand]=track->eta(); 
          TrackCand_phi[nTrackCand]=track->phi(); 
//          TrackCand_charge[nTrackCand]=track->charge();
	  TrackCand_nhits[nTrackCand]=track->numberOfValidHits();
	  TrackCand_vtxdxyz[nTrackCand] = sqrt(((track->vertex().x() - MuMu_vtxx)*(track->vertex().x() - MuMu_vtxx)) + 
					     ((track->vertex().y() - MuMu_vtxy)*(track->vertex().y() - MuMu_vtxy)) +
					     ((track->vertex().z() - MuMu_vtxz)*(track->vertex().z() - MuMu_vtxz)));
	  TrackCand_vtxT[nTrackCand] = sqrt(((track->vertex().x() - MuMu_vtxx)*(track->vertex().x() - MuMu_vtxx)) +
                                             ((track->vertex().y() - MuMu_vtxy)*(track->vertex().y() - MuMu_vtxy)));
	  TrackCand_vtxZ[nTrackCand] = sqrt(((track->vertex().z() - MuMu_vtxz)*(track->vertex().z() - MuMu_vtxz)));

	  if((TrackCand_purity[nTrackCand] == 2) && (TrackCand_nhits[nTrackCand] >= 3))
	    nQualityTrackCand++;


          if(TrackCand_vtxdxyz[nTrackCand] < 0.1)
            MuMu_extratracks1mm++;
          if(TrackCand_vtxdxyz[nTrackCand] < 0.3)
            MuMu_extratracks3mm++;
	  if(TrackCand_vtxdxyz[nTrackCand] < 0.5) 
            MuMu_extratracks5mm++; 
	  if(TrackCand_vtxdxyz[nTrackCand] < 1) 
            MuMu_extratracks1cm++; 
	  if(TrackCand_vtxdxyz[nTrackCand] < 3)
	    MuMu_extratracks3cm++;
          if(TrackCand_vtxdxyz[nTrackCand] < 5) 
            MuMu_extratracks5cm++; 
          if(TrackCand_vtxdxyz[nTrackCand] < 10) 
            MuMu_extratracks10cm++; 
	  if(TrackCand_vtxdxyz[nTrackCand] < closesttrkdxyz)
	    closesttrkdxyz = TrackCand_vtxdxyz[nTrackCand];
          nTrackCand++;  
      } 
      ClosestExtraTrack_vtxdxyz = closesttrkdxyz;
      //      if(ClosestExtraTrack_vtxdxyz < 0.1){cout << "   --> closest < 1 mm !!! : "<< ClosestExtraTrack_vtxdxyz << endl;}
    } 
  else 
    { 
      MuMu_vtxx = 99; 
      MuMu_vtxy = 99; 
      MuMu_vtxz = 99; 
      MuMu_vtxT = 99;
      MuMu_vtxchi2dof = 9999; 
      MuMu_vtxisvalid = 0; 
    } 

  // Check for di-objects with valid vertex
  if(nMuonCand != 2 || MuMu_vtxisvalid != 1)
    passed = false;
  else
    {
      if(MuMu_dphi < mudphimin) 
	passed = false;
      if(fabs(MuonCand_pt[0]-MuonCand_pt[1]) > mudptmax)
	passed = false;      
      if((keepsamesign == false) && (MuonCand_charge[0]*MuonCand_charge[1] > 0))
	passed = false;
    }

  // "Exclusivity" cuts
  if(passed == true){
    thetree->Fill();
    //    cout << "   --> SAVED" << endl;
  }
}

void
GammaGammaMuMu::fillDescriptions(ConfigurationDescriptions & descriptions) {
  
  descriptions.setComment("Exclusive dimuon EDAnalyzer.");
  
  edm::ParameterSetDescription iDesc;  

  iDesc.add<edm::InputTag>("GlobalMuonCollectionLabel", edm::InputTag("selectedLayer1Muons"))->setComment("input muon collection");
  iDesc.add<edm::InputTag>("CaloTowerLabel", edm::InputTag("towerMaker"))->setComment("input calo tower collection"); 
  iDesc.add<edm::InputTag>("RecoTrackLabel", edm::InputTag("generalTracks"))->setComment("input track collection"); 
  iDesc.add<edm::InputTag>("RecoVertexLabel", edm::InputTag("offlinePrimaryVertices"))->setComment("input vertex collection"); 
  iDesc.add<edm::InputTag>("CastorTowerLabel", edm::InputTag("CastorFastTowerReco"))->setComment("input CASTOR tower collection"); 
  iDesc.add<edm::InputTag>("ZDCRecHitsLabel", edm::InputTag("zdchits"))->setComment("input ZDC rechit collection");
  iDesc.add<edm::InputTag>("CastorRecHitsLabel", edm::InputTag("castorreco"))->setComment("input CASTOR rechit collection");
  iDesc.add<edm::InputTag>("JetCollectionLabel", edm::InputTag("selectedLayer1Jets"))->setComment("input jet collection"); 
  iDesc.add<edm::InputTag>("ElectronCollectionLabel", edm::InputTag("selectedLayer1Electrons"))->setComment("input electron collection"); 
  iDesc.add<edm::InputTag>("PhotonCollectionLabel", edm::InputTag("selectedLayer1Photons"))->setComment("input photon collection"); 
  iDesc.add<edm::InputTag>("MetLabel", edm::InputTag("selectedLayer1METs"))->setComment("input MET collection");   
  iDesc.add<double>("CaloTowerdR", 0.3)->setComment("Minimum delta-R to use for finding extra towers");  
  iDesc.add<double>("DimuonMindphi", 0.0)->setComment("Minimum delta-phi of dimuon pair");  
  iDesc.add<double>("DimuonMaxdpt", 2000.0)->setComment("Maximum delta-pT of dimuon pair");  
  iDesc.add<bool>("KeepSameSignDimuons", false)->setComment("Set to true to keep same-sign dimuon combinations");
  iDesc.addOptionalUntracked<std::string>("outfilename", ("mumu.pat.root"))->setComment("output flat ntuple file name");  
  iDesc.add<std::string>("HLTMenuLabel", ("HLT8E29"))->setComment("HLT AOD trigger summary label");
  iDesc.add< vector<std::string> >("AlgoNames")->setComment("Tag-and-probe algorithm names");

  descriptions.add("ParameterDescriptionsForGammaGammaMuMu", iDesc);
}

// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaMuMu::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaMuMu::endJob() {
  thefile->Write();
  thefile->Close();
}
  
