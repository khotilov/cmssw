// -*- C++ -*-
//
// Package:    GammaGammaMuE
// Class:      GammaGammaMuE
// 
/**\class GammaGammaMuE GammaGammaMuE.cc GammaGammaLeptonLepton/GammaGammaMuE/src/GammaGammaMuE.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jonathan Hollar
//         Created:  Wed Sep 20 10:08:38 BST 2006
// $Id: GammaGammaMuE.cc,v 1.6 2011/07/29 14:27:30 jjhollar Exp $
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

#include "FWCore/Framework/interface/Frameworkfwd.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterDescriptionNode.h"
//#include "FWCore/Common/interface/TriggerNames.h" 
#include "FWCore/Common/interface/TriggerNames.h" 
   
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "DataFormats/RecoCandidate/interface/CaloRecHitCandidate.h"
#include "DataFormats/HcalRecHit/interface/ZDCRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"

#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

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
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h" 
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

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/GammaGammaMuE.h"

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
GammaGammaMuE::GammaGammaMuE(const edm::ParameterSet& pset)
{
   //now do what ever initialization is needed
  recTrackLabel      = pset.getParameter<edm::InputTag>("RecoTrackLabel");
  recVertexLabel     = pset.getParameter<edm::InputTag>("RecoVertexLabel");
  theGLBMuonLabel    = pset.getParameter<edm::InputTag>("GlobalMuonCollectionLabel");
  thePixelGsfELabel  = pset.getParameter<edm::InputTag>("ElectronCollectionLabel");
  theJetLabel        = pset.getParameter<edm::InputTag>("JetCollectionLabel");
  theMetLabel        = pset.getParameter<edm::InputTag>("MetLabel");
  theCaloTowLabel    = pset.getParameter<edm::InputTag>("CaloTowerLabel");
  recCastorTowerLabel = pset.getParameter<edm::InputTag>("CastorTowerLabel"); 
  recZDCRecHitsLabel = pset.getParameter<edm::InputTag>("ZDCRecHitsLabel");
  recCastorRecHitsLabel = pset.getParameter<edm::InputTag>("CastorRecHitsLabel");
  hltMenuLabel       = pset.getParameter<std::string>("HLTMenuLabel");

  drisocalo          = pset.getParameter<double>("CaloTowerdR");
  keepsamesign       = pset.getParameter<bool>("KeepSameSign");
  minmuevtxd        = pset.getParameter<double>("MinMuEVertexSeparation"); 

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
  ELEMAX=10;
  JETMAX=30;
  TRACKMAX=500;
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
  thetree->Branch("MuonCand_vtxx",MuonCand_vtxx,"MuonCand_vtxx[nMuonCand]/D"); 
  thetree->Branch("MuonCand_vtxy",MuonCand_vtxy,"MuonCand_vtxy[nMuonCand]/D"); 
  thetree->Branch("MuonCand_vtxz",MuonCand_vtxz,"MuonCand_vtxz[nMuonCand]/D"); 
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
  thetree->Branch("MuonCand_gmPromptTight", MuonCand_gmPromptTight, "MuonCand_gmPromptTight[nMuonCand]/I");
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
  thetree->Branch("MuonCand_validpixelhits", MuonCand_validpixelhits, "MuonCand_validpixelhits[nMuonCand]/I");
  thetree->Branch("MuonCand_validmuonhits", MuonCand_validmuonhits, "MuonCand_validmuonhits[nMuonCand]/I");
  thetree->Branch("MuonCand_matches", MuonCand_matches, "MuonCand_matches[nMuonCand]/I"); 
  thetree->Branch("MuonCand_normchi2", MuonCand_normchi2, "MuonCand_normchi2[nMuonCand]/D");        
  thetree->Branch("MuonCand_normtrackchi2", MuonCand_normtrackchi2, "MuonCand_normtrackchi2[nMuonCand]/D");         
  thetree->Branch("MuonCand_dB", MuonCand_dB, "MuonCand_dB[nMuonCand]/D");
  thetree->Branch("MuonCand_tightID", MuonCand_tightID, "MuonCand_tightID[nMuonCand]/I");  
  thetree->Branch("MuEPairCand",MuEPairCand,"MuEPairCand[2]/I");

  thetree->Branch("nEleCand",&nEleCand,"nEleCand/I");
  thetree->Branch("EleCand_px",EleCand_px,"EleCand_px[nEleCand]/D");
  thetree->Branch("EleCand_py",EleCand_py,"EleCand_py[nEleCand]/D");
  thetree->Branch("EleCand_pz",EleCand_pz,"EleCand_pz[nEleCand]/D");
  thetree->Branch("EleCand_p",EleCand_p,"EleCand_p[nEleCand]/D");
  thetree->Branch("EleCand_e",EleCand_e,"EleCand_e[nEleCand]/D");
  thetree->Branch("EleCand_et",EleCand_et,"EleCand_et[nEleCand]/D");
  thetree->Branch("EleCand_eta",EleCand_eta,"EleCand_eta[nEleCand]/D");
  thetree->Branch("EleCand_phi",EleCand_phi,"EleCand_phi[nEleCand]/D");
  thetree->Branch("EleCand_charge",EleCand_charge,"EleCand_charge[nEleCand]/I");
  thetree->Branch("EleCand_looseid",EleCand_looseid,"EleCand_looseid[nEleCand]/I");
  thetree->Branch("EleCand_likelihoodid",EleCand_likelihoodid,"EleCand_likelihoodid[nEleCand]/D");  
  thetree->Branch("EleCand_robustid",EleCand_robustid,"EleCand_robustid[nEleCand]/I"); 
  thetree->Branch("EleCandTrack_p",EleCandTrack_p,"EleCandTrack_p[nEleCand]/D");
  thetree->Branch("EleCandTrack_pt",EleCandTrack_pt,"EleCandTrack_pt[nEleCand]/D"); 
  thetree->Branch("EleCandTrack_eta",EleCandTrack_eta,"EleCandTrack_eta[nEleCand]/D"); 
  thetree->Branch("EleCandTrack_phi",EleCandTrack_phi,"EleCandTrack_phi[nEleCand]/D"); 
  thetree->Branch("EleCand_vtxx",EleCand_vtxx,"EleCand_vtxx[nEleCand]/D");  
  thetree->Branch("EleCand_vtxy",EleCand_vtxy,"EleCand_vtxy[nEleCand]/D");  
  thetree->Branch("EleCand_vtxz",EleCand_vtxz,"EleCand_vtxz[nEleCand]/D");  
  thetree->Branch("EleCand_deltaPhi",EleCand_deltaPhi,"EleCand_deltaPhi[nEleCand]/D"); 
  thetree->Branch("EleCand_deltaEta",EleCand_deltaEta,"EleCand_deltaEta[nEleCand]/D"); 
  thetree->Branch("EleCand_HoverE",EleCand_HoverE,"EleCand_HoverE[nEleCand]/D"); 
  thetree->Branch("EleCand_trackiso",EleCand_trackiso,"EleCand_trackiso[nEleCand]/D"); 
  thetree->Branch("EleCand_ecaliso",EleCand_ecaliso,"EleCand_ecaliso[nEleCand]/D"); 
  thetree->Branch("EleCand_hcaliso",EleCand_hcaliso,"EleCand_hcaliso[nEleCand]/D"); 
  thetree->Branch("EleCand_sigmaIetaIeta",EleCand_sigmaIetaIeta,"EleCand_sigmaIetaIeta[nEleCand]/D"); 
  thetree->Branch("EleCand_convDist",EleCand_convDist,"EleCand_convDist[nEleCand]/D"); 
  thetree->Branch("EleCand_convDcot",EleCand_convDcot,"EleCand_convDcot[nEleCand]/D"); 
  thetree->Branch("EleCand_ecalDriven",EleCand_ecalDriven,"EleCand_ecalDriven[nEleCand]/D");  
  thetree->Branch("EleCand_wp80", EleCand_wp80, "EleCand_wp80[nEleCand]/I");   

  thetree->Branch("nHLTMu10Ele10MuonCand",&nHLTMu10Ele10MuonCand,"nHLTMu10Ele10MuonCand/I"); 
  thetree->Branch("HLT_Mu10Ele10_MuonCand_pt",&HLT_Mu10Ele10_MuonCand_pt,"HLT_Mu10Ele10_MuonCand_pt[nHLTMu10Ele10MuonCand]/D"); 
  thetree->Branch("HLT_Mu10Ele10_MuonCand_eta",&HLT_Mu10Ele10_MuonCand_eta,"HLT_Mu10Ele10_MuonCand_eta[nHLTMu10Ele10MuonCand]/D"); 
  thetree->Branch("HLT_Mu10Ele10_MuonCand_phi",&HLT_Mu10Ele10_MuonCand_phi,"HLT_Mu10Ele10_MuonCand_phi[nHLTMu10Ele10MuonCand]/D"); 
  thetree->Branch("HLT_Mu10Ele10_MuonCand_charge",&HLT_Mu10Ele10_MuonCand_charge,"HLT_Mu10Ele10_MuonCand_charge[nHLTMu10Ele10MuonCand]/I");  

  thetree->Branch("nHLTMu8Ele17MuonCand",&nHLTMu8Ele17MuonCand,"nHLTMu8Ele17MuonCand/I");  
  thetree->Branch("HLT_Mu8Ele17_MuonCand_pt",&HLT_Mu8Ele17_MuonCand_pt,"HLT_Mu8Ele17_MuonCand_pt[nHLTMu8Ele17MuonCand]/D");  
  thetree->Branch("HLT_Mu8Ele17_MuonCand_eta",&HLT_Mu8Ele17_MuonCand_eta,"HLT_Mu8Ele17_MuonCand_eta[nHLTMu8Ele17MuonCand]/D");  
  thetree->Branch("HLT_Mu8Ele17_MuonCand_phi",&HLT_Mu8Ele17_MuonCand_phi,"HLT_Mu8Ele17_MuonCand_phi[nHLTMu8Ele17MuonCand]/D");  
  thetree->Branch("HLT_Mu8Ele17_MuonCand_charge",&HLT_Mu8Ele17_MuonCand_charge,"HLT_Mu8Ele17_MuonCand_charge[nHLTMu8Ele17MuonCand]/I");   

  thetree->Branch("nHLTMu17Ele8MuonCand",&nHLTMu17Ele8MuonCand,"nHLTMu17Ele8MuonCand/I");  
  thetree->Branch("HLT_Mu17Ele8_MuonCand_pt",&HLT_Mu17Ele8_MuonCand_pt,"HLT_Mu17Ele8_MuonCand_pt[nHLTMu17Ele8MuonCand]/D");  
  thetree->Branch("HLT_Mu17Ele8_MuonCand_eta",&HLT_Mu17Ele8_MuonCand_eta,"HLT_Mu17Ele8_MuonCand_eta[nHLTMu17Ele8MuonCand]/D");  
  thetree->Branch("HLT_Mu17Ele8_MuonCand_phi",&HLT_Mu17Ele8_MuonCand_phi,"HLT_Mu17Ele8_MuonCand_phi[nHLTMu17Ele8MuonCand]/D");  
  thetree->Branch("HLT_Mu17Ele8_MuonCand_charge",&HLT_Mu17Ele8_MuonCand_charge,"HLT_Mu17Ele8_MuonCand_charge[nHLTMu17Ele8MuonCand]/I");   

  thetree->Branch("nCaloCand",&nCaloCand,"nCaloCand/I");
  thetree->Branch("CaloTower_e",CaloTower_e,"CaloTower_e[nCaloCand]/D");
  thetree->Branch("CaloTower_et",CaloTower_et,"CaloTower_et[nCaloCand]/D");
  thetree->Branch("CaloTower_eta",CaloTower_eta,"CaloTower_eta[nCaloCand]/D"); 
  thetree->Branch("CaloTower_phi",CaloTower_phi,"CaloTower_phi[nCaloCand]/D"); 
  thetree->Branch("CaloTower_dr",CaloTower_dr,"CaloTower_dr[nCaloCand]/D");
  thetree->Branch("CaloTower_emE",CaloTower_emE,"CaloTower_emE[nCaloCand]/D");
  thetree->Branch("CaloTower_hadE",CaloTower_hadE,"CaloTower_hadE[nCaloCand]/D");
  thetree->Branch("CaloTower_outE",CaloTower_outE,"CaloTower_outE[nCaloCand]/D");
  thetree->Branch("CaloTower_x",CaloTower_x,"CaloTower_x[nCaloCand]/D");
  thetree->Branch("CaloTower_y",CaloTower_y,"CaloTower_y[nCaloCand]/D");
  thetree->Branch("CaloTower_z",CaloTower_z,"CaloTower_z[nCaloCand]/D");
  thetree->Branch("CaloTower_t",CaloTower_t,"CaloTower_t[nCaloCand]/D");
  thetree->Branch("CaloTower_badhcalcells",CaloTower_badhcalcells,"CaloTower_badhcalcells[nCaloCand]/I"); 
  thetree->Branch("CaloTower_problemhcalcells",CaloTower_problemhcalcells,"CaloTower_problemhcalcells[nCaloCand]/I"); 
  thetree->Branch("CaloTower_badecalcells",CaloTower_badecalcells,"CaloTower_badecalcells[nCaloCand]/I");  
  thetree->Branch("CaloTower_problemecalcells",CaloTower_problemecalcells,"CaloTower_problemecalcells[nCaloCand]/I");  
  thetree->Branch("SumCalo_e",&SumCalo_e,"SumCalo_e/D");

  thetree->Branch("nCastorTowerCand",&nCastorTowerCand,"nCastorTowerCand/I");  
  thetree->Branch("nCastorTowerCandE3",&nCastorTowerCandE3,"nCastorTowerCandE3/I");   
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
  thetree->Branch("TrackCand_charge",TrackCand_charge,"TrackCand_charge[nTrackCand]/D"); 
  thetree->Branch("TrackCand_purity",TrackCand_purity,"TrackCand_purity[nTrackCand]/D");
  thetree->Branch("TrackCand_nhits",TrackCand_nhits,"TrackCand_nhits[nTrackCand]/I");
  thetree->Branch("TrackCand_chi2",TrackCand_chi2,"TrackCand_chi2[nTrackCand]/D");
  thetree->Branch("TrackCand_ndof",TrackCand_ndof,"TrackCand_ndof[nTrackCand]/D");

  thetree->Branch("TrackCand_vtxZ",TrackCand_vtxZ,"TrackCand_vtxZ[nTrackCand]/D");
  thetree->Branch("TrackCand_vtxT",TrackCand_vtxT,"TrackCand_vtxT[nTrackCand]/D");
  thetree->Branch("TrackCand_X",TrackCand_X,"TrackCand_X[nTrackCand]/D");
  thetree->Branch("TrackCand_Y",TrackCand_Y,"TrackCand_Y[nTrackCand]/D");
  thetree->Branch("TrackCand_Z",TrackCand_Z,"TrackCand_Z[nTrackCand]/D");
  thetree->Branch("ClosestExtraTrack_vtxdxyz",&ClosestExtraTrack_vtxdxyz,"ClosestExtraTrack_vtxdxyz/D");
  thetree->Branch("ClosestHighPurityExtraTrack_vtxdxyz",&ClosestHighPurityExtraTrack_vtxdxyz,"ClosestHighPurityExtraTrack_vtxdxyz/D");

  thetree->Branch("MuE_mass",&MuE_mass,"MuE_mass/D");
  thetree->Branch("MuE_dphi",&MuE_dphi,"MuE_dphi/D");
  thetree->Branch("MuE_dpt",&MuE_dpt,"MuE_dpt/D");
  thetree->Branch("MuE_pt",&MuE_pt,"MuE_pt/D"); 
  thetree->Branch("MuE_Kalmanvtxx",&MuE_Kalmanvtxx,"MuE_Kalmanvtxx/D");
  thetree->Branch("MuE_Kalmanvtxy",&MuE_Kalmanvtxy,"MuE_Kalmanvtxy/D"); 
  thetree->Branch("MuE_Kalmanvtxz",&MuE_Kalmanvtxz,"MuE_Kalmanvtxz/D");
  thetree->Branch("MuE_KalmanvtxT",&MuE_KalmanvtxT,"MuE_KalmanvtxT/D"); 
  thetree->Branch("MuE_Kalmanvtxchi2dof",&MuE_Kalmanvtxchi2dof,"MuE_Kalmanvtxchi2dof/D");
  thetree->Branch("MuE_Kalmanvtxisvalid",&MuE_Kalmanvtxisvalid,"MuE_Kalmanvtxisvalid/I");
  thetree->Branch("MuE_extratracks1mm",&MuE_extratracks1mm,"MuE_extratracks1mm/I");
  thetree->Branch("MuE_extratracks2mm",&MuE_extratracks2mm,"MuE_extratracks2mm/I"); 
  thetree->Branch("MuE_extratracks3mm",&MuE_extratracks3mm,"MuE_extratracks3mm/I");
  thetree->Branch("MuE_extratracks4mm",&MuE_extratracks4mm,"MuE_extratracks4mm/I"); 
  thetree->Branch("MuE_extratracks5mm",&MuE_extratracks5mm,"MuE_extratracks5mm/I");
  thetree->Branch("MuE_extratracks1cm",&MuE_extratracks1cm,"MuE_extratracks1cm/I");
  thetree->Branch("MuE_extratracks3cm",&MuE_extratracks3cm,"MuE_extratracks3cm/I");
  thetree->Branch("MuE_extratracks5cm",&MuE_extratracks5cm,"MuE_extratracks5cm/I"); 
  thetree->Branch("MuE_extratracks10cm",&MuE_extratracks10cm,"MuE_extratracks10cm/I"); 

  thetree->Branch("HitInZDC",&HitInZDC,"HitInZDC/I");
  thetree->Branch("HitInCastor",&HitInCastor,"HitInCastor/I");

  thetree->Branch("nGenMuonCand",&nGenMuonCand,"nGenMuonCand/I");  
  thetree->Branch("GenMuonCand_px",GenMuonCand_px,"GenMuonCand_px[nGenMuonCand]/D");  
  thetree->Branch("GenMuonCand_py",GenMuonCand_py,"GenMuonCand_py[nGenMuonCand]/D");   
  thetree->Branch("GenMuonCand_pz",GenMuonCand_pz,"GenMuonCand_pz[nGenMuonCand]/D");   

  thetree->Branch("GenMuE_eta",&GenMuE_eta,"GenMuE_eta/D");    
  thetree->Branch("GenMuE_pt",&GenMuE_pt,"GenMuE_pt/D");     
  
  thetree->Branch("Etmiss",&Etmiss,"Etmiss/D");
  thetree->Branch("Etmiss_phi",&Etmiss_phi,"Etmiss_phi/D"); 
  thetree->Branch("Etmiss_x",&Etmiss_x,"Etmiss_x/D"); 
  thetree->Branch("Etmiss_y",&Etmiss_y,"Etmiss_y/D"); 
  thetree->Branch("Etmiss_z",&Etmiss_z,"Etmiss_z/D"); 
  thetree->Branch("Etmiss_significance",&Etmiss_significance,"Etmiss_significance/D"); 

  thetree->Branch("HLT_Mu17Ele8",&HLT_Mu17Ele8,"HLT_Mu17Ele8/I");
  thetree->Branch("HLT_Mu8Ele17",&HLT_Mu8Ele17,"HLT_Mu8Ele17/I");
  thetree->Branch("HLT_Mu10Ele10",&HLT_Mu10Ele10,"HLT_Mu10Ele10/I");
  thetree->Branch("HLT_Mu17Ele8_Prescl",&HLT_Mu17Ele8_Prescl,"HLT_Mu17Ele8_Prescl/I"); 
  thetree->Branch("HLT_Mu8Ele17_Prescl",&HLT_Mu8Ele17_Prescl,"HLT_Mu8Ele17_Prescl/I"); 
  thetree->Branch("HLT_Mu10Ele10_Prescl",&HLT_Mu10Ele10_Prescl,"HLT_Mu10Ele10_Prescl/I"); 

  thetree->Branch("Run",&Run,"Run/I");
  thetree->Branch("LumiSection",&LumiSection,"LumiSection/I");
  thetree->Branch("BX",&BX,"BX/I");
  thetree->Branch("AvgInstDelLumi",&AvgInstDelLumi,"AvgInstDelLumi/D");
  thetree->Branch("BunchInstLumi",&BunchInstLumi,"BunchInstLumi[3]/D");
  thetree->Branch("EventNum",&EventNum,"EventNum/I");
  thetree->Branch("L1TechnicalTriggers",L1TechnicalTriggers,"L1TechnicalTriggers[128]/I"); 

  thetree->Branch("nPrimVertexCand",&nPrimVertexCand,"nPrimVertexCand/I");
  thetree->Branch("PrimVertexCand_x",&PrimVertexCand_x,"PrimVertexCand_x[nPrimVertexCand]/D");
  thetree->Branch("PrimVertexCand_y",&PrimVertexCand_y,"PrimVertexCand_y[nPrimVertexCand]/D");
  thetree->Branch("PrimVertexCand_z",&PrimVertexCand_z,"PrimVertexCand_z[nPrimVertexCand]/D");
  thetree->Branch("PrimVertexCand_tracks",&PrimVertexCand_tracks,"PrimVertexCand_tracks[nPrimVertexCand]/I");
  thetree->Branch("PrimVertexCand_chi2",&PrimVertexCand_chi2,"PrimVertexCand_chi2[nPrimVertexCand]/D");
  thetree->Branch("PrimVertexCand_ndof",&PrimVertexCand_ndof,"PrimVertexCand_ndof[nPrimVertexCand]/D");
  thetree->Branch("PrimVertexCand_mueTwoTracks",&PrimVertexCand_mueTwoTracks,"PrimVertexCand_mueTwoTracks[nPrimVertexCand]/I"); 
  thetree->Branch("PrimVertexCand_mueExactlyTwoTracks",&PrimVertexCand_mueExactlyTwoTracks,"PrimVertexCand_mueExactlyTwoTracks[nPrimVertexCand]/I");  
  thetree->Branch("PrimVertexCand_mueTwoTracksMuIndex",&PrimVertexCand_mueTwoTracksMuIndex,"PrimVertexCand_mueTwoTracksMuIndex[nPrimVertexCand]/I");  
  thetree->Branch("PrimVertexCand_mueTwoTracksEleIndex",&PrimVertexCand_mueTwoTracksEleIndex,"PrimVertexCand_mueTwoTracksEleIndex[nPrimVertexCand]/I");   


  thetree->Branch("LowPt_pt",LowPt_pt,"LowPt_pt[nMuonCand]/D");
  thetree->Branch("LowPt_eta",LowPt_eta,"LowPt_eta[nMuonCand]/D");

  thetree->Branch("nPU",&nPU,"nPU/I");
  thetree->Branch("PUWeight",&PUWeight,"PUWeight/D");

//  thetree->Branch("evweight",&evweight,"evweight/D");
}


GammaGammaMuE::~GammaGammaMuE()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GammaGammaMuE::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
  nPrimVertexCand=0;
  nMuonCand=0;
  nEleCand=0;
  nHLTMu10Ele10MuonCand=0;
  nHLTMu17Ele8MuonCand=0;  
  nHLTMu8Ele17MuonCand=0;
  nJetCand=0;
  nCaloCand=0;
  nTrackCand=0;
  nQualityTrackCand=0;
  nCastorTowerCand=0;
  nCastorTowerCandE3=0;
  nZDChitCand=0;
  ZDCsumHADminus=0;
  ZDCsumEMminus=0;
  ZDCsumHADplus=0;
  ZDCsumEMplus=0;

  HitInZDC=0;
  HitInCastor=0;
  nGenMuonCand=0;

  MuE_mass = -1;
  MuE_dphi = -1;
  MuE_dpt = -1;
  MuE_pt = -1;
  MuE_extratracks1mm = 0; 
  MuE_extratracks2mm = 0;  
  MuE_extratracks3mm = 0; 
  MuE_extratracks4mm = 0;  
  MuE_extratracks5mm = 0; 
  MuE_extratracks1cm = 0; 
  MuE_extratracks3cm = 0;
  MuE_extratracks5cm = 0;
  MuE_extratracks10cm = 0;
  ClosestExtraTrack_vtxdxyz = 999.;
  double mueprimvtxx = 0.0; 
  double mueprimvtxy = 0.0; 
  double mueprimvtxz = 0.0; 

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

  const edm::LuminosityBlock& iLumi = event.getLuminosityBlock();
  // get LumiSummary
  edm::Handle<LumiSummary> lumiSummary;
  iLumi.getByLabel("lumiProducer", lumiSummary);
  if(lumiSummary->isValid())
    AvgInstDelLumi = lumiSummary->avgInsDelLumi();
  else
    AvgInstDelLumi = -999.;

  BunchInstLumi[0] = -999.;
  BunchInstLumi[1] = -999.;
  BunchInstLumi[2] = -999.;

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
//  trigNames.init(*hltResults) ;
  const edm::TriggerNames & trigNames = event.triggerNames(*hltResults);

  for (unsigned int i=0; i<trigNames.size(); i++)  
    { 
      if ( trigNames.triggerNames().at(i).find("HLT_Mu10_Ele10_CaloIdL_") != string::npos)
        {  
	  HLT_Mu10Ele10_Prescl = hltConfig_.prescaleValue(event, iSetup, "HLT_Mu10Ele10"); 

          if ( hltResults->accept(i) )  
            {HLT_Mu10Ele10 = 1;}
	  else
	    HLT_Mu10Ele10 = 0;
        }  
      if ( trigNames.triggerNames().at(i).find("HLT_Mu8_Ele17_") != string::npos) 
        {   
          HLT_Mu8Ele17_Prescl = hltConfig_.prescaleValue(event, iSetup, "HLT_Mu8Ele17");  

          if ( hltResults->accept(i) )  
	    {HLT_Mu8Ele17 = 1;}
	  else
	    HLT_Mu8Ele17 = 0;
        }  
      if ( trigNames.triggerNames().at(i).find("HLT_Mu17_Ele8_") != string::npos) 
        {   
          HLT_Mu17Ele8_Prescl = hltConfig_.prescaleValue(event, iSetup, "HLT_Mu17Ele8");   

          if ( hltResults->accept(i) )  
	    {HLT_Mu17Ele8 = 1;}
	  else
	    HLT_Mu17Ele8 = 0;
        }  
    }

  Handle<TriggerEvent> hltObjects;
  event.getByLabel(InputTag("hltTriggerSummaryAOD","",hltMenuLabel),hltObjects);
  if (hltObjects.isValid()) 
    {
      size_type mu10ele10index = hltObjects->filterIndex(InputTag("hltL1Mu3EG5L3Filtered10::"+hltMenuLabel)); 
      size_type mu8ele17index = hltObjects->filterIndex(InputTag("hltL1MuOpenEG5L3Filtered8::"+hltMenuLabel)); 
      size_type mu17ele8index = hltObjects->filterIndex(InputTag("hltL1MuOpenEG5L3Filtered17::"+hltMenuLabel)); 

      if( mu8ele17index < hltObjects->sizeFilters() )
	{
	  const trigger::Keys& MU8ELE17KEYS(hltObjects->filterKeys(mu8ele17index)); 
	  const size_type nK(MU8ELE17KEYS.size()); 
	  const TriggerObjectCollection& TOC(hltObjects->getObjects());

	  for(int ipart = 0; ipart != nK; ++ipart)   
	    {
	      const TriggerObject& TO = TOC[MU8ELE17KEYS[ipart]];  
	      
	      if(fabs(TO.id()) == 13)
		{
		  HLT_Mu8Ele17_MuonCand_pt[nHLTMu8Ele17MuonCand] = TO.pt();
		  HLT_Mu8Ele17_MuonCand_eta[nHLTMu8Ele17MuonCand] = TO.eta(); 
		  HLT_Mu8Ele17_MuonCand_phi[nHLTMu8Ele17MuonCand] = TO.phi(); 
		  if(TO.id() > 0)	  
		    HLT_Mu8Ele17_MuonCand_charge[nHLTMu8Ele17MuonCand] = 1; 
		  else
		    HLT_Mu8Ele17_MuonCand_charge[nHLTMu8Ele17MuonCand] = -1;  
		  
		  nHLTMu8Ele17MuonCand++;
		}
	    }
	}
      if( mu8ele17index < hltObjects->sizeFilters() )
	{
	  const trigger::Keys& MU8ELE17KEYS(hltObjects->filterKeys(mu8ele17index)); 
	  const size_type nK(MU8ELE17KEYS.size()); 
	  const TriggerObjectCollection& TOC(hltObjects->getObjects());

	  for(int ipart = 0; ipart != nK; ++ipart)   
	    {
	      const TriggerObject& TO = TOC[MU8ELE17KEYS[ipart]];  
	      
	      if(fabs(TO.id()) == 13)
		{
		  HLT_Mu17Ele8_MuonCand_pt[nHLTMu17Ele8MuonCand] = TO.pt();
		  HLT_Mu17Ele8_MuonCand_eta[nHLTMu17Ele8MuonCand] = TO.eta(); 
		  HLT_Mu17Ele8_MuonCand_phi[nHLTMu17Ele8MuonCand] = TO.phi(); 
		  if(TO.id() > 0)	  
		    HLT_Mu17Ele8_MuonCand_charge[nHLTMu17Ele8MuonCand] = 1; 
		  else
		    HLT_Mu17Ele8_MuonCand_charge[nHLTMu17Ele8MuonCand] = -1;  
		  
		  nHLTMu17Ele8MuonCand++;
		}
	    }
	}
      if( mu10ele10index < hltObjects->sizeFilters() ) 
	{
          const trigger::Keys& MU10ELE10KEYS(hltObjects->filterKeys(mu10ele10index));  
          const size_type nK(MU10ELE10KEYS.size());  
          const TriggerObjectCollection& TOC(hltObjects->getObjects()); 
	  
          for(int ipart = 0; ipart != nK; ++ipart)    
            { 
              const TriggerObject& TO = TOC[MU10ELE10KEYS[ipart]];   
               
              if(fabs(TO.id()) == 13) 
                { 
                  HLT_Mu10Ele10_MuonCand_pt[nHLTMu10Ele10MuonCand] = TO.pt(); 
                  HLT_Mu10Ele10_MuonCand_eta[nHLTMu10Ele10MuonCand] = TO.eta();  
                  HLT_Mu10Ele10_MuonCand_phi[nHLTMu10Ele10MuonCand] = TO.phi();  
                  if(TO.id() > 0)          
                    HLT_Mu10Ele10_MuonCand_charge[nHLTMu10Ele10MuonCand] = 1;  
                  else 
                    HLT_Mu10Ele10_MuonCand_charge[nHLTMu10Ele10MuonCand] = -1;   
                   
                  nHLTMu10Ele10MuonCand++; 
                } 
            } 
	  
	}
    }


  // Get the #PU information
  PUWeight = 0.0;
  Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

  std::vector<PileupSummaryInfo>::const_iterator PVI;

  float sum_nvtx = 0.0;
  int npv = -1;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

    int BX = PVI->getBunchCrossing();

    cout << "PU rewighting - BX = " << BX << endl;

    npv = PVI->getPU_NumInteractions();

    sum_nvtx += float(npv);
    cout << "\tnpv = " << npv << ", sum = " << sum_nvtx << endl;
  }
  PUWeight = sum_nvtx/3.;

  nPU=0;
  edm::Handle<CrossingFrame<edm::HepMCProduct> > crossingFrameHepMCH;
  event.getByLabel("mix","source",crossingFrameHepMCH);
  
  if((crossingFrameHepMCH.isValid())) 
    {  
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

  for (muon = muons->begin(); muon != muons->end() && nMuonCand<MUONMAX; ++muon)
    {
      if((!muon->isTrackerMuon()) && (!muon->isGlobalMuon()))continue;
      if(nMuonCand>0 &&	muon->pt()==MuonCand_pt[nMuonCand-1] && muon->eta()==MuonCand_eta[nMuonCand-1]) continue;     
 
      MuonCand_p[nMuonCand]=muon->p();
      MuonCand_px[nMuonCand]=muon->px();
      MuonCand_py[nMuonCand]=muon->py();
      MuonCand_pz[nMuonCand]=muon->pz();
      MuonCand_pt[nMuonCand]=muon->pt();
      MuonCand_eta[nMuonCand]=muon->eta();
      MuonCand_phi[nMuonCand]=muon->phi();
      MuonCand_charge[nMuonCand]=muon->charge();
      MuonCand_vtxx[nMuonCand]=muon->vertex().x();
      MuonCand_vtxy[nMuonCand]=muon->vertex().y(); 
      MuonCand_vtxz[nMuonCand]=muon->vertex().z();
      
      // Muon ID - 31X compatible
      MuonCand_tmlsloosemuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMLastStationLoose);
      MuonCand_tmlsOptLowPtloosemuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMLastStationOptimizedLowPtLoose);
      MuonCand_tm2dloosemuid[nMuonCand]=muon::isGoodMuon(*muon, muon::TM2DCompatibilityLoose);
      MuonCand_tmlsAngloosemuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMLastStationAngLoose); 
      MuonCand_tmlsAngtightmuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMLastStationAngTight);  
      MuonCand_tmosAngloosemuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMOneStationAngLoose);  
      MuonCand_tmosAngtightmuonid[nMuonCand]=muon::isGoodMuon(*muon, muon::TMOneStationAngTight);  
      MuonCand_arbmuid[nMuonCand]=muon::isGoodMuon(*muon, muon::AllArbitrated);
      MuonCand_gmPromptTight[nMuonCand]=muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight);
      MuonCand_isglobal[nMuonCand]=muon->isGlobalMuon();
      MuonCand_istracker[nMuonCand]=muon->isTrackerMuon(); 
      MuonCand_isstandalone[nMuonCand]=muon->isStandAloneMuon();
      
      if(!(muon::isGoodMuon(*muon, muon::TMLastStationLoose)) && (muon::isGoodMuon(*muon, muon::TMLastStationOptimizedLowPtLoose))){	  
	LowPt_eta[nMuonCand]=muon->eta();
	LowPt_pt[nMuonCand]=muon->pt();
      }else{  LowPt_eta[nMuonCand]=10.;
	LowPt_pt[nMuonCand]=-1.;}
      
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
	  MuonCand_validpixelhits[nMuonCand] = muon->innerTrack()->hitPattern().numberOfValidPixelHits();	  
	  MuonCand_matches[nMuonCand] = muon->numberOfMatches();
	  MuonCand_dB[nMuonCand] = muon->dB();
	  MuonCand_tightID[nMuonCand] = 0;

	  if(muon->isGlobalMuon())
	    {
	      MuonCand_normchi2[nMuonCand]=muon->normChi2(); 
	      MuonCand_validmuonhits[nMuonCand] = muon->outerTrack()->hitPattern().numberOfValidMuonHits(); 

	      if((MuonCand_istracker[nMuonCand] == 1) && 
		 (MuonCand_validpixelhits[nMuonCand] >= 1) && 
		 (MuonCand_normchi2[nMuonCand] < 10) && 
		 (MuonCand_validtrackhits[nMuonCand] > 10) && 
		 (MuonCand_matches[nMuonCand] >= 2) && 
		 (MuonCand_validmuonhits[nMuonCand] >= 1))
		MuonCand_tightID[nMuonCand] = 1;
		 
	    }
	}
      
      nMuonCand++;
    }

  // Get the electron collection from the event

  // PAT
  edm::Handle<edm::View<pat::Electron> > electrons; 
  event.getByLabel(thePixelGsfELabel,electrons); 
  edm::View<pat::Electron>::const_iterator electron;

  for ( electron = electrons->begin(); electron != electrons->end() && nEleCand<ELEMAX; ++electron )
    {
	  EleCand_e[nEleCand]=electron->energy();
	  EleCand_et[nEleCand]=electron->et();
	  EleCand_px[nEleCand]=electron->px();
	  EleCand_py[nEleCand]=electron->py();
	  EleCand_pz[nEleCand]=electron->pz();
	  EleCand_p[nEleCand]=electron->p();
	  EleCand_phi[nEleCand]=electron->phi();
	  EleCand_eta[nEleCand]=electron->eta();
	  EleCand_charge[nEleCand]=electron->charge(); 
	  EleCand_vtxx[nEleCand]=electron->vertex().x(); 
	  EleCand_vtxy[nEleCand]=electron->vertex().y();  
	  EleCand_vtxz[nEleCand]=electron->vertex().z(); 
	  if(electron->closestCtfTrackRef().isNonnull())
	    {
	      EleCandTrack_p[nEleCand]=electron->closestCtfTrackRef()->p(); 
              EleCandTrack_pt[nEleCand]=electron->closestCtfTrackRef()->pt();  
              EleCandTrack_eta[nEleCand]=electron->closestCtfTrackRef()->eta();  
              EleCandTrack_phi[nEleCand]=electron->closestCtfTrackRef()->phi();  
	    }
	  else
	    {
	      EleCandTrack_p[nEleCand]=-999.;
              EleCandTrack_pt[nEleCand]=-999.; 
              EleCandTrack_eta[nEleCand]=-999.; 
              EleCandTrack_phi[nEleCand]=-999.; 
	    }
	  EleCand_deltaPhi[nEleCand] = electron->deltaPhiSuperClusterTrackAtVtx();
	  EleCand_deltaEta[nEleCand] = electron->deltaEtaSuperClusterTrackAtVtx();
	  EleCand_HoverE[nEleCand] = electron->hcalOverEcal();
	  EleCand_trackiso[nEleCand] = electron->dr03TkSumPt() / electron->et();
	  EleCand_ecaliso[nEleCand] = electron->dr03EcalRecHitSumEt() / electron->et();
	  EleCand_hcaliso[nEleCand] = electron->dr03HcalTowerSumEt() / electron->et();
	  EleCand_sigmaIetaIeta[nEleCand] = electron->sigmaIetaIeta();
	  EleCand_convDist[nEleCand] = fabs(electron->convDist()); 
	  EleCand_convDcot[nEleCand] = fabs(electron->convDcot()); 
	  EleCand_ecalDriven[nEleCand] = electron->ecalDrivenSeed();

	  // WP80 cuts
	  const float MAX_MissingHits      = 0.0;
	  const float MIN_Dist             = 0.02;
	  const float MIN_Dcot             = 0.02;
	  const float cut_EB_trackRel03    = 0.09;
	  const float cut_EB_ecalRel03     = 0.07;
	  const float cut_EB_hcalRel03     = 0.10;
	  const float cut_EB_sigmaIetaIeta = 0.01;
	  const float cut_EB_deltaPhi      = 0.06;
	  const float cut_EB_deltaEta      = 0.004;
	  const float cut_EB_HoverE        = 0.04;
	  const float cut_EE_trackRel03    = 0.04;
	  const float cut_EE_ecalRel03     = 0.05;
	  const float cut_EE_hcalRel03     = 0.025;
	  const float cut_EE_sigmaIetaIeta = 0.03;
	  const float cut_EE_deltaPhi      = 0.03;
	  const float cut_EE_deltaEta      = 0.007; 
	  const float cut_EE_HoverE        = 0.025;

	  EleCand_wp80[nEleCand] = 0;
	  bool select = true;
	  if( !(electron->ecalDrivenSeed()==1) )  select = false; 
	  if(EleCand_convDist[nEleCand]<MIN_Dist && EleCand_convDcot[nEleCand]<MIN_Dcot) select = false;
	  if(electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()>MAX_MissingHits) select = false;
	  if(select == true)
	  if( electron->isEB() ) {
	    if( fabs(EleCand_deltaPhi[nEleCand])>cut_EB_deltaPhi ) select = false;
	    if( fabs(EleCand_deltaEta[nEleCand])>cut_EB_deltaEta ) select = false;
	    if( EleCand_HoverE[nEleCand]>cut_EB_HoverE ) select = false;
	    if( EleCand_trackiso[nEleCand]>cut_EB_trackRel03 ) select = false;
	    if( EleCand_ecaliso[nEleCand]>cut_EB_ecalRel03 ) select = false;
	    if( EleCand_hcaliso[nEleCand]>cut_EB_hcalRel03 ) select = false;
	    if( EleCand_sigmaIetaIeta[nEleCand] > cut_EB_sigmaIetaIeta ) select = false;
	  }    
	  else if( electron->isEE() ) {
	    if( fabs(EleCand_deltaPhi[nEleCand])>cut_EE_deltaPhi ) select = false;
	    if( fabs(EleCand_deltaEta[nEleCand])>cut_EE_deltaEta ) select = false;
	    if( EleCand_HoverE[nEleCand]>cut_EE_HoverE ) select = false;
	    if( EleCand_trackiso[nEleCand]>cut_EE_trackRel03 ) select = false;
	    if( EleCand_ecaliso[nEleCand]>cut_EE_ecalRel03 ) select = false;
	    if( EleCand_hcaliso[nEleCand]>cut_EE_hcalRel03 ) select = false;
	    if( EleCand_sigmaIetaIeta[nEleCand] > cut_EE_sigmaIetaIeta ) select = false;
	  }
	  else select = false;

	  if(select == true)
	    EleCand_wp80[nEleCand] = 1; 

	  nEleCand++;
    }	

  // Calculate invariant mass, delta-phi and delta-pT
  bool found_pair(false);
  bool found_muevertex(false);
  MuEPairCand[0]=0; MuEPairCand[1]=0;
  
  if(nMuonCand == 1 && nEleCand == 1)
    {
      if((MuonCand_charge[0]*EleCand_charge[0]<0) || (keepsamesign == true)) found_pair=true;
    }

  if(nMuonCand+nEleCand>2)
    {
    double minimal_distance(999); 
    for(int k=0; k<nMuonCand; k++)
	{
	for(int l=k+1; l<nEleCand; l++)
		{
		if((MuonCand_charge[k]*EleCand_charge[l]<0) || (keepsamesign == true))
			{
			found_pair=true;
			double muonsDist=sqrt(pow(MuonCand_vtxx[k]-EleCand_vtxx[l],2)
        	                             +pow(MuonCand_vtxy[k]-EleCand_vtxy[l],2)
                	                     +pow(MuonCand_vtxz[k]-EleCand_vtxz[l],2));
			if(muonsDist<minimal_distance){minimal_distance=muonsDist; MuEPairCand[0]=k; MuEPairCand[1]=l;}
			}
		}
	}
    }

    if(found_pair){
      TLorentzVector recomuvec; 
      TLorentzVector recoevec; 
      TLorentzVector recomuevec; 
 
      recomuvec.SetXYZM(MuonCand_px[MuEPairCand[0]],MuonCand_py[MuEPairCand[0]],MuonCand_pz[MuEPairCand[0]],0.1057); 
      recoevec.SetXYZM(EleCand_px[MuEPairCand[1]],EleCand_py[MuEPairCand[1]],EleCand_pz[MuEPairCand[1]],0.000511);  
      recomuevec = recomuvec + recoevec; 
      MuE_mass = recomuevec.M();

      MuE_pt = recomuevec.Pt();

      MuE_dpt = fabs(MuonCand_pt[MuEPairCand[0]]-EleCand_et[MuEPairCand[1]]);

      double dphi = fabs(MuonCand_phi[MuEPairCand[0]]-EleCand_phi[MuEPairCand[1]]);
      if(dphi < 3.14159265359)
            MuE_dphi = dphi;
      else
        MuE_dphi = (2.0*3.14159265359)-dphi;
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
  //  edm::Handle<reco::CaloMETCollection> pMET;
  //  const reco::CaloMETCollection* mets = pMET.product();
  //  reco::CaloMETCollection::const_iterator met;
  edm::Handle<reco::PFMETCollection> pMET; 
  event.getByLabel(theMetLabel,pMET); 
  const reco::PFMETCollection* mets = pMET.product(); 
  reco::PFMETCollection::const_iterator met; 

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

    // Now check if a primary vertex is consistent with having exactly 2 muons 
    // and no other tracks

    int track_match_lepton=0;
    PrimVertexCand_mueTwoTracks[nPrimVertexCand] = 0;
    PrimVertexCand_mueExactlyTwoTracks[nPrimVertexCand] = 0; 
    PrimVertexCand_mueTwoTracksMuIndex[nPrimVertexCand] = -1;
    PrimVertexCand_mueTwoTracksEleIndex[nPrimVertexCand] = -1; 

    // Uncomment for keeping background events
    if((PrimVertexCand_tracks[nPrimVertexCand] >= 2) && found_pair)
    //    if((PrimVertexCand_tracks[nPrimVertexCand] == 2) && found_pair) 
      {
        for (reco::Vertex::trackRef_iterator vertex_curTrack = vertex_i->tracks_begin(); vertex_curTrack!=vertex_i->tracks_end(); vertex_curTrack++) {
	  //		if( (fabs((*vertex_curTrack)->pt()-MuonCand_pt[MuEPairCand[0]])<1.e-6   || fabs((*vertex_curTrack)->pt()-EleCandTrack_pt[MuEPairCand[1]])<1.e2) &&
	  if( (fabs((*vertex_curTrack)->p()-MuonCandTrack_p[MuEPairCand[0]])<1.e-2   || fabs((*vertex_curTrack)->pt()-EleCandTrack_pt[MuEPairCand[1]])<1.e2) &&
	      (fabs((*vertex_curTrack)->eta()-MuonCand_eta[MuEPairCand[0]])<1.e-2 || fabs((*vertex_curTrack)->eta()-EleCandTrack_eta[MuEPairCand[1]])<1.e-1) &&
	      (fabs((*vertex_curTrack)->phi()-MuonCand_phi[MuEPairCand[0]])<1.e-2 || fabs((*vertex_curTrack)->phi()-EleCandTrack_phi[MuEPairCand[1]])<1.e-1)
	      ) track_match_lepton++;
	}
    }
    
    // Uncomment for keeping background events
    if((PrimVertexCand_tracks[nPrimVertexCand] >= 2) && found_pair && track_match_lepton>=2)
    //    if((PrimVertexCand_tracks[nPrimVertexCand] == 2) && found_pair && track_match_lepton==2)
    {
	PrimVertexCand_mueTwoTracks[nPrimVertexCand] = 1;
        PrimVertexCand_mueTwoTracksMuIndex[nPrimVertexCand] = MuEPairCand[0]; 
        PrimVertexCand_mueTwoTracksEleIndex[nPrimVertexCand] = MuEPairCand[1]; 

	mueprimvtxx = PrimVertexCand_x[nPrimVertexCand];
	mueprimvtxy = PrimVertexCand_y[nPrimVertexCand]; 
        mueprimvtxz = PrimVertexCand_z[nPrimVertexCand]; 
	found_muevertex = true;

	if((PrimVertexCand_tracks[nPrimVertexCand] == 2) && track_match_lepton==2)
	  PrimVertexCand_mueExactlyTwoTracks[nPrimVertexCand] = 1; 

    }
    nPrimVertexCand++;
  }


  // Get the CASTOR towers collection from the event 
  edm::Handle<reco::CastorTowerCollection> recoCastorTowers;  
  event.getByLabel(recCastorTowerLabel, recoCastorTowers);  

  // Get the ZDC rechits collection from the event

  // Get the PFlow collection from the event
  edm::Handle<reco::PFCandidateCollection> pflows;
  event.getByLabel("particleFlow",pflows);
  reco::PFCandidateCollection::const_iterator pflow;

  double highestejet = -1.0;
  double highestejeteta = -999.0;
  double highestejetphi = -999.0;
  double totalejet = -1.0;
  double highestcastortowerfwd = -999.0; 
  double highestcastortowerbwd = -999.0; 
  double totalecastorfwd = 0.0;
  double totalecastorbwd = 0.0;
  double totalecalo = -1.0; 
  double closesttrkdxyz = 999.0;
  double closesthighpuritytrkdxyz = 999.0;

  // If this event contains a di-mu/e/gamma candidate, look at Jets & MET & CaloTowers & Tracks
  if(nMuonCand >= 1 && nEleCand >= 1)
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
      float e_met = met->et();
      Etmiss = e_met;
      Etmiss_phi = met->phi();
      Etmiss_x = met->px();
      Etmiss_y = met->py();
      Etmiss_z = met->pz();
      Etmiss_significance = met->significance();
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
	  CaloTowerDetId idTower=calo->id();
	  
	  size_t numRecHits = calo->constituentsSize();
	  bool isEB(false),isEE(false),isHB(false),isHE(false),isHF(false),isHO(false);
	  for (size_t j = 0; j < numRecHits; j++) {
	    DetId RecHitDetID=calo->constituent(j);
	    DetId::Detector DetNum=RecHitDetID.det();
	    if( DetNum == DetId::Hcal){
	      HcalDetId HcalID = RecHitDetID;
	      int HcalNum =  HcalID.subdetId();
	      if(HcalNum == HcalForward ){isHF=true;}
	      if(HcalNum == HcalBarrel ) {isHB=true;}
	      if(HcalNum == HcalEndcap ) {isHE=true;}
	      if(HcalNum == HcalOuter ) {isHO=true;}
	    }
	    if( DetNum == DetId::Ecal){
	      int EcalNum = RecHitDetID.subdetId();
	      if(EcalNum==1){isEB=true;}
	      if(EcalNum==2){isEE=true;}
	    }
	  }

	  CaloTower_badhcalcells[nCaloCand]=calo->numBadHcalCells();
	  CaloTower_problemhcalcells[nCaloCand]=calo->numProblematicHcalCells(); 
          CaloTower_badecalcells[nCaloCand]=calo->numBadEcalCells(); 
          CaloTower_problemecalcells[nCaloCand]=calo->numProblematicEcalCells();  
	  
	  float calodr1 = sqrt(((CaloTower_eta[nCaloCand]-MuonCand_eta[MuEPairCand[0]])*(CaloTower_eta[nCaloCand]-MuonCand_eta[MuEPairCand[0]])) + 
			       ((CaloTower_phi[nCaloCand]-MuonCand_phi[MuEPairCand[0]])*(CaloTower_phi[nCaloCand]-MuonCand_phi[MuEPairCand[0]])));
	  float calodr2 = sqrt(((CaloTower_eta[nCaloCand]-MuonCand_eta[MuEPairCand[1]])*(CaloTower_eta[nCaloCand]-MuonCand_eta[MuEPairCand[1]])) + 
			       ((CaloTower_phi[nCaloCand]-MuonCand_phi[MuEPairCand[1]])*(CaloTower_phi[nCaloCand]-MuonCand_phi[MuEPairCand[1]])));
	  
	  if(calodr1 < calodr2)
	    CaloTower_dr[nCaloCand] = calodr1;
	  else
	    CaloTower_dr[nCaloCand] = calodr2;

	  totalecalo = totalecalo + CaloTower_e[nCaloCand]; 

	  nCaloCand++;
	}

      SumCalo_e = totalecalo;

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
	      
	      if(CastorTower_e[nCastorTowerCand] > 3.0)
		nCastorTowerCandE3++;

	      nCastorTowerCand++;  
	    } 
	}

      HighestCastorTowerFwd_e = highestcastortowerfwd; 
      HighestCastorTowerBwd_e = highestcastortowerbwd; 
      SumCastorFwd_e = totalecastorfwd;
      SumCastorBwd_e = totalecastorbwd; 
    }

  // Now ZDC rechits

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

  GenMuE_eta = 0.0;
  GenMuE_pt = 0.0;
  if(nGenMuonCand == 2)
    {
      TLorentzVector muvec1;
      TLorentzVector muvec2;
      TLorentzVector muevec;

      muvec1.SetXYZM(GenMuonCand_px[0],GenMuonCand_py[0],GenMuonCand_pz[0],0.1057);
      muvec2.SetXYZM(GenMuonCand_px[1],GenMuonCand_py[1],GenMuonCand_pz[1],0.1057); 
      muevec = muvec1 + muvec2;
      GenMuE_eta = muevec.Eta();
      GenMuE_pt = muevec.Pt();
    }

  // Now do vertexing and track counting
  edm::ESHandle<TransientTrackBuilder> theVtx;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theVtx);
  //  vector < reco::TransientTrack > mutrks;
  vector<TransientTrack> transmutrks; 
  reco::TrackCollection * mutrks = new reco::TrackCollection;

  // First get "muon" tracks
  bool isLepton = false;
  for( track = tracks->begin(); track != tracks->end(); ++ track ) 
    { 
      isLepton = false;
      for(int j = 0;j < 2; j++)
	{
	  if(MuonCandTrack_p[MuEPairCand[j]] == track->p())
	    {
	      isLepton = true;
	      mutrks->push_back( *track );
	      TransientTrack tmptrk = (*theVtx).build( *track );
	      transmutrks.push_back( tmptrk );
	    }
	  if(EleCandTrack_p[MuEPairCand[j]] == track->p()) 
	    {
              isLepton = true; 
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
      TransientVertex mueVertex = fitter.vertex(transmutrks); 
      if(mueVertex.isValid())
	{
	  MuE_Kalmanvtxx = mueVertex.position().x(); 
	  MuE_Kalmanvtxy = mueVertex.position().y(); 
	  MuE_Kalmanvtxz = mueVertex.position().z(); 
	  MuE_Kalmanvtxchi2dof = mueVertex.normalisedChiSquared();
	  MuE_KalmanvtxT = sqrt(mueVertex.position().x()*mueVertex.position().x() + mueVertex.position().y()*mueVertex.position().y() ); 
	  MuE_Kalmanvtxisvalid = 1;
	  //	  found_muevertex = 1;
	}
      else
	{
	  MuE_Kalmanvtxx = 99;  
	  MuE_Kalmanvtxy = 99;  
	  MuE_Kalmanvtxz = 99;  
	  MuE_KalmanvtxT = 99;
	  MuE_Kalmanvtxchi2dof = 9999;
	  MuE_Kalmanvtxisvalid = 0;
	}
   }

  // OK, now go back and count "extra" tracks on the dimuon vertex
  // Loop2 = compute "track" quantities
  for(track = tracks->begin(); track != tracks->end() && nTrackCand<TRACKMAX; ++ track)
    {
      if((track->p() == MuonCandTrack_p[MuEPairCand[0]]) || (track->p() == MuonCandTrack_p[MuEPairCand[1]])) continue;
      if((track->p() == EleCandTrack_p[MuEPairCand[0]]) || (track->p() == EleCandTrack_p[MuEPairCand[1]])) continue; 

      TrackCand_purity[nTrackCand]=track->quality(TrackBase::highPurity); 
      TrackCand_p[nTrackCand]=track->p(); 
      TrackCand_px[nTrackCand]=track->px(); 
      TrackCand_py[nTrackCand]=track->py(); 
      TrackCand_pz[nTrackCand]=track->pz(); 
      TrackCand_pt[nTrackCand]=track->pt(); 
      TrackCand_eta[nTrackCand]=track->eta(); 
      TrackCand_phi[nTrackCand]=track->phi(); 
      TrackCand_charge[nTrackCand]=track->charge();
      TrackCand_nhits[nTrackCand]=track->numberOfValidHits();
      TrackCand_chi2[nTrackCand]=track->chi2();
      TrackCand_ndof[nTrackCand]=track->ndof();
      TrackCand_vtxdxyz[nTrackCand] = sqrt(((track->vertex().x() - mueprimvtxx)*(track->vertex().x() - mueprimvtxx)) + 
					   ((track->vertex().y() - mueprimvtxy)*(track->vertex().y() - mueprimvtxy)) +
					   ((track->vertex().z() - mueprimvtxz)*(track->vertex().z() - mueprimvtxz)));
      TrackCand_vtxT[nTrackCand] = sqrt(((track->vertex().x() - mueprimvtxx)*(track->vertex().x() - mueprimvtxx)) +
					((track->vertex().y() - mueprimvtxy)*(track->vertex().y() - mueprimvtxy)));
      TrackCand_vtxZ[nTrackCand] = sqrt(((track->vertex().z() - mueprimvtxz)*(track->vertex().z() - mueprimvtxz)));
      TrackCand_X[nTrackCand] = track->vertex().x();
      TrackCand_Y[nTrackCand] = track->vertex().y();
      TrackCand_Z[nTrackCand] = track->vertex().z();

      if((TrackCand_purity[nTrackCand] == 1) && (TrackCand_nhits[nTrackCand] >= 3))
	nQualityTrackCand++;
      
      if(TrackCand_vtxdxyz[nTrackCand] < 0.1)
	MuE_extratracks1mm++;
      if(TrackCand_vtxdxyz[nTrackCand] < 0.2) 
	MuE_extratracks2mm++; 
      if(TrackCand_vtxdxyz[nTrackCand] < 0.3)
	MuE_extratracks3mm++;
      if(TrackCand_vtxdxyz[nTrackCand] < 0.4) 
	MuE_extratracks4mm++; 
      if(TrackCand_vtxdxyz[nTrackCand] < 0.5) 
	MuE_extratracks5mm++; 
      if(TrackCand_vtxdxyz[nTrackCand] < 1) 
	MuE_extratracks1cm++; 
      if(TrackCand_vtxdxyz[nTrackCand] < 3)
	MuE_extratracks3cm++;
      if(TrackCand_vtxdxyz[nTrackCand] < 5) 
	MuE_extratracks5cm++; 
      if(TrackCand_vtxdxyz[nTrackCand] < 10) 
	MuE_extratracks10cm++; 
      if(TrackCand_vtxdxyz[nTrackCand] < closesttrkdxyz)
	closesttrkdxyz = TrackCand_vtxdxyz[nTrackCand];
      if((TrackCand_vtxdxyz[nTrackCand] < closesthighpuritytrkdxyz) && 
	 (TrackCand_purity[nTrackCand] == 1) && 
	 (TrackCand_nhits[nTrackCand] >= 3))
	closesthighpuritytrkdxyz = TrackCand_vtxdxyz[nTrackCand];
      
      nTrackCand++;  
    } 
  ClosestExtraTrack_vtxdxyz = closesttrkdxyz;
  ClosestHighPurityExtraTrack_vtxdxyz = closesthighpuritytrkdxyz;

  // Check for di-objects with valid vertex
  if(nMuonCand < 1 || nEleCand < 1 || !(found_pair)) {passed = false;}

  // Comment for keeping background events
  if(!(found_muevertex)) {passed = false;}
  //  if(ClosestHighPurityExtraTrack_vtxdxyz < minmuevtxd) {passed = false;}

  // "Exclusivity" cuts
  if(passed == true){
    thetree->Fill();
  }
}

void
GammaGammaMuE::fillDescriptions(ConfigurationDescriptions & descriptions) {
  
  descriptions.setComment("Exclusive MuE EDAnalyzer.");
  
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
  iDesc.add<edm::InputTag>("MetLabel", edm::InputTag("selectedLayer1METs"))->setComment("input MET collection");   
  iDesc.add<double>("CaloTowerdR", 0.3)->setComment("Minimum delta-R to use for finding extra towers");  
  iDesc.add<bool>("KeepSameSign", false)->setComment("Set to true to keep same-sign emu combinations");
  iDesc.addOptionalUntracked<std::string>("outfilename", ("mue.pat.root"))->setComment("output flat ntuple file name");  
  iDesc.add<std::string>("HLTMenuLabel", ("HLT8E29"))->setComment("HLT AOD trigger summary label");

  iDesc.add<double>("MinMuEVertexSeparation",0.1)->setComment("Minimum distance in cm between the dimuon vertex and any other track");

  descriptions.add("ParameterDescriptionsForGammaGammaMuE", iDesc);
}

// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaMuE::beginJob()
{
}

void
GammaGammaMuE::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,hltMenuLabel,changed)) {
    if (changed) {
      // check if trigger name in (new) config
      std::string   triggerName_ = "HLT_DoubleMu0";
      if (triggerName_!="@") { // "@" means: analyze all triggers in config
	const unsigned int n(hltConfig_.size());
	const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
	if (triggerIndex>=n) {
	  cout << "GammaGammaMuE::analyze:"
	       << " TriggerName " << triggerName_ 
	       << " not available in (new) config!" << endl;
	  cout << "Available TriggerNames are: " << endl;
//	  hltConfig_.dump("Triggers");
	}
      }
//      hltConfig_.dump("Streams");
//      hltConfig_.dump("Datasets");
//      hltConfig_.dump("PrescaleTable");
//      hltConfig_.dump("ProcessPSet");
    }
  } else {
    cout << "GammaGammaMuE::beginRun:"
	 << " config extraction failure with process name "
	 << hltMenuLabel << endl;
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaMuE::endJob() {
  const edm::ParameterSet &thepset = edm::getProcessParameterSet();
  TList *list = thetree->GetUserInfo();
  list->Add(new TObjString(thepset.dump().c_str()));
  thefile->Write();
  thefile->Close();
}
  
