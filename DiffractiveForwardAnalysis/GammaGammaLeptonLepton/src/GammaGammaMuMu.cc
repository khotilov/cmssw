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
// $Id: GammaGammaMuMu.cc,v 1.107 2011/11/28 16:05:10 jjhollar Exp $
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
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"
 
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
 #include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"

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
//#include "DataFormats/L1GlobalCaloTrigger/src/L1GctJetCounts.cc"
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
  minmumuvtxd        = pset.getParameter<double>("MinMuMuVertexSeparation"); 

  readmcPileup         = pset.getParameter<bool>("ReadMCPileup");
  readmcEffCorrections = pset.getParameter<bool>("ReadMCEffCorrections");
  readmcEffCorrectionsByCharge = pset.getParameter<bool>("ReadMCEffCorrectionsByCharge"); 
  readmcEffCorrectionsBySignedEta = pset.getParameter<bool>("ReadmcEffCorrectionsBySignedEta");
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
  thetree->Branch("MuonPairCand",MuonPairCand,"MuonPairCand[2]/I");

  thetree->Branch("nHLTDiMu7MuonCand",&nHLTDiMu7MuonCand,"nHLTDiMu7MuonCand/I"); 
  thetree->Branch("HLT_DoubleMu7_MuonCand_pt",&HLT_DoubleMu7_MuonCand_pt,"HLT_DoubleMu7_MuonCand_pt[nHLTDiMu7MuonCand]/D"); 
  thetree->Branch("HLT_DoubleMu7_MuonCand_eta",&HLT_DoubleMu7_MuonCand_eta,"HLT_DoubleMu7_MuonCand_eta[nHLTDiMu7MuonCand]/D"); 
  thetree->Branch("HLT_DoubleMu7_MuonCand_phi",&HLT_DoubleMu7_MuonCand_phi,"HLT_DoubleMu7_MuonCand_phi[nHLTDiMu7MuonCand]/D"); 
  thetree->Branch("HLT_DoubleMu7_MuonCand_charge",&HLT_DoubleMu7_MuonCand_charge,"HLT_DoubleMu7_MuonCand_charge[nHLTDiMu7MuonCand]/I");  

  thetree->Branch("nHLTMu13Mu8MuonCand",&nHLTMu13Mu8MuonCand,"nHLTMu13Mu8MuonCand/I");  
  thetree->Branch("HLT_Mu13Mu8_MuonCand_pt",&HLT_Mu13Mu8_MuonCand_pt,"HLT_Mu13Mu8_MuonCand_pt[nHLTMu13Mu8MuonCand]/D");  
  thetree->Branch("HLT_Mu13Mu8_MuonCand_eta",&HLT_Mu13Mu8_MuonCand_eta,"HLT_Mu13Mu8_MuonCand_eta[nHLTMu13Mu8MuonCand]/D");  
  thetree->Branch("HLT_Mu13Mu8_MuonCand_phi",&HLT_Mu13Mu8_MuonCand_phi,"HLT_Mu13Mu8_MuonCand_phi[nHLTMu13Mu8MuonCand]/D");  
  thetree->Branch("HLT_Mu13Mu8_MuonCand_charge",&HLT_Mu13Mu8_MuonCand_charge,"HLT_Mu13Mu8_MuonCand_charge[nHLTMu13Mu8MuonCand]/I");   

  thetree->Branch("nHLTDiMu4AcopMuonCand",&nHLTDiMu4AcopMuonCand,"nHLTDiMu4AcopMuonCand/I");  
  thetree->Branch("HLT_DoubleMu4Acoplanarity_MuonCand_pt",&HLT_DoubleMu4Acoplanarity_MuonCand_pt,"HLT_DoubleMu4Acoplanarity_MuonCand_pt[nHLTDiMu4AcopMuonCand]/D");  
  thetree->Branch("HLT_DoubleMu4Acoplanarity_MuonCand_eta",&HLT_DoubleMu4Acoplanarity_MuonCand_eta,"HLT_DoubleMu4Acoplanarity_MuonCand_eta[nHLTDiMu4AcopMuonCand]/D");  
  thetree->Branch("HLT_DoubleMu4Acoplanarity_MuonCand_phi",&HLT_DoubleMu4Acoplanarity_MuonCand_phi,"HLT_DoubleMu4Acoplanarity_MuonCand_phi[nHLTDiMu4AcopMuonCand]/D");  
  thetree->Branch("HLT_DoubleMu4Acoplanarity_MuonCand_charge",&HLT_DoubleMu4Acoplanarity_MuonCand_charge,"HLT_DoubleMu4Acoplanarity_MuonCand_charge[nHLTDiMu4AcopMuonCand]/I");   

  thetree->Branch("nHLTDiMu5AcopMuonCand",&nHLTDiMu5AcopMuonCand,"nHLTDiMu5AcopMuonCand/I");  
  thetree->Branch("HLT_DoubleMu5Acoplanarity_MuonCand_pt",&HLT_DoubleMu5Acoplanarity_MuonCand_pt,"HLT_DoubleMu5Acoplanarity_MuonCand_pt[nHLTDiMu5AcopMuonCand]/D");  
  thetree->Branch("HLT_DoubleMu5Acoplanarity_MuonCand_eta",&HLT_DoubleMu5Acoplanarity_MuonCand_eta,"HLT_DoubleMu5Acoplanarity_MuonCand_eta[nHLTDiMu5AcopMuonCand]/D");  
  thetree->Branch("HLT_DoubleMu5Acoplanarity_MuonCand_phi",&HLT_DoubleMu5Acoplanarity_MuonCand_phi,"HLT_DoubleMu5Acoplanarity_MuonCand_phi[nHLTDiMu5AcopMuonCand]/D");  
  thetree->Branch("HLT_DoubleMu5Acoplanarity_MuonCand_charge",&HLT_DoubleMu5Acoplanarity_MuonCand_charge,"HLT_DoubleMu5Acoplanarity_MuonCand_charge[nHLTDiMu5AcopMuonCand]/I");   

  thetree->Branch("nHLTDiMu6AcopMuonCand",&nHLTDiMu6AcopMuonCand,"nHLTDiMu6AcopMuonCand/I");
  thetree->Branch("HLT_DoubleMu6Acoplanarity_MuonCand_pt",&HLT_DoubleMu6Acoplanarity_MuonCand_pt,"HLT_DoubleMu6Acoplanarity_MuonCand_pt[nHLTDiMu6AcopMuonCand]/D");
  thetree->Branch("HLT_DoubleMu6Acoplanarity_MuonCand_eta",&HLT_DoubleMu6Acoplanarity_MuonCand_eta,"HLT_DoubleMu6Acoplanarity_MuonCand_eta[nHLTDiMu6AcopMuonCand]/D");
  thetree->Branch("HLT_DoubleMu6Acoplanarity_MuonCand_phi",&HLT_DoubleMu6Acoplanarity_MuonCand_phi,"HLT_DoubleMu6Acoplanarity_MuonCand_phi[nHLTDiMu6AcopMuonCand]/D");
  thetree->Branch("HLT_DoubleMu6Acoplanarity_MuonCand_charge",&HLT_DoubleMu6Acoplanarity_MuonCand_charge,"HLT_DoubleMu6Acoplanarity_MuonCand_charge[nHLTDiMu6AcopMuonCand]/I");

  thetree->Branch("nHLTDiMu7AcopMuonCand",&nHLTDiMu7AcopMuonCand,"nHLTDiMu7AcopMuonCand/I");
  thetree->Branch("HLT_DoubleMu7Acoplanarity_MuonCand_pt",&HLT_DoubleMu7Acoplanarity_MuonCand_pt,"HLT_DoubleMu7Acoplanarity_MuonCand_pt[nHLTDiMu7AcopMuonCand]/D");
  thetree->Branch("HLT_DoubleMu7Acoplanarity_MuonCand_eta",&HLT_DoubleMu7Acoplanarity_MuonCand_eta,"HLT_DoubleMu7Acoplanarity_MuonCand_eta[nHLTDiMu7AcopMuonCand]/D");
  thetree->Branch("HLT_DoubleMu7Acoplanarity_MuonCand_phi",&HLT_DoubleMu7Acoplanarity_MuonCand_phi,"HLT_DoubleMu7Acoplanarity_MuonCand_phi[nHLTDiMu7AcopMuonCand]/D");
  thetree->Branch("HLT_DoubleMu7Acoplanarity_MuonCand_charge",&HLT_DoubleMu7Acoplanarity_MuonCand_charge,"HLT_DoubleMu7Acoplanarity_MuonCand_charge[nHLTDiMu7AcopMuonCand]/I");


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

  thetree->Branch("nExtraCaloTowersE1",&nExtraCaloTowersE1,"nExtraCaloTowersE1/I"); 
  thetree->Branch("nExtraCaloTowersE2",&nExtraCaloTowersE2,"nExtraCaloTowersE2/I"); 
  thetree->Branch("nExtraCaloTowersE3",&nExtraCaloTowersE3,"nExtraCaloTowersE3/I");  
  thetree->Branch("nExtraCaloTowersE4",&nExtraCaloTowersE4,"nExtraCaloTowersE4/I");  
  thetree->Branch("nExtraCaloTowersE5",&nExtraCaloTowersE5,"nExtraCaloTowersE5/I");  
  thetree->Branch("nExtraCaloTowersE6",&nExtraCaloTowersE6,"nExtraCaloTowersE6/I");   
  thetree->Branch("nExtraCaloTowersE7",&nExtraCaloTowersE7,"nExtraCaloTowersE7/I");   
  thetree->Branch("nExtraCaloTowersE8",&nExtraCaloTowersE8,"nExtraCaloTowersE8/I");   
  thetree->Branch("nExtraCaloTowersE9",&nExtraCaloTowersE9,"nExtraCaloTowersE9/I");   

  thetree->Branch("nExtraCaloTowersEt0pt1",&nExtraCaloTowersEt0pt1,"nExtraCaloTowersEt0pt1/I"); 
  thetree->Branch("nExtraCaloTowersEt0pt2",&nExtraCaloTowersEt0pt2,"nExtraCaloTowersEt0pt2/I"); 
  thetree->Branch("nExtraCaloTowersEt0pt5",&nExtraCaloTowersEt0pt5,"nExtraCaloTowersEt0pt5/I"); 
  thetree->Branch("nExtraCaloTowersEt1",&nExtraCaloTowersEt1,"nExtraCaloTowersEt1/I"); 
  thetree->Branch("nExtraCaloTowersEt2",&nExtraCaloTowersEt2,"nExtraCaloTowersEt2/I"); 
  thetree->Branch("nExtraCaloTowersEt3",&nExtraCaloTowersEt3,"nExtraCaloTowersEt3/I");  
  thetree->Branch("nExtraCaloTowersEt4",&nExtraCaloTowersEt4,"nExtraCaloTowersEt4/I");  

  thetree->Branch("nExtraCaloTowersE0pt6eb",&nExtraCaloTowersE0pt6eb, "nExtraCaloTowersE0pt6eb/I");
  thetree->Branch("nExtraCaloTowersE2pt45ee", &nExtraCaloTowersE2pt45ee, "nExtraCaloTowersE2pt45ee/I");
  thetree->Branch("nExtraCaloTowersE1pt25hb", &nExtraCaloTowersE1pt25hb, "nExtraCaloTowersE1pt25hb/I");
  thetree->Branch("nExtraCaloTowersE1pt9he", &nExtraCaloTowersE1pt9he, "nExtraCaloTowersE1pt9he/I");
  thetree->Branch("nExtraCaloTowersE4pt5hfp", &nExtraCaloTowersE4pt5hfp, "nExtraCaloTowersE4pt5hfp/I");
  thetree->Branch("nExtraCaloTowersE4pt0hfm", &nExtraCaloTowersE4pt0hfm, "nExtraCaloTowersE4pt0hfm/I"); 

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

  thetree->Branch("nPFPhotonCand",&nPFPhotonCand,"nPFPhotonCand/I");
  thetree->Branch("PFPhotonCand_pt",PFPhotonCand_pt,"PFPhotonCand_pt[nPFPhotonCand]/D");
  thetree->Branch("PFPhotonCand_eta",PFPhotonCand_eta,"PFPhotonCand_eta[nPFPhotonCand]/D");
  thetree->Branch("PFPhotonCand_phi",PFPhotonCand_phi,"PFPhotonCand_phi[nPFPhotonCand]/D");
  thetree->Branch("PFPhotonCand_drtrue",PFPhotonCand_drtrue,"PFPhotonCand_drtrue[nPFPhotonCand]/D"); 

  thetree->Branch("MuMu_mass",&MuMu_mass,"MuMu_mass/D");
  thetree->Branch("MuMu_dphi",&MuMu_dphi,"MuMu_dphi/D");
  thetree->Branch("MuMu_dpt",&MuMu_dpt,"MuMu_dpt/D");
  thetree->Branch("MuMu_pt",&MuMu_pt,"MuMu_pt/D"); 
  thetree->Branch("MuMu_phi",&MuMu_phi,"MuMu_phi/D");   
  thetree->Branch("MuMu_3Dangle",&MuMu_3Dangle,"MuMu_3Dangle/D");   
  thetree->Branch("MuMu_Kalmanvtxx",&MuMu_Kalmanvtxx,"MuMu_Kalmanvtxx/D");
  thetree->Branch("MuMu_Kalmanvtxy",&MuMu_Kalmanvtxy,"MuMu_Kalmanvtxy/D"); 
  thetree->Branch("MuMu_Kalmanvtxz",&MuMu_Kalmanvtxz,"MuMu_Kalmanvtxz/D");
  thetree->Branch("MuMu_KalmanvtxT",&MuMu_KalmanvtxT,"MuMu_KalmanvtxT/D"); 
  thetree->Branch("MuMu_Kalmanvtxchi2dof",&MuMu_Kalmanvtxchi2dof,"MuMu_Kalmanvtxchi2dof/D");
  thetree->Branch("MuMu_Kalmanvtxisvalid",&MuMu_Kalmanvtxisvalid,"MuMu_Kalmanvtxisvalid/I");
  thetree->Branch("MuMu_extratracks1mm",&MuMu_extratracks1mm,"MuMu_extratracks1mm/I");
  thetree->Branch("MuMu_extratracks2mm",&MuMu_extratracks2mm,"MuMu_extratracks2mm/I"); 
  thetree->Branch("MuMu_extratracks3mm",&MuMu_extratracks3mm,"MuMu_extratracks3mm/I");
  thetree->Branch("MuMu_extratracks4mm",&MuMu_extratracks4mm,"MuMu_extratracks4mm/I"); 
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
  thetree->Branch("Etmiss_phi",&Etmiss_phi,"Etmiss_phi/D");  
  thetree->Branch("Etmiss_x",&Etmiss_x,"Etmiss_x/D");  
  thetree->Branch("Etmiss_y",&Etmiss_y,"Etmiss_y/D");  
  thetree->Branch("Etmiss_z",&Etmiss_z,"Etmiss_z/D");  
  thetree->Branch("Etmiss_significance",&Etmiss_significance,"Etmiss_significance/D");  

  thetree->Branch("HLT_DoubleMu5Acoplanarity",&HLT_DoubleMu5Acoplanarity,"HLT_DoubleMu5Acoplanarity/I");
  thetree->Branch("HLT_DoubleMu4Acoplanarity",&HLT_DoubleMu4Acoplanarity,"HLT_DoubleMu4Acoplanarity/I");
  thetree->Branch("HLT_DoubleMu6Acoplanarity",&HLT_DoubleMu6Acoplanarity,"HLT_DoubleMu6Acoplanarity/I");
  thetree->Branch("HLT_DoubleMu7Acoplanarity",&HLT_DoubleMu7Acoplanarity,"HLT_DoubleMu7Acoplanarity/I");
  thetree->Branch("HLT_DoubleMu7",&HLT_DoubleMu7,"HLT_DoubleMu7/I");
  thetree->Branch("HLT_Mu13Mu8", &HLT_Mu13Mu8, "HLT_Mu13Mu8/I"); 
  thetree->Branch("HLT_DoubleMu5Acoplanarity_Prescl",&HLT_DoubleMu5Acoplanarity_Prescl,"HLT_DoubleMu5Acoplanarity_Prescl/I"); 
  thetree->Branch("HLT_DoubleMu4Acoplanarity_Prescl",&HLT_DoubleMu4Acoplanarity_Prescl,"HLT_DoubleMu4Acoplanarity_Prescl/I"); 
  thetree->Branch("HLT_DoubleMu6Acoplanarity_Prescl",&HLT_DoubleMu6Acoplanarity_Prescl,"HLT_DoubleMu6Acoplanarity_Prescl/I");
  thetree->Branch("HLT_DoubleMu7Acoplanarity_Prescl",&HLT_DoubleMu7Acoplanarity_Prescl,"HLT_DoubleMu7Acoplanarity_Prescl/I");
  thetree->Branch("HLT_DoubleMu7_Prescl",&HLT_DoubleMu7_Prescl,"HLT_DoubleMu7_Prescl/I"); 
  thetree->Branch("HLT_Mu13Mu8_Prescl", &HLT_Mu13Mu8_Prescl, "HLT_Mu13Mu8_Prescl/I");  

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
  thetree->Branch("PrimVertexCand_mumuTwoTracks",&PrimVertexCand_mumuTwoTracks,"PrimVertexCand_mumuTwoTracks[nPrimVertexCand]/I"); 
  thetree->Branch("PrimVertexCand_mumuExactlyTwoTracks",&PrimVertexCand_mumuExactlyTwoTracks,"PrimVertexCand_mumuExactlyTwoTracks[nPrimVertexCand]/I");  
  thetree->Branch("PrimVertexCand_mumuTwoTracksMap",&PrimVertexCand_mumuTwoTracksMap,"PrimVertexCand_mumuTwoTracksMap/I");  

  thetree->Branch("LowPt_pt",LowPt_pt,"LowPt_pt[nMuonCand]/D");
  thetree->Branch("LowPt_eta",LowPt_eta,"LowPt_eta[nMuonCand]/D");

  thetree->Branch("nTruePUforPUWeight",&nTruePUforPUWeight,"nTruePUforPUWeight/I");
  thetree->Branch("nTruePUforPUWeightBXM1", &nTruePUforPUWeightBXM1, "nTruePUforPUWeightBXM1/I");
  thetree->Branch("nTruePUforPUWeightBXP1", &nTruePUforPUWeightBXP1, "nTruePUforPUWeightBXP1/I");
  thetree->Branch("nTruePUforPUWeightBX0", &nTruePUforPUWeightBX0, "nTruePUforPUWeightBX0/I");
  thetree->Branch("PUWeightTrue",&PUWeightTrue,"PUWeightTrue/D");

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
  PUWeightTrue = 0.0;
  nTruePUforPUWeight = 0;
  nTruePUforPUWeightBXM1 = 0;
  nTruePUforPUWeightBXP1 = 0;
  nTruePUforPUWeightBX0 = 0;

  nMuonCand=0;
  nHLTDiMu7MuonCand=0;
  nHLTMu13Mu8MuonCand=0; 
  nHLTDiMu7AcopMuonCand=0;
  nHLTDiMu6AcopMuonCand=0;
  nHLTDiMu5AcopMuonCand=0;  
  nHLTDiMu4AcopMuonCand=0;
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

  nExtraCaloTowersE0pt6eb=0;
  nExtraCaloTowersE2pt45ee=0;
  nExtraCaloTowersE1pt25hb=0; 
  nExtraCaloTowersE1pt9he=0; 
  nExtraCaloTowersE4pt5hfp=0;
  nExtraCaloTowersE4pt0hfm=0; 

  HitInZDC=0;
  HitInCastor=0;
  nGenPhotCand=0;
  nGenMuonCand=0;

  nPFPhotonCand=0;

  MuMu_mass = -1;
  MuMu_dphi = -1;
  MuMu_dpt = -1;
  MuMu_pt = -1;
  MuMu_phi = -999.;
  MuMu_3Dangle = -999.; 
  MuMu_extratracks1mm = 0; 
  MuMu_extratracks2mm = 0;  
  MuMu_extratracks3mm = 0; 
  MuMu_extratracks4mm = 0;  
  MuMu_extratracks5mm = 0; 
  MuMu_extratracks1cm = 0; 
  MuMu_extratracks3cm = 0;
  MuMu_extratracks5cm = 0;
  MuMu_extratracks10cm = 0;
  ClosestExtraTrack_vtxdxyz = 999.;
  double mumuprimvtxx = 0.0; 
  double mumuprimvtxy = 0.0; 
  double mumuprimvtxz = 0.0; 

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

  /*
  iLumi.getByLabel("lumiProducer", lumiSummary);

  if(lumiSummary->isValid())
    AvgInstDelLumi = lumiSummary->avgInsDelLumi();
  else
    AvgInstDelLumi = -999.;
  */

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
      if ( trigNames.triggerNames().at(i).find("HLT_DoubleMu4_Acoplanarity03_v") != string::npos)
	{  
	  HLT_DoubleMu4Acoplanarity_Prescl = hltConfig_.prescaleValue(event, iSetup, trigNames.triggerNames().at(i)); 
	  
          if ( hltResults->accept(i) )  
            {HLT_DoubleMu4Acoplanarity = 1;}
	  else
	    HLT_DoubleMu4Acoplanarity = 0;
        }  
      if ( trigNames.triggerNames().at(i).find("HLT_DoubleMu5_Acoplanarity03_v") != string::npos)
	{   
          HLT_DoubleMu5Acoplanarity_Prescl = hltConfig_.prescaleValue(event, iSetup, trigNames.triggerNames().at(i));

          if ( hltResults->accept(i) )  
	    {HLT_DoubleMu5Acoplanarity = 1;}
	  else
	    HLT_DoubleMu5Acoplanarity = 0;
        }  
      if ( trigNames.triggerNames().at(i).find("HLT_DoubleMu6_Acoplanarity03_v") != string::npos)
        {
          HLT_DoubleMu6Acoplanarity_Prescl = hltConfig_.prescaleValue(event, iSetup, trigNames.triggerNames().at(i));

          if ( hltResults->accept(i) )
            {HLT_DoubleMu6Acoplanarity = 1;}
          else
            HLT_DoubleMu6Acoplanarity = 0;
        }
      if ( trigNames.triggerNames().at(i).find("HLT_DoubleMu7_Acoplanarity03_v") != string::npos)
        {
          HLT_DoubleMu7Acoplanarity_Prescl = hltConfig_.prescaleValue(event, iSetup, trigNames.triggerNames().at(i));

          if ( hltResults->accept(i) )
            {HLT_DoubleMu7Acoplanarity = 1;}
          else
            HLT_DoubleMu7Acoplanarity = 0;
        }
      if ( trigNames.triggerNames().at(i).find("HLT_Mu13_Mu8") != string::npos ) 
        {   
          HLT_Mu13Mu8_Prescl = hltConfig_.prescaleValue(event, iSetup, trigNames.triggerNames().at(i)); 

          if ( hltResults->accept(i) )  
	    {HLT_Mu13Mu8 = 1;}
	  else
	    HLT_Mu13Mu8 = 0;
        }  
      if ( trigNames.triggerNames().at(i).find("HLT_DoubleMu7") != string::npos )
        {   
          HLT_DoubleMu7_Prescl = hltConfig_.prescaleValue(event, iSetup, trigNames.triggerNames().at(i)); 
 
          if ( hltResults->accept(i) )   
            {HLT_DoubleMu7 = 1;}
          else 
            HLT_DoubleMu7 = 0; 
        }   
    }

  Handle<TriggerEvent> hltObjects;
  event.getByLabel(InputTag("hltTriggerSummaryAOD","",hltMenuLabel),hltObjects);
  if (hltObjects.isValid()) 
    {
      size_type dimu4acopindex = hltObjects->filterIndex(InputTag("hltDoubleMu4ExclL3PreFiltered::"+hltMenuLabel)); 
      size_type dimu5acopindex = hltObjects->filterIndex(InputTag("hltDoubleMu5ExclL3PreFiltered::"+hltMenuLabel)); 
      size_type dimu6acopindex = hltObjects->filterIndex(InputTag("hltDoubleMu6ExclL3PreFiltered::"+hltMenuLabel));
      size_type dimu7acopindex = hltObjects->filterIndex(InputTag("hltDoubleMu7ExclL3PreFiltered::"+hltMenuLabel));
      size_type dimu7index = hltObjects->filterIndex(InputTag("hltDiMuonL3PreFiltered7::"+hltMenuLabel));
      size_type dimu8index = hltObjects->filterIndex(InputTag("hltDiMuonL3PreFiltered8::"+hltMenuLabel)); 

      if( dimu4acopindex < hltObjects->sizeFilters() )
	{
	  const trigger::Keys& DIMU4ACOPKEYS(hltObjects->filterKeys(dimu4acopindex)); 
	  const size_type nK(DIMU4ACOPKEYS.size()); 
	  const TriggerObjectCollection& TOC(hltObjects->getObjects());

	  for(int ipart = 0; ipart != nK; ++ipart)   
	    {
	      const TriggerObject& TO = TOC[DIMU4ACOPKEYS[ipart]];  
	      
	      if(fabs(TO.id()) == 13)
		{
		  HLT_DoubleMu4Acoplanarity_MuonCand_pt[nHLTDiMu4AcopMuonCand] = TO.pt();
		  HLT_DoubleMu4Acoplanarity_MuonCand_eta[nHLTDiMu4AcopMuonCand] = TO.eta(); 
		  HLT_DoubleMu4Acoplanarity_MuonCand_phi[nHLTDiMu4AcopMuonCand] = TO.phi(); 
		  if(TO.id() > 0)	  
		    HLT_DoubleMu4Acoplanarity_MuonCand_charge[nHLTDiMu4AcopMuonCand] = 1; 
		  else
		    HLT_DoubleMu4Acoplanarity_MuonCand_charge[nHLTDiMu4AcopMuonCand] = -1;  
		  
		  nHLTDiMu4AcopMuonCand++;
		}
	    }
	}
      if( dimu5acopindex < hltObjects->sizeFilters() ) 
        { 
          const trigger::Keys& DIMU5ACOPKEYS(hltObjects->filterKeys(dimu5acopindex));  
          const size_type nK(DIMU5ACOPKEYS.size());  
          const TriggerObjectCollection& TOC(hltObjects->getObjects()); 
 
          for(int ipart = 0; ipart != nK; ++ipart)    
            { 
              const TriggerObject& TO = TOC[DIMU5ACOPKEYS[ipart]];   
               
              if(fabs(TO.id()) == 13) 
                { 
                  HLT_DoubleMu5Acoplanarity_MuonCand_pt[nHLTDiMu5AcopMuonCand] = TO.pt(); 
                  HLT_DoubleMu5Acoplanarity_MuonCand_eta[nHLTDiMu5AcopMuonCand] = TO.eta();  
                  HLT_DoubleMu5Acoplanarity_MuonCand_phi[nHLTDiMu5AcopMuonCand] = TO.phi();  
                  if(TO.id() > 0)          
                    HLT_DoubleMu5Acoplanarity_MuonCand_charge[nHLTDiMu5AcopMuonCand] = 1;  
                  else 
                    HLT_DoubleMu5Acoplanarity_MuonCand_charge[nHLTDiMu5AcopMuonCand] = -1;   
                   
                  nHLTDiMu5AcopMuonCand++; 
                } 
            } 
        } 
      if( dimu6acopindex < hltObjects->sizeFilters() )
        {
          const trigger::Keys& DIMU6ACOPKEYS(hltObjects->filterKeys(dimu6acopindex));
          const size_type nK(DIMU6ACOPKEYS.size());
          const TriggerObjectCollection& TOC(hltObjects->getObjects());

          for(int ipart = 0; ipart != nK; ++ipart)
            {
              const TriggerObject& TO = TOC[DIMU6ACOPKEYS[ipart]];

              if(fabs(TO.id()) == 13)
                {
                  HLT_DoubleMu6Acoplanarity_MuonCand_pt[nHLTDiMu6AcopMuonCand] = TO.pt();
                  HLT_DoubleMu6Acoplanarity_MuonCand_eta[nHLTDiMu6AcopMuonCand] = TO.eta();
                  HLT_DoubleMu6Acoplanarity_MuonCand_phi[nHLTDiMu6AcopMuonCand] = TO.phi();
                  if(TO.id() > 0)
                    HLT_DoubleMu6Acoplanarity_MuonCand_charge[nHLTDiMu6AcopMuonCand] = 1;
                  else
                    HLT_DoubleMu6Acoplanarity_MuonCand_charge[nHLTDiMu6AcopMuonCand] = -1;

                  nHLTDiMu6AcopMuonCand++;
                }
            }
        }
      if( dimu7acopindex < hltObjects->sizeFilters() )
        {
          const trigger::Keys& DIMU7ACOPKEYS(hltObjects->filterKeys(dimu7acopindex));
          const size_type nK(DIMU7ACOPKEYS.size());
          const TriggerObjectCollection& TOC(hltObjects->getObjects());

          for(int ipart = 0; ipart != nK; ++ipart)
            {
              const TriggerObject& TO = TOC[DIMU7ACOPKEYS[ipart]];

              if(fabs(TO.id()) == 13)
                {
                  HLT_DoubleMu7Acoplanarity_MuonCand_pt[nHLTDiMu7AcopMuonCand] = TO.pt();
                  HLT_DoubleMu7Acoplanarity_MuonCand_eta[nHLTDiMu7AcopMuonCand] = TO.eta();
                  HLT_DoubleMu7Acoplanarity_MuonCand_phi[nHLTDiMu7AcopMuonCand] = TO.phi();
                  if(TO.id() > 0)
                    HLT_DoubleMu7Acoplanarity_MuonCand_charge[nHLTDiMu7AcopMuonCand] = 1;
                  else
                    HLT_DoubleMu7Acoplanarity_MuonCand_charge[nHLTDiMu7AcopMuonCand] = -1;

                  nHLTDiMu7AcopMuonCand++;
                }
            }
        }

      if( dimu7index < hltObjects->sizeFilters() )  
        {  
          const trigger::Keys& DIMU7KEYS(hltObjects->filterKeys(dimu7index));   
          const size_type nK(DIMU7KEYS.size());   
          const TriggerObjectCollection& TOC(hltObjects->getObjects());  
  
          for(int ipart = 0; ipart != nK; ++ipart)     
            {  
              const TriggerObject& TO = TOC[DIMU7KEYS[ipart]];    
                
              if(fabs(TO.id()) == 13)  
                {  
                  HLT_DoubleMu7_MuonCand_pt[nHLTDiMu7MuonCand] = TO.pt();  
                  HLT_DoubleMu7_MuonCand_eta[nHLTDiMu7MuonCand] = TO.eta();   
                  HLT_DoubleMu7_MuonCand_phi[nHLTDiMu7MuonCand] = TO.phi();   
                  if(TO.id() > 0)           
                    HLT_DoubleMu7_MuonCand_charge[nHLTDiMu7MuonCand] = 1;   
                  else  
                    HLT_DoubleMu7_MuonCand_charge[nHLTDiMu7MuonCand] = -1;    
                    
                  nHLTDiMu7MuonCand++;  
                }  
            }  
        }  
      if( dimu8index < hltObjects->sizeFilters() )   
        {   
          const trigger::Keys& DIMU8KEYS(hltObjects->filterKeys(dimu8index));    
          const size_type nK(DIMU8KEYS.size());    
          const TriggerObjectCollection& TOC(hltObjects->getObjects());   
   
          for(int ipart = 0; ipart != nK; ++ipart)      
            {   
              const TriggerObject& TO = TOC[DIMU8KEYS[ipart]];     
                 
              if(fabs(TO.id()) == 13)   
                {   
                  HLT_Mu13Mu8_MuonCand_pt[nHLTMu13Mu8MuonCand] = TO.pt();   
                  HLT_Mu13Mu8_MuonCand_eta[nHLTMu13Mu8MuonCand] = TO.eta();    
                  HLT_Mu13Mu8_MuonCand_phi[nHLTMu13Mu8MuonCand] = TO.phi();    
                  if(TO.id() > 0)            
                    HLT_Mu13Mu8_MuonCand_charge[nHLTMu13Mu8MuonCand] = 1;    
                  else   
                    HLT_Mu13Mu8_MuonCand_charge[nHLTMu13Mu8MuonCand] = -1;     
                     
                  nHLTMu13Mu8MuonCand++;   
                }   
            }   
        }   
    }


  // Get the #PU information
  if(readmcPileup)
    { 
      edm::Lumi3DReWeighting *LumiWeights;
      LumiWeights = new edm::Lumi3DReWeighting("PUMC_dist.root", "PUData_dist.root", "pileup", "pileup");
      LumiWeights->weight3D_init( 1.0 );
      const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&event);
      PUWeightTrue = LumiWeights->weight3D( (*iEventB) );
      
      Handle<std::vector< PileupSummaryInfo > >  PupInfo;
      event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
      
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      
      float sum_nvtx = 0.0;
      int npv = -1;
      int npvtrue = -1;
      int npvm1true = -1;
      int npvp1true = -1;
      int npv0true = -1;
      
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
	
	int BX = PVI->getBunchCrossing();
	
	//    cout << "PU rewighting - BX = " << BX << endl;
	if(BX == -1)
	  npvm1true++;
	if(BX == 0)
	  npv0true++;
	if(BX == 1)
	  npvp1true++;
	
	npv = PVI->getPU_NumInteractions();
	npvtrue = PVI->getTrueNumInteractions();
	
	sum_nvtx += float(npvtrue);
	
      }

      nTruePUforPUWeight = sum_nvtx;
      nTruePUforPUWeightBXM1 = npvm1true;
      nTruePUforPUWeightBXP1 = npvp1true;
      nTruePUforPUWeightBX0 = npv0true;
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

  //  edm::Handle<edm::ValueMap<reco::MuonCosmicCompatibility> > muonCosmicCompatibilityValueMapH;
  //  event.getByLabel(InputTag(“cosmicsVeto”), muonCosmicCompatibilityValueMapH);
  //  unsigned int muonIdx = 0;
  
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
      
      std::string algoname;
      double totalmuoneff = 1.0;
      double totalmuonmceff = 1.0;

      MuonCand_efficiency[nMuonCand] = totalmuoneff/totalmuonmceff;

      // Cosmics compatibility
      //     MuonRef muonRef(muonCollectionH, muonIdx);
      //      MuonCosmicCompatibility muonCosmicCompatibility = (*muonCosmicCompatibilityValueMapH)[muonRef];
      //      float combinedCompat = muonCosmicCompatibility.cosmicCompatibility;
      //      float B2BCosmicCompat = muonCosmicCompatibility.backToBackCompatibility;
      
      nMuonCand++;
    }

  // Calculate invariant mass, delta-phi and delta-pT
  bool found_pair(false);
  bool found_mumuvertex(false);
  MuonPairCand[0]=0; MuonPairCand[1]=1;
  
  if(nMuonCand == 2)
    {
      if((MuonCand_charge[0]*MuonCand_charge[1]<0) || (keepsamesign == true)) found_pair=true;
    }

  if(nMuonCand>2)
    {
    double minimal_distance(999); 
    for(int k=0; k<nMuonCand; k++)
	{
	for(int l=k+1; l<nMuonCand; l++)
		{
		if((MuonCand_charge[k]*MuonCand_charge[l]<0) || (keepsamesign == true))
			{
			found_pair=true;
			double muonsDist=sqrt(pow(MuonCand_vtxx[k]-MuonCand_vtxx[l],2)
        	                             +pow(MuonCand_vtxy[k]-MuonCand_vtxy[l],2)
                	                     +pow(MuonCand_vtxz[k]-MuonCand_vtxz[l],2));
			if(muonsDist<minimal_distance){minimal_distance=muonsDist; MuonPairCand[0]=k; MuonPairCand[1]=l;}
			}
		}
	}
    }
    if(found_pair){
      TLorentzVector recomuvec1; 
      TLorentzVector recomuvec2; 
      TLorentzVector recomumuvec; 
 
      recomuvec1.SetXYZM(MuonCand_px[MuonPairCand[0]],MuonCand_py[MuonPairCand[0]],MuonCand_pz[MuonPairCand[0]],0.1057); 
      recomuvec2.SetXYZM(MuonCand_px[MuonPairCand[1]],MuonCand_py[MuonPairCand[1]],MuonCand_pz[MuonPairCand[1]],0.1057);  
      recomumuvec = recomuvec1 + recomuvec2; 
      MuMu_mass = recomumuvec.M();
      
      MuMu_pt = recomumuvec.Pt();
      MuMu_phi = recomumuvec.Phi(); 
      MuMu_3Dangle = (recomuvec1.Angle(recomuvec2.Vect()))/3.14159265359; 

      MuMu_dpt = fabs(MuonCand_pt[MuonPairCand[0]]-MuonCand_pt[MuonPairCand[1]]);

      double dphi = fabs(MuonCand_phi[MuonPairCand[0]]-MuonCand_phi[MuonPairCand[1]]);
      if(dphi < 3.14159265359)
            MuMu_dphi = dphi;
      else
        MuMu_dphi = (2.0*3.14159265359)-dphi;
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

    int track_match_muon=0;
    PrimVertexCand_mumuTwoTracks[nPrimVertexCand] = 0;
    PrimVertexCand_mumuExactlyTwoTracks[nPrimVertexCand] = 0; 
    PrimVertexCand_mumuTwoTracksMap = -1;

    // Uncomment for keeping background events
    if((PrimVertexCand_tracks[nPrimVertexCand] >= 2) && found_pair)
    //    if((PrimVertexCand_tracks[nPrimVertexCand] == 2) && found_pair) 
      {
        for (reco::Vertex::trackRef_iterator vertex_curTrack = vertex_i->tracks_begin(); vertex_curTrack!=vertex_i->tracks_end(); vertex_curTrack++) {
		if( (fabs((*vertex_curTrack)->p()-MuonCandTrack_p[MuonPairCand[0]])<1.e-2   || fabs((*vertex_curTrack)->pt()-MuonCand_pt[MuonPairCand[1]])<1.e-2) &&
		    (fabs((*vertex_curTrack)->eta()-MuonCand_eta[MuonPairCand[0]])<1.e-2 || fabs((*vertex_curTrack)->eta()-MuonCand_eta[MuonPairCand[1]])<1.e-2) &&
                    (fabs((*vertex_curTrack)->phi()-MuonCand_phi[MuonPairCand[0]])<1.e-2 || fabs((*vertex_curTrack)->phi()-MuonCand_phi[MuonPairCand[1]])<1.e-2)
		) track_match_muon++;
	}
    }
    
    // Uncomment for keeping background events
    if((PrimVertexCand_tracks[nPrimVertexCand] >= 2) && found_pair && track_match_muon>=2)
    //    if((PrimVertexCand_tracks[nPrimVertexCand] == 2) && found_pair && track_match_muon==2)
    {
	PrimVertexCand_mumuTwoTracks[nPrimVertexCand] = 1;
	mumuprimvtxx = PrimVertexCand_x[nPrimVertexCand];
	mumuprimvtxy = PrimVertexCand_y[nPrimVertexCand]; 
        mumuprimvtxz = PrimVertexCand_z[nPrimVertexCand]; 
	if(PrimVertexCand_tracks[nPrimVertexCand] == 2)
	  {
	    PrimVertexCand_mumuExactlyTwoTracks[nPrimVertexCand] = 1;
	    PrimVertexCand_mumuTwoTracksMap = nPrimVertexCand;
	  }
	found_mumuvertex = true;

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
  double closesthighpuritytrkdxyz = 999.0;

  // If this event contains a di-mu/e/gamma candidate, look at Jets & MET & CaloTowers & Tracks
  if(nMuonCand >= 2)
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
	  if(isHF&&!isHB&&!isHE&&!isEB&&!isEE){CaloTower_ID[nCaloCand]=1;} //HF
	  else if(!isHF&&!isHB&&!isHE&&isEB&&!isEE){CaloTower_ID[nCaloCand]=2;} //EB
	  else if(!isHF&&!isHB&&!isHE&&!isEB&&isEE){CaloTower_ID[nCaloCand]=3;} //EE
	  else if(!isHF&&isHB&&!isHE&&!isEB&&!isEE){CaloTower_ID[nCaloCand]=4;} //HB
	  else if(!isHF&&!isHB&&isHE&&!isEB&&!isEE){CaloTower_ID[nCaloCand]=5;} //HE
	  else if(!isHF&&isHB&&!isHE&&isEB&&!isEE){CaloTower_ID[nCaloCand]=6;} // HB+EB
	  else if(!isHF&&!isHB&&isHE&&isEB&&!isEE){CaloTower_ID[nCaloCand]=7;} // HE+EB
	  else if(!isHF&&!isHB&&isHE&&!isEB&&isEE){CaloTower_ID[nCaloCand]=8;} // HE+EE
	  else if(!isHF&&isHB&&isHE&&!isEB&&!isEE){CaloTower_ID[nCaloCand]=9;} // HB+HE

	  CaloTower_badhcalcells[nCaloCand]=calo->numBadHcalCells();
	  CaloTower_problemhcalcells[nCaloCand]=calo->numProblematicHcalCells(); 
          CaloTower_badecalcells[nCaloCand]=calo->numBadEcalCells(); 
          CaloTower_problemecalcells[nCaloCand]=calo->numProblematicEcalCells();  
	  
	  float calodr1 = sqrt(((CaloTower_eta[nCaloCand]-MuonCand_eta[MuonPairCand[0]])*(CaloTower_eta[nCaloCand]-MuonCand_eta[MuonPairCand[0]])) + 
			       ((CaloTower_phi[nCaloCand]-MuonCand_phi[MuonPairCand[0]])*(CaloTower_phi[nCaloCand]-MuonCand_phi[MuonPairCand[0]])));
	  float calodr2 = sqrt(((CaloTower_eta[nCaloCand]-MuonCand_eta[MuonPairCand[1]])*(CaloTower_eta[nCaloCand]-MuonCand_eta[MuonPairCand[1]])) + 
			       ((CaloTower_phi[nCaloCand]-MuonCand_phi[MuonPairCand[1]])*(CaloTower_phi[nCaloCand]-MuonCand_phi[MuonPairCand[1]])));
	  
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

	      // Thresholds tuned on ZeroBias data
	      if(CaloTower_emE[nCaloCand] > 0.6 && CaloTower_ID[nCaloCand] == 1) 
		nExtraCaloTowersE0pt6eb=0; 
	      if(CaloTower_emE[nCaloCand] > 2.45 && CaloTower_ID[nCaloCand] == 2)
		nExtraCaloTowersE2pt45ee=0; 
	      if(CaloTower_hadE[nCaloCand] > 1.25 && CaloTower_ID[nCaloCand] == 4)
		nExtraCaloTowersE1pt25hb=0;  
	      if(CaloTower_hadE[nCaloCand] > 1.9 && CaloTower_ID[nCaloCand] == 5) 
		nExtraCaloTowersE1pt9he=0;  
	      if(CaloTower_e[nCaloCand] > 4.5 && CaloTower_eta[nCaloCand] > 2.95)
		  nExtraCaloTowersE4pt5hfp=0; 
	      if(CaloTower_e[nCaloCand] > 4.0 && CaloTower_eta[nCaloCand] < -2.95)      
		  nExtraCaloTowersE4pt0hfm=0;  
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

	  if(nMuonCand >= 2)  
	    {  
	      double mmgmass = pow(MuonCand_p[MuonPairCand[0]]+MuonCand_p[MuonPairCand[1]]+pflow->p(),2);   
	      mmgmass-=pow(MuonCand_px[MuonPairCand[0]]+MuonCand_px[MuonPairCand[1]]+pflow->px(),2);   
	      mmgmass-=pow(MuonCand_py[MuonPairCand[0]]+MuonCand_py[MuonPairCand[1]]+pflow->py(),2);   
	      mmgmass-=pow(MuonCand_pz[MuonPairCand[0]]+MuonCand_pz[MuonPairCand[1]]+pflow->pz(),2);   
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
      for(int j = 0;j < 2; j++)
	{
	  if(MuonCandTrack_p[MuonPairCand[j]] == track->p())
	    {
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
	  MuMu_Kalmanvtxx = mumuVertex.position().x(); 
	  MuMu_Kalmanvtxy = mumuVertex.position().y(); 
	  MuMu_Kalmanvtxz = mumuVertex.position().z(); 
	  MuMu_Kalmanvtxchi2dof = mumuVertex.normalisedChiSquared();
	  MuMu_KalmanvtxT = sqrt(mumuVertex.position().x()*mumuVertex.position().x() + mumuVertex.position().y()*mumuVertex.position().y() ); 
	  MuMu_Kalmanvtxisvalid = 1;
	  //	  found_mumuvertex = 1;
	}
      else
	{
	  MuMu_Kalmanvtxx = 99;  
	  MuMu_Kalmanvtxy = 99;  
	  MuMu_Kalmanvtxz = 99;  
	  MuMu_KalmanvtxT = 99;
	  MuMu_Kalmanvtxchi2dof = 9999;
	  MuMu_Kalmanvtxisvalid = 0;
	}
   }
  
  // OK, now go back and count "extra" tracks on the dimuon vertex
  // Loop2 = compute "track" quantities
  for(track = tracks->begin(); track != tracks->end() && nTrackCand<TRACKMAX; ++ track)
    {
      if(track->p() == MuonCandTrack_p[MuonPairCand[0]] || track->p() == MuonCandTrack_p[MuonPairCand[1]]) continue;
      
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
      TrackCand_vtxdxyz[nTrackCand] = sqrt(((track->vertex().x() - mumuprimvtxx)*(track->vertex().x() - mumuprimvtxx)) + 
					   ((track->vertex().y() - mumuprimvtxy)*(track->vertex().y() - mumuprimvtxy)) +
					   ((track->vertex().z() - mumuprimvtxz)*(track->vertex().z() - mumuprimvtxz)));
      TrackCand_vtxT[nTrackCand] = sqrt(((track->vertex().x() - mumuprimvtxx)*(track->vertex().x() - mumuprimvtxx)) +
					((track->vertex().y() - mumuprimvtxy)*(track->vertex().y() - mumuprimvtxy)));
      TrackCand_vtxZ[nTrackCand] = sqrt(((track->vertex().z() - mumuprimvtxz)*(track->vertex().z() - mumuprimvtxz)));
      TrackCand_X[nTrackCand] = track->vertex().x();
      TrackCand_Y[nTrackCand] = track->vertex().y();
      TrackCand_Z[nTrackCand] = track->vertex().z();
      
      if((TrackCand_purity[nTrackCand] == 1) && (TrackCand_nhits[nTrackCand] >= 3))
	nQualityTrackCand++;
      
      if(TrackCand_vtxdxyz[nTrackCand] < 0.1)
	MuMu_extratracks1mm++;
      if(TrackCand_vtxdxyz[nTrackCand] < 0.2) 
	MuMu_extratracks2mm++; 
      if(TrackCand_vtxdxyz[nTrackCand] < 0.3)
	MuMu_extratracks3mm++;
      if(TrackCand_vtxdxyz[nTrackCand] < 0.4) 
	MuMu_extratracks4mm++; 
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
      if((TrackCand_vtxdxyz[nTrackCand] < closesthighpuritytrkdxyz) && 
	 (TrackCand_purity[nTrackCand] == 1) && 
	 (TrackCand_nhits[nTrackCand] >= 3))
	closesthighpuritytrkdxyz = TrackCand_vtxdxyz[nTrackCand];
      
      nTrackCand++;  
    } 
  ClosestExtraTrack_vtxdxyz = closesttrkdxyz;
  ClosestHighPurityExtraTrack_vtxdxyz = closesthighpuritytrkdxyz;

  // Check for di-objects with valid vertex
  if(nMuonCand < 2 || !(found_pair)) {passed = false;}

  // Comment for keeping background events
  if(!(found_mumuvertex)) {passed = false;}
  //  if(ClosestHighPurityExtraTrack_vtxdxyz < minmumuvtxd) {passed = false;}

  // "Exclusivity" cuts
  if(passed == true){
    thetree->Fill();
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
  iDesc.add<bool>("ReadMCEffCorrections", false)->setComment("Flag to read Tag-&-Probe eff. corrections when running on MC"); 
  iDesc.add<bool>("ReadMCEffCorrectionsByCharge", false)->setComment("Flag to read Tag-&-Probe eff. corrections vs. Charge when running on MC");
  iDesc.add<bool>("ReadmcEffCorrectionsBySignedEta", false)->setComment("Flag to read Tag-&-Probe eff. corrections vs. signed eta when running on MC");
  iDesc.add<bool>("ReadMCPileup", false)->setComment("Flag to read pileup information when running on MC");

  iDesc.add<double>("MinMuMuVertexSeparation",0.1)->setComment("Minimum distance in cm between the dimuon vertex and any other track");

  descriptions.add("ParameterDescriptionsForGammaGammaMuMu", iDesc);
}

// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaMuMu::beginJob()
{
}

void
GammaGammaMuMu::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
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
	  cout << "GammaGammaMuMu::analyze:"
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
    cout << "GammaGammaMuMu::beginRun:"
	 << " config extraction failure with process name "
	 << hltMenuLabel << endl;
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaMuMu::endJob() {
  const edm::ParameterSet &thepset = edm::getProcessParameterSet();
  TList *list = thetree->GetUserInfo();
  list->Add(new TObjString(thepset.dump().c_str()));
  thefile->Write();
  thefile->Close();
}
  
