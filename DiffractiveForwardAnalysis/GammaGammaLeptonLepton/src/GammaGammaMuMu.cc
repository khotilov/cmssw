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
// $Id: GammaGammaMuMu.cc,v 1.26 2008/09/29 12:16:25 jjhollar Exp $
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
 
#include "FWCore/ParameterSet/interface/ParameterSet.h"  
#include "DataFormats/Common/interface/Ref.h"  
 
#include "DataFormats/Common/interface/TriggerResults.h"  
#include "FWCore/Framework/interface/TriggerNames.h"  
   
#include "FWCore/Framework/interface/ESHandle.h" 
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h" 
#include "DataFormats/JetReco/interface/CaloJetCollection.h" 
#include "DataFormats/JetReco/interface/CaloJet.h"  
#include "DataFormats/EgammaCandidates/interface/Electron.h"  
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"   
#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h"  
#include "DataFormats/TrackReco/interface/Track.h"  
#include "DataFormats/TrackReco/interface/TrackFwd.h"  
#include "DataFormats/MuonReco/interface/Muon.h"  
#include "DataFormats/MuonReco/interface/MuonFwd.h"   
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
  theGLBMuonLabel    = pset.getParameter<edm::InputTag>("GlobalMuonCollectionLabel");
  thePixelGsfELabel  = pset.getParameter<edm::InputTag>("ElectronCollectionLabel");
  theJetLabel        = pset.getParameter<edm::InputTag>("JetCollectionLabel");
  theMetLabel        = pset.getParameter<edm::InputTag>("MetLabel");
  thePhotonLabel     = pset.getParameter<edm::InputTag>("PhotonCollectionLabel");
  theCaloTowLabel    = pset.getParameter<edm::InputTag>("CaloTowerLabel");

  mudptmax           = pset.getParameter<double>("DimuonMaxdpt");
  mudphimin          = pset.getParameter<double>("DimuonMindphi");
  drisocalo          = pset.getParameter<double>("CaloTowerdR");

  rootfilename       = pset.getUntrackedParameter<std::string>("outfilename","test.root");

  edm::FileInPath myDataFile("FastSimulation/ProtonTaggers/data/acceptance_420_220.root");  
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

  //  nEvt=0;
  MUONMAX=10;
  JETMAX=30;
  TRACKMAX=100;

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
  thetree->Branch("MuonCand_tm2dloosemuid",MuonCand_tm2dloosemuid,"MuonCand_tm2dloosemuid[nMuonCand]/I");
  thetree->Branch("MuonCand_arbmuid",MuonCand_arbmuid,"MuonCand_arbmuid[nMuonCand]/I");
  thetree->Branch("MuonCand_isglobal",MuonCand_isglobal,"MuonCand_isglobal[nMuonCand]/I");
  thetree->Branch("MuonCand_istracker",MuonCand_istracker,"MuonCand_istracker[nMuonCand]/I"); 
  thetree->Branch("MuonCand_isstandalone",MuonCand_isstandalone,"MuonCand_isstandalone[nMuonCand]/I"); 
  thetree->Branch("MuonCand_hcalisor3",MuonCand_hcalisor3,"MuonCand_hcalisor3[nMuonCand]/D"); 
  thetree->Branch("MuonCand_ecalisor3",MuonCand_ecalisor3,"MuonCand_ecalisor3[nMuonCand]/D");  
  thetree->Branch("MuonCand_trkisor3",MuonCand_trkisor3,"MuonCand_trkisor3[nMuonCand]/D");  
  thetree->Branch("MuonCand_hcalisor5",MuonCand_hcalisor5,"MuonCand_hcalisor5[nMuonCand]/D");  
  thetree->Branch("MuonCand_ecalisor5",MuonCand_ecalisor5,"MuonCand_ecalisor5[nMuonCand]/D");  
  thetree->Branch("MuonCand_trkisor5",MuonCand_trkisor5,"MuonCand_trkisor5[nMuonCand]/D"); 

  thetree->Branch("nCaloCand",&nCaloCand,"nCaloCand/I");
  thetree->Branch("CaloTower_e",CaloTower_e,"CaloTower_e[nCaloCand]/D");
  thetree->Branch("CaloTower_et",CaloTower_et,"CaloTower_et[nCaloCand]/D");
  thetree->Branch("CaloTower_eta",CaloTower_eta,"CaloTower_eta[nCaloCand]/D"); 
  thetree->Branch("CaloTower_phi",CaloTower_phi,"CaloTower_phi[nCaloCand]/D"); 
  thetree->Branch("CaloTower_dr",CaloTower_dr,"CaloTower_dr[nCaloCand]/D");
  thetree->Branch("HighestCaloTower_e",&HighestCaloTower_e,"HighestCaloTower_e/D");
  thetree->Branch("HighestCaloTower_eta",&HighestCaloTower_eta,"HighestCaloTower_eta/D");
  thetree->Branch("HighestCaloTower_phi",&HighestCaloTower_phi,"HighestCaloTower_phi/D"); 
  thetree->Branch("HighestCaloTower_dr",&HighestCaloTower_dr,"HighestCaloTower_dr/D");
  thetree->Branch("HighestEtCaloTower_et",&HighestEtCaloTower_et,"HighestEtCaloTower_et/D");
  thetree->Branch("HighestEtCaloTower_eta",&HighestEtCaloTower_eta,"HighestEtCaloTower_eta/D");
  thetree->Branch("HighestEtCaloTower_phi",&HighestEtCaloTower_phi,"HighestEtCaloTower_phi/D"); 
  thetree->Branch("HighestEtCaloTower_dr",&HighestEtCaloTower_dr,"HighestEtCaloTower_dr/D");
  thetree->Branch("SumCalo_e",&SumCalo_e,"SumCalo_e/D");

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
  thetree->Branch("nExtraCaloTowersE2hf", &nExtraCaloTowersE2hf, "nExtraCaloTowersE12hf/I"); 
  thetree->Branch("nExtraCaloTowersE1he", &nExtraCaloTowersE1he, "nExtraCaloTowersE1he/I"); 
  thetree->Branch("nExtraCaloTowersE2he", &nExtraCaloTowersE2he, "nExtraCaloTowersE2he/I");  
  thetree->Branch("nExtraCaloTowersE3he", &nExtraCaloTowersE3he, "nExtraCaloTowersE3he/I");  
  thetree->Branch("nExtraCaloTowersE2hb", &nExtraCaloTowersE2hb, "nExtraCaloTowersE2hb/I"); 
  thetree->Branch("nExtraCaloTowersE3hb", &nExtraCaloTowersE3hb, "nExtraCaloTowersE3hb/I");  
  thetree->Branch("nExtraCaloTowersE4hb", &nExtraCaloTowersE4hb, "nExtraCaloTowersE4hb/I");  

  thetree->Branch("nExtraCaloTowersEt0pt1",&nExtraCaloTowersEt0pt1,"nExtraCaloTowersEt0pt1/I"); 
  thetree->Branch("nExtraCaloTowersEt0pt2",&nExtraCaloTowersEt0pt2,"nExtraCaloTowersEt0pt2/I"); 
  thetree->Branch("nExtraCaloTowersEt0pt5",&nExtraCaloTowersEt0pt5,"nExtraCaloTowersEt0pt5/I"); 
  thetree->Branch("nExtraCaloTowersEt1",&nExtraCaloTowersEt1,"nExtraCaloTowersEt1/I"); 
  thetree->Branch("nExtraCaloTowersEt2",&nExtraCaloTowersEt2,"nExtraCaloTowersEt2/I"); 
  thetree->Branch("nExtraCaloTowersEt3",&nExtraCaloTowersEt3,"nExtraCaloTowersEt3/I");  
  thetree->Branch("nExtraCaloTowersEt4",&nExtraCaloTowersEt4,"nExtraCaloTowersEt4/I");  

  thetree->Branch("nTrackCand",&nTrackCand,"nTrackCand/I");
  thetree->Branch("nExtraTrackCand",&nExtraTrackCand,"nExtraTrackCand/I"); 
  thetree->Branch("TrackCand_px",TrackCand_px,"TrackCand_px[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_py",TrackCand_py,"TrackCand_py[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_pz",TrackCand_pz,"TrackCand_pz[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_p",TrackCand_p,"TrackCand_p[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_pt",TrackCand_pt,"TrackCand_pt[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_eta",TrackCand_eta,"TrackCand_eta[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_phi",TrackCand_phi,"TrackCand_phi[nExtraTrackCand]/D");
  thetree->Branch("TrackCand_vtxdxyz",TrackCand_vtxdxyz,"TrackCand_vtxdxyz[nExtraTrackCand]/D");
  thetree->Branch("ClosestExtraTrack_vtxdxyz",&ClosestExtraTrack_vtxdxyz,"ClosestExtraTrack_vtxdxyz/D");
  
  thetree->Branch("MuMu_mass",&MuMu_mass,"MuMu_mass/D");
  thetree->Branch("MuMu_dphi",&MuMu_dphi,"MuMu_dphi/D");
  thetree->Branch("MuMu_vtxx",&MuMu_vtxx,"MuMu_vtxx/D");
  thetree->Branch("MuMu_vtxy",&MuMu_vtxy,"MuMu_vtxy/D"); 
  thetree->Branch("MuMu_vtxz",&MuMu_vtxz,"MuMu_vtxz/D"); 
  thetree->Branch("MuMu_vtxchi2dof",&MuMu_vtxchi2dof,"MuMu_vtxchi2dof/D");
  thetree->Branch("MuMu_vtxisvalid",&MuMu_vtxisvalid,"MuMu_vtxisvalid/I");

  thetree->Branch("MuMu_extratracks1mm",&MuMu_extratracks1mm,"MuMu_extratracks1mm/I");
  thetree->Branch("MuMu_extratracks2mm",&MuMu_extratracks2mm,"MuMu_extratracks2mm/I");
  thetree->Branch("MuMu_extratracks3mm",&MuMu_extratracks3mm,"MuMu_extratracks3mm/I");
  thetree->Branch("MuMu_extratracks4mm",&MuMu_extratracks4mm,"MuMu_extratracks4mm/I"); 
  thetree->Branch("MuMu_extratracks5mm",&MuMu_extratracks5mm,"MuMu_extratracks5mm/I");  
  thetree->Branch("MuMu_extratracks1cm",&MuMu_extratracks1cm,"MuMu_extratracks1cm/I"); 
  thetree->Branch("MuMu_extratracks2cm",&MuMu_extratracks2cm,"MuMu_extratracks2cm/I");
  thetree->Branch("MuMu_extratracks5cm",&MuMu_extratracks5cm,"MuMu_extratracks5cm/I"); 
  thetree->Branch("MuMu_extratracks10cm",&MuMu_extratracks10cm,"MuMu_extratracks10cm/I"); 

  thetree->Branch("nPFlowCand",&nPFlowCand,"nPFlowCand/I");
  thetree->Branch("PFowCandIds",PFlowCandIds,"PFlowCandIds[nPFlowCand]/I");

  thetree->Branch("HitInZDC",&HitInZDC,"HitInZDC/I");
  thetree->Branch("HitInCastor",&HitInCastor,"HitInCastor/I");
  
  thetree->Branch("Etmiss",&Etmiss,"Etmiss/D");

  thetree->Branch("HLT2MuonNonIso",&HLT2MuonNonIso,"HLT2MuonNonIso/I");
  thetree->Branch("HLT1MuonPrescalePt3",&HLT1MuonPrescalePt3,"HLT1MuonPrescalePt3/I"); 

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
  nMuonCand=0;
  nJetCand=0;
  nCaloCand=0;
  nTrackCand=0;
  nExtraTrackCand=0;

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
  nExtraCaloTowersE1he=0; 
  nExtraCaloTowersE2he=0; 
  nExtraCaloTowersE3he=0;  
  nExtraCaloTowersE2hb=0; 
  nExtraCaloTowersE3hb=0; 
  nExtraCaloTowersE4hb=0;  

  nPFlowCand=0;

  HitInZDC=0;
  HitInCastor=0;

  MuMu_mass = -1;
  MuMu_dphi = -1;

  MuMu_extratracks1mm = 0;
  MuMu_extratracks2mm = 0;
  MuMu_extratracks3mm = 0;
  MuMu_extratracks4mm = 0;
  MuMu_extratracks5mm = 0;
  MuMu_extratracks1cm = 0; 
  MuMu_extratracks2cm = 0;
  MuMu_extratracks5cm = 0;
  MuMu_extratracks10cm = 0;

  bool passed = true;

 //using namespace edm;
  using reco::TrackCollection;

  //  Handle< double> weightHandle;
  //  event.getByLabel ("weight", weightHandle);
  //  evweight = * weightHandle;


  // Get the trigger information from the event
  edm::Handle<edm::TriggerResults> hltResults ; 
  event.getByLabel(InputTag("TriggerResults::HLT"),hltResults) ; 
  trigNames.init(*hltResults) ; 
  for (unsigned int i=0; i<trigNames.size(); i++)  
    { 
      // This is for CMSSW_2_1_X!!! 
      if ( trigNames.triggerNames().at(i) == "HLT_Mu3" )       

      // This is for CMSSW_2_0_X!!!
      //      if ( trigNames.triggerNames().at(i) == "HLT1MuonPrescalePt3" )
        {  
          if ( hltResults->accept(i) )  
            HLT1MuonPrescalePt3 = 1;
	  else
	    HLT1MuonPrescalePt3 = 0;
        }  
	
      // This is for CMSSW_2_1_X!!! 	
      if ( trigNames.triggerNames().at(i) == "HLT_DoubleMu3" ) 

      // This is for CMSSW_2_0_X!!! 
	//      if ( trigNames.triggerNames().at(i) == "HLT2MuonNonIso" )
        {   
          if ( hltResults->accept(i) )  
	    HLT2MuonNonIso = 1;
	  else
	    HLT2MuonNonIso = 0;
        }   
    }

  // Get the muon collection from the event
  // PAT
  edm::Handle<edm::View<pat::Muon> > muons; 
  event.getByLabel(theGLBMuonLabel,muons); 
  edm::View<pat::Muon>::const_iterator muon;

  // AOD
  //  Handle<reco::MuonCollection> muons;
  //  event.getByLabel(theGLBMuonLabel, muons);
  //  reco::MuonCollection::const_iterator muon;

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

// 	  // Muon ID for CMSSW_2_0_X
// 	  MuonCand_tmlsloosemuonid[nMuonCand]=muonid::isGoodMuon(*muon,muonid::TMLastStationLoose);
// 	  MuonCand_tm2dloosemuid[nMuonCand]=muonid::isGoodMuon(*muon,muonid::TM2DCompatibilityLoose);
// 	  MuonCand_arbmuid[nMuonCand]=-1; // This doesn't exist in 2_0_X!!!
// 	  MuonCand_isglobal[nMuonCand]=muon->isGlobalMuon();
//           MuonCand_istracker[nMuonCand]=muon->isTrackerMuon(); 
//           MuonCand_isstandalone[nMuonCand]=muon->isStandAloneMuon(); 

// 	  // Isolation for CMSSW_2_0_X
// 	  MuonCand_hcalisor3[nMuonCand]=muon->isolationR03().hadEt; 
// 	  MuonCand_ecalisor3[nMuonCand]=muon->isolationR03().emEt;  
// 	  MuonCand_trkisor3[nMuonCand]=muon->isolationR03().nTracks;  
// 	  MuonCand_hcalisor5[nMuonCand]=muon->isolationR05().hadEt;  
// 	  MuonCand_ecalisor5[nMuonCand]=muon->isolationR05().emEt;   
// 	  MuonCand_trkisor5[nMuonCand]=muon->isolationR05().nTracks;   

	  // Muon ID for CMSSW_2_1_X
	  MuonCand_tmlsloosemuonid[nMuonCand]=muon->isGood(reco::Muon::TMLastStationLoose);
	  MuonCand_tm2dloosemuid[nMuonCand]=muon->isGood(reco::Muon::TM2DCompatibilityLoose);
	  MuonCand_arbmuid[nMuonCand]=muon->isGood(reco::Muon::AllArbitrated);
	  MuonCand_isglobal[nMuonCand]=muon->isGlobalMuon();
	  MuonCand_istracker[nMuonCand]=muon->isTrackerMuon(); 
	  MuonCand_isstandalone[nMuonCand]=muon->isStandAloneMuon(); 

	  // Isolation for CMSSW_2_1_X
	  MuonCand_hcalisor3[nMuonCand]=muon->isolationR03().hadEt; 
	  MuonCand_ecalisor3[nMuonCand]=muon->isolationR03().emEt;  
	  MuonCand_trkisor3[nMuonCand]=muon->isolationR03().nTracks;  
	  MuonCand_hcalisor5[nMuonCand]=muon->isolationR05().hadEt;  
	  MuonCand_ecalisor5[nMuonCand]=muon->isolationR05().emEt;   
	  MuonCand_trkisor5[nMuonCand]=muon->isolationR05().nTracks;   

	  if(MuonCand_isglobal[nMuonCand] || MuonCand_istracker[nMuonCand])
	    MuonCandTrack_p[nMuonCand] = muon->track()->p(); 
	  else
	    MuonCandTrack_p[nMuonCand] = muon->p();

	  nMuonCand++;
	}  

      // Calculate invariant mass and delta-phi
      if(nMuonCand == 2)
	{
	  if(MuonCand_charge[0]*MuonCand_charge[1]<0)
	    {
	      double mass = pow(MuonCand_p[0]+MuonCand_p[1],2);
	      mass-=pow(MuonCand_px[0]+MuonCand_px[1],2);
	      mass-=pow(MuonCand_py[0]+MuonCand_py[1],2);
	      mass-=pow(MuonCand_pz[0]+MuonCand_pz[1],2);
	      MuMu_mass = sqrt(mass);
	      
	      double dphi = fabs(MuonCand_phi[0]-MuonCand_phi[1]);
	      if(dphi < 3.14159)
		MuMu_dphi = dphi;
	      else
		MuMu_dphi = (2.0*3.14159)-dphi;
	      
	    }
	}
    }


  // Get the Jet collection from the event
  // PAT
  edm::Handle<edm::View<pat::Jet> > jets; 
  event.getByLabel(theJetLabel,jets); 
  edm::View<pat::Jet>::const_iterator jet;

  // AOD
  //  edm::Handle<reco::CaloJetCollection> pJets;
  //  event.getByLabel(theJetLabel,pJets);
  //  const reco::CaloJetCollection* jets = pJets.product();
  //  reco::CaloJetCollection::const_iterator jet;

  // Get the MET collection from the event

  // PAT
  edm::Handle<edm::View<pat::MET> > mets; 
  event.getByLabel(theMetLabel,mets); 
  edm::View<pat::MET>::const_iterator met;

  // AOD
  //  edm::Handle<reco::CaloMETCollection> pMET;
  //  event.getByLabel(theMetLabel,pMET);
  //  const reco::CaloMETCollection* mets = pMET.product();
  //  reco::CaloMETCollection::const_iterator met;

  // Count PFlow  objects
  edm::Handle<reco::PFCandidateCollection> pflows;
  event.getByLabel("particleFlow",pflows);
  reco::PFCandidateCollection::const_iterator pflow;

  for(pflow = pflows->begin(); pflow != pflows->end(); ++pflow)
    {
      int parttype = PFCandidate::ParticleType (pflow->particleId());
      if(parttype != 2)
	{
	  PFlowCandIds[nPFlowCand] = parttype;
	  nPFlowCand++;
	}
    }

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
              if(CaloTower_e[nCaloCand] > 1.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5)      
		nExtraCaloTowersE1he++;  
              if(CaloTower_e[nCaloCand] > 2.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5)   
		nExtraCaloTowersE2he++;  
	      if(CaloTower_e[nCaloCand] > 3.0 && abs(CaloTower_eta[nCaloCand]) < 3.0 && abs(CaloTower_eta[nCaloCand]) > 1.5)   
		nExtraCaloTowersE3he++;   
              if(CaloTower_e[nCaloCand] > 2.0 && abs(CaloTower_eta[nCaloCand]) < 1.5)   
		nExtraCaloTowersE2hb++;  
              if(CaloTower_e[nCaloCand] > 3.0 && abs(CaloTower_eta[nCaloCand]) < 1.5)   
		nExtraCaloTowersE3hb++;  
              if(CaloTower_e[nCaloCand] > 4.0 && abs(CaloTower_eta[nCaloCand]) < 1.5)   
		nExtraCaloTowersE4hb++;   


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
    }

  // Check for particles in ZDC/Castor acceptance. 
  // Use MC truth for now, replace with real RECO when available
  double MCPar_px,MCPar_py,MCPar_pz,MCPar_e,MCPar_eta,MCPar_mass;
  int MCPar_pdgid;
  
  Handle<GenParticleCollection> genParticles;
  event.getByLabel( "genParticles", genParticles );
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

  // Now do vertexing and track counting
  edm::ESHandle<TransientTrackBuilder> theVtx;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theVtx);
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
	    {
	      isMuon = true;
	      mutrks->push_back( *track );
	      TransientTrack tmptrk = (*theVtx).build( *track );
	      transmutrks.push_back( tmptrk );
	    }
	}
      if(isMuon == false)
	nTrackCand++;
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
	  MuMu_vtxisvalid = 1;
	}
      else
	{
	  MuMu_vtxx = 0;  
	  MuMu_vtxy = 0;  
	  MuMu_vtxz = 0;  
	  MuMu_vtxchi2dof = 0;
	  MuMu_vtxisvalid = 0;
	}

      // OK, now go back and count "extra" tracks on the dimuon vertex
      for(track = tracks->begin(); track != tracks->end() && nExtraTrackCand<TRACKMAX; ++ track) 
        { 
	  if(track->p() == MuonCandTrack_p[0] || track->p() == MuonCandTrack_p[1])
	    continue;
	  
          TrackCand_p[nExtraTrackCand]=track->p(); 
          TrackCand_px[nExtraTrackCand]=track->px(); 
          TrackCand_py[nExtraTrackCand]=track->py(); 
          TrackCand_pz[nExtraTrackCand]=track->pz(); 
          TrackCand_pt[nExtraTrackCand]=track->pt(); 
          TrackCand_eta[nExtraTrackCand]=track->eta(); 
          TrackCand_phi[nExtraTrackCand]=track->phi(); 
          TrackCand_charge[nExtraTrackCand]=track->charge(); 
	  TrackCand_vtxdxyz[nExtraTrackCand] = sqrt(((track->vertex().x() - MuMu_vtxx)*(track->vertex().x() - MuMu_vtxx)) + 
					     ((track->vertex().y() - MuMu_vtxy)*(track->vertex().y() - MuMu_vtxy)) +
					     ((track->vertex().z() - MuMu_vtxz)*(track->vertex().z() - MuMu_vtxz)));
	  
	  if(TrackCand_vtxdxyz[nExtraTrackCand] < 0.1) 
            MuMu_extratracks1mm++; 
          if(TrackCand_vtxdxyz[nExtraTrackCand] < 0.2)  
            MuMu_extratracks2mm++;  
          if(TrackCand_vtxdxyz[nExtraTrackCand] < 0.3)  
            MuMu_extratracks3mm++;  
          if(TrackCand_vtxdxyz[nExtraTrackCand] < 0.4)   
            MuMu_extratracks4mm++;   
	  if(TrackCand_vtxdxyz[nExtraTrackCand] < 0.5)    
            MuMu_extratracks5mm++;
	  if(TrackCand_vtxdxyz[nExtraTrackCand] < 1) 
            MuMu_extratracks1cm++; 
	  if(TrackCand_vtxdxyz[nExtraTrackCand] < 2)
	    MuMu_extratracks2cm++;
          if(TrackCand_vtxdxyz[nExtraTrackCand] < 5) 
            MuMu_extratracks5cm++; 
          if(TrackCand_vtxdxyz[nExtraTrackCand] < 10) 
            MuMu_extratracks10cm++; 
	  if(TrackCand_vtxdxyz[nExtraTrackCand] < closesttrkdxyz)
	    closesttrkdxyz = TrackCand_vtxdxyz[nExtraTrackCand];

          nExtraTrackCand++;  
        } 
      ClosestExtraTrack_vtxdxyz = closesttrkdxyz;
    } 
  else 
    { 
      MuMu_vtxx = 0; 
      MuMu_vtxy = 0; 
      MuMu_vtxz = 0; 
      MuMu_vtxchi2dof = 0; 
      MuMu_vtxisvalid = 0; 
      nExtraTrackCand = 0;
    } 

  // Check for di-objects
  if(nMuonCand != 2)
    passed = false;
  else
    {
      if(MuMu_dphi < mudphimin) 
	passed = false;
      if(fabs(MuonCand_pt[0]-MuonCand_pt[1]) > mudptmax)
	passed = false;      
    }

  // "Exclusivity" cuts

  if(passed == true)
    thetree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaMuMu::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaMuMu::endJob() {
  thefile->Write();
  thefile->Close();
}
  
