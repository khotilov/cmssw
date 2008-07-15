
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2008/07/12 10:47:19 $
 *  $Revision: 1.11 $
 *  \author G. Mila - INFN Torino
 */

#include "DQMOffline/Muon/src/MuonRecoAnalyzer.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h" 
#include "DataFormats/MuonReco/interface/MuonEnergy.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>
using namespace std;
using namespace edm;



MuonRecoAnalyzer::MuonRecoAnalyzer(const edm::ParameterSet& pSet, MuonServiceProxy *theService):MuonAnalyzerBase(theService) {

 cout<<"[MuonRecoAnalyzer] Constructor called!"<<endl;
  parameters = pSet;

}


MuonRecoAnalyzer::~MuonRecoAnalyzer() { }


void MuonRecoAnalyzer::beginJob(edm::EventSetup const& iSetup,DQMStore * dbe) {

  metname = "muRecoAnalyzer";

  LogTrace(metname)<<"[MuonRecoAnalyzer] Parameters initialization";
  dbe->setCurrentFolder("Muons/MuonRecoAnalyzer");

  muReco = dbe->book1D("muReco", "muReco", 6, 1, 7);
  muReco->setBinLabel(1,"glb+tk+sta");
  muReco->setBinLabel(2,"glb+sta");
  muReco->setBinLabel(3,"tk+sta");
  muReco->setBinLabel(4,"tk");
  muReco->setBinLabel(5,"sta");
  muReco->setBinLabel(6,"calo");

  int binFactor = 4;

  // monitoring of eta parameter
  etaBin = parameters.getParameter<int>("etaBin");
  etaMin = parameters.getParameter<double>("etaMin");
  etaMax = parameters.getParameter<double>("etaMax");
  std::string histname = "GlbMuon_";
  etaGlbTrack.push_back(dbe->book1D(histname+"Glb_eta", histname+"Glb_eta", etaBin, etaMin, etaMax));
  etaGlbTrack.push_back(dbe->book1D(histname+"Tk_eta", histname+"Tk_eta", etaBin, etaMin, etaMax));
  etaGlbTrack.push_back(dbe->book1D(histname+"Sta_eta", histname+"Sta_eta", etaBin, etaMin, etaMax));
  etaResolution.push_back(dbe->book1D("Res_TkGlb_eta", "Res_TkGlb_eta", etaBin*binFactor, etaMin/3000, etaMax/3000));
  etaResolution.push_back(dbe->book1D("Res_GlbSta_eta", "Res_GlbSta_eta", etaBin*binFactor, etaMin/100, etaMax/100));
  etaResolution.push_back(dbe->book1D("Res_TkSta_eta", "Res_TkSta_eta", etaBin*binFactor, etaMin/100, etaMax/100));
  etaResolution.push_back(dbe->book2D("ResVsEta_TkGlb_eta", "ResVsEta_TkGlb_eta", etaBin, etaMin, etaMax, etaBin*binFactor, etaMin/3000, etaMax/3000));
  etaResolution.push_back(dbe->book2D("ResVsEta_GlbSta_eta", "ResVsEta_GlbSta_eta", etaBin, etaMin, etaMax, etaBin*binFactor, etaMin/100, etaMax/100));
  etaResolution.push_back(dbe->book2D("ResVsEta_TkSta_eta", "ResVsTkEta_TkSta_eta", etaBin, etaMin, etaMax, etaBin*binFactor, etaMin/100, etaMax/100));
  etaTrack = dbe->book1D("TkMuon_eta", "TkMuon_eta", etaBin, etaMin, etaMax);
  etaStaTrack = dbe->book1D("StaMuon_eta", "StaMuon_eta", etaBin, etaMin, etaMax);
  etaEfficiency.push_back(dbe->book1D("StaEta", "StaEta", etaBin, etaMin, etaMax));
  etaEfficiency.push_back(dbe->book1D("StaEta_ifCombinedAlso", "StaEta_ifCombinedAlso", etaBin, etaMin, etaMax));

  // monitoring of theta parameter
  thetaBin = parameters.getParameter<int>("thetaBin");
  thetaMin = parameters.getParameter<double>("thetaMin");
  thetaMax = parameters.getParameter<double>("thetaMax");
  thetaGlbTrack.push_back(dbe->book1D(histname+"Glb_theta", histname+"Glb_theta", thetaBin, thetaMin, thetaMax));
  thetaGlbTrack.push_back(dbe->book1D(histname+"Tk_theta", histname+"Tk_theta", thetaBin, thetaMin, thetaMax));
  thetaGlbTrack.push_back(dbe->book1D(histname+"Sta_theta", histname+"Sta_theta", thetaBin, thetaMin, thetaMax));
  thetaResolution.push_back(dbe->book1D("Res_TkGlb_theta", "Res_TkGlb_theta", thetaBin*binFactor, -(thetaMax/3000), thetaMax/3000));
  thetaResolution.push_back(dbe->book1D("Res_GlbSta_theta", "Res_GlbSta_theta", thetaBin*binFactor,-(thetaMax/100), thetaMax/100));
  thetaResolution.push_back(dbe->book1D("Res_TkSta_theta", "Res_TkSta_theta", thetaBin*binFactor, -(thetaMax/100), thetaMax/100));
  thetaResolution.push_back(dbe->book2D("ResVsTheta_TkGlb_theta", "ResVsTheta_TkGlb_theta", thetaBin, thetaMin, thetaMax, thetaBin*binFactor, -(thetaMax/3000), thetaMax/3000));
  thetaResolution.push_back(dbe->book2D("ResVsTheta_GlbSta_theta", "ResVsTheta_GlbSta_theta", thetaBin, thetaMin, thetaMax, thetaBin*binFactor, -(thetaMax/100), thetaMax/100));
  thetaResolution.push_back(dbe->book2D("ResVsTheta_TkSta_theta", "ResVsTheta_TkSta_theta", thetaBin, thetaMin, thetaMax, thetaBin*binFactor, -(thetaMax/100), thetaMax/100));
  thetaTrack = dbe->book1D("TkMuon_theta", "TkMuon_theta", thetaBin, thetaMin, thetaMax);
  thetaStaTrack = dbe->book1D("StaMuon_theta", "StaMuon_theta", thetaBin, thetaMin, thetaMax);

  // monitoring of phi paramater
  phiBin = parameters.getParameter<int>("phiBin");
  phiMin = parameters.getParameter<double>("phiMin");
  phiMax = parameters.getParameter<double>("phiMax");
  phiGlbTrack.push_back(dbe->book1D(histname+"Glb_phi", histname+"Glb_phi", phiBin, phiMin, phiMax));
  phiGlbTrack.push_back(dbe->book1D(histname+"Tk_phi", histname+"Tk_phi", phiBin, phiMin, phiMax));
  phiGlbTrack.push_back(dbe->book1D(histname+"Sta_phi", histname+"Sta_phi", phiBin, phiMin, phiMax));
  phiResolution.push_back(dbe->book1D("Res_TkGlb_phi", "Res_TkGlb_phi", phiBin*binFactor, phiMin/3000, phiMax/3000));
  phiResolution.push_back(dbe->book1D("Res_GlbSta_phi", "Res_GlbSta_phi", phiBin*binFactor, phiMin/100, phiMax/100));
  phiResolution.push_back(dbe->book1D("Res_TkSta_phi", "Res_TkSta_phi", phiBin*binFactor, phiMin/100, phiMax/100));
  phiResolution.push_back(dbe->book2D("ResVsPhi_TkGlb_phi", "ResVsPhi_TkGlb_phi", phiBin, phiMin, phiMax, phiBin*binFactor, phiMin/3000, phiMax/3000));
  phiResolution.push_back(dbe->book2D("ResVsPhi_GlbSta_phi", "ResVsPhi_GlbSta_phi", phiBin, phiMin, phiMax, phiBin*binFactor, phiMin/100, phiMax/100));
  phiResolution.push_back(dbe->book2D("ResVsPhi_TkSta_phi", "ResVsTkPhi_TkSta_phi", phiBin, phiMin, phiMax, phiBin*binFactor, phiMin/100, phiMax/100));
  phiTrack = dbe->book1D("TkMuon_phi", "TkMuon_phi", phiBin, phiMin, phiMax);
  phiStaTrack = dbe->book1D("StaMuon_phi", "StaMuon_phi", phiBin, phiMin, phiMax);
  phiEfficiency.push_back(dbe->book1D("StaPhi", "StaPhi", phiBin, phiMin, phiMax));
  phiEfficiency.push_back(dbe->book1D("StaPhi_ifCombinedAlso", "StaPhi_ifCombinedAlso", phiBin, phiMin, phiMax));

  // monitoring of the momentum
  pBin = parameters.getParameter<int>("pBin");
  pMin = parameters.getParameter<double>("pMin");
  pMax = parameters.getParameter<double>("pMax");
  pGlbTrack.push_back(dbe->book1D(histname+"Glb_p", histname+"Glb_p", pBin, pMin, pMax));
  pGlbTrack.push_back(dbe->book1D(histname+"Tk_p", histname+"Tk_p", pBin, pMin, pMax));
  pGlbTrack.push_back(dbe->book1D(histname+"Sta_p", histname+"Sta_p", pBin, pMin, pMax));
  pTrack = dbe->book1D("TkMuon_p", "TkMuon_p", pBin, pMin, pMax);
  pStaTrack = dbe->book1D("StaMuon_p", "StaMuon_p", pBin, pMin, pMax);

  // monitoring of the transverse momentum
  ptBin = parameters.getParameter<int>("ptBin");
  ptMin = parameters.getParameter<double>("ptMin");
  ptMax = parameters.getParameter<double>("ptMax");
  ptGlbTrack.push_back(dbe->book1D(histname+"Glb_pt", histname+"Glb_pt", ptBin, ptMin, ptMax));
  ptGlbTrack.push_back(dbe->book1D(histname+"Tk_pt", histname+"Tk_pt", ptBin, ptMin, ptMax));
  ptGlbTrack.push_back(dbe->book1D(histname+"Sta_pt", histname+"Sta_pt", ptBin, ptMin, ptMax));
  ptTrack = dbe->book1D("TkMuon_pt", "TkMuon_pt", ptBin, ptMin, ptMax);
  ptStaTrack = dbe->book1D("StaMuon_pt", "StaMuon_pt", ptBin, ptMin, pMax);

  // monitoring of the muon charge
  qGlbTrack.push_back(dbe->book1D(histname+"Glb_q", histname+"Glb_q", 5, -2.5, 2.5));
  qGlbTrack.push_back(dbe->book1D(histname+"Tk_q", histname+"Tk_q", 5, -2.5, 2.5));
  qGlbTrack.push_back(dbe->book1D(histname+"Sta_q", histname+"Sta_q", 5, -2.5, 2.5));
  qGlbTrack.push_back(dbe->book1D(histname+"qComparison", histname+"qComparison", 8, 0.5, 8.5));
  qGlbTrack[3]->setBinLabel(1,"qGlb=qSta");
  qGlbTrack[3]->setBinLabel(2,"qGlb!=qSta");
  qGlbTrack[3]->setBinLabel(3,"qGlb=qTk");
  qGlbTrack[3]->setBinLabel(4,"qGlb!=qTk");
  qGlbTrack[3]->setBinLabel(5,"qSta=qTk");
  qGlbTrack[3]->setBinLabel(6,"qSta!=qTk");
  qGlbTrack[3]->setBinLabel(7,"qGlb!=qSta,qGlb!=Tk");
  qGlbTrack[3]->setBinLabel(8,"qGlb=qSta,qGlb=Tk");
  qTrack = dbe->book1D("TkMuon_q", "TkMuon_q", 5, -2.5, 2.5);
  qStaTrack = dbe->book1D("StaMuon_q", "StaMuon_q", 5, -2.5, 2.5);

  // monitoring of the momentum resolution
  pResBin = parameters.getParameter<int>("pResBin");
  pResMin = parameters.getParameter<double>("pResMin");
  pResMax = parameters.getParameter<double>("pResMax");
  qOverpResolution.push_back(dbe->book1D("Res_TkGlb_qOverp", "Res_TkGlb_qOverp", pResBin*binFactor*2, pResMin/10, pResMax/10));
  qOverpResolution.push_back(dbe->book1D("Res_GlbSta_qOverp", "Res_GlbSta_qOverp", pResBin*binFactor, pResMin, pResMax));
  qOverpResolution.push_back(dbe->book1D("Res_TkSta_qOverp", "Res_TkSta_qOverp", pResBin*binFactor, pResMin, pResMax));
  oneOverpResolution.push_back(dbe->book1D("Res_TkGlb_oneOverp", "Res_TkGlb_oneOverp", pResBin*binFactor*2, pResMin/10, pResMax/10));
  oneOverpResolution.push_back(dbe->book1D("Res_GlbSta_oneOverp", "Res_GlbSta_oneOverp", pResBin*binFactor, pResMin, pResMax));
  oneOverpResolution.push_back(dbe->book1D("Res_TkSta_oneOverp", "Res_TkSta_oneOverp", pResBin*binFactor, pResMin, pResMax));
  qOverptResolution.push_back(dbe->book1D("Res_TkGlb_qOverpt", "Res_TkGlb_qOverpt", pResBin*binFactor*2, pResMin/10, pResMax/10));
  qOverptResolution.push_back(dbe->book1D("Res_GlbSta_qOverpt", "Res_GlbSta_qOverpt", pResBin*binFactor, pResMin, pResMax));
  qOverptResolution.push_back(dbe->book1D("Res_TkSta_qOverpt", "Res_TkSta_qOverpt", pResBin*binFactor, pResMin, pResMax));
  oneOverptResolution.push_back(dbe->book1D("Res_TkGlb_oneOverpt", "Res_TkGlb_oneOverpt", pResBin*binFactor*2, pResMin/10, pResMax/10));
  oneOverptResolution.push_back(dbe->book1D("Res_GlbSta_oneOverpt", "Res_GlbSta_oneOverpt", pResBin*binFactor, pResMin, pResMax));
  oneOverptResolution.push_back(dbe->book1D("Res_TkSta_oneOverpt", "Res_TkSta_oneOverpt", pResBin*binFactor, pResMin, pResMax));
  oneOverptResolution.push_back(dbe->book2D("ResVsEta_TkGlb_oneOverpt", "ResVsEta_TkGlb_oneOverpt", etaBin, etaMin, etaMax, pResBin*binFactor*2, pResMin/10, pResMax/10));
  oneOverptResolution.push_back(dbe->book2D("ResVsEta_GlbSta_oneOverpt", "ResVsEta_GlbSta_oneOverpt", etaBin, etaMin, etaMax, pResBin*binFactor, pResMin, pResMax));
  oneOverptResolution.push_back(dbe->book2D("ResVsEta_TkSta_oneOverpt", "ResVsEta_TkSta_oneOverpt", etaBin, etaMin, etaMax, pResBin*binFactor, pResMin, pResMax));
  oneOverptResolution.push_back(dbe->book2D("ResVsPhi_TkGlb_oneOverpt", "ResVsPhi_TkGlb_oneOverpt", phiBin, phiMin, phiMax, pResBin*binFactor*2, pResMin/10, pResMax/10));
  oneOverptResolution.push_back(dbe->book2D("ResVsPhi_GlbSta_oneOverpt", "ResVsPhi_GlbSta_oneOverpt", phiBin, phiMin, phiMax, pResBin*binFactor, pResMin, pResMax));
  oneOverptResolution.push_back(dbe->book2D("ResVsPhi_TkSta_oneOverpt", "ResVsPhi_TkSta_oneOverpt", phiBin, phiMin, phiMax, pResBin*binFactor, pResMin, pResMax));
  oneOverptResolution.push_back(dbe->book2D("ResVsPt_TkGlb_oneOverpt", "ResVsPt_TkGlb_oneOverpt", ptBin, ptMin, ptMax, pResBin*binFactor*2, pResMin/10, pResMax/10));
  oneOverptResolution.push_back(dbe->book2D("ResVsPt_GlbSta_oneOverpt", "ResVsPt_GlbSta_oneOverpt", ptBin, ptMin, ptMax, pResBin*binFactor, pResMin, pResMax));
  oneOverptResolution.push_back(dbe->book2D("ResVsPt_TkSta_oneOverpt", "ResVsPt_TkSta_oneOverpt", ptBin, ptMin, ptMax, pResBin*binFactor, pResMin, pResMax));

  // monitoring of the recHits provenance
  rhBin=parameters.getParameter<int>("rhBin");
  rhMin=parameters.getParameter<double>("rhMin");
  rhMax=parameters.getParameter<double>("rhMax");
  rhAnalysis.push_back(dbe->book1D("glb_hit_percentual_from_sta", "glb_hit_percentual_from_sta", rhBin, rhMin, rhMax));
  rhAnalysis.push_back(dbe->book1D("glb_hit_percentual_from_tk", "glb_hit_percentual_from_tk", rhBin, rhMin, rhMax));
  rhAnalysis.push_back(dbe->book1D("glb_hit_percentual_sta", "glb_hit_percentual_sta", rhBin, rhMin, rhMax));
  rhAnalysis.push_back(dbe->book1D("glb_hit_percentual_tk", "glb_hit_percentual_tk", rhBin, rhMin, rhMax));
  rhAnalysis.push_back(dbe->book1D("used_hit_percentual", "used_hit_percentual", rhBin, rhMin, rhMax));
  rhAnalysis.push_back(dbe->book1D("invalid_hit_percentual", "invalid_hit_percentual", rhBin, rhMin, rhMax));

}


void MuonRecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::Muon& recoMu) {

  LogTrace(metname)<<"[MuonRecoAnalyzer] Analyze the mu";

  if(recoMu.isGlobalMuon()) {

    LogTrace(metname)<<"[MuonRecoAnalyzer] The mu is global - filling the histos";
    if(recoMu.isTrackerMuon() && recoMu.isStandAloneMuon())
      muReco->Fill(1);
    if(!(recoMu.isTrackerMuon()) && recoMu.isStandAloneMuon())
      muReco->Fill(2);
    if(!recoMu.isStandAloneMuon())
      LogTrace(metname)<<"[MuonRecoAnalyzer] ERROR: the mu is global but not standalone!";

    // get the track combinig the information from both the Tracker and the Spectrometer
    reco::TrackRef recoCombinedGlbTrack = recoMu.combinedMuon();
    // get the track using only the tracker data
    reco::TrackRef recoGlbTrack = recoMu.track();
    // get the track using only the mu spectrometer data
    reco::TrackRef recoStaGlbTrack = recoMu.standAloneMuon();
  
    etaGlbTrack[0]->Fill(recoCombinedGlbTrack->eta());
    etaGlbTrack[1]->Fill(recoGlbTrack->eta());
    etaGlbTrack[2]->Fill(recoStaGlbTrack->eta());
    etaResolution[0]->Fill(recoGlbTrack->eta()-recoCombinedGlbTrack->eta());
    etaResolution[1]->Fill(-recoStaGlbTrack->eta()+recoCombinedGlbTrack->eta());
    etaResolution[2]->Fill(recoGlbTrack->eta()-recoStaGlbTrack->eta());
    etaResolution[3]->Fill(recoCombinedGlbTrack->eta(), recoGlbTrack->eta()-recoCombinedGlbTrack->eta());
    etaResolution[4]->Fill(recoCombinedGlbTrack->eta(), -recoStaGlbTrack->eta()+recoCombinedGlbTrack->eta());
    etaResolution[5]->Fill(recoCombinedGlbTrack->eta(), recoGlbTrack->eta()-recoStaGlbTrack->eta());

    thetaGlbTrack[0]->Fill(recoCombinedGlbTrack->theta());
    thetaGlbTrack[1]->Fill(recoGlbTrack->theta());
    thetaGlbTrack[2]->Fill(recoStaGlbTrack->theta());
    thetaResolution[0]->Fill(recoGlbTrack->theta()-recoCombinedGlbTrack->theta());
    thetaResolution[1]->Fill(-recoStaGlbTrack->theta()+recoCombinedGlbTrack->theta());
    thetaResolution[2]->Fill(recoGlbTrack->theta()-recoStaGlbTrack->theta());
    thetaResolution[3]->Fill(recoCombinedGlbTrack->theta(), recoGlbTrack->theta()-recoCombinedGlbTrack->theta());
    thetaResolution[4]->Fill(recoCombinedGlbTrack->theta(), -recoStaGlbTrack->theta()+recoCombinedGlbTrack->theta());
    thetaResolution[5]->Fill(recoCombinedGlbTrack->theta(), recoGlbTrack->theta()-recoStaGlbTrack->theta());
     
    phiGlbTrack[0]->Fill(recoCombinedGlbTrack->phi());
    phiGlbTrack[1]->Fill(recoGlbTrack->phi());
    phiGlbTrack[2]->Fill(recoStaGlbTrack->phi());
    phiResolution[0]->Fill(recoGlbTrack->phi()-recoCombinedGlbTrack->phi());
    phiResolution[1]->Fill(-recoStaGlbTrack->phi()+recoCombinedGlbTrack->phi());
    phiResolution[2]->Fill(recoGlbTrack->phi()-recoStaGlbTrack->phi());
    phiResolution[3]->Fill(recoCombinedGlbTrack->phi(), recoGlbTrack->phi()-recoCombinedGlbTrack->phi());
    phiResolution[4]->Fill(recoCombinedGlbTrack->phi(), -recoStaGlbTrack->phi()+recoCombinedGlbTrack->phi());
    phiResolution[5]->Fill(recoCombinedGlbTrack->phi(), recoGlbTrack->phi()-recoStaGlbTrack->phi());
    
    pGlbTrack[0]->Fill(recoCombinedGlbTrack->p());
    pGlbTrack[1]->Fill(recoGlbTrack->p());
    pGlbTrack[2]->Fill(recoStaGlbTrack->p());

    ptGlbTrack[0]->Fill(recoCombinedGlbTrack->pt());
    ptGlbTrack[1]->Fill(recoGlbTrack->pt());
    ptGlbTrack[2]->Fill(recoStaGlbTrack->pt());

    qGlbTrack[0]->Fill(recoCombinedGlbTrack->charge());
    qGlbTrack[1]->Fill(recoGlbTrack->charge());
    qGlbTrack[2]->Fill(recoStaGlbTrack->charge());
    if(recoCombinedGlbTrack->charge()==recoStaGlbTrack->charge()) qGlbTrack[3]->Fill(1);
    else qGlbTrack[3]->Fill(2);
    if(recoCombinedGlbTrack->charge()==recoGlbTrack->charge()) qGlbTrack[3]->Fill(3);
    else qGlbTrack[3]->Fill(4);
    if(recoStaGlbTrack->charge()==recoGlbTrack->charge()) qGlbTrack[3]->Fill(5);
    else qGlbTrack[3]->Fill(6);
    if(recoCombinedGlbTrack->charge()!=recoStaGlbTrack->charge() && recoCombinedGlbTrack->charge()!=recoGlbTrack->charge()) qGlbTrack[3]->Fill(7);
    if(recoCombinedGlbTrack->charge()==recoStaGlbTrack->charge() && recoCombinedGlbTrack->charge()==recoGlbTrack->charge()) qGlbTrack[3]->Fill(8);
    
    qOverpResolution[0]->Fill((recoGlbTrack->charge()/recoGlbTrack->p())-(recoCombinedGlbTrack->charge()/recoCombinedGlbTrack->p()));
    qOverpResolution[1]->Fill(-(recoStaGlbTrack->charge()/recoStaGlbTrack->p())+(recoCombinedGlbTrack->charge()/recoCombinedGlbTrack->p()));
    qOverpResolution[2]->Fill((recoGlbTrack->charge()/recoGlbTrack->p())-(recoStaGlbTrack->charge()/recoStaGlbTrack->p()));
    oneOverpResolution[0]->Fill((1/recoGlbTrack->p())-(1/recoCombinedGlbTrack->p()));
    oneOverpResolution[1]->Fill(-(1/recoStaGlbTrack->p())+(1/recoCombinedGlbTrack->p()));
    oneOverpResolution[2]->Fill((1/recoGlbTrack->p())-(1/recoStaGlbTrack->p()));
    qOverptResolution[0]->Fill((recoGlbTrack->charge()/recoGlbTrack->pt())-(recoCombinedGlbTrack->charge()/recoCombinedGlbTrack->pt()));
    qOverptResolution[1]->Fill(-(recoStaGlbTrack->charge()/recoStaGlbTrack->pt())+(recoCombinedGlbTrack->charge()/recoCombinedGlbTrack->pt()));
    qOverptResolution[2]->Fill((recoGlbTrack->charge()/recoGlbTrack->pt())-(recoStaGlbTrack->charge()/recoStaGlbTrack->pt()));
    oneOverptResolution[0]->Fill((1/recoGlbTrack->pt())-(1/recoCombinedGlbTrack->pt()));
    oneOverptResolution[1]->Fill(-(1/recoStaGlbTrack->pt())+(1/recoCombinedGlbTrack->pt()));
    oneOverptResolution[2]->Fill((1/recoGlbTrack->pt())-(1/recoStaGlbTrack->pt()));
    oneOverptResolution[3]->Fill(recoCombinedGlbTrack->eta(),(1/recoGlbTrack->pt())-(1/recoCombinedGlbTrack->pt()));
    oneOverptResolution[4]->Fill(recoCombinedGlbTrack->eta(),-(1/recoStaGlbTrack->pt())+(1/recoCombinedGlbTrack->pt()));
    oneOverptResolution[5]->Fill(recoCombinedGlbTrack->eta(),(1/recoGlbTrack->pt())-(1/recoStaGlbTrack->pt()));
    oneOverptResolution[6]->Fill(recoCombinedGlbTrack->phi(),(1/recoGlbTrack->pt())-(1/recoCombinedGlbTrack->pt()));
    oneOverptResolution[7]->Fill(recoCombinedGlbTrack->phi(),-(1/recoStaGlbTrack->pt())+(1/recoCombinedGlbTrack->pt()));
    oneOverptResolution[8]->Fill(recoCombinedGlbTrack->phi(),(1/recoGlbTrack->pt())-(1/recoStaGlbTrack->pt()));
    oneOverptResolution[9]->Fill(recoCombinedGlbTrack->pt(),(1/recoGlbTrack->pt())-(1/recoCombinedGlbTrack->pt()));
    oneOverptResolution[10]->Fill(recoCombinedGlbTrack->pt(),-(1/recoStaGlbTrack->pt())+(1/recoCombinedGlbTrack->pt()));
    oneOverptResolution[11]->Fill(recoCombinedGlbTrack->pt(),(1/recoGlbTrack->pt())-(1/recoStaGlbTrack->pt()));
    
    // valid hits Glb track
    double rhGlb = recoCombinedGlbTrack->found();
    // invalid hits Glb track
    double rhGlb_notValid = recoCombinedGlbTrack->recHitsSize()-recoCombinedGlbTrack->found();
    // valid hits Glb track from Tracker
    double rhGlb_StaProvenance=0;
    // valid hits Glb track from Sta system
    double rhGlb_TkProvenance=0;
    for (trackingRecHit_iterator recHit = recoCombinedGlbTrack->recHitsBegin();
	 recHit!=recoCombinedGlbTrack->recHitsEnd(); ++recHit){
      if((*recHit)->isValid()){
	DetId id = (*recHit)->geographicalId();
	if (id.det() == DetId::Muon)
	  rhGlb_StaProvenance++;
	if (id.det() == DetId::Tracker)
	  rhGlb_TkProvenance++;
      }
    }
    // valid hits Sta track associated to Glb track
    double rhStaGlb = recoStaGlbTrack->found();
    // valid hits Traker track associated to Glb track
    double rhTkGlb = recoGlbTrack->found();
    // invalid hits Traker track associated to Glb track
    double rhTkGlb_notValid = recoGlbTrack->recHitsSize()-recoGlbTrack->found();
   
    // fill the histos
    rhAnalysis[0]->Fill(rhGlb_StaProvenance/rhGlb);
    rhAnalysis[1]->Fill(rhGlb_TkProvenance/rhGlb);
    rhAnalysis[2]->Fill(rhGlb_StaProvenance/rhStaGlb);
    rhAnalysis[3]->Fill(rhGlb_TkProvenance/rhTkGlb);
    rhAnalysis[4]->Fill(rhGlb/(rhStaGlb+rhTkGlb));
    if(rhGlb_notValid!=0)
      rhAnalysis[5]->Fill(rhGlb_notValid/rhGlb);
    else{
      if(rhTkGlb_notValid!=0)
	rhAnalysis[5]->Fill(rhTkGlb_notValid/rhGlb);
    }
    
  }


  if(recoMu.isTrackerMuon() && !(recoMu.isGlobalMuon())) {
    LogTrace(metname)<<"[MuonRecoAnalyzer] The mu is tracker only - filling the histos";
     if(recoMu.isStandAloneMuon())
      muReco->Fill(3);
    if(!(recoMu.isStandAloneMuon()))
      muReco->Fill(4);

    // get the track using only the tracker data
    reco::TrackRef recoTrack = recoMu.track();

    etaTrack->Fill(recoTrack->eta());
    thetaTrack->Fill(recoTrack->theta());
    phiTrack->Fill(recoTrack->phi());
    pTrack->Fill(recoTrack->p());
    ptTrack->Fill(recoTrack->pt());
    qTrack->Fill(recoTrack->charge());
    
  }

  if(recoMu.isStandAloneMuon() && !(recoMu.isGlobalMuon())) {
    LogTrace(metname)<<"[MuonRecoAnalyzer] The mu is STA only - filling the histos";
    if(!(recoMu.isTrackerMuon()))
      muReco->Fill(5);
     
    // get the track using only the mu spectrometer data
    reco::TrackRef recoStaTrack = recoMu.standAloneMuon();

    etaStaTrack->Fill(recoStaTrack->eta());
    thetaStaTrack->Fill(recoStaTrack->theta());
    phiStaTrack->Fill(recoStaTrack->phi());
    pStaTrack->Fill(recoStaTrack->p());
    ptStaTrack->Fill(recoStaTrack->pt());
    qStaTrack->Fill(recoStaTrack->charge());

  }
    
  if(recoMu.isCaloMuon() && !(recoMu.isGlobalMuon()) && !(recoMu.isTrackerMuon()) && !(recoMu.isStandAloneMuon()))
    muReco->Fill(6);
  
  //efficiency plots
  
  // get the track using only the mu spectrometer data
  reco::TrackRef recoStaGlbTrack = recoMu.standAloneMuon();
  
  if(recoMu.isStandAloneMuon()){
    etaEfficiency[0]->Fill(recoStaGlbTrack->eta());
    phiEfficiency[0]->Fill(recoStaGlbTrack->phi());
  }
  if(recoMu.isStandAloneMuon() && recoMu.isGlobalMuon()){
    etaEfficiency[1]->Fill(recoStaGlbTrack->eta());
    phiEfficiency[1]->Fill(recoStaGlbTrack->phi());
  }




}
