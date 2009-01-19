
// -*- C++ -*-
//
// Package:    SusyDiJetAnalysis
// Class:      SusyDiJetAnalysis
// 
/**\class SusyDiJetAnalysis SusyDiJetAnalysis.cc SusyAnalysis/AnalysisSkeleton/src/SusyDiJetAnalysis.cc

Description: Skeleton analysis for SUSY search with Jets + MET

Implementation:Uses the EventSelector interface for event selection and TFileService for plotting.

*/
//
// Original Author:  Markus Stoye
//         Created:  Mon Feb 18 15:40:44 CET 2008
// $Id: SusyDiJetAnalysis.cpp,v 1.9 2009/01/06 09:51:35 trommers Exp $
//
//
//#include "SusyAnalysis/EventSelector/interface/BJetEventSelector.h"
#include "SusyAnalysis/AnalysisSkeleton/test/SusyDiJetAnalysis.h"

using namespace std;
using namespace reco;
using namespace edm;

//________________________________________________________________________________________
SusyDiJetAnalysis::SusyDiJetAnalysis(const edm::ParameterSet& iConfig):
  sequence_( iConfig.getParameter<edm::ParameterSet>("selections") ),
  plotSelection_( iConfig.getParameter<std::vector<std::string> >("plotSelection") ),
  eventWeight_( iConfig.getParameter<double>("eventWeight") ),
  nrEventTotalRaw_(0), nrEventTotalWeighted_(0.0),
  pathNames_(0), nEvents_(0), nWasRun_(0), nAccept_(0), nErrors_(0), hlWasRun_(0), hlAccept_(0), hlErrors_(0), init_(false), //georgia
  genTag_(iConfig.getParameter<edm::InputTag>("genTag"))
{ 
  // Translate plotSelection strings to indices
   plotSelectionIndices_.reserve(plotSelection_.size());
  for ( size_t i=0; i<plotSelection_.size(); ++i )  plotSelectionIndices_.push_back(sequence_.selectorIndex(plotSelection_[i]));
   

  // List all selectors and selection variables
  edm::LogVerbatim("SusyDiJet") << "Selectors are:" << std::endl;
 
  for ( std::vector<const SusyEventSelector*>::const_iterator it = sequence_.selectors().begin();
        it != sequence_.selectors().end(); ++it )
    {  edm::LogVerbatim("SusyDiJet") << " * " << (*it)->name()
                                          << " selects on following " 
                                          << (*it)->numberOfVariables() << " variable(s):";
      for ( unsigned int i=0; i<(*it)->numberOfVariables(); ++i ){
	edm::LogVerbatim("SusyDiJet") << "    - " << (*it)->variableNames()[i];
	}
      edm::LogVerbatim("SusyDiJet") << std::endl;
    }

  mSelectorResults = new unsigned int[sequence_.size()];
  
  // Say something about event weights
 
    edm::LogInfo("SusyDiJet") << "Global event weight set to " << eventWeight_;

  // get the data tags
  jetTag_ = iConfig.getParameter<edm::InputTag>("jetTag");
  metTag_ = iConfig.getParameter<edm::InputTag>("metTag");
  photTag_ = iConfig.getParameter<edm::InputTag>("photTag");
  elecTag_ = iConfig.getParameter<edm::InputTag>("elecTag");
  muonTag_ = iConfig.getParameter<edm::InputTag>("muonTag");
  tauTag_ = iConfig.getParameter<edm::InputTag>("tauTag");
  vtxTag_ = iConfig.getParameter<edm::InputTag>("vtxTag"); 
  ccjetTag_ = iConfig.getParameter<edm::InputTag>("ccjetTag");
  ccmetTag_ = iConfig.getParameter<edm::InputTag>("ccmetTag");
  ccelecTag_ = iConfig.getParameter<edm::InputTag>("ccelecTag"); 
  ccmuonTag_ = iConfig.getParameter<edm::InputTag>("ccmuonTag");

  // trigger stuff
  triggerResults_ = iConfig.getParameter<edm::InputTag>("triggerResults");
  // trigger path names
  pathNames_ = iConfig.getParameter< std::vector<std::string> >("pathNames");

  
  // Initialise counters
  nrEventSelected_.resize( sequence_.size(), 0.0 );
  nrEventAllButOne_.resize( sequence_.size(), 0.0 );
  nrEventCumulative_.resize( sequence_.size(), 0.0 );

  localPi = acos(-1.0);


  // Initialise plots [should improve in the future]
    initPlots();

}


//________________________________________________________________________________________
SusyDiJetAnalysis::~SusyDiJetAnalysis() {}


//filter---------------------------------------------------------
bool
SusyDiJetAnalysis::filter(const edm::Event& iEvent,const edm::EventSetup& iSetup ){

  run_   = iEvent.id().run();
  event_ = iEvent.id().event();
  processId_   = 0;

 // Retrieve the decision of each selector module
  SelectorDecisions decisions = sequence_.decisions(iEvent);

 // Count all events
  nrEventTotalRaw_++;
  nrEventTotalWeighted_ += eventWeight_;

   // Fill plots with all variables
  bool dec(true);
  for ( size_t i=0; i<sequence_.size(); ++i ) {
    dec = dec && decisions.decision(i);
 
    // Add the decision to the tree
    mSelectorResults[i] = (decisions.decision(i)?1:0);
    
    // Update counters
    if ( decisions.decision(i) ) nrEventSelected_[i] += eventWeight_;
    if ( decisions.complementaryDecision(i) ) nrEventAllButOne_[i] += eventWeight_;
    if ( decisions.cumulativeDecision(i) ) nrEventCumulative_[i] += eventWeight_;
    
  }

  // Fill event with variables we computed
  if(dec)fillPlots( iEvent, decisions );

  // Print summary so far (every 10 till 100, every 100 till 1000, etc.)
  for ( unsigned int i=10; i<nrEventTotalRaw_; i*=10 )
    if ( nrEventTotalRaw_<=10*i && (nrEventTotalRaw_%i)==0 )
      printSummary();

  return dec;
}



//________________________________________________________________________________________
// Method called to for each event
void
SusyDiJetAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  
  //std::cout<< " in analyze " << std::endl;

   edm::LogVerbatim("SusyDiJetEvent") << " Start  " << std::endl;

  std::ostringstream dbg;

  ALPGENParticleId myALPGENParticleId;
 

  // markus
  ///
  if ( !filter(iEvent,iSetup)) {
    // Just fill the failure
    mGlobalDecision = 0;
    mSelectorData->Fill();
    return;
  }
  

  // Just fill the success
  mGlobalDecision = 1;
  mSelectorData->Fill();
  
 edm::LogVerbatim("SusyDiJetEvent") << " Trigger decision  " << std::endl;
 // std::cout << " trigger decision " << std::endl;

  //get the trigger decision
  mTempTreeHLT1JET=false;
  mTempTreeHLT2JET=false;
  mTempTreeHLT1MET1HT=false;

  // Get the trigger results and check validity
  edm::Handle<edm::TriggerResults> hltHandle;
  iEvent.getByLabel(triggerResults_, hltHandle);
  if ( !hltHandle.isValid() ) {
    edm::LogWarning("HLTEventSelector") << "No trigger results for InputTag " << triggerResults_;
    return;
  }

  //  std::cout << " get results " << std::endl;

  // Get results
  edm::TriggerNames trgNames;
  trgNames.init(*hltHandle);
  unsigned int trgSize = trgNames.size();

  
  // Example for OR of all specified triggers

  edm::LogWarning("HLTEventSelector") << " triggers " << trgNames.size() << std::endl;

  /*   for (unsigned int itrig = 0; itrig != trgSize; ++itrig) {
        TString trigName = trgNames.triggerName(itrig);
	  	std::cout << " trigName " << trigName <<  " accept " << hltHandle->accept(trgNames.triggerIndex(trgNames.triggerName(itrig))) <<std::endl;
	edm::LogWarning("HLTEventSelector") << " trigName " << trigName << " accept " << hltHandle->accept(trgNames.triggerIndex(trgNames.triggerName(itrig))) << std::endl;
	}*/
  
  //  std::cout << " accept " << hltHandle->size()<<std::endl;
  
  
    // GEORGIA
  if (!hltHandle.isValid()) {
    // triggerExists = false;
    std::cout << "HLTriggerResult Not Valid!" << endl;
  }
  else {  
    if (hltHandle->wasrun()) nWasRun_++;
    const bool accept(hltHandle->accept());
    LogDebug("") << "HL TriggerResults decision: " << accept;
    if (accept) ++nAccept_;
    if (hltHandle->error() ) nErrors_++;
  }
  if (!init_) {
    init_=true;
    triggerNames_.init(*hltHandle);
    pathNames_=triggerNames_.triggerNames();
    const unsigned int n(pathNames_.size());
    hlWasRun_.resize(n);
    hlAccept_.resize(n);
    hlErrors_.resize(n);
    for (unsigned int i=0; i!=n; ++i) {
      hlWasRun_[i]=0;
      hlAccept_[i]=0;
      hlErrors_[i]=0;
    }
  }

  // decision for each HL algorithm
  const unsigned int n(pathNames_.size());
  for (unsigned int i=0; i!=n; ++i) {
    if (hltHandle->wasrun(i)) hlWasRun_[i]++;
    if (hltHandle->accept(i)) hlAccept_[i]++;
    if (hltHandle->error(i) ) hlErrors_[i]++;
  }
  
  nHLT=static_cast<int>(n);
  for(unsigned int i=0; i!=n; ++i) {
    HLTArray[i]=hltHandle->accept(i);
  }



  //looping over list of trig path names
  for ( std::vector<std::string>::const_iterator i=pathNames_.begin();
	i!=pathNames_.end(); ++i ) {
    // Get index
 
    //std::cout << " accept " << hltHandle->accept(trgNames.triggerIndex(*i)) <<  std::endl;

    unsigned int index = trgNames.triggerIndex(*i);
    if ( index==trgSize ) {
      edm::LogWarning("HLTEventSelector") << "Unknown trigger name " << *i;
      continue;
    }
    if ( hltHandle->accept(index) ) {
      LogDebug("HLTEventSelector") << "Event selected by " << *i ;
      std::string trigName = *i;
      if (trigName == "HLT_Jet180") mTempTreeHLT1JET=true;
      if (trigName == "HLT_DiJetAve130") mTempTreeHLT2JET=true;
      if (trigName == "HLT_MET50") mTempTreeHLT1MET1HT=true;
      if (trigName == "HLT_Mu9") mTempTreeHLT1Muon=true; 
      
    } 
  }

  //std::cout << " after
 
  mTempTreeEventWeight =eventWeight_;
 

  mTempTreeRun = run_;
  mTempTreeEvent = event_;

 

  // GEORGIA 
  // get the Vertex collection
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vtxTag_, vertices);
  if ( !vertices.isValid() ) {
    LogDebug("") << "No Vertex results for InputTag" << vtxTag_;
    return;
  } 

  // Should assume that 1st index in VertexCollection corresponds to the primary vertex ?? ( verices are ordered?)

  mTempTreenVtx = (*vertices).size();

  for (int i=0; i< (int) (*vertices).size(); i++){  
    //    int indPrim=0;
    const reco::Vertex* pVertex = &(*vertices)[i];
    mTempTreeVtxChi2[i] = pVertex->chi2();
    mTempTreeVtxNdof[i] = pVertex->ndof();
    mTempTreeVtxNormalizedChi2[i] = pVertex->normalizedChi2();
    mTempTreeVtxX[i] = pVertex->x();
    mTempTreeVtxY[i] = pVertex->y();
    mTempTreeVtxZ[i] = pVertex->z();
    mTempTreeVtxdX[i] = pVertex->xError();
    mTempTreeVtxdY[i] = pVertex->yError();
    mTempTreeVtxdZ[i] = pVertex->zError();
  } 
  // end GEORGIA

  
  // get the photons
  edm::Handle< std::vector<pat::Photon> > photHandle;
  iEvent.getByLabel(photTag_, photHandle);
  if ( !photHandle.isValid() ) {
    edm::LogWarning("SusySelectorExample") << "No Photon results for InputTag " << photTag_;
    return;
  }


  
  // Add the photons
  mTempTreeNphot = photHandle->size();
  if ( mTempTreeNphot > 50 ) mTempTreeNphot = 50;
  for (int i=0;i<mTempTreeNphot;i++){
    mTempTreePhotPt[i]  = (*photHandle)[i].pt();
    mTempTreePhotE[i]   = (*photHandle)[i].energy();
    mTempTreePhotEt[i]  = (*photHandle)[i].et();
    mTempTreePhotPx[i]  = (*photHandle)[i].momentum().X();
    mTempTreePhotPy[i]  = (*photHandle)[i].momentum().Y();
    mTempTreePhotPz[i]  = (*photHandle)[i].momentum().Z();
    mTempTreePhotEta[i] = (*photHandle)[i].eta();
    mTempTreePhotPhi[i] = (*photHandle)[i].phi();
    mTempTreePhotTrkIso[i] = (*photHandle)[i].trackIso();
    mTempTreePhotECalIso[i] = (*photHandle)[i].ecalIso();
    mTempTreePhotHCalIso[i] = (*photHandle)[i].hcalIso();
    mTempTreePhotAllIso[i] = (*photHandle)[i].caloIso();
  }

  
  
  // get the electrons
  edm::Handle< std::vector<pat::Electron> > elecHandle;
  iEvent.getByLabel(elecTag_, elecHandle);
  if ( !elecHandle.isValid() ) {
    edm::LogWarning("SusySelectorExample") << "No Electron results for InputTag " << elecTag_;
    return;
  }
  
  
  // Add the electrons
  mTempTreeNelec= elecHandle->size();
  if ( mTempTreeNelec > 50 ) mTempTreeNelec = 50;
  for (int i=0;i<mTempTreeNelec;i++){
    mTempTreeElecPt[i] = (*elecHandle)[i].pt();
    mTempTreeElecE[i] = (*elecHandle)[i].energy();
    mTempTreeElecEt[i] = (*elecHandle)[i].et();
    mTempTreeElecPx[i] = (*elecHandle)[i].momentum().X();
    mTempTreeElecPy[i] = (*elecHandle)[i].momentum().Y();
    mTempTreeElecPz[i] = (*elecHandle)[i].momentum().Z();
    mTempTreeElecEta[i] = (*elecHandle)[i].eta();
    mTempTreeElecPhi[i] = (*elecHandle)[i].phi();
    mTempTreeElecTrkIso[i] = ( (*elecHandle)[i].trackIso()+(*elecHandle)[i].et() )/(*elecHandle)[i].et();

    mTempTreeElecECalIso[i] = (*elecHandle)[i].ecalIso();
    mTempTreeElecHCalIso[i] = (*elecHandle)[i].hcalIso() ;
    mTempTreeElecAllIso[i] = (*elecHandle)[i].caloIso() ;
    mTempTreeElecCharge[i] = (*elecHandle)[i].charge();

    //MICHELE
    mTempTreeElecIdLoose[i] = (*elecHandle)[i].electronID("eidLoose");
    mTempTreeElecIdTight[i] = (*elecHandle)[i].electronID("eidTight");
    mTempTreeElecIdRobLoose[i] = (*elecHandle)[i].electronID("eidRobustLoose");
    mTempTreeElecIdRobTight[i] = (*elecHandle)[i].electronID("eidRobustTight"); 



    mTempTreeElecCaloEnergy[i] = (*elecHandle)[i].caloEnergy();
    mTempTreeElecHOverE[i] = (*elecHandle)[i].hadronicOverEm();
    mTempTreeElecVx[i] = (*elecHandle)[i].vx();
    mTempTreeElecVy[i] = (*elecHandle)[i].vy();
    mTempTreeElecVz[i] = (*elecHandle)[i].vz();
    mTempTreeElecD0[i] = (*elecHandle)[i].gsfTrack()->d0();
    mTempTreeElecDz[i] = (*elecHandle)[i].gsfTrack()->dz();
    mTempTreeElecChargeMode[i] = (*elecHandle)[i].gsfTrack()->chargeMode();	
    mTempTreeElecPtTrkMode[i] = (*elecHandle)[i].gsfTrack()->ptMode();
    mTempTreeElecQOverPErrTrkMode[i] = (*elecHandle)[i].gsfTrack()->qoverpModeError();
    mTempTreeElecCharge[i] = (*elecHandle)[i].gsfTrack()->charge();
    mTempTreeElecPtTrk[i] = (*elecHandle)[i].gsfTrack()->pt();
    mTempTreeElecQOverPErrTrk[i] = (*elecHandle)[i].gsfTrack()->qoverpError();
    mTempTreeElecLostHits[i] = (*elecHandle)[i].gsfTrack()->lost();
    mTempTreeElecValidHits[i] = (*elecHandle)[i].gsfTrack()->found();
    mTempTreeElecNCluster[i] = (*elecHandle)[i].numberOfClusters();
    mTempTreeElecEtaTrk[i] = (*elecHandle)[i].trackMomentumAtVtx().Eta();
    mTempTreeElecPhiTrk[i] = (*elecHandle)[i].trackMomentumAtVtx().Phi();
    mTempTreeElecWidthClusterEta[i] = (*elecHandle)[i].superCluster()->etaWidth();
    mTempTreeElecWidthClusterPhi[i] = (*elecHandle)[i].superCluster()->phiWidth();
    mTempTreeElecPinTrk[i] = sqrt((*elecHandle)[i].trackMomentumAtVtx().Mag2());
    mTempTreeElecPoutTrk[i] = sqrt((*elecHandle)[i].trackMomentumOut().Mag2());

    if (&(*(*elecHandle)[i].genLepton())!=0){
      mTempTreeGenElecPdgId[i] = (*elecHandle)[i].genLepton()->pdgId();
      mTempTreeGenElecPx[i] = (*elecHandle)[i].genLepton()->px();
      mTempTreeGenElecPy[i] = (*elecHandle)[i].genLepton()->py();
      mTempTreeGenElecPz[i] = (*elecHandle)[i].genLepton()->pz();
      if(&(*(*elecHandle)[i].genLepton()->mother())!=0){
	mTempTreeGenElecMother[i]= (*elecHandle)[i].genLepton()->mother()->pdgId();
      }
    }
    else {
      mTempTreeGenElecPdgId[i] = 999.;
      mTempTreeGenElecPx[i]=999.;
      mTempTreeGenElecPy[i]=999.;
      mTempTreeGenElecPz[i]=999.;
      mTempTreeGenElecMother[i] = 999.;
    }
    //PIOPPI
  }


   // get the ccelectrons
  edm::Handle< std::vector<pat::Electron> > ccelecHandle;
  iEvent.getByLabel(ccelecTag_, ccelecHandle);
  if ( !ccelecHandle.isValid() ) {
    edm::LogWarning("SusySelectorExample") << "No ccElectron results for InputTag " << ccelecTag_;
    return;
  }
  //std::cout << " add the ccelectrons " << std::endl;
  
  // Add the electrons
  mTempTreeNccelec= ccelecHandle->size();
  if ( mTempTreeNccelec > 50 ) mTempTreeNccelec = 50;
  for (int i=0;i<mTempTreeNccelec;i++){
    mTempTreeccElecPt[i] = (*ccelecHandle)[i].pt();
    mTempTreeccElecE[i] = (*ccelecHandle)[i].energy();
    mTempTreeccElecEt[i] = (*ccelecHandle)[i].et();
    mTempTreeccElecPx[i] = (*ccelecHandle)[i].momentum().X();
    mTempTreeccElecPy[i] = (*ccelecHandle)[i].momentum().Y();
    mTempTreeccElecPz[i] = (*ccelecHandle)[i].momentum().Z();
    mTempTreeccElecEta[i] = (*ccelecHandle)[i].eta();
    mTempTreeccElecPhi[i] = (*ccelecHandle)[i].phi();

  }
  
  // get the muons
  edm::Handle< std::vector<pat::Muon> > muonHandle;
  iEvent.getByLabel(muonTag_, muonHandle);
  if ( !muonHandle.isValid() ) {
    edm::LogWarning("SusySelectorExample") << "No Muon results for InputTag " << muonTag_;
    return;
  }
  

  // Add the muons
  mTempTreeNmuon= muonHandle->size();
  if ( mTempTreeNmuon > 50 ) mTempTreeNmuon = 50;
  for (int i=0;i<mTempTreeNmuon;i++){
   
    mTempTreeMuonPt[i] = (*muonHandle)[i].pt();
    mTempTreeMuonE[i] = (*muonHandle)[i].energy();
    mTempTreeMuonEt[i] = (*muonHandle)[i].et();
    mTempTreeMuonPx[i] = (*muonHandle)[i].momentum().X();
    mTempTreeMuonPy[i] = (*muonHandle)[i].momentum().Y();
    mTempTreeMuonPz[i] = (*muonHandle)[i].momentum().Z();
    mTempTreeMuonEta[i] = (*muonHandle)[i].eta();
    mTempTreeMuonPhi[i] = (*muonHandle)[i].phi();
    mTempTreeMuonTrkIso[i] = ( (*muonHandle)[i].trackIso()+(*muonHandle)[i].et() )/(*muonHandle)[i].et();
    mTempTreeMuonCharge[i] = (*muonHandle)[i].charge();
    mTempTreeMuonECalIso[i] = (*muonHandle)[i].ecalIso();
    mTempTreeMuonHCalIso[i] = (*muonHandle)[i].hcalIso() ;
    mTempTreeMuonAllIso[i] = (*muonHandle)[i].caloIso() ;
    mTempTreeMuonIsGlobal[i] = (*muonHandle)[i].isGlobalMuon();
    mTempTreeMuonIsStandAlone[i] = (*muonHandle)[i].isStandAloneMuon();
    mTempTreeMuonIsTracker[i] = (*muonHandle)[i].isTrackerMuon();
    mTempTreeMuonIsGlobalTight[i] = (*muonHandle)[i].isGood(pat::Muon::SelectionType(6));
    mTempTreeMuonIsTMLastStationLoose[i] = (*muonHandle)[i].isGood(pat::Muon::SelectionType(7));

    
    //MICHELE    
    //    mTempTreeMuonId[i]=2.; //FIX NEEDED

    mTempTreeMuonCombVx[i]=(*muonHandle)[i].globalTrack()->vx();
    mTempTreeMuonCombVy[i]=(*muonHandle)[i].globalTrack()->vy();
    mTempTreeMuonCombVz[i]=(*muonHandle)[i].globalTrack()->vz();
    mTempTreeMuonCombD0[i]=(*muonHandle)[i].globalTrack()->d0();
    mTempTreeMuonCombDz[i]=(*muonHandle)[i].globalTrack()->dz();
    if((*muonHandle)[i].isStandAloneMuon() && (*muonHandle)[i].standAloneMuon().isNonnull()){
      mTempTreeMuonStandValidHits[i]=(*muonHandle)[i].standAloneMuon()->found();
      mTempTreeMuonStandLostHits[i]=(*muonHandle)[i].standAloneMuon()->lost();
      mTempTreeMuonStandPt[i]=(*muonHandle)[i].standAloneMuon()->pt();
      mTempTreeMuonStandPz[i]=(*muonHandle)[i].standAloneMuon()->pz();
      mTempTreeMuonStandP[i]=(*muonHandle)[i].standAloneMuon()->p();
      mTempTreeMuonStandEta[i]=(*muonHandle)[i].standAloneMuon()->eta();
      mTempTreeMuonStandPhi[i]=(*muonHandle)[i].standAloneMuon()->phi();
      mTempTreeMuonStandChi[i]=(*muonHandle)[i].standAloneMuon()->chi2();
      mTempTreeMuonStandCharge[i]=(*muonHandle)[i].standAloneMuon()->charge();
      mTempTreeMuonStandQOverPError[i]=(*muonHandle)[i].standAloneMuon()->qoverpError();
    } 
    else{
      mTempTreeMuonStandValidHits[i]=999.;
      mTempTreeMuonStandLostHits[i]=999.;
      mTempTreeMuonStandPt[i]=999.;
      mTempTreeMuonStandPz[i]=999.;
      mTempTreeMuonStandP[i]=999.;
      mTempTreeMuonStandEta[i]=999.;
      mTempTreeMuonStandPhi[i]=999.;
      mTempTreeMuonStandChi[i]=999.;
      mTempTreeMuonStandCharge[i]=999.;
      mTempTreeMuonStandQOverPError[i]=999.;
    }
    if((*muonHandle)[i].isTrackerMuon() && (*muonHandle)[i].track().isNonnull()){
      mTempTreeMuonTrkChiNorm[i] = (*muonHandle)[i].track()->normalizedChi2();
      mTempTreeMuonTrkValidHits[i]=(*muonHandle)[i].track()->found();
      mTempTreeMuonTrkLostHits[i]=(*muonHandle)[i].track()->lost();
      mTempTreeMuonTrkPt[i]=(*muonHandle)[i].track()->pt();
      mTempTreeMuonTrkPz[i]=(*muonHandle)[i].track()->pz();
      mTempTreeMuonTrkP[i]=(*muonHandle)[i].track()->p();
      mTempTreeMuonTrkEta[i]=(*muonHandle)[i].track()->eta();
      mTempTreeMuonTrkPhi[i]=(*muonHandle)[i].track()->phi();
      mTempTreeMuonTrkChi[i]=(*muonHandle)[i].track()->chi2();
      mTempTreeMuonTrkCharge[i]=(*muonHandle)[i].track()->charge();
      mTempTreeMuonTrkQOverPError[i]=(*muonHandle)[i].track()->qoverpError();

    }
    else{
      mTempTreeMuonTrkChiNorm[i] = 999.;
      mTempTreeMuonTrkValidHits[i]=999.;
      mTempTreeMuonTrkLostHits[i]=999.;
      mTempTreeMuonTrkPt[i]=999.;
      mTempTreeMuonTrkPz[i]=999.;
      mTempTreeMuonTrkP[i]=999.;
      mTempTreeMuonTrkEta[i]=999.;
      mTempTreeMuonTrkPhi[i]=999.;
      mTempTreeMuonTrkChi[i]=999.;
      mTempTreeMuonTrkCharge[i]=999.;
      mTempTreeMuonTrkQOverPError[i]=999.;
    }
    //PIOPPI
    if (&(*(*muonHandle)[i].genLepton())!=0){
      mTempTreeGenMuonPdgId[i]=(*muonHandle)[i].genLepton()->pdgId();
      mTempTreeGenMuonPx[i]=(*muonHandle)[i].genLepton()->px();
      mTempTreeGenMuonPy[i]=(*muonHandle)[i].genLepton()->py();
      mTempTreeGenMuonPz[i]=(*muonHandle)[i].genLepton()->pz();
      if (&(*(*muonHandle)[i].genLepton()->mother())!=0)
	mTempTreeGenMuonMother[i]=(*muonHandle)[i].genLepton()->mother()->pdgId();
      else mTempTreeGenMuonMother[i]=999.;
    }
    else{
      mTempTreeGenMuonPdgId[i]=999.;
      mTempTreeGenMuonMother[i]=999.;
      mTempTreeGenMuonPx[i]=999.;
      mTempTreeGenMuonPy[i]=999.;
      mTempTreeGenMuonPz[i]=999.;
    }

  }
  //std::cout << " add the ccmuons " << std::endl;
  

  // get the ccmuons
  edm::Handle< std::vector<pat::Muon> > ccmuonHandle;
  iEvent.getByLabel(ccmuonTag_, ccmuonHandle);
  if ( !ccmuonHandle.isValid() ) {
    edm::LogWarning("SusySelectorExample") << "No ccMuon results for InputTag " << ccmuonTag_;
    return;
  }
  

  // Add the ccmuons
  mTempTreeNccmuon= ccmuonHandle->size();
  if ( mTempTreeNccmuon > 50 ) mTempTreeNccmuon = 50;
  for (int i=0;i<mTempTreeNccmuon;i++){
    mTempTreeccMuonPt[i] = (*ccmuonHandle)[i].pt();
    mTempTreeccMuonE[i] = (*ccmuonHandle)[i].energy();
    mTempTreeccMuonEt[i] = (*ccmuonHandle)[i].et();
    mTempTreeccMuonPx[i] = (*ccmuonHandle)[i].momentum().X();
    mTempTreeccMuonPy[i] = (*ccmuonHandle)[i].momentum().Y();
    mTempTreeccMuonPz[i] = (*ccmuonHandle)[i].momentum().Z();
    mTempTreeccMuonEta[i] = (*ccmuonHandle)[i].eta();
    mTempTreeccMuonPhi[i] = (*ccmuonHandle)[i].phi();
 
  }

  // get the taus
  edm::Handle< std::vector<pat::Tau> > tauHandle;
  iEvent.getByLabel(tauTag_, tauHandle);
  if ( !tauHandle.isValid() ) {
    edm::LogWarning("SusySelectorExample") << "No Tau results for InputTag " << tauTag_;
    return;
  }
 

  // Add the taus
  mTempTreeNtau= tauHandle->size();
  if ( mTempTreeNtau > 50 ) mTempTreeNtau = 50;
  for (int i=0;i<mTempTreeNtau;i++){
    mTempTreeTauPt[i] = (*tauHandle)[i].pt();
    mTempTreeTauE[i] = (*tauHandle)[i].energy();
    mTempTreeTauEt[i] = (*tauHandle)[i].et();
    mTempTreeTauPx[i] = (*tauHandle)[i].momentum().X();
    mTempTreeTauPy[i] = (*tauHandle)[i].momentum().Y();
    mTempTreeTauPz[i] = (*tauHandle)[i].momentum().Z();
    mTempTreeTauEta[i] = (*tauHandle)[i].eta();
    mTempTreeTauPhi[i] = (*tauHandle)[i].phi();
    mTempTreeTauTrkIso[i] = ( (*tauHandle)[i].trackIso()+(*tauHandle)[i].et() )/(*tauHandle)[i].et();

 edm::LogVerbatim("SusyDiJetEvent") << "Taus " << i << " iso " <<mTempTreeTauTrkIso[i]  << std::endl;
    mTempTreeTauECalIso[i] = (*tauHandle)[i].ecalIso();
    mTempTreeTauHCalIso[i] = (*tauHandle)[i].hcalIso() ;
    mTempTreeTauAllIso[i] = (*tauHandle)[i].caloIso() ;

//     //MICHELE
    mTempTreeTauVx[i] =(*tauHandle)[i].vx();
    mTempTreeTauVy[i] =(*tauHandle)[i].vy();
    mTempTreeTauVz[i] =(*tauHandle)[i].vz();
    mTempTreeTauNTks[i] =(*tauHandle)[i].isolationTracks().size();
    //NEUTRAL
    if ((*tauHandle)[i].isPFTau()){
      int ntnsize=(*tauHandle)[i].signalPFNeutrHadrCands().size();
      mTempTreeTauNNeutrals[i] =ntnsize;
      float hentau=0.;
      float eentau=0.;
      for (int itneu=0;itneu<ntnsize;itneu++){
         eentau+=(*(*tauHandle)[i].signalPFNeutrHadrCands()[itneu]).ecalEnergy();
         hentau+=(*(*tauHandle)[i].signalPFNeutrHadrCands()[itneu]).hcalEnergy();
      }
      mTempTreeTauNeutralE[i]=hentau+eentau;
      mTempTreeTauNeutralHOverHPlusE[i] = (hentau+eentau>0.) ? hentau/(hentau+eentau) : 999.;
   }
    else {
      mTempTreeTauNeutralE[i]=999.;
      mTempTreeTauNeutralHOverHPlusE[i]=999.;
      mTempTreeTauNNeutrals[i] =999.;
    }
    if ((*tauHandle)[i].isolationTracks().size()>0){
      //TK1
      mTempTreeTauTk1Vx[i]=(*tauHandle)[i].isolationTracks()[0]->vx();
      mTempTreeTauTk1Vy[i]=(*tauHandle)[i].isolationTracks()[0]->vy();
      mTempTreeTauTk1Vz[i]=(*tauHandle)[i].isolationTracks()[0]->vz();
      mTempTreeTauTk1D0[i]=(*tauHandle)[i].isolationTracks()[0]->d0();
      mTempTreeTauTk1Dz[i]=(*tauHandle)[i].isolationTracks()[0]->dz();
      mTempTreeTauTk1Pt[i]=(*tauHandle)[i].isolationTracks()[0]->pt(); 
      mTempTreeTauTk1Pz[i]=(*tauHandle)[i].isolationTracks()[0]->pz();
      mTempTreeTauTk1Eta[i]=(*tauHandle)[i].isolationTracks()[0]->eta();
      mTempTreeTauTk1Phi[i]=(*tauHandle)[i].isolationTracks()[0]->phi();
      mTempTreeTauTk1Chi[i]=(*tauHandle)[i].isolationTracks()[0]->chi2();
      mTempTreeTauTk1Charge[i]=(*tauHandle)[i].isolationTracks()[0]->charge();
      mTempTreeTauTk1QOverPError[i]=(*tauHandle)[i].isolationTracks()[0]->qoverpError();
      mTempTreeTauTk1ValidHits[i]=(*tauHandle)[i].isolationTracks()[0]->found();
      mTempTreeTauTk1LostHits[i]=(*tauHandle)[i].isolationTracks()[0]->lost();
      mTempTreeTauTk1CaloE[i]=2.;
    }
    else {
      //TK1
      mTempTreeTauTk1Vx[i]=999.;
      mTempTreeTauTk1Vy[i]=999.;
      mTempTreeTauTk1Vz[i]=999.;
      mTempTreeTauTk1D0[i]=999.;
      mTempTreeTauTk1Dz[i]=999.;
      mTempTreeTauTk1Pt[i]=999.; 
      mTempTreeTauTk1Pz[i]=999.;
      mTempTreeTauTk1Eta[i]=999.;
      mTempTreeTauTk1Phi[i]=999.;
      mTempTreeTauTk1Chi[i]=999.;
      mTempTreeTauTk1Charge[i]=999.;
      mTempTreeTauTk1QOverPError[i]=999.;
      mTempTreeTauTk1ValidHits[i]=999.;
      mTempTreeTauTk1LostHits[i]=999.;
      mTempTreeTauTk1CaloE[i]=999.;
    }
    //TK2
    if ((*tauHandle)[i].isolationTracks().size()>1){
      //TK2
      mTempTreeTauTk2Vx[i]=(*tauHandle)[i].isolationTracks()[1]->vx();
      mTempTreeTauTk2Vy[i]=(*tauHandle)[i].isolationTracks()[1]->vy();
      mTempTreeTauTk2Vz[i]=(*tauHandle)[i].isolationTracks()[1]->vz();
      mTempTreeTauTk2D0[i]=(*tauHandle)[i].isolationTracks()[1]->d0();
      mTempTreeTauTk2Dz[i]=(*tauHandle)[i].isolationTracks()[1]->dz();
      mTempTreeTauTk2Pt[i]=(*tauHandle)[i].isolationTracks()[1]->pt(); 
      mTempTreeTauTk2Pz[i]=(*tauHandle)[i].isolationTracks()[1]->pz();
      mTempTreeTauTk2Eta[i]=(*tauHandle)[i].isolationTracks()[1]->eta();
      mTempTreeTauTk2Phi[i]=(*tauHandle)[i].isolationTracks()[1]->phi();
      mTempTreeTauTk2Chi[i]=(*tauHandle)[i].isolationTracks()[1]->chi2();
      mTempTreeTauTk2Charge[i]=(*tauHandle)[i].isolationTracks()[1]->charge();
      mTempTreeTauTk2QOverPError[i]=(*tauHandle)[i].isolationTracks()[1]->qoverpError();
      mTempTreeTauTk2ValidHits[i]=(*tauHandle)[i].isolationTracks()[1]->found();
      mTempTreeTauTk2LostHits[i]=(*tauHandle)[i].isolationTracks()[1]->lost();
      mTempTreeTauTk2CaloE[i]=2.;
    }
    else {
      //TK2
      mTempTreeTauTk2Vx[i]=999.;
      mTempTreeTauTk2Vy[i]=999.;
      mTempTreeTauTk2Vz[i]=999.;
      mTempTreeTauTk2D0[i]=999.;
      mTempTreeTauTk2Dz[i]=999.;
      mTempTreeTauTk2Pt[i]=999.; 
      mTempTreeTauTk2Pz[i]=999.;
      mTempTreeTauTk2Eta[i]=999.;
      mTempTreeTauTk2Phi[i]=999.;
      mTempTreeTauTk2Chi[i]=999.;
      mTempTreeTauTk2Charge[i]=999.;
      mTempTreeTauTk2QOverPError[i]=999.;
      mTempTreeTauTk2ValidHits[i]=999.;
      mTempTreeTauTk2LostHits[i]=999.;
      mTempTreeTauTk2CaloE[i]=999.;
    }
    //TK3
    if ((*tauHandle)[i].isolationTracks().size()>2){
      //TK2
      mTempTreeTauTk3Vx[i]=(*tauHandle)[i].isolationTracks()[2]->vx();
      mTempTreeTauTk3Vy[i]=(*tauHandle)[i].isolationTracks()[2]->vy();
      mTempTreeTauTk3Vz[i]=(*tauHandle)[i].isolationTracks()[2]->vz();
      mTempTreeTauTk3D0[i]=(*tauHandle)[i].isolationTracks()[2]->d0();
      mTempTreeTauTk3Dz[i]=(*tauHandle)[i].isolationTracks()[2]->dz();
      mTempTreeTauTk3Pt[i]=(*tauHandle)[i].isolationTracks()[2]->pt(); 
      mTempTreeTauTk3Pz[i]=(*tauHandle)[i].isolationTracks()[2]->pz();
      mTempTreeTauTk3Eta[i]=(*tauHandle)[i].isolationTracks()[2]->eta();
      mTempTreeTauTk3Phi[i]=(*tauHandle)[i].isolationTracks()[2]->phi();
      mTempTreeTauTk3Chi[i]=(*tauHandle)[i].isolationTracks()[2]->chi2();
      mTempTreeTauTk3Charge[i]=(*tauHandle)[i].isolationTracks()[2]->charge();
      mTempTreeTauTk3QOverPError[i]=(*tauHandle)[i].isolationTracks()[2]->qoverpError();
      mTempTreeTauTk3ValidHits[i]=(*tauHandle)[i].isolationTracks()[2]->found();
      mTempTreeTauTk3LostHits[i]=(*tauHandle)[i].isolationTracks()[2]->lost();
      mTempTreeTauTk3CaloE[i]=2.;
    }
    else {
      //TK2
      mTempTreeTauTk3Vx[i]=999.;
      mTempTreeTauTk3Vy[i]=999.;
      mTempTreeTauTk3Vz[i]=999.;
      mTempTreeTauTk3D0[i]=999.;
      mTempTreeTauTk3Dz[i]=999.;
      mTempTreeTauTk3Pt[i]=999.; 
      mTempTreeTauTk3Pz[i]=999.;
      mTempTreeTauTk3Eta[i]=999.;
      mTempTreeTauTk3Phi[i]=999.;
      mTempTreeTauTk3Chi[i]=999.;
      mTempTreeTauTk3Charge[i]=999.;
      mTempTreeTauTk3QOverPError[i]=999.;
      mTempTreeTauTk3ValidHits[i]=999.;
      mTempTreeTauTk3LostHits[i]=999.;
      mTempTreeTauTk3CaloE[i]=999.;
    }

    if (&(*(*tauHandle)[i].genLepton())!=0){
      mTempTreeGenTauPdgId[i]=(*tauHandle)[i].genLepton()->pdgId();
      mTempTreeGenTauPx[i]=(*tauHandle)[i].genLepton()->px();
      mTempTreeGenTauPy[i]=(*tauHandle)[i].genLepton()->py();
      mTempTreeGenTauPz[i]=(*tauHandle)[i].genLepton()->pz();
      if (&(*(*tauHandle)[i].genLepton()->mother())!=0)
	mTempTreeGenTauMother[i]=(*tauHandle)[i].genLepton()->mother()->pdgId();
      else mTempTreeGenTauMother[i]=999.;
    }
    else{
      mTempTreeGenTauPdgId[i]=999.;
      mTempTreeGenTauMother[i]=999.;
      mTempTreeGenTauPx[i]=999.;
      mTempTreeGenTauPy[i]=999.;
      mTempTreeGenTauPz[i]=999.;
    }
    
    //PIOPPI
    
  }

   
  // get the jets
  edm::Handle< std::vector<pat::Jet> > jetHandle;
  iEvent.getByLabel(jetTag_, jetHandle);
  if ( !jetHandle.isValid() ) {
    edm::LogWarning("SusySelectorExample") << "No Jet results for InputTag " << jetTag_;
    return;
  }


  //get number of jets
  mTempTreeNjets = jetHandle->size();



  // Add the jets
  int i=0;
  if ( mTempTreeNjets >50 ) mTempTreeNjets = 50;
  for (int k=0;k<mTempTreeNjets;k++){
    if ( (*jetHandle)[k].et() > 20. ){
      mTempTreeJetsPt[i] = (*jetHandle)[i].pt();
      mTempTreeJetsE[i] = (*jetHandle)[i].energy();
      mTempTreeJetsEt[i] = (*jetHandle)[i].et();
      mTempTreeJetsPx[i] = (*jetHandle)[i].momentum().X();
      mTempTreeJetsPy[i] = (*jetHandle)[i].momentum().Y();
      mTempTreeJetsPz[i] = (*jetHandle)[i].momentum().Z();
      mTempTreeJetsEta[i] = (*jetHandle)[i].eta();
      mTempTreeJetsPhi[i] = (*jetHandle)[i].phi();
      if ((*jetHandle)[i].isCaloJet())
	mTempTreeJetsFem[i] = (*jetHandle)[i].emEnergyFraction();
      if ((*jetHandle)[i].isPFJet())
	mTempTreeJetsFem[i] = (*jetHandle)[i].neutralEmEnergyFraction()+
	  (*jetHandle)[i].chargedEmEnergyFraction();  


       mTempTreeJetsBTag_TkCountHighEff[i] = (*jetHandle)[i].bDiscriminator("trackCountingHighEffBJetTags");
       // std::cout << " discri " << (*jetHandle)[i].bDiscriminator("trackCountingHighEffBJetTags") << std::endl;
      mTempTreeJetsBTag_SimpleSecVtx[i] = (*jetHandle)[i].bDiscriminator("simpleSecondaryVertexBJetTags");
      mTempTreeJetsBTag_CombSecVtx[i] = (*jetHandle)[i].bDiscriminator("combinedSecondaryVertexBJetTags");
      mTempTreeJetPartonFlavour[i] = (*jetHandle)[i].partonFlavour();
    
      if((*jetHandle)[i].genJet()!= 0) {
	mTempTreeGenJetsPt[i]=(*jetHandle)[i].genJet()->pt();
	mTempTreeGenJetsE[i]=(*jetHandle)[i].genJet()->energy();
	mTempTreeGenJetsEt[i]=(*jetHandle)[i].genJet()->et();
	mTempTreeGenJetsPx[i]=(*jetHandle)[i].genJet()->momentum().X();
	mTempTreeGenJetsPy[i]=(*jetHandle)[i].genJet()->momentum().Y();
	mTempTreeGenJetsPz[i]=(*jetHandle)[i].genJet()->momentum().z();
	mTempTreeGenJetsEta[i]=(*jetHandle)[i].genJet()->eta();
	mTempTreeGenJetsPhi[i]=(*jetHandle)[i].genJet()->phi();
      }
      else {
	mTempTreeGenJetsPt[i]  =-999;
	mTempTreeGenJetsE[i]   =-999;
	mTempTreeGenJetsEt[i]  =-999;
	mTempTreeGenJetsPx[i]  =-999;
	mTempTreeGenJetsPy[i]  =-999;
	mTempTreeGenJetsPz[i]  =-999;
	mTempTreeGenJetsEta[i] =-999;
	mTempTreeGenJetsPhi[i] =-999;
      }

      if((*jetHandle)[i].genParton() != 0){
	mTempTreeJetPartonId[i] = (*jetHandle)[i].genParton()->pdgId();
	mTempTreeJetPartonMother[i] = (*jetHandle)[i].genParton()->mother()->pdgId();
	mTempTreeJetPartonPx[i] = (*jetHandle)[i].genParton()->px();
	mTempTreeJetPartonPy[i] = (*jetHandle)[i].genParton()->py();
	mTempTreeJetPartonPz[i] = (*jetHandle)[i].genParton()->pz();
	mTempTreeJetPartonEt[i] = (*jetHandle)[i].genParton()->et();
	mTempTreeJetPartonEnergy[i] = (*jetHandle)[i].genParton()->energy();
	mTempTreeJetPartonPhi[i] = (*jetHandle)[i].genParton()->phi();
	mTempTreeJetPartonEta[i] = (*jetHandle)[i].genParton()->eta();
      }
      else{
	mTempTreeJetPartonId[i] = -999;
	mTempTreeJetPartonMother[i] = -999;
	mTempTreeJetPartonPx[i] = -999;
	mTempTreeJetPartonPy[i] = -999;
	mTempTreeJetPartonPz[i] = -999;
	mTempTreeJetPartonEt[i] = -999;
	mTempTreeJetPartonEnergy[i] = -999;
	mTempTreeJetPartonPhi[i] = -999;
	mTempTreeJetPartonEta[i] = -999;
      }
     
      //  const reco::TrackRefVector & associatedTracks();      
      /*   mTempTreeNJetsT[i] = ((*jetHandle)[i].associatedTracks()).size();
      if(((*jetHandle)[i].associatedTracks()).isAvailable()) { 
	
	float maxPt = 0;
	int maxPtIndex = 0;
	for(int pttest=0; pttest<mTempTreeNJetsT[i];pttest++) { 
	  if(maxPt<float(((*jetHandle)[i].associatedTracks())[pttest]->pt())) {
	    maxPt = float(((*jetHandle)[i].associatedTracks())[pttest]->pt());
	    maxPtIndex = pttest;
	  }
	  
	}
	mTempTreeJetsTPt[i] = float(((*jetHandle)[i].associatedTracks())[maxPtIndex]->pt());
	mTempTreeJetsTEta[i] = float(((*jetHandle)[i].associatedTracks())[maxPtIndex]->eta());
	mTempTreeJetsTPhi[i] = float(((*jetHandle)[i].associatedTracks())[maxPtIndex]->phi());
      }
      else {
	mTempTreeJetsTPt[i] = -1 ;
	mTempTreeJetsTEta[i] = -10 ;
	mTempTreeJetsTPhi[i] = -10 ;
	}   */
      i++;
    }
  }
  
  mTempTreeNjets = i;

  //std::cout << " add the ccjets " << std::endl;
  
  //Get the cross-cleaned Jets
  edm::Handle< std::vector<pat::Jet> > ccjetHandle;
  iEvent.getByLabel(ccjetTag_, ccjetHandle);
  if ( !ccjetHandle.isValid() ) {
    edm::LogWarning("SusySelectorExample") << "No ccJet results for InputTag " << ccjetTag_;
    return;
  }

  //get number of cc-jets
  mTempTreeNccjets = ccjetHandle->size();

  // Add the cc jets
  i=0;
  if ( mTempTreeNccjets >50 ) mTempTreeNccjets = 50;
  for (int k=0;k<mTempTreeNccjets;k++){
    if ( (*ccjetHandle)[k].et() > 20. ){
      mTempTreeccJetsPt[i] = (*ccjetHandle)[i].pt();
      mTempTreeccJetsE[i] = (*ccjetHandle)[i].energy();
      mTempTreeccJetsEt[i] = (*ccjetHandle)[i].et();
      mTempTreeccJetsPx[i] = (*ccjetHandle)[i].momentum().X();
      mTempTreeccJetsPy[i] = (*ccjetHandle)[i].momentum().Y();
      mTempTreeccJetsPz[i] = (*ccjetHandle)[i].momentum().Z();
      mTempTreeccJetsEta[i] = (*ccjetHandle)[i].eta();
      mTempTreeccJetsPhi[i] = (*ccjetHandle)[i].phi();
      // std::cout << " cc jet " << (*ccjetHandle)[i].pt() << std::endl;
      i++;
    }
  }
  
  mTempTreeNccjets = i;

// Get the hemispheres
  Handle< edm::View<pat::Hemisphere> > hemisphereHandle;
  iEvent.getByLabel("selectedLayer1Hemispheres", hemisphereHandle);
  if ( !hemisphereHandle.isValid() ) {
    edm::LogWarning("SusySelectorExample") << "No Hemisphere results for InputTag ";
    return;
  }
  const edm::View<pat::Hemisphere>& hemispheres = (*hemisphereHandle); // For simplicity...
  
  mTempTreeNhemispheres = 2;
  for (unsigned int i=0;i <  hemispheres.size() ;i++){
    mTempTreeHemispheresPt[i] = hemispheres[i].pt();
    
    mTempTreeHemispheresE[i] = hemispheres[i].energy();
    mTempTreeHemispheresEt[i] = hemispheres[i].et();
    mTempTreeHemispheresPx[i] = hemispheres[i].momentum().X();
    mTempTreeHemispheresPy[i] = hemispheres[i].momentum().Y();
    mTempTreeHemispheresPz[i] = hemispheres[i].momentum().Z();
    mTempTreeHemispheresEta[i] = hemispheres[i].eta();
    mTempTreeHemispheresPhi[i] = hemispheres[i].phi();

    for(unsigned int dau = 0; dau < hemispheres[i].numberOfDaughters();dau++){
      for (int k=0;k<mTempTreeNjets;k++){
	//	mTempTreeJetsHemi[k]= -1;
	if(  hemispheres[i].daughter(dau)->phi() >= mTempTreeJetsPhi[k]-0.0001 &  hemispheres[i].daughter(dau)->phi() <= mTempTreeJetsPhi[k]+0.0001 )  mTempTreeJetsHemi[k] = i; 
	//	std::cout << " k " << k << " hemiphi " << hemispheres[i].daughter(dau)->phi() << " jetphi " << mTempTreeJetsPhi[k] << " hemi " << mTempTreeJetsHemi[k] << std::endl;
	}
     
    }


  }   

 
  // for (int k=0;k<mTempTreeNjets;k++)  std::cout << " after " << " k " << k << " hemiphi " <<  " jetphi " << mTempTreeJetsPhi[k] <<" hemi " << mTempTreeJetsHemi[k] << std::endl;


  //
  // get the MET result
  //
  edm::Handle< std::vector<pat::MET> > metHandle;
  iEvent.getByLabel(metTag_, metHandle);
  if ( !metHandle.isValid() ) {
    edm::LogWarning("METEventSelector") << "No Met results for InputTag " << metTag_;
    return;
  }
  
  //
  // sanity check on collection
  //
  if ( metHandle->size()!=1 ) {
    edm::LogWarning("METEventSelector") << "MET collection size is "
					<< metHandle->size() << " instead of 1";
    return;
  }
  
  // Do the MET save
  mTempTreeMET = metHandle->front().et();
  mTempTreeMEX = metHandle->front().momentum().X();
  mTempTreeMEY = metHandle->front().momentum().Y();
  mTempTreeSumET = metHandle->front().sumEt();
  mTempTreeMETphi = metHandle->front().phi();
  mTempTreeSumETSignif = metHandle->front().mEtSig();
  
  mTempTreeMETuncor = metHandle->front().uncorrectedPt(pat::MET::UncorrectionType(2));
  mTempTreeMETphiuncor = metHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(2));

 //
  // get the cc MET result
  //

  //std::cout << " add the ccmet " << std::endl;
  
  edm::Handle< std::vector<pat::MET> > ccmetHandle;
  iEvent.getByLabel(ccmetTag_, ccmetHandle);
  if ( !ccmetHandle.isValid() ) {
    edm::LogWarning("METEventSelector") << "No Met results for InputTag " << metTag_;
    return;
  }
  
  //
  // sanity check on collection
  //
  if ( ccmetHandle->size()!=1 ) {
    edm::LogWarning("METEventSelector") << "ccMET collection size is "
					<< ccmetHandle->size() << " instead of 1";
    return;
  }
  
  // Do the MET save
  mTempTreeccMET = ccmetHandle->front().et();
  mTempTreeccMEX = ccmetHandle->front().momentum().X();
  mTempTreeccMEY = ccmetHandle->front().momentum().Y();
  mTempTreeccSumET = ccmetHandle->front().sumEt();
  mTempTreeccMETphi = ccmetHandle->front().phi();

  
  mTempAlpIdTest = myALPGENParticleId.AplGenParID(iEvent,genTag_);
  mTempAlpPtScale = myALPGENParticleId.getPt();
 
 
  //set information of event is affected by b-bug
  // is_ok = true;
   length = 1000;
  length =  myALPGENParticleId.AplGenParID(iEvent,genTag_,  ids , refs ,genE, genPx, genPy, genPz ,genPhi ,genEta ,genStatus, length);
  /*
bool BadFlag = false;
 
   for (int i=0;i<jetHandle->size();i++) 
   {
    
    if( (*jetHandle)[i].genParton() != 0) //<---Check if jet is matched
      {  
	if(fabs((*jetHandle)[i].genParton()->pdgId()) == 5 && (*jetHandle)[i].genParton()->pt() > 100 ) {
	  if(fabs((*jetHandle)[i].genParton()->pt()-(*jetHandle)[i].pt( ))/(*jetHandle)[i].pt() > 5 ) BadFlag == true;
	}

      }
   }//<---End of i loop
  */
   // std::cout << " bad flag " << BadFlag << std::endl;
  
  /*  mTempSimuCheck = -1 ;


  float min_dR;
  int matched_jet = 0;
  int matched_b = 0;
  for(int i=0 ; i < length;i++) {
    min_dR = 100;
    matched_jet = -1;
    matched_b = -1;
    if(fabs(ids[i]) == 5 && genPt[i] > 100){
      for (int jeti=0;jeti < mTempTreeNjets;jeti++){
	
	if(sqrt(reco::deltaPhi((*jetHandle)[jeti].phi(),genPhi[i])
		*reco::deltaPhi((*jetHandle)[jeti].phi(),genPhi[i])
		+((*jetHandle)[jeti].eta()-genEta[i])*((*jetHandle)[jeti].eta()-genEta[i])) < min_dR) {
	  min_dR = sqrt(reco::deltaPhi((*jetHandle)[jeti].phi(),genPhi[i]) 
			* reco::deltaPhi((*jetHandle)[jeti].phi(),genPhi[i])
			+((*jetHandle)[jeti].eta()-genEta[i]) * ((*jetHandle)[jeti].eta()-genEta[i]));
	  matched_jet = jeti;
	  matched_b = i; 
	}
      }
      if(fabs((*jetHandle)[matched_jet].pt()- genPt[matched_b])/(*jetHandle)[matched_jet].pt() > 5 && min_dR < 0.5)  is_ok = false;
    }
  }

  */
 
   
  
  // Fill the tree
  mAllData->Fill();

 
}


//________________________________________________________________________________________
void 
SusyDiJetAnalysis::beginJob(const edm::EventSetup&) {}

//________________________________________________________________________________________
void 
SusyDiJetAnalysis::endJob() {

  printSummary();
  printHLTreport();

}


//________________________________________________________________________________________
void
SusyDiJetAnalysis::printSummary( void ) {

  edm::LogWarning("SusyDiJet|SummaryCount") << "*** Summary of counters: ";
  edm::LogWarning("SusyDiJet|SummaryCount") 
    << "Total number of events = " << nrEventTotalWeighted_ 
    << " (" << nrEventTotalRaw_ << " unweighted)"
    << " ; selected = " << nrEventCumulative_.back();

  std::ostringstream summary;
  summary << std::setw(21) << std::left << "Name"
          << std::setw(21) << "Selected"
          << std::setw(21) << "AllButOne"
          << std::setw(21) << "Cumulative" << std::endl;
         
  for ( size_t i=0; i<sequence_.size(); ++i ) {
    summary << std::setw(20) << std::left << sequence_.selectorName(i) << std::right;
    summary << std::setw(10) << std::setprecision(2)  << std::fixed
            << nrEventSelected_[i] 
            << "[" << std::setw(6) 
            << (nrEventSelected_[i]/nrEventTotalWeighted_)*100. << "%]  ";
    summary << std::setw(10) << nrEventAllButOne_[i] 
            << "[" << std::setw(6) 
            << (nrEventAllButOne_[i]/nrEventTotalWeighted_)*100. << "%]  ";
    summary << std::setw(10) << nrEventCumulative_[i] 
            << "[" << std::setw(6) 
            << (nrEventCumulative_[i]/nrEventTotalWeighted_)*100. << "%]  ";
    summary << std::endl; 
  }
  edm::LogWarning("SusyDiJet|SummaryCount") << summary.str();

}


// GEORGIA
void
SusyDiJetAnalysis::printHLTreport( void ) {

  // georgia :  prints an HLT report -- associates trigger bits with trigger names (prints #events fired the trigger etc)
 const unsigned int n(pathNames_.size());
  std::cout << "\n";
  std::cout << "HLT-Report " << "---------- Event  Summary ------------\n";
  std::cout << "HLT-Report"
	    << " Events total = " << nEvents_
	    << " wasrun = " << nWasRun_
	    << " passed = " << nAccept_
	    << " errors = " << nErrors_
	    << "\n";

  std::cout << endl;
  std::cout << "HLT-Report " << "---------- HLTrig Summary ------------\n";
  std::cout << "HLT-Report "
	    << right << setw(10) << "HLT  Bit#" << " "
	    << right << setw(10) << "WasRun" << " "
	    << right << setw(10) << "Passed" << " "
	    << right << setw(10) << "Errors" << " "
	    << "Name" << "\n";

  if (init_) {
    for (unsigned int i=0; i!=n; ++i) {
      std::cout << "HLT-Report "
		<< right << setw(10) << i << " "
		<< right << setw(10) << hlWasRun_[i] << " "
		<< right << setw(10) << hlAccept_[i] << " "
		<< right << setw(10) << hlErrors_[i] << " "
		<< pathNames_[i] << "\n";
    }
  } else {
    std::cout << "HLT-Report - No HL TriggerResults found!" << endl;
  }
  
  std::cout << endl;
  std::cout << "HLT-Report end!" << endl;
  std::cout << endl;

}
// end GEORGIA


//________________________________________________________________________________________
void
SusyDiJetAnalysis::initPlots() {

  std::ostringstream variables; // Container for all variables

  // 1. Event variables
  variables << "weight:process";
 
  // 2. Decision from all selectors
  for ( std::vector<const SusyEventSelector*>::const_iterator it = sequence_.selectors().begin();
        it != sequence_.selectors().end(); ++it ) {
    std::string var( (*it)->name() );
    var += "_result";
    // Push to list of variables
    variables << ":" << var;
  }
  variables << ":all_result"; // Also store global decision
  
  // 3. All variables from sequence
  for ( std::vector<const SusyEventSelector*>::const_iterator it = sequence_.selectors().begin();
        it != sequence_.selectors().end(); ++it ) {
    for ( unsigned int i=0; i<(*it)->numberOfVariables(); ++i ) {
      std::string var( (*it)->name() ); // prefix variable with selector name
      var += "." + (*it)->variableNames()[i];
      // Push to list of variables
      variables << ":" << var;
    }
  }


  // Register this ntuple
  edm::Service<TFileService> fs;
  ntuple_ = fs->make<TNtuple>( "ntuple","SusyDiJetAnalysis variables",
                               variables.str().c_str() );

  // Now we add some additional ones for the dijet analysis
  mAllData = fs->make<TTree>( "allData", "data after cuts" );
  mSelectorData = fs->make<TTree>( "selectorData" , "Bit results for selectors");

  mAllData->SetAutoSave(10000);
  mSelectorData->SetAutoSave(10000);

  std::vector<std::string> names = sequence_.selectorNames();
  for ( size_t i = 0 ; i < sequence_.size() ; ++i ) {
    std::string tempName = names[i] + "/i";
    mSelectorData->Branch(names[i].c_str(),&mSelectorResults[i],tempName.c_str());
  }
  mSelectorData->Branch("globalDecision",&mGlobalDecision,"globalDecision/i");
  

  // Add the branches
  mAllData->Branch("run",&mTempTreeRun,"run/int");
  mAllData->Branch("event",&mTempTreeEvent,"event/int");

  // GEORGIA
  mAllData->Branch("nHLT",&nHLT,"nHLT/I");
  mAllData->Branch("HLTArray",HLTArray,"HLTArray[nHLT]/I");
  // end GEORGIA 

  mAllData->Branch("HLT1JET",&mTempTreeHLT1JET,"HLT1JET/bool");
  mAllData->Branch("HLT2JET",&mTempTreeHLT2JET,"HLT2JET/bool");
  mAllData->Branch("HLT1MET1HT",&mTempTreeHLT1MET1HT,"HLT1MET1HT/bool");
  mAllData->Branch("HLT1MUON",&mTempTreeHLT1Muon,"HLT1MUON/bool");
  

  mAllData->Branch("met",&mTempTreeMET,"met/double");
  mAllData->Branch("mex",&mTempTreeMEY,"mex/double");
  mAllData->Branch("mey",&mTempTreeMEX,"mey/double");
  mAllData->Branch("metphi",&mTempTreeMETphi,"metphi/double");
  mAllData->Branch("metuncor",&mTempTreeMETuncor,"metuncor/double");
  mAllData->Branch("metphiuncor",&mTempTreeMETphiuncor,"metphiuncor/double");


  mAllData->Branch("evtWeight",&mTempTreeEventWeight,"evtWeight/double");
  mAllData->Branch("procID",&mTempTreeProcID,"procID/int");
  mAllData->Branch("pthat",&mTempTreePthat,"pthat/double");


  // GEORGIA
  mAllData->Branch("nVtx",&mTempTreenVtx,"nVtx/int");
  mAllData->Branch("VertexChi2",mTempTreeVtxChi2,"VertexChi2[nVtx]/double");
  mAllData->Branch("VertexNdof",mTempTreeVtxNdof,"VertexNdof[nVtx]/double");
  mAllData->Branch("VertexNormalizedChi2",mTempTreeVtxNormalizedChi2,"VertexNormalizedChi2[nVtx]/double");
  mAllData->Branch("VertexX",mTempTreeVtxX,"VertexX[nVtx]/double");
  mAllData->Branch("VertexY",mTempTreeVtxY,"VertexY[nVtx]/double");
  mAllData->Branch("VertexZ",mTempTreeVtxZ,"VertexZ[nVtx]/double");
  mAllData->Branch("VertexdX",mTempTreeVtxdX,"VertexdX[nVtx]/double");
  mAllData->Branch("VertexdY",mTempTreeVtxdY,"VertexdY[nVtx]/double");
  mAllData->Branch("VertexdZ",mTempTreeVtxdZ,"VertexdZ[nVtx]/double");
  // end GEORGIA

 
 //add hemispheres
  mAllData->Branch("Nhemispheres" ,&mTempTreeNhemispheres ,"Nhemispheres/int");  
  mAllData->Branch("HemisphereE" ,mTempTreeHemispheresE ,"HemisphereE[Nhemispheres]/double");
  mAllData->Branch("HemisphereEt",mTempTreeHemispheresEt,"HemisphereEt[Nhemispheres]/double");
  mAllData->Branch("Hemispherept",mTempTreeHemispheresPt,"Hemispherept[Nhemispheres]/double");
  mAllData->Branch("Hemispherepx",mTempTreeHemispheresPx,"Hemispherepx[Nhemispheres]/double");
  mAllData->Branch("Hemispherepy",mTempTreeHemispheresPy,"Hemispherepy[Nhemispheres]/double");
  mAllData->Branch("Hemispherepz",mTempTreeHemispheresPz,"Hemispherepz[Nhemispheres]/double");
  mAllData->Branch("Hemisphereeta",mTempTreeHemispheresEta,"Hemisphereeta[Nhemispheres]/double");
  mAllData->Branch("Hemispherephi",mTempTreeHemispheresPhi,"Hemispherephi[Nhemispheres]/double");

  //  mAllData->Branch("is_ok",&is_ok,"is_ok/bool");

  //  mAllData->Branch("SimuCheck",&mTempSimuCheck,"mTempSimuCheck/int");

  //add jets
  mAllData->Branch("Njets" ,&mTempTreeNjets ,"Njets/int");  
  mAllData->Branch("JetE" ,mTempTreeJetsE ,"JetE[Njets]/double");
  mAllData->Branch("JetEt",mTempTreeJetsEt,"JetEt[Njets]/double");
  mAllData->Branch("Jetpt",mTempTreeJetsPt,"Jetpt[Njets]/double");
  mAllData->Branch("Jetpx",mTempTreeJetsPx,"Jetpx[Njets]/double");
  mAllData->Branch("Jetpy",mTempTreeJetsPy,"Jetpy[Njets]/double");
  mAllData->Branch("Jetpz",mTempTreeJetsPz,"Jetpz[Njets]/double");
  mAllData->Branch("Jeteta",mTempTreeJetsEta,"Jeteta[Njets]/double");
  mAllData->Branch("Jetphi",mTempTreeJetsPhi,"Jetphi[Njets]/double");
  mAllData->Branch("JetFem",mTempTreeJetsFem,"JetFem[Njets]/double");
  /*  mAllData->Branch("NJetsT",mTempTreeNJetsT,"NJetsT[Njets]/int");
  mAllData->Branch("JetsTPt",mTempTreeJetsTPt,"JetsTPt[Njets]/float");
  mAllData->Branch("JetsTEta",mTempTreeJetsTEta,"JetsTEta[Njets]/float");
  mAllData->Branch("JetsTPhi",mTempTreeJetsTPhi,"JetsTPhi[Njets]/float");*/
  mAllData->Branch("JetHemi",mTempTreeJetsHemi,"JetHemi[Njets]/int");
 
  mAllData->Branch("JetPartonId" ,mTempTreeJetPartonId ,"JetPartonId[Njets]/int"); 
  mAllData->Branch("JetPartonMother" ,mTempTreeJetPartonMother ,"JetPartonMother[Njets]/int"); 
  mAllData->Branch("JetPartonPx", mTempTreeJetPartonPx ,"JetPartonPx[Njets]/double"); 
  mAllData->Branch("JetPartonPy" ,mTempTreeJetPartonPy ,"JetPartonPy[Njets]/double"); 
  mAllData->Branch("JetPartonPz" ,mTempTreeJetPartonPz ,"JetPartonPz[Njets]/double"); 
  mAllData->Branch("JetPartonEt" ,mTempTreeJetPartonEt ,"JetPartonEt[Njets]/double"); 
  mAllData->Branch("JetPartonE" ,mTempTreeJetPartonEnergy ,"JetPartonE[Njets]/double"); 
  mAllData->Branch("JetPartonPhi" ,mTempTreeJetPartonPhi ,"JetPartonPhi[Njets]/double"); 
  mAllData->Branch("JetPartonEta" ,mTempTreeJetPartonEta ,"JetPartonEta[Njets]/double"); 
  
 

  mAllData->Branch("JetPartonFlavour",mTempTreeJetPartonFlavour,"JetPartonFlavour[Njets]/int");

 

  mAllData->Branch("JetBTag_TkCountHighEff",mTempTreeJetsBTag_TkCountHighEff,"JetBTag_TkCountHighEff[Njets]/float");
  mAllData->Branch("JetBTag_SimpleSecVtx",mTempTreeJetsBTag_SimpleSecVtx,"JetBTag_SimpleSecVtx[Njets]/float");
  mAllData->Branch("JetBTag_CombSecVtx",mTempTreeJetsBTag_CombSecVtx,"JetBTag_CombSecVtx[Njets]/float");



  mAllData->Branch("GenJetE" ,mTempTreeGenJetsE ,"GenJetE[Njets]/double");
  mAllData->Branch("GenJetEt",mTempTreeGenJetsEt,"GenJetEt[Njets]/double");
  mAllData->Branch("GenJetpt",mTempTreeGenJetsPt,"GenJetpt[Njets]/double");
  mAllData->Branch("GenJetpx",mTempTreeGenJetsPx,"GenJetpx[Njets]/double");
  mAllData->Branch("GenJetpy",mTempTreeGenJetsPy,"GenJetpy[Njets]/double");
  mAllData->Branch("GenJetpz",mTempTreeGenJetsPz,"GenJetpz[Njets]/double");
  mAllData->Branch("GenJeteta",mTempTreeGenJetsEta,"GenJeteta[Njets]/double");
  mAllData->Branch("GenJetphi",mTempTreeGenJetsPhi,"GenJetphi[Njets]/double");

  

  //ccJets
  mAllData->Branch("Nccjets" ,&mTempTreeNccjets ,"Nccjets/int");  
  mAllData->Branch("ccJetE" ,mTempTreeccJetsE ,"ccJetE[Njets]/double");
  mAllData->Branch("ccJetEt",mTempTreeccJetsEt,"ccJetEt[Njets]/double");
  mAllData->Branch("ccJetpt",mTempTreeccJetsPt,"ccJetpt[Njets]/double");
  mAllData->Branch("ccJetpx",mTempTreeccJetsPx,"ccJetpx[Njets]/double");
  mAllData->Branch("ccJetpy",mTempTreeccJetsPy,"ccJetpy[Njets]/double");
  mAllData->Branch("ccJetpz",mTempTreeccJetsPz,"ccJetpz[Njets]/double");
  mAllData->Branch("ccJeteta",mTempTreeccJetsEta,"ccJeteta[Njets]/double");
  mAllData->Branch("ccJetphi",mTempTreeccJetsPhi,"ccJetphi[Njets]/double");
 
  //add photons
  mAllData->Branch("Nphot" ,&mTempTreeNphot ,"Nphot/int");  
  mAllData->Branch("PhotE" ,mTempTreePhotE ,"PhotE[Nphot]/double");
  mAllData->Branch("PhotEt",mTempTreePhotEt,"PhotEt[Nphot]/double");
  mAllData->Branch("Photpt",mTempTreePhotPt,"Photpt[Nphot]/double");
  mAllData->Branch("Photpx",mTempTreePhotPx,"Photpx[Nphot]/double");
  mAllData->Branch("Photpy",mTempTreePhotPy,"Photpy[Nphot]/double");
  mAllData->Branch("Photpz",mTempTreePhotPz,"Photpz[Nphot]/double");
  mAllData->Branch("Photeta",mTempTreePhotEta,"Photeta[Nphot]/double");
  mAllData->Branch("Photphi",mTempTreePhotPhi,"Photphi[Nphot]/double");
  mAllData->Branch("PhotTrkIso",mTempTreePhotTrkIso,"mTempTreePhotTrkIso[Nphot]/double");
  mAllData->Branch("PhotECalIso",mTempTreePhotECalIso,"mTempTreePhotECalIso[Nphot]/double");
  mAllData->Branch("PhotHCalIso",mTempTreePhotHCalIso,"mTempTreePhotHCalIso[Nphot]/double");
  mAllData->Branch("PhotAllIso",mTempTreePhotAllIso,"mTempTreePhotAllIso[Nphot]/double");

  //add electrons
  mAllData->Branch("Nelec" ,&mTempTreeNelec ,"Nelec/int");  
  mAllData->Branch("ElecE" ,mTempTreeElecE ,"ElecE[Nelec]/double");
  mAllData->Branch("ElecEt",mTempTreeElecEt,"ElecEt[Nelec]/double");
  mAllData->Branch("Elecpt",mTempTreeElecPt,"Elecpt[Nelec]/double");
  mAllData->Branch("Elecpx",mTempTreeElecPx,"Elecpx[Nelec]/double");
  mAllData->Branch("Elecpy",mTempTreeElecPy,"Elecpy[Nelec]/double");
  mAllData->Branch("Elecpz",mTempTreeElecPz,"Elecpz[Nelec]/double");
  mAllData->Branch("Eleceta",mTempTreeElecEta,"Eleceta[Nelec]/double");
  mAllData->Branch("Elecphi",mTempTreeElecPhi,"Elecphi[Nelec]/double");
  mAllData->Branch("ElecCharge",mTempTreeElecCharge,"ElecCharge[Nelec]/double");
  mAllData->Branch("ElecTrkIso",mTempTreeElecTrkIso,"ElecTrkIso[Nelec]/double");
  mAllData->Branch("ElecECalIso", mTempTreeElecECalIso,"ElecECalIso[Nelec]/double");
  mAllData->Branch("ElecHCalIso", mTempTreeElecHCalIso ,"ElecHCalIso[Nelec]/double");
  mAllData->Branch("ElecAllIso",  mTempTreeElecAllIso ,"ElecAllIso[Nelec]/double");
  mAllData->Branch("ElecTrkChiNorm",mTempTreeElecTrkChiNorm  ,"ElecTrkChiNorm[Nelec]/double");
  //MICHELE
  mAllData->Branch("ElecIdLoose",mTempTreeElecIdLoose,"ElecIdLoose [Nelec]/double");
  mAllData->Branch("ElecIdTight",mTempTreeElecIdTight,"ElecIdTight [Nelec]/double");
  mAllData->Branch("ElecIdRobLoose",mTempTreeElecIdRobLoose,"ElecIdRobLoose [Nelec]/double");
  mAllData->Branch("ElecIdRobTight",mTempTreeElecIdRobTight,"ElecIdRobTight [Nelec]/double");

  mAllData->Branch("ElecChargeMode",mTempTreeElecChargeMode,"ElecChargeMode [Nelec]/double");
  mAllData->Branch("ElecPtMode",mTempTreeElecPtTrkMode,"ElecPtMode [Nelec]/double");
  mAllData->Branch("ElecQOverPErrTrkMode",mTempTreeElecQOverPErrTrkMode,"ElecQOverPErrTrkMode [Nelec]/double");
  mAllData->Branch("ElecCaloEnergy",mTempTreeElecCaloEnergy,"ElecCaloEnergy[Nelec]/double");
  mAllData->Branch("ElecHOverE", mTempTreeElecHOverE,"ElecHOverE[Nelec]/double");
  mAllData->Branch("ElecVx",  mTempTreeElecVx,"ElecVx[Nelec]/double");
  mAllData->Branch("ElecVy",mTempTreeElecVy  ,"ElecVy[Nelec]/double");
  mAllData->Branch("ElecVz",mTempTreeElecVz,"ElecVz[Nelec]/double");
  mAllData->Branch("ElecD0",mTempTreeElecD0,"ElecD0[Nelec]/double");
  mAllData->Branch("ElecDz", mTempTreeElecDz,"ElecDz[Nelec]/double");
  mAllData->Branch("ElecPtTrk", mTempTreeElecPtTrk ,"ElecPtTrk[Nelec]/double");
  mAllData->Branch("ElecQOverPErrTrk",  mTempTreeElecQOverPErrTrk ,"ElecQOverPErrTrk[Nelec]/double");
  mAllData->Branch("ElecPinTrk",mTempTreeElecPinTrk  ,"ElecPinTrk[Nelec]/double");
  mAllData->Branch("ElecPoutTrk",mTempTreeElecPoutTrk  ,"ElecPoutTrk[Nelec]/double"); 
  mAllData->Branch("ElecLostHits",mTempTreeElecLostHits  ,"ElecLostHits[Nelec]/double"); 
  mAllData->Branch("ElecValidHits",mTempTreeElecValidHits  ,"ElecValidHits[Nelec]/double"); 
  mAllData->Branch("ElecNCluster",mTempTreeElecNCluster  ,"ElecNCluster[Nelec]/double"); 
  mAllData->Branch("ElecEtaTrk",mTempTreeElecEtaTrk,"ElecEtaTrk[Nelec]/double"); 
  mAllData->Branch("ElecPhiTrk",mTempTreeElecPhiTrk ,"ElecPhiTrk[Nelec]/double"); 
  mAllData->Branch("ElecWidthClusterEta",mTempTreeElecWidthClusterEta ,"ElecWidthClusterEta[Nelec]/double"); 
  mAllData->Branch("ElecWidthClusterPhi",mTempTreeElecWidthClusterPhi ,"ElecWidthClusterPhi[Nelec]/double"); 
  mAllData->Branch("ElecGenPdgId",mTempTreeGenElecPdgId,"ElecGenPdgId[Nelec]/double");
  mAllData->Branch("ElecGenMother",mTempTreeGenElecMother,"ElecGenMother[Nelec]/double");
  mAllData->Branch("ElecGenPx",mTempTreeGenElecPx,"ElecGenPx[Nelec]/double");
  mAllData->Branch("ElecGenPy",mTempTreeGenElecPy,"ElecGenPy[Nelec]/double");
  mAllData->Branch("ElecGenPz",mTempTreeGenElecPz,"ElecGenPz[Nelec]/double");
  //PIOPPI


  mAllData->Branch("Nccelec" ,&mTempTreeNccelec ,"Nccelec/int");  
  mAllData->Branch("ccElecE" ,mTempTreeccElecE ,"ccElecE[Nelec]/double");
  mAllData->Branch("ccElecEt",mTempTreeccElecEt,"ccElecEt[Nelec]/double");
  mAllData->Branch("ccElecpt",mTempTreeccElecPt,"ccElecpt[Nelec]/double");
  mAllData->Branch("ccElecpx",mTempTreeccElecPx,"ccElecpx[Nelec]/double");
  mAllData->Branch("ccElecpy",mTempTreeccElecPy,"ccElecpy[Nelec]/double");
  mAllData->Branch("ccElecpz",mTempTreeccElecPz,"ccElecpz[Nelec]/double");
  mAllData->Branch("ccEleceta",mTempTreeccElecEta,"ccEleceta[Nelec]/double");
  mAllData->Branch("ccElecphi",mTempTreeccElecPhi,"ccElecphi[Nelec]/double");




  //add muons
  mAllData->Branch("Nmuon" ,&mTempTreeNmuon ,"Nmuon/int");  
  mAllData->Branch("MuonE" ,mTempTreeMuonE ,"MuonE[Nmuon]/double");
  mAllData->Branch("MuonEt",mTempTreeMuonEt,"MuonEt[Nmuon]/double");
  mAllData->Branch("Muonpt",mTempTreeMuonPt,"Muonpt[Nmuon]/double");
  mAllData->Branch("Muonpx",mTempTreeMuonPx,"Muonpx[Nmuon]/double");
  mAllData->Branch("Muonpy",mTempTreeMuonPy,"Muonpy[Nmuon]/double");
  mAllData->Branch("Muonpz",mTempTreeMuonPz,"Muonpz[Nmuon]/double");
  mAllData->Branch("Muoneta",mTempTreeMuonEta,"Muoneta[Nmuon]/double");
  mAllData->Branch("Muonphi",mTempTreeMuonPhi,"Muonphi[Nmuon]/double");
  mAllData->Branch("MuonCharge",mTempTreeMuonCharge,"MuonCharge[Nmuon]/double");
  mAllData->Branch("MuonTrkIso",mTempTreeMuonTrkIso,"MuonTrkIso[Nmuon]/double");
  mAllData->Branch("MuonECalIso", mTempTreeMuonECalIso,"MuonECalIso[Nmuon]/double");
  mAllData->Branch("MuonHCalIso", mTempTreeMuonHCalIso ,"MuonHCalIso[Nmuon]/double");
  mAllData->Branch("MuonAllIso",  mTempTreeMuonAllIso ,"MuonAllIso[Nmuon]/double");
  mAllData->Branch("MuonTrkChiNorm",mTempTreeMuonTrkChiNorm  ,"MuonTrkChiNorm[Nmuon]/double");

  mAllData->Branch("MuonIsGlobal",mTempTreeMuonIsGlobal,"mTempTreeMuonIsGlobal[Nmuon]/bool");
  mAllData->Branch("MuonIsStandAlone",mTempTreeMuonIsStandAlone,"mTempTreeMuonIsStandAlone[Nmuon]/bool");
  mAllData->Branch("MuonIsGlobalTight",mTempTreeMuonIsGlobalTight,"mTempTreeMuonIsGlobalTight[Nmuon]/bool");
  mAllData->Branch("MuonIsTMLastStationLoose",mTempTreeMuonIsTMLastStationLoose,"mTempTreeMuonIsTMLastStationLoose[Nmuon]/bool");
  mAllData->Branch("MuonIsTracker",mTempTreeMuonIsTracker,"mTempTreeMuonIsTracker[Nmuon]/bool");
  //MICHELE
  //  mAllData->Branch("MuonId",mTempTreeMuonId,"mTempTreeMuonId[Nmuon]/double");
  mAllData->Branch("MuonCombVx",mTempTreeMuonCombVx,"mTempTreeMuonCombVx[Nmuon]/double");
  mAllData->Branch("MuonCombVy",mTempTreeMuonCombVy,"mTempTreeMuonCombVy[Nmuon]/double");
  mAllData->Branch("MuonCombVz",mTempTreeMuonCombVz,"mTempTreeMuonCombVz[Nmuon]/double");
  mAllData->Branch("MuonCombD0",mTempTreeMuonCombD0,"mTempTreeMuonCombD0[Nmuon]/double");
  mAllData->Branch("MuonCombDz",mTempTreeMuonCombDz,"mTempTreeMuonCombDz[Nmuon]/double");

  mAllData->Branch("MuonStandValidHits",mTempTreeMuonStandValidHits,"mTempTreeMuonStandValidHits[Nmuon]/double");
  mAllData->Branch("MuonStandLostHits",mTempTreeMuonStandLostHits,"mTempTreeMuonStandLostHits[Nmuon]/double");
  mAllData->Branch("MuonStandPt",mTempTreeMuonStandPt,"mTempTreeMuonStandPt[Nmuon]/double");
  mAllData->Branch("MuonStandPz",mTempTreeMuonStandPz,"mTempTreeMuonStandPz[Nmuon]/double");
  mAllData->Branch("MuonStandP",mTempTreeMuonStandP,"mTempTreeMuonStandP[Nmuon]/double");
  mAllData->Branch("MuonStandEta",mTempTreeMuonStandEta,"mTempTreeMuonStandEta[Nmuon]/double");
  mAllData->Branch("MuonStandPhi",mTempTreeMuonStandPhi,"mTempTreeMuonStandPhi[Nmuon]/double");
  mAllData->Branch("MuonStandCharge",mTempTreeMuonStandCharge,"mTempTreeMuonStandCharge[Nmuon]/double");
  mAllData->Branch("MuonStandChi",mTempTreeMuonStandChi,"mTempTreeMuonStandChi[Nmuon]/double");
  mAllData->Branch("MuonStandQOverPError",mTempTreeMuonStandQOverPError,"mTempTreeMuonStandQOverPError[Nmuon]/double");

  mAllData->Branch("MuonTrkValidHits",mTempTreeMuonTrkValidHits,"mTempTreeMuonTrkValidHits[Nmuon]/double");
  mAllData->Branch("MuonTrkLostHits",mTempTreeMuonTrkLostHits,"mTempTreeMuonTrkLostHits[Nmuon]/double");
  mAllData->Branch("MuonTrkPt",mTempTreeMuonTrkPt,"mTempTreeMuonTrkPt[Nmuon]/double");
  mAllData->Branch("MuonTrkPz",mTempTreeMuonTrkPz,"mTempTreeMuonTrkPz[Nmuon]/double");
  mAllData->Branch("MuonTrkP",mTempTreeMuonTrkP,"mTempTreeMuonTrkP[Nmuon]/double");
  mAllData->Branch("MuonTrkEta",mTempTreeMuonTrkEta,"mTempTreeMuonTrkEta[Nmuon]/double");
  mAllData->Branch("MuonTrkPhi",mTempTreeMuonTrkPhi,"mTempTreeMuonTrkPhi[Nmuon]/double");
  mAllData->Branch("MuonTrkCharge",mTempTreeMuonTrkCharge,"mTempTreeMuonTrkCharge[Nmuon]/double");
  mAllData->Branch("MuonTrkChi",mTempTreeMuonTrkChi,"mTempTreeMuonTrkChi[Nmuon]/double");
  mAllData->Branch("MuonTrkQOverPError",mTempTreeMuonTrkQOverPError,"mTempTreeMuonTrkQOverPError[Nmuon]/double"); 
  mAllData->Branch("MuonGenMother",mTempTreeGenMuonMother,"MuonGenMother[Nmuon]/double");
  mAllData->Branch("MuonGenPx",mTempTreeGenMuonPx,"MuonGenPx[Nmuon]/double");
  mAllData->Branch("MuonGenPy",mTempTreeGenMuonPy,"MuonGenPy[Nmuon]/double");
  mAllData->Branch("MuonGenPz",mTempTreeGenMuonPz,"MuonGenPz[Nmuon]/double");

  //PIOPPI


  mAllData->Branch("Nccmuon" ,&mTempTreeNccmuon ,"Nccmuon/int");  
  mAllData->Branch("ccMuonE" ,mTempTreeccMuonE ,"ccMuonE[Nmuon]/double");
  mAllData->Branch("ccMuonEt",mTempTreeccMuonEt,"ccMuonEt[Nmuon]/double");
  mAllData->Branch("ccMuonpt",mTempTreeccMuonPt,"ccMuonpt[Nmuon]/double");
  mAllData->Branch("ccMuonpx",mTempTreeccMuonPx,"ccMuonpx[Nmuon]/double");
  mAllData->Branch("ccMuonpy",mTempTreeccMuonPy,"ccMuonpy[Nmuon]/double");
  mAllData->Branch("ccMuonpz",mTempTreeccMuonPz,"ccMuonpz[Nmuon]/double");
  mAllData->Branch("ccMuoneta",mTempTreeccMuonEta,"ccMuoneta[Nmuon]/double");
  mAllData->Branch("ccMuonphi",mTempTreeccMuonPhi,"ccMuonphi[Nmuon]/double");


  //add taus
  mAllData->Branch("Ntau" ,&mTempTreeNtau ,"Ntau/int");  
  mAllData->Branch("TauE" ,mTempTreeTauE ,"TauE[Ntau]/double");
  mAllData->Branch("TauEt",mTempTreeTauEt,"TauEt[Ntau]/double");
  mAllData->Branch("Taupt",mTempTreeTauPt,"Taupt[Ntau]/double");
  mAllData->Branch("Taupx",mTempTreeTauPx,"Taupx[Ntau]/double");
  mAllData->Branch("Taupy",mTempTreeTauPy,"Taupy[Ntau]/double");
  mAllData->Branch("Taupz",mTempTreeTauPz,"Taupz[Ntau]/double");
  mAllData->Branch("Taueta",mTempTreeTauEta,"Taueta[Ntau]/double");
  mAllData->Branch("Tauphi",mTempTreeTauPhi,"Tauphi[Ntau]/double");
  mAllData->Branch("TauTrkIso",mTempTreeTauTrkIso,"TauTrkIso[Ntau]/double");
  mAllData->Branch("TauECalIso", mTempTreeTauECalIso,"TauECalIso[Ntau]/double");
  mAllData->Branch("TauHCalIso", mTempTreeTauHCalIso ,"TauHCalIso[Ntau]/double");
  mAllData->Branch("TauAllIso",  mTempTreeTauAllIso ,"TauAllIso[Ntau]/double");
  //MICHELE
  mAllData->Branch("TauVx",mTempTreeTauVx,"TauVx[Ntau]/double");
  mAllData->Branch("TauVy",mTempTreeTauVy,"TauVy[Ntau]/double");
  mAllData->Branch("TauVz",mTempTreeTauVz,"TauVz[Ntau]/double");
  mAllData->Branch("TauNTks",mTempTreeTauNTks,"TauNTks[Ntau]/double");
  mAllData->Branch("TauNNeutrals",mTempTreeTauNNeutrals,"TauNNeutrals[Ntau]/double");
  mAllData->Branch("TauNeutralE",mTempTreeTauNeutralE,"TauNeutralE[Ntau]/double");
  mAllData->Branch("TauNeutralHOverHPlusE",mTempTreeTauNeutralHOverHPlusE,"TauNeutralHOverHPlusE[Ntau]/double");

  //TK1
  mAllData->Branch("TauTk1Vx",mTempTreeTauTk1Vx,"TauTk1Vx[Ntau]/double");
  mAllData->Branch("TauTk1Vy",mTempTreeTauTk1Vy,"TauTk1Vy[Ntau]/double");
  mAllData->Branch("TauTk1Vz",mTempTreeTauTk1Vz,"TauTk1Vz[Ntau]/double");
  mAllData->Branch("TauTk1D0",mTempTreeTauTk1D0,"TauTk1D0[Ntau]/double");
  mAllData->Branch("TauTk1Dz",mTempTreeTauTk1Dz,"TauTk1Dz[Ntau]/double");
  mAllData->Branch("TauTk1Pt",mTempTreeTauTk1Pt,"TauTk1Pt[Ntau]/double");
  mAllData->Branch("TauTk1Pz",mTempTreeTauTk1Pz,"TauTk1Pz[Ntau]/double");
  mAllData->Branch("TauTk1Eta",mTempTreeTauTk1Eta,"TauTk1Eta[Ntau]/double");
  mAllData->Branch("TauTk1Phi",mTempTreeTauTk1Phi,"TauTk1Phi[Ntau]/double");
  mAllData->Branch("TauTk1Chi",mTempTreeTauTk1Chi,"TauTk1Chi[Ntau]/double");
  mAllData->Branch("TauTk1Charge",mTempTreeTauTk1Charge,"TauTk1Charge[Ntau]/double");
  mAllData->Branch("TauTk1QOverPError",mTempTreeTauTk1QOverPError,"TauTk1QOverPError[Ntau]/double");
  mAllData->Branch("TauTk1ValidHits",mTempTreeTauTk1ValidHits,"TauTk1ValidHits[Ntau]/double");
  mAllData->Branch("TauTk1LostHits",mTempTreeTauTk1LostHits,"TauTk1LostHits[Ntau]/double");
  mAllData->Branch("TauTk1CaloE",mTempTreeTauTk1CaloE,"TauTk1CaloE[Ntau]/double");
  //TK2
  mAllData->Branch("TauTk2Vx",mTempTreeTauTk2Vx,"TauTk2Vx[Ntau]/double");
  mAllData->Branch("TauTk2Vy",mTempTreeTauTk2Vy,"TauTk2Vy[Ntau]/double");
  mAllData->Branch("TauTk2Vz",mTempTreeTauTk2Vz,"TauTk2Vz[Ntau]/double");
  mAllData->Branch("TauTk2D0",mTempTreeTauTk2D0,"TauTk2D0[Ntau]/double");
  mAllData->Branch("TauTk2Dz",mTempTreeTauTk2Dz,"TauTk2Dz[Ntau]/double");
  mAllData->Branch("TauTk2Pt",mTempTreeTauTk2Pt,"TauTk2Pt[Ntau]/double");
  mAllData->Branch("TauTk2Pz",mTempTreeTauTk2Pz,"TauTk2Pz[Ntau]/double");
  mAllData->Branch("TauTk2Eta",mTempTreeTauTk2Eta,"TauTk2Eta[Ntau]/double");
  mAllData->Branch("TauTk2Phi",mTempTreeTauTk2Phi,"TauTk2Phi[Ntau]/double");
  mAllData->Branch("TauTk2Chi",mTempTreeTauTk2Chi,"TauTk2Chi[Ntau]/double");
  mAllData->Branch("TauTk2Charge",mTempTreeTauTk2Charge,"TauTk2Charge[Ntau]/double");
  mAllData->Branch("TauTk2QOverPError",mTempTreeTauTk2QOverPError,"TauTk2QOverPError[Ntau]/double");
  mAllData->Branch("TauTk2ValidHits",mTempTreeTauTk2ValidHits,"TauTk2ValidHits[Ntau]/double");
  mAllData->Branch("TauTk2LostHits",mTempTreeTauTk2LostHits,"TauTk2LostHits[Ntau]/double");
  mAllData->Branch("TauTk2CaloE",mTempTreeTauTk2CaloE,"TauTk2CaloE[Ntau]/double");
  //TK3
  mAllData->Branch("TauTk3Vx",mTempTreeTauTk3Vx,"TauTk3Vx[Ntau]/double");
  mAllData->Branch("TauTk3Vy",mTempTreeTauTk3Vy,"TauTk3Vy[Ntau]/double");
  mAllData->Branch("TauTk3Vz",mTempTreeTauTk3Vz,"TauTk3Vz[Ntau]/double");
  mAllData->Branch("TauTk3D0",mTempTreeTauTk3D0,"TauTk3D0[Ntau]/double");
  mAllData->Branch("TauTk3Dz",mTempTreeTauTk3Dz,"TauTk3Dz[Ntau]/double");
  mAllData->Branch("TauTk3Pt",mTempTreeTauTk3Pt,"TauTk3Pt[Ntau]/double");
  mAllData->Branch("TauTk3Pz",mTempTreeTauTk3Pz,"TauTk3Pz[Ntau]/double");
  mAllData->Branch("TauTk3Eta",mTempTreeTauTk3Eta,"TauTk3Eta[Ntau]/double");
  mAllData->Branch("TauTk3Phi",mTempTreeTauTk3Phi,"TauTk3Phi[Ntau]/double");
  mAllData->Branch("TauTk3Chi",mTempTreeTauTk3Chi,"TauTk3Chi[Ntau]/double");
  mAllData->Branch("TauTk3Charge",mTempTreeTauTk3Charge,"TauTk3Charge[Ntau]/double");
  mAllData->Branch("TauTk3QOverPError",mTempTreeTauTk3QOverPError,"TauTk3QOverPError[Ntau]/double");
  mAllData->Branch("TauTk3ValidHits",mTempTreeTauTk3ValidHits,"TauTk3ValidHits[Ntau]/double");
  mAllData->Branch("TauTk3LostHits",mTempTreeTauTk3LostHits,"TauTk3LostHits[Ntau]/double");
  mAllData->Branch("TauTk3CaloE",mTempTreeTauTk3CaloE,"TauTk3CaloE[Ntau]/double");

  mAllData->Branch("TauGenPdgId",mTempTreeGenTauPdgId,"TauGenPdgId[Ntau]/double");
  mAllData->Branch("TauGenMother",mTempTreeGenTauMother,"TauGenMother[Ntau]/double");
  mAllData->Branch("TauGenPx",mTempTreeGenTauPx,"TauGenPx[Ntau]/double");
  mAllData->Branch("TauGenPy",mTempTreeGenTauPy,"TauGenPy[Ntau]/double");
  mAllData->Branch("TauGenPz",mTempTreeGenTauPz,"TauGenPz[Ntau]/double");
  //PIOPPI
 
    mAllData->Branch("genN",&length,"genN/int");
    mAllData->Branch("genid",ids,"ids[genN]/int");
    mAllData->Branch("genMother",refs,"refs[genN]/int");
    mAllData->Branch("genPhi",genPhi,"genPhi[genN]/float");
    mAllData->Branch("genE",genE,"genE[genN]/float");
    mAllData->Branch("genPx",genPx,"genPx[genN]/float");
    mAllData->Branch("genPy",genPy,"genPy[genN]/float");
    mAllData->Branch("genPz",genPz,"genPz[genN]/float");
    mAllData->Branch("genEta",genEta,"genEta[genN]/float");
    mAllData->Branch("genStatus",genStatus,"genStatus[genN]/int");
    

  // add test stuff
 
  mAllData->Branch("AlpPtScale" ,&mTempAlpPtScale,"mTempAlpPtScale/double");
  mAllData->Branch("AlpIdTest" ,&mTempAlpIdTest ,"AlpIdTest/int");  

 
  edm::LogInfo("SusyDiJet") << "Ntuple variables " << variables.str();
  
}


//________________________________________________________________________________________
void
SusyDiJetAnalysis::fillPlots( const edm::Event& iEvent, 
    const SelectorDecisions& decisions ) {
  
  // Container array
  float* x = new float[ntuple_->GetNbranches()];
  int ivar = 0; 

  // 1. Event variables
  x[ivar++] = eventWeight_;
  x[ivar++] = processId_;

  // 2. Decision from all selectors
  for ( size_t i=0; i<sequence_.size(); ++i ) x[ivar++] = decisions.decision(i);
  x[ivar++] = decisions.globalDecision();

  // 3. All variables from sequence
  std::vector<double> values = sequence_.values();
  for ( size_t i=0; i<values.size(); ++i ) x[ivar++] = values[i];

  if ( ntuple_->Fill( x ) < 0 ) { // Fill returns number of bytes committed, -1 on error
    edm::LogWarning("SusyDiJet") << "@SUB=fillPlots()" << "Problem filling ntuple";
  }

  delete [] x; // Important! otherwise we'll leak...

}

//________________________________________________________________________________________
// Define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SusyDiJetAnalysis);


