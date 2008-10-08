#include "FastSimDataFormats/External/interface/FastL1BitInfo.h"

#include "ElectroWeakAnalysis/TauTriggerEfficiency/interface/L1TauEfficiencyAnalyzer.h"
#include "PhysicsTools/Utilities/interface/deltaR.h"

// Default constructor
L1TauEfficiencyAnalyzer::L1TauEfficiencyAnalyzer()
{}


L1TauEfficiencyAnalyzer::~L1TauEfficiencyAnalyzer(){
  /*
	cout << endl;
	cout << "Events analyzed " << nEvents << endl;
	cout << endl;
  */
}

void L1TauEfficiencyAnalyzer::Setup(const edm::ParameterSet& iConfig,TTree *trigtree)
{
  L1extraTauJetSource = iConfig.getParameter<edm::InputTag>("L1extraTauJetSource");
  L1bitInfoSource = iConfig.getParameter<edm::InputTag>("L1bitInfoSource");
  jetMatchingCone = iConfig.getParameter<double>("JetMatchingCone");
  //rootFile_ = iConfig.getParameter<std::string>("outputFileName");
  nEvents = 0; nSelectedEvents = 0;

  l1tree = trigtree;

  // Setup branches
  l1tree->Branch("L1JetPt", &jetPt, "L1JetPt/F");
  l1tree->Branch("L1JetEt", &jetPt, "L1JetEt/F");
  l1tree->Branch("L1JetEta", &jetEta, "L1JetEta/F");
  l1tree->Branch("L1JetPhi", &jetPhi, "L1JetPhi/F");
  l1tree->Branch("L1TauVeto", &hasTauVeto, "L1TauVeto/B");
  l1tree->Branch("L1EmTauVeto", &hasEmTauVeto, "L1EmTauVeto/B");
  l1tree->Branch("L1HadTauVeto", &hasHadTauVeto, "L1HadTauVeto/B");
  l1tree->Branch("L1IsolationVeto", &hasIsolationVeto, "L1IsolationVeto/B");
  l1tree->Branch("L1SumEtBelowThreshold", &hasSumEtBelowThres, "L1SumEtBelowThrehold/B");
  l1tree->Branch("L1MaxEt", &hasMaxEt, "L1MaxEt/B");
  l1tree->Branch("L1Soft", &hasSoft, "L1Soft/B");
  l1tree->Branch("L1Hard", &hasHard, "L1Hard/B");
  l1tree->Branch("hasMatchedL1Jet", &hasL1Jet, "hasMatchedL1Jet/B");
}

void L1TauEfficiencyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup ){
        Handle<L1JetParticleCollection> l1TauHandle;
	nEvents++;

	try{
          iEvent.getByLabel(L1extraTauJetSource,l1TauHandle);
	}catch(...) {;}


	Handle<CaloTauCollection> caloTauHandle;
	try{
          iEvent.getByLabel("IdentifiedTaus",caloTauHandle);
        }catch(...) {;}

	if(caloTauHandle.isValid()){
	  const CaloTauCollection & caloTaus = *(caloTauHandle.product());

	  LogDebug("L1TauEfficiency") << "calotau collection size " << caloTaus.size() << endl;

	  CaloTauCollection::const_iterator iTau;
          for(iTau = caloTaus.begin(); iTau != caloTaus.end(); ++iTau){
            //if(L1TauFound(iTau->p4())){
            //}
	  }
	}

        Handle<PFTauCollection> pfTauHandle;
        try{
          iEvent.getByLabel("IdentifiedTaus",pfTauHandle);
        }catch(...) {;}

        if(pfTauHandle.isValid()){
          const PFTauCollection & pfTaus = *(pfTauHandle.product());

          LogDebug("L1TauEfficiency") << "pftau collection size " << pfTaus.size() << endl;

          PFTauCollection::const_iterator iTau;
          for(iTau = pfTaus.begin(); iTau != pfTaus.end(); ++iTau){
            jetPt = 0;
            jetEta = 0;
            jetPhi = 0;
            hasL1Jet = 0;

            //if(L1TauFound(iTau->p4())){
            //}
            //l1tree->Fill();
          }
        }
}

void L1TauEfficiencyAnalyzer::fill(const edm::Event& iEvent, const reco::PFTau& tau) {
  // Reset variables
  jetPt = 0.0;
  jetEta = 0.0;
  jetPhi = 0.0;
  hasL1Jet = 0;
  hasTauVeto = 0;
  hasIsolationVeto = 0;

  // Get data from event 
  Handle<FastL1BitInfoCollection> bitInfos;
  iEvent.getByLabel(L1bitInfoSource, bitInfos);

  Handle<L1JetParticleCollection> l1TauHandle;
  iEvent.getByLabel(L1extraTauJetSource, l1TauHandle);

  // Process L1 triggered taus
  if(l1TauHandle.isValid()) {  
    const L1JetParticleCollection & l1Taus = *(l1TauHandle.product());
    L1JetParticleCollection::const_iterator iTau;

    float minDR = 99999999.;
    for(iTau = l1Taus.begin(); iTau != l1Taus.end(); ++iTau){
      double DR = deltaR(iTau->eta(), iTau->phi(), tau.eta(), tau.phi());
      if(DR < jetMatchingCone && DR < minDR) {
        minDR = DR;
        jetPt = iTau->pt();
        jetEta = iTau->eta();
        jetPhi = iTau->phi();
        hasL1Jet = 1;

      }
    }
  }

  // Process bit info
  if(bitInfos.isValid()) {
    float minDR = 99999999.;
    for(FastL1BitInfoCollection::const_iterator bitInfo = bitInfos->begin(); bitInfo != bitInfos->end(); ++bitInfo) {
      double DR = deltaR(bitInfo->getEta(), bitInfo->getPhi(), tau.eta(), tau.phi());
      if(DR < jetMatchingCone && DR < minDR) {
        minDR = DR;

        //std::cout << "Foo dr" << DR << std::endl;

        hasTauVeto = bitInfo->getTauVeto() ? 1 : 0;
        hasEmTauVeto = bitInfo->getEmTauVeto() ? 1 : 0;
        hasHadTauVeto = bitInfo->getHadTauVeto() ? 1 : 0;
        hasIsolationVeto = bitInfo->getIsolationVeto() ? 1 : 0;
        hasSumEtBelowThres = bitInfo->getSumEtBelowThres() ? 1 : 0;
        hasMaxEt = bitInfo->getMaxEt() ? 1 : 0;
        hasSoft = bitInfo->getSoft() ? 1 : 0;
        hasHard = bitInfo->getHard() ? 1 : 0;
      }
    }
  }
} 

void L1TauEfficiencyAnalyzer::beginJob(const edm::EventSetup& iSetup){}

void L1TauEfficiencyAnalyzer::endJob(){
        LogInfo("L1TauEfficiency") << "Events analyzed: " << nEvents << endl;
        //l1file->Write();
}

#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(L1TauEfficiencyAnalyzer);


