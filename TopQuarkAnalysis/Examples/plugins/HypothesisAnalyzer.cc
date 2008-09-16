#include "DataFormats/Candidate/interface/Candidate.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TopQuarkAnalysis/TopTools/interface/TtSemiLepEvtPartons.h"
#include "TopQuarkAnalysis/Examples/plugins/HypothesisAnalyzer.h"


HypothesisAnalyzer::HypothesisAnalyzer(const edm::ParameterSet& cfg):
  semiLepEvt_ (cfg.getParameter<edm::InputTag>("semiLepEvent")),
  hypoKey_ (cfg.getParameter<edm::InputTag>("hypoKey"  ))
{
}

void
HypothesisAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<TtSemiLeptonicEvent> semiLepEvt;
  evt.getByLabel(semiLepEvt_, semiLepEvt);

  edm::Handle<int> hypoKeyHandle;
  evt.getByLabel(hypoKey_, hypoKeyHandle);
  TtSemiLeptonicEvent::HypoKey& hypoKey = (TtSemiLeptonicEvent::HypoKey&) *hypoKeyHandle;

  if( !semiLepEvt->isHypoAvailable(hypoKey) ){
    edm::LogInfo ( "NonValidHyp" ) << "Hypothesis not available for this event";
    return;
  }
  if( !semiLepEvt->isHypoValid(hypoKey) ){
    edm::LogInfo ( "NonValidHyp" ) << "Hypothesis not valid for this event";
    return;
  }
  
  const reco::Candidate* hadTop = semiLepEvt->hadronicTop(hypoKey);
  const reco::Candidate* hadW   = semiLepEvt->hadronicW  (hypoKey);
  const reco::Candidate* lepTop = semiLepEvt->leptonicTop(hypoKey);
  const reco::Candidate* lepW   = semiLepEvt->leptonicW  (hypoKey);

  if(hadTop && hadW && lepTop && lepW){
    hadWPt_    ->Fill( hadW->pt()    );
    hadWMass_  ->Fill( hadW->mass()  );
    hadTopPt_  ->Fill( hadTop->pt()  );
    hadTopMass_->Fill( hadTop->mass());
    
    lepWPt_    ->Fill( lepW->pt()    );
    lepWMass_  ->Fill( lepW->mass()  );
    lepTopPt_  ->Fill( lepTop->pt()  );
    lepTopMass_->Fill( lepTop->mass());
  }
}

void 
HypothesisAnalyzer::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  hadWPt_     = fs->make<TH1F>("hadWPt",     "p_{t} (W_{had}) [GeV]", 100,  0.,  500.);
  hadWMass_   = fs->make<TH1F>("hadWMass",   "M (W_{had}) [GeV]"    ,  50,  0. , 150.);
  hadTopPt_   = fs->make<TH1F>("hadTopPt",   "p_{t} (t_{had}) [GeV]", 100,  0. , 500.);
  hadTopMass_ = fs->make<TH1F>("hadTopMass", "M (t_{had}) [GeV]",      50, 50. , 250.);

  lepWPt_     = fs->make<TH1F>("lepWPt",     "p_{t} (W_{lep}) [GeV]", 100,  0.,  500.);
  lepWMass_   = fs->make<TH1F>("lepWMass",   "M (W_{lep}) [GeV]"    ,  50,  0. , 150.);
  lepTopPt_   = fs->make<TH1F>("lepTopPt",   "p_{t} (t_{lep}) [GeV]", 100,  0. , 500.);
  lepTopMass_ = fs->make<TH1F>("lepTopMass", "M (t_{lep}) [GeV]",      50, 50. , 250.);
}

void
HypothesisAnalyzer::endJob() 
{
}
