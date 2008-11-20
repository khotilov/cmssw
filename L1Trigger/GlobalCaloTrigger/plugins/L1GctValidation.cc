#include "L1Trigger/GlobalCaloTrigger/plugins/L1GctValidation.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/L1TObjects/interface/L1CaloEtScale.h"
#include "CondFormats/L1TObjects/interface/L1GctJetEtCalibrationFunction.h"
#include "CondFormats/DataRecord/interface/L1JetEtScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1GctJetCalibFunRcd.h"

#include <math.h>

L1GctValidation::L1GctValidation(const edm::ParameterSet& iConfig) :
  m_energy_tag(iConfig.getUntrackedParameter<edm::InputTag>("inputTag", edm::InputTag("gctDigis"))),
  m_missHt_tag(iConfig.getUntrackedParameter<edm::InputTag>("missHtTag",edm::InputTag("gctDigis:missingHt")))
{
}


L1GctValidation::~L1GctValidation()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
L1GctValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  // Get the scales from the event setup
  ESHandle< L1GctJetEtCalibrationFunction > calibFun ;
  iSetup.get< L1GctJetCalibFunRcd >().get( calibFun ) ; // which record?
  ESHandle< L1CaloEtScale > etScale ;
  iSetup.get< L1JetEtScaleRcd >().get( etScale ) ; // which record?

  double lsbForEt = etScale.product()->linearLsb();
  double lsbForHt = calibFun.product()->getHtScaleLSB();

  // Get the Gct energy sums from the event
  Handle< L1GctEtTotalCollection > sumEtColl ;
  iEvent.getByLabel( m_energy_tag, sumEtColl ) ;
  Handle< L1GctEtHadCollection >   sumHtColl ;
  iEvent.getByLabel( m_energy_tag, sumHtColl ) ;
  Handle< L1GctEtMissCollection >  missEtColl ;
  iEvent.getByLabel( m_energy_tag, missEtColl ) ;
  Handle< L1GctEtMissCollection >  missHtColl ;
  iEvent.getByLabel( m_missHt_tag, missHtColl ) ;

  // Get the internal jet data from the event (for checking Ht)
  Handle < L1GctInternJetDataCollection > internalJetsColl;
  iEvent.getByLabel( m_energy_tag, internalJetsColl ) ;

  double etTot = 0.0;
  for (L1GctEtTotalCollection::const_iterator jbx=sumEtColl->begin(); jbx!=sumEtColl->end(); jbx++) {
    if (jbx->bx()==0) { etTot  = static_cast<double>(jbx->et()); }
  } 

  double etHad  = 0.0;
  for (L1GctEtHadCollection::const_iterator jbx=sumHtColl->begin(); jbx!=sumHtColl->end(); jbx++) {
    if (jbx->bx()==0) { etHad  = static_cast<double>(jbx->et()); }
  }

  double etMiss = 0.0;
  double etMAng = 0.0;
  for (L1GctEtMissCollection::const_iterator jbx=missEtColl->begin(); jbx!=missEtColl->end(); jbx++) {
    if (jbx->bx()==0) {
      etMiss = static_cast<double>(jbx->et());
      int phibin = jbx->phi();
      if (phibin>=36) phibin -= 72;
      double etMPhi = static_cast<double>(phibin);

      etMAng = (etMPhi+0.5)*M_PI/36.;
    }
  }

  double htMiss = 0.0;
  double htMAng = 0.0;
  for (L1GctEtMissCollection::const_iterator jbx=missHtColl->begin(); jbx!=missHtColl->end(); jbx++) {
    if (jbx->bx()==0) {
      htMiss = static_cast<double>(jbx->et());
      int phibin = jbx->phi();
      if (phibin>=36) phibin -= 72;
      double htMPhi = static_cast<double>(phibin);

      htMAng = (htMPhi+0.5)*M_PI/36.;
    }
  }

  double htFromJets = 0.0;
  double hxFromJets = 0.0;
  double hyFromJets = 0.0;
  for (L1GctInternJetDataCollection::const_iterator jet=internalJetsColl->begin(); jet!=internalJetsColl->end(); jet++) {
    if (jet->bx()==0 && !jet->empty()) {
      double jetEt = static_cast<double>(jet->et());
      int phibin = jet->regionId().iphi();
      if (phibin>=9) phibin -= 18;
      double jetAng = (static_cast<double>(phibin)+0.5)*M_PI/9.;
      htFromJets += jetEt;
      hxFromJets += jetEt*cos(jetAng);
      hyFromJets += jetEt*sin(jetAng);
    }
  }

  theSumEtInLsb->Fill(etTot);
  theSumHtInLsb->Fill(etHad);
  theMissEtInLsb->Fill(etMiss);
  theMissHtInLsb->Fill(htMiss);
  theSumEtInGeV->Fill(etTot*lsbForEt);
  theSumHtInGeV->Fill(etHad*lsbForHt);
  theMissEtInGeV->Fill(etMiss*lsbForEt);
  theMissEtAngle->Fill(etMAng);
  theMissEtVector->Fill(etMiss*lsbForEt*cos(etMAng),etMiss*lsbForEt*sin(etMAng));
  theMissHtInGeV->Fill(htMiss*lsbForHt*8.);
  theMissHtAngle->Fill(htMAng);
  theMissHtVector->Fill(htMiss*lsbForHt*cos(htMAng)*8.,htMiss*lsbForHt*sin(htMAng)*8.);

  theMissEtVsMissHt->Fill(etMiss*lsbForEt, htMiss*lsbForHt*8);
  theMissEtVsMissHtAngle->Fill(etMAng, htMAng);

  theHtVsInternalJetsSum->Fill(etHad*lsbForHt, htFromJets*lsbForHt);
  theMissHtVsInternalJetsSum->Fill(htMiss*lsbForHt*8, sqrt(hxFromJets*hxFromJets + hyFromJets*hyFromJets)*lsbForHt);
  theMissHxVsInternalJetsSum->Fill(htMiss*lsbForHt*cos(htMAng)*8, hxFromJets*lsbForHt);
  theMissHyVsInternalJetsSum->Fill(htMiss*lsbForHt*sin(htMAng)*8, hyFromJets*lsbForHt);

  // Get jet counts from the event
  Handle< L1GctJetCountsCollection > jetCountColl ;
  iEvent.getByLabel( m_energy_tag, jetCountColl ) ;

  for (L1GctJetCountsCollection::const_iterator jbx=jetCountColl->begin(); jbx!=jetCountColl->end(); jbx++) {
    if (jbx->bx()==0) {
      for (unsigned jc=0; jc<L1GctJetCounts::MAX_TOTAL_COUNTS; jc++) {
	theJetCounts.at(jc)->Fill(jbx->count(jc));
      }
    }

  }

  // Get minbias trigger quantities from HF
  Handle<L1GctHFRingEtSumsCollection> HFEtSumsColl;
  Handle<L1GctHFBitCountsCollection>  HFCountsColl;
  iEvent.getByLabel( m_energy_tag, HFEtSumsColl ) ;
  iEvent.getByLabel( m_energy_tag, HFCountsColl ) ;

  for (L1GctHFRingEtSumsCollection::const_iterator es = HFEtSumsColl->begin(); es != HFEtSumsColl->end(); es++) {
    if (es->bx()==0) {
      theHfRing0EtSumPositiveEta->Fill(es->etSum(0));  
      theHfRing0EtSumNegativeEta->Fill(es->etSum(1));  
      theHfRing1EtSumPositiveEta->Fill(es->etSum(2));  
      theHfRing1EtSumNegativeEta->Fill(es->etSum(3));
    }
  }  

  for (L1GctHFBitCountsCollection::const_iterator bc = HFCountsColl->begin(); bc != HFCountsColl->end(); bc++) {
    if (bc->bx()==0) {
      theHfRing0CountPositiveEta->Fill(bc->bitCount(0));  
      theHfRing0CountNegativeEta->Fill(bc->bitCount(1));  
      theHfRing1CountPositiveEta->Fill(bc->bitCount(2));  
      theHfRing1CountNegativeEta->Fill(bc->bitCount(3));
    }
  }

}

// ------------ method called once each job just before starting event loop  ------------
void 
L1GctValidation::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;

  TFileDirectory dir0 = fs->mkdir("L1GctEtSums");

  theSumEtInLsb   = dir0.make<TH1F>("SumEtInLsb",   "Total Et (GCT units)",
                                    128, 0., 2048.); 
  theSumHtInLsb   = dir0.make<TH1F>("SumHtInLsb",   "Total Ht (GCT units)",
                                    128, 0., 2048.); 
  theMissEtInLsb  = dir0.make<TH1F>("MissEtInLsb",  "Missing Et magnitude (GCT units)",
                                    128, 0., 1024.); 
  theMissHtInLsb  = dir0.make<TH1F>("MissHtInLsb",  "Missing Ht magnitude (GCT units)",
                                    64, 0., 63.); 
  theSumEtInGeV   = dir0.make<TH1F>("SumEtInGeV",   "Total Et (in GeV)",
                                    100, 0., 1000.); 
  theSumHtInGeV   = dir0.make<TH1F>("SumHtInGeV",   "Total Ht (in GeV)",
                                    100, 0., 1000.); 
  theMissEtInGeV  = dir0.make<TH1F>("MissEtInGeV",  "Missing Et magnitude (in GeV)",
                                    100, 0., 500.); 
  theMissEtAngle  = dir0.make<TH1F>("MissEtAngle",  "Missing Et angle",
                                    72, -M_PI, M_PI);
  theMissEtVector = dir0.make<TH2F>("MissEtVector", "Missing Ex vs Missing Ey",
                                    100, -100., 100., 100, -100., 100.); 
  theMissHtInGeV  = dir0.make<TH1F>("MissHtInGeV",  "Missing Ht magnitude (in GeV)",
                                    100, 0., 500.); 
  theMissHtAngle  = dir0.make<TH1F>("MissHtAngle",  "Missing Ht angle",
                                    72, -M_PI, M_PI);
  theMissHtVector = dir0.make<TH2F>("MissHtVector", "Missing Hx vs Missing Hy",
                                    100, -100., 100., 100, -100., 100.); 
  theMissEtVsMissHt = dir0.make<TH2F>("MissEtVsMissHt", "Missing Et vs Missing Ht",
				      100, 0., 500., 100, 0., 500.);
  theMissEtVsMissHtAngle = dir0.make<TH2F>("MissEtVsMissHtAngle", "Angle correlation Missing Et vs Missing Ht",
					   72, -M_PI, M_PI, 72, -M_PI, M_PI);

  theHtVsInternalJetsSum     = dir0.make<TH2F>("HtVsInternalJetsSum", "Ht vs scalar sum of jet Et values (in GCT units)",
					       128, 0., 2048., 128, 0., 2048.);
  theMissHtVsInternalJetsSum = dir0.make<TH2F>("MissHtVsInternalJetsSum", "Missing Ht vs vector sum of jet Et values (in GeV)",
					       128, 0., 512., 128, 0., 512.);
  theMissHxVsInternalJetsSum = dir0.make<TH2F>("MissHxVsInternalJetsSum", "Missing Ht x component vs sum of jet Et values (in GeV)",
					       128, -256., 256., 128, -256., 256.);
  theMissHyVsInternalJetsSum = dir0.make<TH2F>("MissHyVsInternalJetsSum", "Missing Ht y component vs sum of jet Et values (in GeV)",
					       128, -256., 256., 128, -256., 256.);


  TFileDirectory dir1 = fs->mkdir("L1GctHfSumsAndJetCounts");

  for (unsigned jc=0; jc<L1GctJetCounts::MAX_TOTAL_COUNTS; jc++) {
    std::stringstream ss;
    std::string title;
    std::string header;
    ss << "JetCount#" << jc;
    ss >> title;
    ss << "Jet Count number " << jc;
    if (jc== 6 || jc== 7) { ss << " (Hf tower count)"; }
    if (jc== 8 || jc== 9) { ss << " (Hf Et0 sum MSB)"; }
    if (jc==10 || jc==11) { ss << " (Hf Et1 sum MSB)"; }
    ss >> header;
    theJetCounts.push_back(dir1.make<TH1F>(title.c_str(), header.c_str(), 32, 0., 32.));
  }

  // Minimum bias triggers from Hf inner rings
  theHfRing0EtSumPositiveEta = dir1.make<TH1F>("HfRing0EtSumPositiveEta", "Hf Inner Ring0 Et eta+",
                                               60, 0., 30.);
  theHfRing0EtSumNegativeEta = dir1.make<TH1F>("HfRing0EtSumNegativeEta", "Hf Inner Ring0 Et eta-",
                                               60, 0., 30.);
  theHfRing1EtSumPositiveEta = dir1.make<TH1F>("HfRing1EtSumPositiveEta", "Hf Inner Ring1 Et eta+",
                                               60, 0., 30.);
  theHfRing1EtSumNegativeEta = dir1.make<TH1F>("HfRing1EtSumNegativeEta", "Hf Inner Ring1 Et eta-",
                                               60, 0., 30.);
  theHfRing0CountPositiveEta = dir1.make<TH1F>("HfRing0CountPositiveEta", "Hf Threshold bits Ring0 eta+",
                                               20, 0., 20.);
  theHfRing0CountNegativeEta = dir1.make<TH1F>("HfRing0CountNegativeEta", "Hf Threshold bits Ring0 eta-",
                                               20, 0., 20.);
  theHfRing1CountPositiveEta = dir1.make<TH1F>("HfRing1CountPositiveEta", "Hf Threshold bits Ring1 eta+",
                                               20, 0., 20.);
  theHfRing1CountNegativeEta = dir1.make<TH1F>("HfRing1CountNegativeEta", "Hf Threshold bits Ring1 eta-",
                                               20, 0., 20.);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1GctValidation::endJob() {
}

DEFINE_ANOTHER_FWK_MODULE(L1GctValidation);

