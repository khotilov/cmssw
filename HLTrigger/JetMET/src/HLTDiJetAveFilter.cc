/** \class HLTDiJetAveFilter
 *
 *
 *  \author Leonard Apanasevich
 *
 */

#include "HLTrigger/JetMET/interface/HLTDiJetAveFilter.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Math/interface/deltaPhi.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

//
// constructors and destructor
//
HLTDiJetAveFilter::HLTDiJetAveFilter(const edm::ParameterSet& iConfig) : HLTFilter(iConfig) 
{
   inputJetTag_ = iConfig.getParameter< edm::InputTag > ("inputJetTag");
   minPtAve_    = iConfig.getParameter<double> ("minPtAve"); 
   minPtJet3_   = iConfig.getParameter<double> ("minPtJet3"); 
   minDphi_     = iConfig.getParameter<double> ("minDphi"); 
}

HLTDiJetAveFilter::~HLTDiJetAveFilter(){}

void
HLTDiJetAveFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("inputJetTag",edm::InputTag("hltIterativeCone5CaloJets"));
  desc.add<bool>("saveTags",false);
  desc.add<double>("minPtAve",100.0);
  desc.add<double>("minPtJet3",99999.0);
  desc.add<double>("minDphi",-1.0);
  descriptions.add("hltDiJetAveFilter",desc);
}

// ------------ method called to produce the data  ------------
bool
HLTDiJetAveFilter::hltFilter(edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs & filterproduct)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  // The filter object
  if (saveTags()) filterproduct.addCollectionTag(inputJetTag_);

  Handle<CaloJetCollection> recocalojets;
  iEvent.getByLabel(inputJetTag_,recocalojets);

  // look at all candidates,  check cuts and add to filter object
  int n(0);

  if(recocalojets->size() > 1){
    // events with two or more jets

    double ptjet1=0., ptjet2=0.,ptjet3=0.;
    double phijet1=0.,phijet2=0;
    int countjets =0;

    int nmax=1;
    if (recocalojets->size() > 2) nmax=2;

    CaloJetRef JetRef1,JetRef2;

    for (CaloJetCollection::const_iterator recocalojet = recocalojets->begin(); 
	 recocalojet<=(recocalojets->begin()+nmax); ++recocalojet) {
      
      if(countjets==0) {
	ptjet1 = recocalojet->pt();
	phijet1 = recocalojet->phi();
	JetRef1 = CaloJetRef(recocalojets,distance(recocalojets->begin(),recocalojet));
      }
      if(countjets==1) {
	ptjet2 = recocalojet->pt();
	phijet2 = recocalojet->phi();
	JetRef2 = CaloJetRef(recocalojets,distance(recocalojets->begin(),recocalojet));
      }
      if(countjets==2) {
	ptjet3 = recocalojet->pt();
      }
      ++countjets;
    }
    
    double PtAve=(ptjet1 + ptjet2) / 2.;
    double Dphi = fabs(deltaPhi(phijet1,phijet2));

    if( PtAve>minPtAve_  && ptjet3<minPtJet3_ && Dphi>minDphi_){
      filterproduct.addObject(TriggerJet,JetRef1);
      filterproduct.addObject(TriggerJet,JetRef2);
      ++n;
    }
    
  } // events with two or more jets
  
  
  
  // filter decision
  bool accept(n>=1);
  
  return accept;
}
