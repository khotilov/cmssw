#include "RecoTauTag/HLTProducers/interface/L1HLTJetsMatching.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h" 
#include "DataFormats/Math/interface/Point3D.h"
//
// class decleration
//
using namespace reco;
using namespace std;
using namespace edm;
using namespace l1extra;

L1HLTJetsMatching::L1HLTJetsMatching(const edm::ParameterSet& iConfig)
{
  jetSrc = iConfig.getParameter<InputTag>("JetSrc");
  tauTrigger = iConfig.getParameter<InputTag>("L1TauTrigger");
  mEt_Min = iConfig.getParameter<double>("EtMin");
  
  produces<std::vector<reco::LeafCandidate> >();
  produces<reco::CaloJetCollection>();
}

L1HLTJetsMatching::~L1HLTJetsMatching(){ }

void L1HLTJetsMatching::produce(edm::Event& iEvent, const edm::EventSetup& iES)
{

 using namespace edm;
 using namespace std;
 using namespace reco;
 using namespace trigger;
 using namespace l1extra;

 typedef vector<LeafCandidate> LeafCandidateCollection;

 auto_ptr<LeafCandidateCollection> tauL2LC(new LeafCandidateCollection);
 auto_ptr<CaloJetCollection> tauL2jets(new CaloJetCollection);
 
 double deltaR = 1.0;
 double matchingR = 0.5;
 //Getting HLT jets to be matched
 edm::Handle<edm::View<Candidate> > tauJets;
 if(iEvent.getByLabel( jetSrc, tauJets ))
   {
     
     Handle<trigger::TriggerFilterObjectWithRefs> l1TriggeredTaus;
     if(iEvent.getByLabel(tauTrigger,l1TriggeredTaus)){
	
       
       tauCandRefVec.clear();
       jetCandRefVec.clear();
	
       l1TriggeredTaus->getObjects( trigger::TriggerL1TauJet,tauCandRefVec);
       l1TriggeredTaus->getObjects( trigger::TriggerL1CenJet,jetCandRefVec);
       math::XYZPoint a(0.,0.,0.);
       CaloJet::Specific f;
	
       for( unsigned int iL1Tau=0; iL1Tau <tauCandRefVec.size();iL1Tau++)
	 {  
	   for(unsigned int iJet=0;iJet<tauJets->size();iJet++)
	     {
	       //Find the relative L2TauJets, to see if it has been reconstructed
	       const Candidate &  myJet = (*tauJets)[iJet];
	       deltaR = ROOT::Math::VectorUtil::DeltaR(myJet.p4().Vect(), (tauCandRefVec[iL1Tau]->p4()).Vect());
	       if(deltaR < matchingR ) {
		 LeafCandidate myLC(myJet);
		 CaloJet myCaloJet(myJet.p4(),a,f);
		if(myJet.pt() > mEt_Min) {
		  tauL2LC->push_back(myLC);
		  tauL2jets->push_back(myCaloJet);
		}
		break;
	       }
	     }
	 }  
       
       for(unsigned int iL1Tau=0; iL1Tau <jetCandRefVec.size();iL1Tau++)
	 {  
	   for(unsigned int iJet=0;iJet<tauJets->size();iJet++)
	     {
	       const Candidate & myJet = (*tauJets)[iJet];
	       //Find the relative L2TauJets, to see if it has been reconstructed
	       deltaR = ROOT::Math::VectorUtil::DeltaR(myJet.p4().Vect(), (jetCandRefVec[iL1Tau]->p4()).Vect());
	       if(deltaR < matchingR ) {
		 LeafCandidate myLC(myJet);
		 CaloJet myCaloJet(myJet.p4(),a,f);
		 if(myJet.pt() > mEt_Min) {
		   tauL2LC->push_back(myLC);
		   tauL2jets->push_back(myCaloJet);
		 }
		break;
	       }
	     }
	 }
     }
    
   }
    
 //  cout <<"Size of L2 jets "<<tauL2jets->size()<<endl;

 iEvent.put(tauL2jets);
 iEvent.put(tauL2LC);

}
