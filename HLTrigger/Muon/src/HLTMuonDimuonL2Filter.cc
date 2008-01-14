/** \class HLTMuonDimuonL2Filter
 *
 * See header file for documentation
 *  
 *  \author J. Alcaraz, P. Garcia, M.Vander Donckt
 *
 */


#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeedCollection.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"

#include "HLTrigger/Muon/interface/HLTMuonDimuonL2Filter.h"

using namespace edm;
using namespace std;
using namespace reco;
using namespace trigger;
using namespace l1extra;
//
// constructors and destructor
//
HLTMuonDimuonL2Filter::HLTMuonDimuonL2Filter(const edm::ParameterSet& iConfig) :
   candTag_     (iConfig.getParameter< edm::InputTag > ("CandTag")),
   previousCandTag_   (iConfig.getParameter<InputTag > ("PreviousCandTag")),
   linksTag_   (iConfig.getParameter<InputTag > ("SeedTag")),
   fast_Accept_ (iConfig.getParameter<bool> ("FastAccept")),
   max_Eta_     (iConfig.getParameter<double> ("MaxEta")),
   min_Nhits_   (iConfig.getParameter<int> ("MinNhits")),
   max_Dr_      (iConfig.getParameter<double> ("MaxDr")),
   max_Dz_      (iConfig.getParameter<double> ("MaxDz")),
   chargeOpt_   (iConfig.getParameter<int> ("ChargeOpt")),
   min_PtPair_  (iConfig.getParameter<double> ("MinPtPair")),
   min_PtMax_   (iConfig.getParameter<double> ("MinPtMax")),
   min_PtMin_   (iConfig.getParameter<double> ("MinPtMin")),
   min_InvMass_ (iConfig.getParameter<double> ("MinInvMass")),
   max_InvMass_ (iConfig.getParameter<double> ("MaxInvMass")),
   min_Acop_    (iConfig.getParameter<double> ("MinAcop")),
   max_Acop_    (iConfig.getParameter<double> ("MaxAcop")),
   min_PtBalance_ (iConfig.getParameter<double> ("MinPtBalance")),
   max_PtBalance_ (iConfig.getParameter<double> ("MaxPtBalance")),
   nsigma_Pt_   (iConfig.getParameter<double> ("NSigmaPt")) 
{

   LogDebug("HLTMuonDimuonL2Filter")
      << " CandTag/MinN/MaxEta/MinNhits/MaxDr/MaxDz/MinPt1/MinPt2/MinInvMass/MaxInvMass/MinAcop/MaxAcop/MinPtBalance/MaxPtBalance/NSigmaPt : " 
      << candTag_.encode()
      << " " << fast_Accept_
      << " " << max_Eta_
      << " " << min_Nhits_
      << " " << max_Dr_
      << " " << max_Dz_
      << " " << chargeOpt_ << " " << min_PtPair_
      << " " << min_PtMax_ << " " << min_PtMin_
      << " " << min_InvMass_ << " " << max_InvMass_
      << " " << min_Acop_ << " " << max_Acop_
      << " " << min_PtBalance_ << " " << max_PtBalance_
      << " " << nsigma_Pt_;

   //register your products
   produces<trigger::TriggerFilterObjectWithRefs>();
}

HLTMuonDimuonL2Filter::~HLTMuonDimuonL2Filter()
{
}

//
// member functions
//

// ------------ method called to produce the data  ------------
bool
HLTMuonDimuonL2Filter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   double const MuMass = 0.106;
   double const MuMass2 = MuMass*MuMass;
   // All HLT filters must create and fill an HLT filter object,
   // recording any reconstructed physics objects satisfying (or not)
   // this HLT filter, and place it in the Event.

   // The filter object
   auto_ptr<TriggerFilterObjectWithRefs>
     filterproduct (new TriggerFilterObjectWithRefs(path(),module()));

   // Ref to Candidate object to be recorded in filter object
   RecoChargedCandidateRef ref1;
   RecoChargedCandidateRef ref2;

   // get hold of trks
   Handle<RecoChargedCandidateCollection> mucands;
   iEvent.getByLabel (candTag_,mucands);
   Handle<TriggerFilterObjectWithRefs> previousLevelCands;
   iEvent.getByLabel (previousCandTag_,previousLevelCands);

   Handle<L2MuonTrajectorySeedCollection> museeds; 
   vector<L1MuonParticleRef> vl1cands;
   iEvent.getByLabel (linksTag_,museeds);
   previousLevelCands->getObjects(TriggerL1Mu,vl1cands);


   // look at all mucands,  check cuts and add to filter object
   int n = 0;
   double e1,e2;
   Particle::LorentzVector p,p1,p2;

   RecoChargedCandidateCollection::const_iterator cand1;
   RecoChargedCandidateCollection::const_iterator cand2;
   for (cand1=mucands->begin(); cand1!=mucands->end(); cand1++) {
      TrackRef tk1 = cand1->get<TrackRef>();
      // eta cut
      LogDebug("HLTMuonDimuonL2Filter") << " 1st muon in loop: q*pt= "
            << tk1->charge()*tk1->pt() << ", eta= " << tk1->eta() << ", hits= " << tk1->numberOfValidHits();
      // find the L1 Particle corresponding to the L2 Track
      if (!triggeredByLevel1(tk1,museeds,vl1cands)) continue;
 
      if (fabs(tk1->eta())>max_Eta_) continue;

      // cut on number of hits
      if (tk1->numberOfValidHits()<min_Nhits_) continue;

      //dr cut
      if (fabs(tk1->d0())>max_Dr_) continue;

      //dz cut
      if (fabs(tk1->dz())>max_Dz_) continue;

      // Pt threshold cut
      double pt1 = tk1->pt();
      double err1 = tk1->error(0);
      double abspar1 = fabs(tk1->parameter(0));
      double ptLx1 = pt1;
      // convert 50% efficiency threshold to 90% efficiency threshold
      if (abspar1>0) ptLx1 += nsigma_Pt_*err1/abspar1*pt1;
      LogDebug("HLTMuonDimuonL2Filter") << " ... 1st muon in loop, pt1= "
            << pt1 << ", ptLx1= " << ptLx1;

      cand2 = cand1; cand2++;
      for (; cand2!=mucands->end(); cand2++) {
            TrackRef tk2 = cand2->get<TrackRef>();

            // eta cut
            LogDebug("HLTMuonDimuonL2Filter") << " 2nd muon in loop: q*pt= " << tk2->charge()*tk2->pt() << ", eta= " << tk2->eta() << ", hits= " << tk2->numberOfValidHits() << ", d0= " << tk2->d0();
	    if (!triggeredByLevel1(tk2,museeds,vl1cands)) continue;

            if (fabs(tk2->eta())>max_Eta_) continue;

            // cut on number of hits
            if (tk2->numberOfValidHits()<min_Nhits_) continue;

            //dr cut
            if (fabs(tk2->d0())>max_Dr_) continue;

            //dz cut
            if (fabs(tk2->dz())>max_Dz_) continue;

            // Pt threshold cut
            double pt2 = tk2->pt();
            double err2 = tk2->error(0);
            double abspar2 = fabs(tk2->parameter(0));
            double ptLx2 = pt2;
            // convert 50% efficiency threshold to 90% efficiency threshold
            if (abspar2>0) ptLx2 += nsigma_Pt_*err2/abspar2*pt2;
            LogDebug("HLTMuonDimuonL2Filter") << " ... 2nd muon in loop, pt2= "
                  << pt2 << ", ptLx2= " << ptLx2;

            if (ptLx1>ptLx2) {
                  if (ptLx1<min_PtMax_) continue;
                  if (ptLx2<min_PtMin_) continue;
            } else {
                  if (ptLx1<min_PtMax_) continue;
                  if (ptLx2<min_PtMin_) continue;
            }

            if (chargeOpt_<0) {
                  if (tk1->charge()*tk2->charge()>0) continue;
            } else if (chargeOpt_>0) {
                  if (tk1->charge()*tk2->charge()<0) continue;
            }

            // Acoplanarity
            double acop = fabs(tk1->phi()-tk2->phi());
            if (acop>M_PI) acop = 2*M_PI - acop;
            acop = M_PI - acop;
            LogDebug("HLTMuonDimuonL2Filter") << " ... 1-2 acop= " << acop;
            if (acop<min_Acop_) continue;
            if (acop>max_Acop_) continue;

            // Pt balance
            double ptbalance = fabs(tk1->pt()-tk2->pt());
            if (ptbalance<min_PtBalance_) continue;
            if (ptbalance>max_PtBalance_) continue;

            // Combined dimuon system
            e1 = sqrt(tk1->momentum().Mag2()+MuMass2);
            e2 = sqrt(tk2->momentum().Mag2()+MuMass2);
            p1 = Particle::LorentzVector(tk1->px(),tk1->py(),tk1->pz(),e1);
            p2 = Particle::LorentzVector(tk2->px(),tk2->py(),tk2->pz(),e2);
            p = p1+p2;

            double pt12 = p.pt();
            LogDebug("HLTMuonDimuonL2Filter") << " ... 1-2 pt12= " << pt12;
            if (pt12<min_PtPair_) continue;

            double invmass = abs(p.mass());
         // if (invmass>0) invmass = sqrt(invmass); else invmass = 0;
            LogDebug("HLTMuonDimuonL2Filter") << " ... 1-2 invmass= " << invmass;
            if (invmass<min_InvMass_) continue;
            if (invmass>max_InvMass_) continue;

            // Add this pair
            n++;
            LogDebug("HLTMuonDimuonL2Filter") << " Track1 passing filter: pt= " << tk1->pt() << ", eta: " << tk1->eta();
            LogDebug("HLTMuonDimuonL2Filter") << " Track2 passing filter: pt= " << tk2->pt() << ", eta: " << tk2->eta();
            LogDebug("HLTMuonDimuonL2Filter") << " Invmass= " << invmass;

            bool i1done = false;
            bool i2done = false;
            vector<RecoChargedCandidateRef> vref;
	    filterproduct->getObjects(TriggerMuon,vref);
            for (unsigned int i=0; i<vref.size(); i++) {
	      RecoChargedCandidateRef candref =  RecoChargedCandidateRef(vref[i]);
	      TrackRef tktmp = candref->get<TrackRef>();
	      if (tktmp==tk1) {
                        i1done = true;
	      } else if (tktmp==tk2) {
		i2done = true;
	      }
	      if (i1done && i2done) break;
            }
	    
            if (!i1done) { 
	      ref1=RecoChargedCandidateRef( Ref<RecoChargedCandidateCollection> (mucands,distance(mucands->begin(), cand1)));
	      filterproduct->addObject(TriggerMuon,ref1);
            }
            if (!i2done) { 
	      ref2=RecoChargedCandidateRef( Ref<RecoChargedCandidateCollection> (mucands,distance(mucands->begin(),cand2 )));
	    filterproduct->addObject(TriggerMuon,ref2);
            }

            if (fast_Accept_) break;
      }

   }

   // filter decision
   const bool accept (n >= 1);

   // put filter object into the Event
   iEvent.put(filterproduct);

   LogDebug("HLTMuonDimuonL2Filter") << " >>>>> Result of HLTMuonDimuonL2Filter is "<< accept << ", number of muon pairs passing thresholds= " << n; 

   return accept;
}


bool
HLTMuonDimuonL2Filter::triggeredByLevel1(TrackRef& tk,Handle<L2MuonTrajectorySeedCollection> &museeds,vector<L1MuonParticleRef>& vcands)
{
  bool ok=false;
  L1MuonParticleRef l1ref;
  RefToBase<TrajectorySeed> seed=tk->seedRef();   
  L2MuonTrajectorySeedCollection::const_iterator l2seed;
  for (l2seed=museeds->begin() ; l2seed != museeds->end();++l2seed){
    if ( (l2seed->startingState()).detId() == (seed->startingState()).detId() ) {
      l1ref = l2seed->l1Particle();
      break;
    }
  }
  for (unsigned int i=0; i<vcands.size(); i++) {
    L1MuonParticleRef candref =  L1MuonParticleRef(vcands[i]);
    if ( l1ref == candref ) {
      ok=true;
      LogTrace("HLTMuonL2PreFilter") << "The L1 stub triggered";
      break;
    }
  } 
  return ok;
}
