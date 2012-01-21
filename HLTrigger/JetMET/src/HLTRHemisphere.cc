#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/JetMET/interface/HLTRHemisphere.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"


#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TVector3.h"
#include "TLorentzVector.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include<vector>

//
// constructors and destructor
//
HLTRHemisphere::HLTRHemisphere(const edm::ParameterSet& iConfig) : HLTFilter(iConfig),
  inputTag_    (iConfig.getParameter<edm::InputTag>("inputTag")),
  min_Jet_Pt_  (iConfig.getParameter<double>       ("minJetPt" )),
  max_Eta_     (iConfig.getParameter<double>       ("maxEta" )),
  max_NJ_      (iConfig.getParameter<int>          ("maxNJ" )),
  accNJJets_   (iConfig.getParameter<bool>         ("acceptNJ" ))
{
   LogDebug("") << "Input/minJetPt/maxEta/maxNJ/acceptNJ : "
		<< inputTag_.encode() << " "
		<< min_Jet_Pt_ << "/"
		<< max_Eta_ << "/"
		<< max_NJ_ << "/"
		<< accNJJets_ << ".";

   //register your products
   produces<std::vector<math::XYZTLorentzVector> >();
}

HLTRHemisphere::~HLTRHemisphere()
{
}

void
HLTRHemisphere::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("inputTag",edm::InputTag("hltMCJetCorJetIcone5HF07"));
  desc.add<double>("minJetPt",30.0);
  desc.add<double>("maxEta",3.0);
  desc.add<int>("maxNJ",7);
  desc.add<bool>("acceptNJ",true);
  descriptions.add("hltRHemisphere",desc);
}

//
// member functions
//

// ------------ method called to produce the data  ------------
bool 
HLTRHemisphere::hltFilter(edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs & filterproduct)
{
   using namespace std;
   using namespace edm;
   using namespace reco;
   using namespace math;

   // get hold of collection of objects
   Handle<CaloJetCollection> jets;
   iEvent.getByLabel (inputTag_,jets);

   // The output Collection
   std::auto_ptr<vector<math::XYZTLorentzVector> > Hemispheres(new vector<math::XYZTLorentzVector> );

   // look at all objects, check cuts and add to filter object
   int n(0);
   reco::CaloJetCollection JETS;
   CaloJetCollection::const_iterator i ( jets->begin() );
   for (unsigned int i=0; i<jets->size(); i++) {
     if(fabs(jets->at(i).eta()) < max_Eta_ && jets->at(i).pt() >= min_Jet_Pt_){
       JETS.push_back(jets->at(i));
       n++;
     }
   }

  if(n<2){
    return false; //need at least 2 jets to build the hemispheres
  }

  if(n>max_NJ_ && max_NJ_!=-1){
    iEvent.put(Hemispheres);
    return accNJJets_; // 
  }
  unsigned int N_comb(1); // compute the number of combinations of jets possible
  for(unsigned int i = 0; i < JETS.size(); i++){
    N_comb *= 2;                
  }
  //Make the hemispheres
  XYZTLorentzVector j1R(0.1, 0., 0., 0.1);
  XYZTLorentzVector j2R(0.1, 0., 0., 0.1);
  double M_minR = 9999999999.0;
  unsigned int j_count;
  for (unsigned int i = 0; i < N_comb; i++) {       
    XYZTLorentzVector j_temp1, j_temp2;
    unsigned int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while (j_count > 0) {
      if (itemp/j_count == 1){
	j_temp1 += JETS.at(count).p4();
      } else {
	j_temp2 += JETS.at(count).p4();
      }
      itemp -= j_count * (itemp/j_count);
      j_count /= 2;
      count++;
    }
    double M_temp = j_temp1.M2() + j_temp2.M2();
    if (M_temp < M_minR) {
      M_minR = M_temp;
      j1R = j_temp1; 
      j2R = j_temp2; 
    }
  }

  Hemispheres->push_back(j1R);
  Hemispheres->push_back(j2R);

  iEvent.put(Hemispheres);
  return true;
}

DEFINE_FWK_MODULE(HLTRHemisphere);
