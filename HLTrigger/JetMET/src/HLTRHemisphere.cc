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

#include<vector>

//
// constructors and destructor
//
HLTRHemisphere::HLTRHemisphere(const edm::ParameterSet& iConfig) :
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

//
// member functions
//

// ------------ method called to produce the data  ------------
bool 
HLTRHemisphere::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
   int N_comb(1); // compute the number of combinations of jets possible
  for(unsigned int i = 0; i < JETS.size(); i++){
    N_comb *= 2;                
  }
  //Make the hemispheres
  XYZTLorentzVector j1,j2;
  double M_min = 9999999999.0;
  double dHT_min = 99999999.0;
  int j_count;
  for(int i=0;i<N_comb;i++){       
    XYZTLorentzVector j_temp1, j_temp2;
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while(j_count > 0){
      if(itemp/j_count == 1){
	j_temp1 += JETS.at(count).p4();
      } else {
	j_temp2 += JETS.at(count).p4();
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }
    double M_temp = j_temp1.M2()+j_temp2.M2();
    if(M_temp < M_min){
      M_min = M_temp;
      j1= j_temp1;
      j2= j_temp2;
    }
    double dHT_temp = fabs(j_temp1.E()-j_temp2.E());
    if(dHT_temp < dHT_min){
      dHT_min = dHT_temp;
    }
  }

  Hemispheres->push_back(j1);
  Hemispheres->push_back(j2);

  iEvent.put(Hemispheres);
  return true;
}

DEFINE_FWK_MODULE(HLTRHemisphere);
