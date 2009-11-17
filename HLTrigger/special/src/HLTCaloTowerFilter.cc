// -*- C++ -*-
//
// Package:    HLTCaloTowerFilter
// Class:      HLTCaloTowerFilter
// 
/**\class HLTCaloTowerFilter HLTCaloTowerFilter.cc Work/HLTCaloTowerFilter/src/HLTCaloTowerFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
// Original Author:  Yen-Jie Lee
//         Created:  Wed Nov 13 16:12:29 CEST 2009


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"

//
// class declaration
//
class HLTCaloTowerFilter : public HLTFilter {
public:
    explicit HLTCaloTowerFilter(const edm::ParameterSet&);
    ~HLTCaloTowerFilter();
    
private:
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    
    // ----------member data ---------------------------
    edm::InputTag inputTag_;    // input tag identifying product
    double        min_Pt_;      // pt threshold in GeV 
    double        max_Eta_;     // eta range (symmetric)
    unsigned int  min_N_;       // number of objects passing cuts required

};

//
// constructors and destructor
//
HLTCaloTowerFilter::HLTCaloTowerFilter(const edm::ParameterSet& iConfig) :
  inputTag_ (iConfig.getParameter<edm::InputTag>("inputTag")),
  min_Pt_   (iConfig.getParameter<double>       ("MinPt"   )),
  max_Eta_  (iConfig.getParameter<double>       ("MaxEta"  )),
  min_N_    (iConfig.getParameter<unsigned int> ("MinN"    ))
{
}


HLTCaloTowerFilter::~HLTCaloTowerFilter()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HLTCaloTowerFilter::filter(edm::Event& event, const edm::EventSetup& setup) {
  using namespace std;
  using namespace edm;
  using namespace reco;

  // All HLT filters must create and fill an HLT filter object,
  // recording any reconstructed physics objects satisfying (or not)
  // this HLT filter, and place it in the Event.

  // get hold of collection of objects
  Handle<CaloTowerCollection> caloTowers;
  event.getByLabel(inputTag_, caloTowers);

  // look at all objects, check cuts and add to filter object
  unsigned int n = 0;
  for (CaloTowerCollection::const_iterator i = caloTowers->begin(); i != caloTowers->end(); ++i) {
    if ( (i->pt() >= min_Pt_) and ( (max_Eta_ < 0.0) or (abs(i->eta()) <= max_Eta_) ) )
      ++n;
  }

  // filter decision
  return (n >= min_N_);
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTCaloTowerFilter);
