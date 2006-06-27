#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalSourcePositionData.h"
#include "FWCore/Framework/interface/Selector.h"
#include <iostream>

using namespace std;

namespace cms {

  /** \class HcalRecHitDump
      
  $Date: 2006/01/25 02:10:24 $
  $Revision: 1.4 $
  \author J. Mans - Minnesota
  */
  class HcalRecHitDump : public edm::EDAnalyzer {
  public:
    explicit HcalRecHitDump(edm::ParameterSet const& conf);
    virtual void analyze(edm::Event const& e, edm::EventSetup const& c);
  };

  HcalRecHitDump::HcalRecHitDump(edm::ParameterSet const& conf) {
  }
  
  void HcalRecHitDump::analyze(edm::Event const& e, edm::EventSetup const& c) {
    std::vector<edm::Handle<HORecHitCollection> > ho;
    std::vector<edm::Handle<HFRecHitCollection> > hf;
    edm::Handle<HcalSourcePositionData> spd;

    try {
      std::vector<edm::Handle<HBHERecHitCollection> > colls;
      e.getManyByType(colls);
      std::vector<edm::Handle<HBHERecHitCollection> >::iterator i;
      for (i=colls.begin(); i!=colls.end(); i++) {
	for (HBHERecHitCollection::const_iterator j=(*i)->begin(); j!=(*i)->end(); j++) 
	  cout << *j << std::endl;
      }
    } catch (...) {
      cout << "No HB/HE RecHits." << endl;
    }
    
    try {
      std::vector<edm::Handle<HFRecHitCollection> > colls;
      e.getManyByType(colls);
      std::vector<edm::Handle<HFRecHitCollection> >::iterator i;
      for (i=colls.begin(); i!=colls.end(); i++) {
	for (HFRecHitCollection::const_iterator j=(*i)->begin(); j!=(*i)->end(); j++) 
	  cout << *j << std::endl;
      }
    } catch (...) {
      cout << "No HF RecHits." << endl;
    }
    
    try {
      std::vector<edm::Handle<HORecHitCollection> > colls;
      e.getManyByType(colls);
      std::vector<edm::Handle<HORecHitCollection> >::iterator i;
      for (i=colls.begin(); i!=colls.end(); i++) {
	for (HORecHitCollection::const_iterator j=(*i)->begin(); j!=(*i)->end(); j++) 
	  cout << *j << std::endl;
      }
    } catch (...) {
      cout << "No HO RecHits." << endl;
    }

    try {
      std::vector<edm::Handle<HcalCalibRecHitCollection> > colls;
      e.getManyByType(colls);
      std::vector<edm::Handle<HcalCalibRecHitCollection> >::iterator i;
      for (i=colls.begin(); i!=colls.end(); i++) {
	for (HcalCalibRecHitCollection::const_iterator j=(*i)->begin(); j!=(*i)->end(); j++) 
	  cout << *j << std::endl;
      }
    } catch (...) {
      //      cout << "No Calib RecHits." << endl;
    }

    try {
      std::vector<edm::Handle<ZDCRecHitCollection> > colls;
      e.getManyByType(colls);
      std::vector<edm::Handle<ZDCRecHitCollection> >::iterator i;
      for (i=colls.begin(); i!=colls.end(); i++) {
	for (ZDCRecHitCollection::const_iterator j=(*i)->begin(); j!=(*i)->end(); j++) 
	  cout << *j << std::endl;
      }
    } catch (...) {
      //      cout << "No ZDC RecHits." << endl;
    }


    try {
      e.getByType(spd);
      cout << *spd << std::endl;
    } catch (...) {
//      cout << "No Source Position Data" << endl;
    }

    cout << endl;    
  }
}

#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace cms;

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(HcalRecHitDump)

