#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "FWCore/Framework/interface/Selector.h"
#include <iostream>

using namespace std;


/** \class HcalDigiDump
      
$Date: 2006/02/10 21:45:20 $
$Revision: 1.7 $
\author J. Mans - Minnesota
*/
class HcalDigiDump : public edm::EDAnalyzer {
public:
  explicit HcalDigiDump(edm::ParameterSet const& conf);
  virtual void analyze(edm::Event const& e, edm::EventSetup const& c);
};


HcalDigiDump::HcalDigiDump(edm::ParameterSet const& conf) {
}

void HcalDigiDump::analyze(edm::Event const& e, edm::EventSetup const& c) {
  std::vector<edm::Handle<HBHEDigiCollection> > hbhe;
  std::vector<edm::Handle<HODigiCollection> > ho;
  std::vector<edm::Handle<HFDigiCollection> > hf;
  std::vector<edm::Handle<ZDCDigiCollection> > zdc;
  std::vector<edm::Handle<HcalCalibDigiCollection> > hc;
  std::vector<edm::Handle<HcalTrigPrimDigiCollection> > htp;
  std::vector<edm::Handle<HcalHistogramDigiCollection> > hh;  

  try {
    e.getManyByType(hbhe);
    std::vector<edm::Handle<HBHEDigiCollection> >::iterator i;
    for (i=hbhe.begin(); i!=hbhe.end(); i++) {
      const HBHEDigiCollection& c=*(*i);
      
      for (HBHEDigiCollection::const_iterator j=c.begin(); j!=c.end(); j++)
	cout << *j << std::endl;
    }
  } catch (...) {
    cout << "No HB/HE Digis." << endl;
  }
  
  try {
    e.getManyByType(hf);
    std::vector<edm::Handle<HFDigiCollection> >::iterator i;
    for (i=hf.begin(); i!=hf.end(); i++) {
      const HFDigiCollection& c=*(*i);
      
      for (HFDigiCollection::const_iterator j=c.begin(); j!=c.end(); j++)
	cout << *j << std::endl;
    }
  } catch (...) {
    cout << "No HF Digis." << endl;
  }
  
  try {
    e.getManyByType(ho);
    std::vector<edm::Handle<HODigiCollection> >::iterator i;
    for (i=ho.begin(); i!=ho.end(); i++) {
      const HODigiCollection& c=*(*i);
      
      for (HODigiCollection::const_iterator j=c.begin(); j!=c.end(); j++)
	cout << *j << std::endl;

    }
  } catch (...) {
    cout << "No HO Digis." << endl;
  }

  try {
    e.getManyByType(htp);
    std::vector<edm::Handle<HcalTrigPrimDigiCollection> >::iterator i;
    for (i=htp.begin(); i!=htp.end(); i++) {
      const HcalTrigPrimDigiCollection& c=*(*i);
      
      for (HcalTrigPrimDigiCollection::const_iterator j=c.begin(); j!=c.end(); j++)
	cout << *j << std::endl;

    }
  } catch (...) {
    cout << "No HCAL Trigger Primitive Digis." << endl;
  }

  try {
    e.getManyByType(hc);
    std::vector<edm::Handle<HcalCalibDigiCollection> >::iterator i;
    for (i=hc.begin(); i!=hc.end(); i++) {
      const HcalCalibDigiCollection& c=*(*i);
      
      for (HcalCalibDigiCollection::const_iterator j=c.begin(); j!=c.end(); j++)
	cout << *j << std::endl;
    }
  } catch (...) {
  }

  try {
    e.getManyByType(zdc);
    std::vector<edm::Handle<ZDCDigiCollection> >::iterator i;
    for (i=zdc.begin(); i!=zdc.end(); i++) {
      const ZDCDigiCollection& c=*(*i);
      
      for (ZDCDigiCollection::const_iterator j=c.begin(); j!=c.end(); j++)
	cout << *j << std::endl;
    }
  } catch (...) {
  }


  try {
    e.getManyByType(hh);
    std::vector<edm::Handle<HcalHistogramDigiCollection> >::iterator i;
    for (i=hh.begin(); i!=hh.end(); i++) {
      const HcalHistogramDigiCollection& c=*(*i);
      
      for (HcalHistogramDigiCollection::const_iterator j=c.begin(); j!=c.end(); j++)
	cout << *j << std::endl;

    }
  } catch (...) {
  }
  
  cout << endl;    
}

#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(HcalDigiDump)

