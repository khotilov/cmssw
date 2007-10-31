#ifndef L1COMPARATOR_H
#define L1COMPARATOR_H

/*\class L1Comparator
 *\description L1 trigger data|emulation comparison and validation
 *\author Nuno Leonardo (CERN)
 *\date 07.02
 */

// common/system includes
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// l1 dataformats, d|e record includes
#include "L1Trigger/HardwareValidation/interface/DEtrait.h"

// comparator template
#include "L1Trigger/HardwareValidation/interface/DEcompare.h"

// extra
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

template <class T> class DEcompare;

using dedefs::DEnsys;

class L1Comparator : public edm::EDProducer {

public:

  explicit L1Comparator(const edm::ParameterSet&);
  ~L1Comparator();
  
private:

  virtual void beginJob(const edm::EventSetup&);
  virtual void produce (edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  template <class T> 
    void process( T const*, T const*, const int);
  template <class T> 
    void process(const edm::Handle<T> data, const edm::Handle<T> emul, 
		 const int sys) {
    if(data.isValid()&&emul.isValid())
      process(data.product(),emul.product(),sys);
  }

  // gt, fedraw, extra
  bool compareCollections(edm::Handle<L1GlobalTriggerReadoutRecord>   data, 
			  edm::Handle<L1GlobalTriggerReadoutRecord>   emul);
  bool compareCollections(edm::Handle<L1GlobalTriggerEvmReadoutRecord>data, 
			  edm::Handle<L1GlobalTriggerEvmReadoutRecord>emul);
  bool compareCollections(edm::Handle<L1GlobalTriggerObjectMapRecord> data, 
			  edm::Handle<L1GlobalTriggerObjectMapRecord> emul);
  bool compareFedRawCollections(edm::Handle<FEDRawDataCollection>     data, 
				edm::Handle<FEDRawDataCollection>     emul, int fedid);
  template <class T> bool CompareCollections(edm::Handle<T> data, edm::Handle<T> emul);
  template <class T> bool dumpCandidate(const T& dt, const T& em, std::ostream& s);

  int verbose() {return verbose_;}

 private:

  int nevt_;
  int evtNum_;
  int runNum_;
  int verbose_;
  bool dumpEvent_;

  edm::InputTag m_DEsource[DEnsys][4];
  bool m_doSys[DEnsys];
  std::string m_dumpFileName;
  std::ofstream m_dumpFile;
  int m_dumpMode;
  bool m_match;
  bool DEmatchEvt[DEnsys]; 
  int DEncand[DEnsys][2];
  L1DEDigiCollection m_dedigis;

  int m_fedId;
  edm::InputTag m_FEDsource[2];

};

#endif
