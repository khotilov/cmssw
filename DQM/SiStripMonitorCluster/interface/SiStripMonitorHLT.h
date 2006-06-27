#ifndef SiStripMonitorCluster_SiStripMonitorHLT_h
#define SiStripMonitorCluster_SiStripMonitorHLT_h
// -*- C++ -*-
//
// Package:     SiStripMonitorCluster
// Class  :     SiStripMonitorHLT



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/MonitorElement.h"

class SiStripMonitorHLT : public edm::EDAnalyzer {
   public:
      explicit SiStripMonitorHLT(const edm::ParameterSet&);
      ~SiStripMonitorHLT(){};

      virtual void analyze(const edm::Event&, const edm::EventSetup&);
       virtual void beginJob(edm::EventSetup const&) ;
       virtual void endJob() ;

   private:
       DaqMonitorBEInterface* dbe_;
       edm::ParameterSet conf_;
       MonitorElement * HLTDecision;
       MonitorElement * ClusterCharge;
};

#endif
