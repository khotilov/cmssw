#ifndef _dqm_sistripcommissioningsources_SiStripCommissioningRunTypeFilter_h_
#define _dqm_sistripcommissioningsources_SiStripCommissioningRunTypeFilter_h_

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/SiStripCommon/interface/SiStripEnumsAndStrings.h"
#include <FWCore/ParameterSet/interface/InputTag.h>
#include <DataFormats/SiStripCommon/interface/SiStripEventSummary.h>


//
// class declaration
//

class SiStripCommissioningRunTypeFilter : public edm::EDFilter {
   public:
      explicit SiStripCommissioningRunTypeFilter(const edm::ParameterSet&);
      ~SiStripCommissioningRunTypeFilter() {}

   private:
      virtual void beginJob(const edm::EventSetup&) {}
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() {}
      
      // ----------member data ---------------------------
      edm::InputTag inputModuleLabel_;
      sistrip::RunType runType_;
};

#endif
