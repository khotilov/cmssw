#ifndef _dqm_sistripcommissioningsources_SiStripCommissioningSeedFilter_h_
#define _dqm_sistripcommissioningsources_SiStripCommissioningSeedFilter_h_

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <FWCore/ParameterSet/interface/InputTag.h>


//
// class declaration
//

class SiStripCommissioningSeedFilter : public edm::EDFilter {

   public:

      explicit SiStripCommissioningSeedFilter(const edm::ParameterSet&);
      ~SiStripCommissioningSeedFilter() {}

   private:

      virtual bool filter(edm::Event&, const edm::EventSetup&);
      
      // ----------member data ---------------------------
      edm::InputTag inputModuleLabel_;

};

#endif
