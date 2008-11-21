#ifndef HLTHcalNoiseFilter_h
#define HLTHcalNoiseFilter_h

#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

class HLTHcalNoiseFilter : public HLTFilter {
   public:
      explicit HLTHcalNoiseFilter(const edm::ParameterSet&);
      ~HLTHcalNoiseFilter();
      virtual bool filter(edm::Event&, const edm::EventSetup&);

   private:
      edm::InputTag JetSource_;
      edm::InputTag MetSource_;
      edm::InputTag TowerSource_;
      bool useMet_;
      bool useJet_;
      double MetCut_;
      double JetMinE_;
      double JetHCALminEnergyFraction_;
      int nAnomalousEvents;
      int nEvents;
};

#endif
