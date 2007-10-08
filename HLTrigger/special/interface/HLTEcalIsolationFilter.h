#ifndef HLTEcalIsolationFilter_h
#define HLTEcalIsolationFilter_h

#include "HLTrigger/HLTcore/interface/HLTFilter.h"

class HLTEcalIsolationFilter : public HLTFilter {

   public:
      explicit HLTEcalIsolationFilter(const edm::ParameterSet&);
      ~HLTEcalIsolationFilter();
      virtual bool filter(edm::Event&, const edm::EventSetup&);

   private:
      edm::InputTag candTag_; 
      double maxennearby; 
      double minen;        
      int maxhitout;
      int maxhitin;
      double maxenin;
      double maxenout;
      double maxetacand;

};

#endif 
