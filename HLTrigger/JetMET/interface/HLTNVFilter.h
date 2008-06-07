#ifndef HLTNVFilter_h
#define HLTNVFilter_h

/** \class HLTNVFilter
 *
 *  \author Dominique J. Mangeol
 *
 */

#include "HLTrigger/HLTcore/interface/HLTFilter.h"

//
// class declaration
//

class HLTNVFilter : public HLTFilter {

   public:
      explicit HLTNVFilter(const edm::ParameterSet&);
      ~HLTNVFilter();
      virtual bool filter(edm::Event&, const edm::EventSetup&);

   private:
      edm::InputTag inputJetTag_; // input tag identifying jets
      edm::InputTag inputMETTag_; // input tag identifying for MET
      bool saveTags_;             // whether to save these tags

      double minEtjet1_;
      double minEtjet2_;
      double minNV_;
};

#endif //HLTNVFilter_h
