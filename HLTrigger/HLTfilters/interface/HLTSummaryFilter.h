#ifndef HLTSummaryFilter_h
#define HLTSummaryFilter_h

/** \class HLTSummaryFilter
 *
 *  
 *  This class is an HLTFilter (-> EDFilter) implementing a smart HLT
 *  trigger cut, specified as a string such as "pt>15 && -3<eta<3",
 *  for objects in the TriggerSummaryAOD product, allowing to cut on
 *  variables relating to their 4-momentum representation
 *
 *  $Date: 2009/07/08 01:13:08 $
 *  $Revision: 1.1 $
 *
 *  \author Martin Grunewald
 *
 */

#include "HLTrigger/HLTcore/interface/HLTFilter.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include<string>

//
// class declaration
//

class HLTSummaryFilter : public HLTFilter {

   public:

      explicit HLTSummaryFilter(const edm::ParameterSet&);
      ~HLTSummaryFilter();
      virtual bool filter(edm::Event&, const edm::EventSetup&);

   private:
      edm::InputTag summaryTag_; // input tag identifying TriggerSummaryAOD
      edm::InputTag memberTag_;  // which packed-up collection or filter
      std::string   cut_;        // smart cut
      int           min_N_;      // number of objects passing cuts required

      StringCutObjectSelector<trigger::TriggerObject> select_; // smart selector
};

#endif //HLTSummaryFilter_h
