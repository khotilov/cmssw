#ifndef TauAnalysis_RecoTools_TriggerResultEventSelector_h
#define TauAnalysis_RecoTools_TriggerResultEventSelector_h

/** \class TriggerResultEventSelector
 *
 * Selects an event if a at least one of the HLT triggers
 * specified by the "triggerPaths" configuration parameter triggered the event
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: TriggerResultEventSelector.h,v 1.3 2010/04/28 14:51:15 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/UtilAlgos/interface/EventSelectorBase.h"

class TriggerResultEventSelector : public EventSelectorBase
{
 public:
  // constructor 
  explicit TriggerResultEventSelector(const edm::ParameterSet&);
  
  bool operator()(edm::Event&, const edm::EventSetup&);
 
 private:
  // source collection label
  edm::InputTag src_;
 
  // list of HLT trigger path to be checked
  typedef std::vector<std::string> vstring;
  vstring triggerPaths_;
};

#endif

