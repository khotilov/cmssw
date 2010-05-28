#ifndef SusyAnalysis_EventSelector_EMFractionSelector_h_
#define SusyAnalysis_EventSelector_EMFractionSelector_h_
///
/// Selection based on EM fraction.
///
/// Calculates the total EM fraction of all jets
///
/// $Id: EMFractionSelector.h,v 1.2 2008/05/23 15:48:21 fronga Exp $

// user include files
#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

class EMFractionSelector: public SusyEventSelector {
public:
   EMFractionSelector(const edm::ParameterSet&);
   virtual bool select(const edm::Event&) const;
   virtual ~EMFractionSelector() {
   }

private:
   edm::InputTag jetTag_; ///< tag for input collection
   float minFraction_; ///< lower EM fraction cut
   float maxEta_; ///< upper |eta| cut on all jets
   float minPt_; ///< lower pt cut on all jets

};
#endif
