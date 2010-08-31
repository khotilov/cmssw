#ifndef SusyAnalysis_MetEventSelector_h_
#define SusyAnalysis_MetEventSelector_h_
/// Trivial example for a MET selector.
/// To be modified for analysis!
///
/// Original author: W. Adam, 10/4/08
///
/// $Id: MetEventSelector.h,v 1.9 2010/05/28 08:01:48 csander Exp $

//if uncorrNumber =
//0: bare uncorrected MET
//1: MET uncorrected for JES
//2: MET uncorrected for Muons
//>2: totally corrected MET
//T. Rommerskirchen /27/06/08

// system include files
#include <memory>

// user include files
#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include <vector>
#include <string>

class MetEventSelector: public SusyEventSelector {
public:
   MetEventSelector(const edm::ParameterSet&);
   virtual bool select(const edm::Event&) const;
   virtual ~MetEventSelector() {
   }

private:
   edm::InputTag metTag_; ///< tag for input collection
   float minMet_; ///< lower cut on MET
   float maxMet_; ///< upper cut on MET

   pat::MET::UncorrectionType uncorrType_;
};
#endif
