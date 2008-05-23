#ifndef SusyAnalysis_EventSelector_MeffSelector_h_
#define SusyAnalysis_EventSelector_MeffSelector_h_
///
/// MeffSelector
///
/// Calculates the "effective mass" and selects on it
///
/// \author Frederic Ronga - Fri May 23 14:21:31 CEST 2008
///
/// $Id: $
#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
class MeffSelector : public SusyEventSelector {
public:
  MeffSelector (const edm::ParameterSet&);
  virtual bool select(const edm::Event&) const;
  virtual ~MeffSelector () {}

 private:
  edm::InputTag jetTag_; ///< tag for input jet collection
  edm::InputTag metTag_; ///< tag for input MET collection

  double minMeff_;   ///< Cut on minimum M_eff

};
#endif
