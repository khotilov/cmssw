#ifndef TauAnalysis_Core_HistManagerAdapter_h
#define TauAnalysis_Core_HistManagerAdapter_h

/** \class HistManagerAdapter
 *
 * Provide classes derrived from HistManagerBase with an EDAnalyzer interface 
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: HistManagerAdapter.h,v 1.1 2009/01/22 16:30:02 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

template<typename T>
class HistManagerAdapter : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit HistManagerAdapter(const edm::ParameterSet& cfg) : 
    histManager_( cfg ) {
  }
    
  // destructor
  virtual ~HistManagerAdapter() {}
    
 private:
  void beginJob(const edm::EventSetup& es) { histManager_.bookHistograms(es); }
  void analyze(const edm::Event& evt, const edm::EventSetup& es) { histManager_.fillHistograms(evt, es); }
  void endJob() {}
  
  T histManager_;
};

#endif   

