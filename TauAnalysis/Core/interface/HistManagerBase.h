#ifndef TauAnalysis_Core_HistManagerBase_h
#define TauAnalysis_Core_HistManagerBase_h

/** \class HistManagerBase
 *
 * Base-class for histogram booking and filling in physics analyses
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: HistManagerBase.h,v 1.1 2009/01/22 16:30:02 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class HistManagerBase
{
 public:
  // constructor 
  explicit HistManagerBase() {}
  
  // destructor
  virtual ~HistManagerBase() {}

  // methods for booking and filling of histograms
  virtual void bookHistograms(const edm::EventSetup&) = 0;
  virtual void fillHistograms(const edm::Event&, const edm::EventSetup&) = 0;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<HistManagerBase* (const edm::ParameterSet&)> HistManagerPluginFactory;

#endif  

