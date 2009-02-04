#ifndef TauAnalysis_Core_EventDumpBase_h
#define TauAnalysis_Core_EventDumpBase_h

/** \class EventDumpBase
 *
 * Base-class for print-out of event level information
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: EventDumpBase.h,v 1.1 2009/01/22 16:30:02 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>
#include <string>
#include <iostream>

class EventDumpBase
{
 public:
  // constructor 
  explicit EventDumpBase(const edm::ParameterSet&);
  
  // destructor
  virtual ~EventDumpBase();

  // methods for print-out of event level information
  typedef std::vector<std::pair<std::string, bool> > filterResults_type;
  virtual void analyze(const edm::Event&, const edm::EventSetup&, 
		       const filterResults_type&, const filterResults_type&, double) const;

 protected:
  virtual void print(const edm::Event&, const edm::EventSetup&, 
		     const filterResults_type&, const filterResults_type&, double) const = 0;

  std::ostream* outputStream_;

  int cfgError_;

 private:
//--- configuration parameters
  std::string output_; 
  bool isOutputFile_;

  typedef std::vector<std::string> vstring;
  std::map<std::string, vstring> triggerConditions_;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<EventDumpBase* (const edm::ParameterSet&)> EventDumpPluginFactory;

#endif       

