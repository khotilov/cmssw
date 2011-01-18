#ifndef HcalBaseDQClient_GUARD_H
#define HcalBaseDQClient_GUARD_H

#include <string>
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DQM/HcalMonitorTasks/interface/HcalEtaPhiHists.h"
#include "DQM/HcalMonitorClient/interface/HcalHistoUtils.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

/*
 * \file HcalBaseDQClient.h
 * 
 * $Date: 2010/03/20 20:55:43 $
 * $Revision: 1.1.2.6 $
 * \author J. Temple
 * \brief Hcal Monitor Client base class
 * based on code in EcalBarrelMonitorClient/interface/EBClient.h
 */


class HcalBaseDQClient
{
 public:
  HcalBaseDQClient(){name_="HcalBaseDQClient";subdir_="HcalInfo";badChannelStatusMask_=0;enoughevents_=true;minerrorrate_=0;minevents_=0;};
  HcalBaseDQClient(std::string s, const edm::ParameterSet& ps);
  ~HcalBaseDQClient(void){}
  
  // Overload these functions with client-specific instructions
  virtual void beginJob(void);
  virtual void beginRun(void)          {}
  virtual void setup(void)             {}

  virtual void analyze(void)           {enoughevents_=true;} // fill new histograms
  virtual void calculateProblems(void) {} // update/fill ProblemCell histograms
  
  virtual void endRun(void)            {}
  virtual void endJob(void)            {}
  virtual void cleanup(void)           {}
  
  virtual bool hasErrors_Temp(void)    {return false;};
  virtual bool hasWarnings_Temp(void)  {return false;};
  virtual bool hasOther_Temp(void)     {return false;};
  virtual bool test_enabled(void)      {return false;};

  virtual void htmlOutput(std::string htmlDir);
  virtual void setStatusMap(std::map<HcalDetId, unsigned int>& map);
  virtual void updateChannelStatus(std::map<HcalDetId, unsigned int>& myqual){};     
  
  virtual bool validHtmlOutput();

  inline void setEventSetup(const edm::EventSetup& es) 
    { c = &(es);  }
  const edm::EventSetup *c;
  std::string name(){return name_;};
  // make these private, with public accessors, at some point?
  std::string name_;
  std::string prefixME_;
  std::string subdir_;
  bool cloneME_;
  bool enableCleanup_;
  int debug_;
  int badChannelStatusMask_; 
  bool validHtmlOutput_;

  bool testenabled_;
  int minevents_; // minimum number of events for test to pass
  double minerrorrate_;

  MonitorElement* ProblemCells;
  EtaPhiHists* ProblemCellsByDepth;
  std::vector<std::string> problemnames_;

  std::map<HcalDetId, unsigned int> badstatusmap;
  DQMStore* dqmStore_;
  bool enoughevents_;
}; // class HcalBaseDQClient



#endif
