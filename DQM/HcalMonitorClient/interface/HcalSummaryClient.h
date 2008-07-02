#ifndef HcalSummaryClient_H
#define HcalSummaryClient_H

/*
 * \file HcalSummaryClient.h
 *
 * Code ported from DQM/EcalBarrelMonitorClient/interface/EBSummaryClient.h
 * $Date: 2008/07/01 20:22:45 $
 * $Revision: 1.9 $
 * \author Jeff Temple
 *
*/

#include <vector>
#include <string>
#include <fstream>
#include <sys/time.h>

#include "TROOT.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQM/HcalMonitorClient/interface/HcalBaseClient.h"

class MonitorElement;
class DQMStore;

class HcalSummaryClient : public HcalBaseClient {

 public:

  // Constructor
   
  HcalSummaryClient(const edm::ParameterSet& ps);
   
  // Destructor
  virtual ~HcalSummaryClient();
     
  // BeginJob
  void beginJob(DQMStore* dqmStore);
    
  // EndJob
  void endJob(void);
  
  // BeginRun
  void beginRun(void);
  
  // EndRun
  void endRun(void);
  
  // Setup
  void setup(void);
  
  // Cleanup
  void cleanup(void);
  

  // Analyze
  void analyze(void);
  float analyze_everything(std::string name, int type, float& subdet);

  float analyze_deadcell(std::string name, float& subdet); 
  float analyze_hotcell(std::string name, float& subdet);  
  float analyze_digi(std::string name, float& subdet);

  void incrementCounters(void);

  // HtmlOutput
  void htmlOutput(int& run, time_t& mytime, int& minlumi, int& maxlumi, std::string& htmlDir, std::string& htmlName);
  
  // Get Functions
  inline int getEvtPerJob() { return ievt_; }
  inline int getEvtPerRun() { return jevt_; }


 private:

  int ievt_;
  int jevt_;
  int lastupdate_;

  bool cloneME_;
  
  bool verbose_;
  bool debug_;
  
  std::string prefixME_;
  
  bool enableCleanup_;

  DQMStore* dqmStore_;

  MonitorElement* meGlobalSummary_;

  bool dataFormatClient_, digiClient_, recHitClient_, pedestalClient_;
  bool ledClient_, hotCellClient_, deadCellClient_, trigPrimClient_, caloTowerClient_;


  std::map<std::string, int> subdetCells_;
  int totalcells_; // stores total possible # of cells being checked

  // Individual status values for each task by subdetector
  float status_digi[4];
  float status_deadcell[4];
  float status_hotcell[4];

  float status_HB_;
  float status_HE_;
  float status_HO_;
  float status_HF_;
  float status_global_;

  double etaMin_, etaMax_, phiMin_, phiMax_;
  int phiBins_, etaBins_;
    
  ofstream htmlFile;

}; // end of class declaration

#endif
