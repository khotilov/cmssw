#ifndef ESIntegrityClient_H
#define ESIntegrityClient_H

#include <vector>
#include <string>

#include "TROOT.h"
#include "TProfile2D.h"
#include "TH1F.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQM/EcalPreshowerMonitorClient/interface/ESClient.h"

class MonitorElement;
class DQMStore;
class ESCondDBInterface;
class RunIOV;
class ESMonRunIOV;

class ESIntegrityClient : public ESClient {
  
  friend class ESSummaryClient;
  
 public:
  
  /// Constructor
  ESIntegrityClient(const edm::ParameterSet& ps);
  
  /// Destructor
  virtual ~ESIntegrityClient();
  
  /// Analyze
  void analyze(void);
  
  /// BeginJob
  // void beginJob(DQMStore* dqmStore);
  void beginJob(void);
  /// EndJob
  void endJob(void);
  
  /// BeginRun
  void beginRun(void);
  
  /// EndRun
  void endRun(void);
  
  /// Setup
  void setup(void);
  
  /// Cleanup
  void cleanup(void);
  
  /// SoftReset
  void softReset(bool flag);

  ///  writeDb
  void  writeDb(ESCondDBInterface* econn, RunIOV* runiov, ESMonRunIOV* moniov, int side);
  
  /// Get Functions
  inline int getEvtPerJob() { return ievt_; }
  inline int getEvtPerRun() { return jevt_; }
  
 private:
  
  int ievt_;
  int jevt_;
  int fed_[2][2][40][40];  
  int kchip_[2][2][40][40];  
  int fedStatus_[56];
  int fiberStatus_[56];
  int syncStatus_[56];
  int slinkCRCStatus_[56];

  bool cloneME_;
  bool verbose_;
  bool debug_;
  bool enableCleanup_;
  
  std::string prefixME_;

  edm::FileInPath lookup_;
  
  DQMStore* dqmStore_;

  MonitorElement* meFED_[2][2];
  MonitorElement* meKCHIP_[2][2];
  MonitorElement* meDIErrors_[2][2];

  TH1F *hFED_;  
  TH2F *hFiber_;
  TH2F *hKF1_;
  TH2F *hKF2_;
  TH1F *hKBC_;
  TH1F *hKEC_; 
  TH1F *hL1ADiff_;
  TH1F *hBXDiff_;
  TH1F *hOrbitNumberDiff_;
  TH1F *hSLinkCRCErr_; 

};

#endif
