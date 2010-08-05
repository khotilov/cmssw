#ifndef RPCDqmClient_H
#define RPCDqmClient_H

#include "DQM/RPCMonitorClient/interface/RPCClient.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include "DQMServices/Core/interface/MonitorElement.h"

#include <memory>
#include <string>
#include <vector>

class RPCDqmClient:public edm::EDAnalyzer{

public:

  /// Constructor
 RPCDqmClient(const edm::ParameterSet& ps);
  
  /// Destructor
  virtual ~ RPCDqmClient();

  /// BeginJob
  void beginJob( );

  //Begin Run
   void beginRun(const edm::Run& , const edm::EventSetup&);
  
  
  /// Begin Lumi block 
  void beginLuminosityBlock(edm::LuminosityBlock const& , edm::EventSetup const& ) ;

  /// Analyze  
  void analyze(const edm::Event& , const edm::EventSetup& );

  /// End Lumi Block
  void endLuminosityBlock(edm::LuminosityBlock const& , edm::EventSetup const& );
 
  //End Run
  void endRun(const edm::Run& , const edm::EventSetup& ); 		
  
  /// Endjob
  void endJob();

 protected:
  void makeClientMap(void);

 private:
  int prescaleGlobalFactor_, minimumEvents_, numLumBlock_;
 
  bool enableDQMClients_ , init_; 
  std::string  prefixDir_;
  std::string  globalFolder_;
  std::vector<std::string>  clientList_;

  std::string  recHitType_;
  std::vector<std::string> clientNames_,clientHisto_; 
  std::vector<RPCClient*> clientModules_;

  std::vector<int> clientTag_;
  //std::map<RPCClient *, std::string> clientMap_;
  edm::ParameterSet parameters_;

  DQMStore* dbe_;
 
  
};
#endif
