#ifndef DQMOFFLINE_TRIGGER_EGHLTOFFLINECLIENT
#define DQMOFFLINE_TRIGGER_EGHLTOFFLINECLIENT

// -*- C++ -*-
//
// Package:    EgammaHLTOfflineClient
// Class:      EgammaHLTOffline
// 
/*
 Description: This is a DQM client meant to plot high-level HLT trigger 
 quantities as stored in the HLT results object TriggerResults for the Egamma triggers

 Notes:
  Currently I would like to plot simple histograms of three seperate types of variables
  1) global event quantities: eg nr of electrons
  2) di-object quanities: transverse mass, di-electron mass
  3) single object kinematic and id variables: eg et,eta,isolation

*/
//
// Original Author:  Sam Harper
//         Created:  June 2008
// 
//
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <vector>
#include <string>

class DQMStore;
class MonitorElement;


class EgHLTOfflineClient : public edm::EDAnalyzer {

 private:
  DQMStore* dbe_; //dbe seems to be the standard name for this, I dont know why. We of course dont own it
  std::string dirName_;
  
  std::vector<std::string> eleHLTFilterNames_;//names of the filters monitored using electrons to make plots for
  std::vector<std::string> eleTightLooseTrigNames_;
  std::vector<std::string> phoHLTFilterNames_;//names of the filters monitored using photons to make plots for
  std::vector<std::string> phoTightLooseTrigNames_;
  

  std::vector<std::string> eleEffVars_;
  std::vector<std::string> phoEffVars_;
  std::vector<std::string> eleTrigTPEffVsVars_;
  std::vector<std::string> phoTrigTPEffVsVars_;
  std::vector<std::string> eleLooseTightTrigEffVsVars_;
  std::vector<std::string> phoLooseTightTrigEffVsVars_;

  

  //disabling copying/assignment (in theory this is copyable but lets not just in case)
  EgHLTOfflineClient(const EgHLTOfflineClient& rhs){}
  EgHLTOfflineClient& operator=(const EgHLTOfflineClient& rhs){return *this;}

 public:
  explicit EgHLTOfflineClient(const edm::ParameterSet& );
  virtual ~EgHLTOfflineClient();
  
  
  virtual void beginJob(const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&); //dummy
  virtual void endJob();
  virtual void beginRun(const edm::Run& run, const edm::EventSetup& c);
  virtual void endRun(const edm::Run& run, const edm::EventSetup& c);
  
  
  virtual void beginLuminosityBlock(const edm::LuminosityBlock& lumiSeg,const edm::EventSetup& context){}
  // DQM Client Diagnostic
  virtual void endLuminosityBlock(const edm::LuminosityBlock& lumiSeg,const edm::EventSetup& c);

  //at somepoint these all may migrate to a helper class
  void createN1EffHists(const std::string& baseName,const std::string& region,const std::vector<std::string>& varNames);
  void createLooseTightTrigEff(const std::vector<std::string>&  tightLooseTrigNames,const std::string& region,const std::vector<std::string>& vsVarNames,const std::string& objName);
  void createTrigTagProbeEffHists(const std::string& filterName,const std::string& region,const std::vector<std::string>& vsVarNames,const std::string& objName);
  
  MonitorElement* makeEffMonElemFromPassAndAll(const std::string& name,const MonitorElement* pass,const MonitorElement* all);
  MonitorElement* makeEffMonElemFromPassAndFail(const std::string& name,const MonitorElement* pass,const MonitorElement* fail);

private:
  void runClient_(); //master function which runs the client
  
};
 


#endif
