// -*-c++-*-
// 
//
// $Id: HLTScalers.h,v 1.12 2008/09/03 13:59:05 wittich Exp $
// Class to collect HLT scaler information 
// for Trigger Cross Section Monitor
// [wittich 11/07] 

// $Log: HLTScalers.h,v $
// Revision 1.12  2008/09/03 13:59:05  wittich
// make HLT DQM path configurable via python parameter,
// which defaults to HLT/HLTScalers_EvF
//
// Revision 1.11  2008/09/03 02:13:47  wittich
// - bug fix in L1Scalers
// - configurable dqm directory in L1SCalers
// - other minor tweaks in HLTScalers
//
// Revision 1.10  2008/09/02 02:37:21  wittich
// - split L1 code from HLTScalers into L1Scalers
// - update cfi file accordingly
// - make sure to cd to correct directory before booking ME's
//
// Revision 1.9  2008/08/22 20:56:55  wittich
// - add client for HLT Scalers
// - Move rate calculation to HLTScalersClient and slim down the
//   filter-farm part of HLTScalers
//
// Revision 1.8  2008/08/15 15:40:57  wteo
// split hltScalers into smaller histos, calculate rates
//
// Revision 1.7  2008/08/01 14:37:33  bjbloom
// Added ability to specify which paths are cross-correlated
//
// Revision 1.6  2008/07/04 15:57:18  wittich
// - move histograms to HLT directory (was in L1T)
// - add counter for number of lumi sections
// - attempt to hlt label histo axes locally; disabled (it's illegible)
//
// Revision 1.5  2008/03/01 00:40:16  lat
// DQM core migration.
//
// Revision 1.4  2007/12/11 17:24:54  wittich
// - add extra monitoring histos (eg hlt exceptions and correlations)
//
// Revision 1.3  2007/12/04 20:24:32  wittich
// - make hlt histograms variable width depending on path
// - add strings for path names
// - add int for nprocessed
// - add L1 scaler locally derived on Kaori's suggestion
//   + updates to cfi file for this, need to include unpacking of GT
//
// Revision 1.2  2007/12/01 19:28:56  wittich
// - fix cfi file (debug -> verbose, HLT -> FU for TriggerResults  label)
// - handle multiple beginRun for same run (don't call reset on DQM )
// - remove PathTimerService from cfg file in test subdir
//
// Revision 1.1  2007/11/26 16:37:50  wittich
// Prototype HLT scaler information.
//

#ifndef HLTSCALERS_H
#define HLTSCALERS_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DQMServices/Core/interface/DQMStore.h"

#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

class HLTScalers: public edm::EDAnalyzer
{
public:
  /// Constructors
  HLTScalers(const edm::ParameterSet& ps);
  
  /// Destructor
  virtual ~HLTScalers() {};
  
  /// BeginJob
  void beginJob(void);

//   /// Endjob
//   void endJob(void);
  
  /// BeginRun
  void beginRun(const edm::Run& run, const edm::EventSetup& c);

  /// EndRun
  void endRun(const edm::Run& run, const edm::EventSetup& c);

  
//   /// Begin LumiBlock
//   void beginLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
//                             const edm::EventSetup& c) ;

  /// End LumiBlock
  /// DQM Client Diagnostic should be performed here
  void endLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
                          const edm::EventSetup& c);

  void analyze(const edm::Event& e, const edm::EventSetup& c) ;


private:
  DQMStore * dbe_;
  MonitorElement *scalers_;
  MonitorElement *scalersException_;
  MonitorElement *hltCorrelations_;
  MonitorElement *detailedScalers_;
  std::string folderName_; // dqm folder name
  MonitorElement *nProc_;
  MonitorElement *nLumiBlock_;
  std::vector<MonitorElement*> hltPathNames_;
  edm::InputTag trigResultsSource_;
  bool resetMe_, monitorDaemon_; 

  int nev_; // Number of events processed
  int nLumi_; // number of lumi blocks
  int currentRun_;

};

#endif // HLTSCALERS_H
