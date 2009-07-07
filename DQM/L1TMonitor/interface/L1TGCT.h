// -*-C++-*-
#ifndef L1TGCT_H
#define L1TGCT_H

/*
 * \file L1TGCT.h
 *
 * $Date: 2008/09/21 14:33:12 $
 * $Revision: 1.15 $
 * \author J. Berryhill
 * $Id: L1TGCT.h,v 1.15 2008/09/21 14:33:12 jad Exp $
 * $Log: L1TGCT.h,v $
 * Revision 1.15  2008/09/21 14:33:12  jad
 * updated HF Sums & Counts and added individual Jet Candidates and differences
 *
 * Revision 1.14  2008/06/09 11:08:05  tapper
 * Removed electron sub-folders with histograms per eta and phi bin.
 *
 * Revision 1.13  2008/06/02 11:08:58  tapper
 * Added HF ring histograms....
 *
 * Revision 1.12  2008/04/28 09:23:07  tapper
 * Added 1D eta and phi histograms for electrons and jets as input to Q tests.
 *
 * Revision 1.11  2008/04/25 15:40:21  tapper
 * Added histograms to EventInfo//errorSummarySegments.
 *
 * Revision 1.10  2008/03/01 00:40:00  lat
 * DQM core migration.
 *
 * Revision 1.9  2008/02/20 19:24:24  tapper
 * Removed noisy include.
 *
 * Revision 1.8  2008/02/20 18:59:29  tapper
 * Ported GCTMonitor histograms into L1TGCT
 *
 * Revision 1.7  2007/09/04 02:54:21  wittich
 * - fix dupe ME in RCT
 * - put in rank>0 req in GCT
 * - various small other fixes
 *
 * Revision 1.6  2007/08/31 18:14:20  wittich
 * update GCT packages to reflect GctRawToDigi, and move to raw plots
 *
 * Revision 1.5  2007/08/31 11:02:55  wittich
 * cerr -> LogInfo
 *
 * Revision 1.4  2007/02/22 19:43:52  berryhil
 *
 *
 *
 * InputTag parameters added for all modules
 *
 * Revision 1.3  2007/02/19 22:49:53  wittich
 * - Add RCT monitor
 *
 * Revision 1.2  2007/02/19 21:11:23  wittich
 * - Updates for integrating GCT monitor.
 *   + Adapted right now only the L1E elements thereof.
 *   + added DataFormats/L1Trigger to build file.
 *
 *
*/

// system include files
#include <memory>
#include <unistd.h>


#include <iostream>
#include <fstream>
#include <vector>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

// DQM
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"





//
// class declaration
//

class L1TGCT : public edm::EDAnalyzer {

public:

// Constructor
  L1TGCT(const edm::ParameterSet& ps);

// Destructor
 virtual ~L1TGCT();

protected:
// Analyze
 void analyze(const edm::Event& e, const edm::EventSetup& c);

// BeginJob
 void beginJob(const edm::EventSetup& c);

// EndJob
void endJob(void);

private:
  // ----------member data ---------------------------
  DQMStore * dbe;

  // GCT stuff
  MonitorElement* l1GctAllJetsEtEtaPhi_; 
  MonitorElement* l1GctCenJetsEtEtaPhi_; 
  MonitorElement* l1GctForJetsEtEtaPhi_;
  MonitorElement* l1GctTauJetsEtEtaPhi_;
  MonitorElement* l1GctIsoEmRankEtaPhi_;
  MonitorElement* l1GctNonIsoEmRankEtaPhi_;

  MonitorElement* l1GctCenJetsOccEtaPhi_;
  MonitorElement* l1GctForJetsOccEtaPhi_;  
  MonitorElement* l1GctTauJetsOccEtaPhi_;  
  MonitorElement* l1GctIsoEmOccEtaPhi_;    
  MonitorElement* l1GctNonIsoEmOccEtaPhi_; 

  MonitorElement* l1GctCenJetsRank_;
  MonitorElement* l1GctForJetsRank_;
  MonitorElement* l1GctTauJetsRank_;
  MonitorElement* l1GctIsoEmRank_;
  MonitorElement* l1GctNonIsoEmRank_;

  MonitorElement* l1GctEtMiss_;
  MonitorElement* l1GctEtMissPhi_;
  MonitorElement* l1GctEtTotal_;
  MonitorElement* l1GctEtHad_;
  
  //HF Rings stuff
  MonitorElement* l1GctHFRing1PosEtaNegEta_;
  MonitorElement* l1GctHFRing2PosEtaNegEta_;
  MonitorElement* l1GctHFRing1TowerCountPosEtaNegEta_;
  MonitorElement* l1GctHFRing2TowerCountPosEtaNegEta_;
  MonitorElement* l1GctHFRing1TowerCountPosEta_;
  MonitorElement* l1GctHFRing1TowerCountNegEta_;
  MonitorElement* l1GctHFRing2TowerCountPosEta_;
  MonitorElement* l1GctHFRing2TowerCountNegEta_;
  MonitorElement* l1GctHFRing1ETSumPosEta_;
  MonitorElement* l1GctHFRing1ETSumNegEta_;
  MonitorElement* l1GctHFRing2ETSumPosEta_;
  MonitorElement* l1GctHFRing2ETSumNegEta_;
  MonitorElement* l1GctHFRingRatioPosEta_;
  MonitorElement* l1GctHFRingRatioNegEta_;

  // GCT electron stuff
  MonitorElement* l1GctIsoEmRankCand0_;
  MonitorElement* l1GctIsoEmRankCand1_;
  MonitorElement* l1GctIsoEmRankCand2_;
  MonitorElement* l1GctIsoEmRankCand3_;

  MonitorElement* l1GctNonIsoEmRankCand0_;
  MonitorElement* l1GctNonIsoEmRankCand1_;
  MonitorElement* l1GctNonIsoEmRankCand2_;
  MonitorElement* l1GctNonIsoEmRankCand3_;

  MonitorElement* l1GctIsoEmRankDiff01_;
  MonitorElement* l1GctIsoEmRankDiff12_;
  MonitorElement* l1GctIsoEmRankDiff23_;
  MonitorElement* l1GctNonIsoEmRankDiff01_;
  MonitorElement* l1GctNonIsoEmRankDiff12_;
  MonitorElement* l1GctNonIsoEmRankDiff23_;

  //GCT jet stuff
  MonitorElement* l1GctCenJetsRankCand0_;
  MonitorElement* l1GctCenJetsRankCand1_;
  MonitorElement* l1GctCenJetsRankCand2_;
  MonitorElement* l1GctCenJetsRankCand3_;
  MonitorElement* l1GctForJetsRankCand0_;
  MonitorElement* l1GctForJetsRankCand1_;
  MonitorElement* l1GctForJetsRankCand2_;
  MonitorElement* l1GctForJetsRankCand3_;
  MonitorElement* l1GctTauJetsRankCand0_;
  MonitorElement* l1GctTauJetsRankCand1_;
  MonitorElement* l1GctTauJetsRankCand2_;
  MonitorElement* l1GctTauJetsRankCand3_;

  MonitorElement* l1GctCenJetsRankDiff01_;
  MonitorElement* l1GctCenJetsRankDiff12_;
  MonitorElement* l1GctCenJetsRankDiff23_;
  MonitorElement* l1GctForJetsRankDiff01_;
  MonitorElement* l1GctForJetsRankDiff12_;
  MonitorElement* l1GctForJetsRankDiff23_;
  MonitorElement* l1GctTauJetsRankDiff01_;
  MonitorElement* l1GctTauJetsRankDiff12_;
  MonitorElement* l1GctTauJetsRankDiff23_;


  int nev_; // Number of events processed
  std::string outputFile_; //file name for ROOT ouput
  bool verbose_;
  bool monitorDaemon_;
  ofstream logFile_;

  edm::InputTag gctCenJetsSource_;
  edm::InputTag gctForJetsSource_;
  edm::InputTag gctTauJetsSource_;
  edm::InputTag gctEnergySumsSource_;
  edm::InputTag gctIsoEmSource_;
  edm::InputTag gctNonIsoEmSource_;

};

#endif
