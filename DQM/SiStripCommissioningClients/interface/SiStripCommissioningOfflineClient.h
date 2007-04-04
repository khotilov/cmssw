// Last commit: $Id: $

#ifndef DQM_SiStripCommissioningClients_SiStripCommissioningOfflineClient_H
#define DQM_SiStripCommissioningClients_SiStripCommissioningOfflineClient_H

#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include "DQM/SiStripCommissioningClients/interface/SiStripTFile.h"
#include "DQM/SiStripCommissioningClients/interface/ConfigParser.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <string>
#include <vector>
#include <map>

class CommissioningHistograms;
class DaqMonitorBEInterface;
class TH1;

/**
   @class SiStripCommissioningOfflineClient 
   @author M.Wingham, R.Bainbridge
   
   @brief Class which reads a root file containing "commissioning
   histograms", analyzes the histograms to extract "monitorables", and
   creates summary histograms.
*/
class SiStripCommissioningOfflineClient : public edm::EDAnalyzer {

 public:
  
  SiStripCommissioningOfflineClient( const edm::ParameterSet& );
  virtual ~SiStripCommissioningOfflineClient();
  
  virtual void beginJob( edm::EventSetup const& );
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  virtual void endJob() {;}
  
 protected:
  
  /** */
  virtual void processHistos( DaqMonitorBEInterface* const,
			      const std::vector<std::string>& );
  
  /** Action "executor" */
  CommissioningHistograms* histos_;
  
  /** Input .root file. */
  std::string rootFile_;
  
  /** Input .xml file. */
  std::string xmlFile_;

  /** Flag. */
  bool createSummaryPlots_;

  /** Flag. */
  bool saveSummaryPlots_;
  
  bool uploadToDb_;

  /** Commissioning runType. */
  sistrip::RunType runType_;
  
  /** Logical view. */
  sistrip::View view_;
  
  /** Run number. */
  uint16_t run_;

  /** */
  typedef std::vector<TH1*> Histos;

  /** */
  typedef std::map<uint32_t,Histos> HistosMap;

  /** Map containing commissioning histograms. */
  HistosMap map_;
  
  /** SummaryPlot objects. */
  std::vector<ConfigParser::SummaryPlot> plots_;
  
};

#endif // DQM_SiStripCommissioningClients_SiStripCommissioningOfflineClient_H

