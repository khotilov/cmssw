#ifndef DQM_SiStripCommissioningSources_SiStripCommissioningSource_H
#define DQM_SiStripCommissioningSources_SiStripCommissioningSource_H

#include "CalibFormats/SiStripObjects/interface/SiStripFecCabling.h"
#include "CondFormats/SiStripObjects/interface/SiStripFedCabling.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include "DataFormats/SiStripCommon/interface/SiStripHistoTitle.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include <boost/cstdint.hpp>
#include <string>
#include <vector>
#include <map>

class DaqMonitorBEInterface;
class CommissioningTask;
class FedChannelConnection;
class SiStripEventSummary;

/**
   @class SiStripCommissioningSource
*/
class SiStripCommissioningSource : public edm::EDAnalyzer {
  
 public: // ---------- Public interface ----------
  
  /** Map of task objects, identified through FedChanelId */
  typedef std::map<unsigned int, CommissioningTask*> TaskMap;
  typedef std::vector<CommissioningTask*> VecOfTasks;
  typedef std::vector<VecOfTasks> VecOfVecOfTasks;
  
  SiStripCommissioningSource( const edm::ParameterSet& );
  ~SiStripCommissioningSource();
  
  void beginJob( edm::EventSetup const& );
  void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob();
  
 private: // ---------- Private methods ----------

  /** Private default constructor. */
  SiStripCommissioningSource();
  
  /** */
  DaqMonitorBEInterface* const dqm( std::string method = "" ) const;
  
  /** */
  void createRunNumber();

  /** */
  void createTask( const SiStripEventSummary* const, const edm::EventSetup& );
  
  /** */
  void createCablingTasks();

  /** */
  void createTasks( sistrip::RunType, const edm::EventSetup& );
  
  /** */
  void clearCablingTasks();

  /** */
  void clearTasks();
  
  /** */
  void fillCablingHistos( const SiStripEventSummary* const,
			  const edm::DetSetVector<SiStripRawDigi>& );

  /** */
  void fillHistos( const SiStripEventSummary* const,
		   const edm::DetSetVector<SiStripRawDigi>& );
  
  /** */
  void remove();
  
  void directory( std::stringstream&, 
		  uint32_t run_number = 0 );
  
  // ---------- DQM fwk and cabling ----------

  /** Interface to Data Quality Monitoring framework. */
  DaqMonitorBEInterface* dqm_;

  /** */
  SiStripFedCabling* fedCabling_;

  /** */
  SiStripFecCabling* fecCabling_;
  
  // ---------- Input / output ----------

  /** Name of digi input module. */
  std::string inputModuleLabel_;
  std::string inputModuleLabelSummary_;

  /** Filename of output root file containing source histos. */
  std::string filename_;

  /** Run number used for naming of root file. */
  uint32_t run_;

  /** Record of time used to calculate event rate. */
  int32_t time_;

  // ---------- Histogram-related ----------

  /** Identifies commissioning task read from cfg file. */
  std::string taskConfigurable_; 

  /** Identifies commissioning task. */
  sistrip::RunType task_; 

  /** Vector of vector of task objects (indexed using FED id.ch. */
  VecOfVecOfTasks tasks_;

  /** Map of cabling task objects (indexed using FEC key). */
  TaskMap cablingTasks_;

  /** Flag to indicate whether histo objects exist or not. */
  bool tasksExist_;

  /** Flag to indicate whether task is FED cabling or not. */
  bool cablingTask_;
  
  /** Update frequency for histograms (ignored for cabling). */
  int updateFreq_;

  /** */
  std::string base_;

};

#endif // DQM_SiStripCommissioningSources_SiStripCommissioningSource_H

