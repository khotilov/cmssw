#ifndef DQM_SiStripCommissioningClients_OptoScanHistograms_H
#define DQM_SiStripCommissioningClients_OptoScanHistograms_H

#include "DQM/SiStripCommissioningClients/interface/CommissioningHistograms.h"
#include "DQM/SiStripCommissioningSummary/interface/OptoScanSummaryFactory.h"
#include "DQM/SiStripCommissioningAnalysis/interface/OptoScanAnalysis.h"

class MonitorUserInterface;
class DaqMonitorBEInterface;

class OptoScanHistograms : public CommissioningHistograms {

 public:

  OptoScanHistograms( MonitorUserInterface* );
  OptoScanHistograms( DaqMonitorBEInterface* );
  virtual ~OptoScanHistograms();

  typedef SummaryPlotFactory<OptoScanAnalysis*> Factory;
  typedef std::map<uint32_t,OptoScanAnalysis*> Analyses;
  
  /** */
  void histoAnalysis( bool debug );

  /** */
  void printAnalyses();

  /** */
  void createSummaryHisto( const sistrip::Monitorable&,
			   const sistrip::Presentation&,
			   const std::string& top_level_dir,
			   const sistrip::Granularity& );

 protected:

  Analyses data_;

  std::auto_ptr<Factory> factory_;
  
};

#endif // DQM_SiStripCommissioningClients_OptoScanHistograms_H


