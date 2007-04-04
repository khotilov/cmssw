#ifndef DQM_SiStripCommissioningClients_VpspScanHistograms_H
#define DQM_SiStripCommissioningClients_VpspScanHistograms_H

#include "DQM/SiStripCommissioningClients/interface/CommissioningHistograms.h"
#include "DQM/SiStripCommissioningSummary/interface/VpspScanSummaryFactory.h"
#include "DQM/SiStripCommissioningAnalysis/interface/VpspScanAnalysis.h"

class MonitorUserInterface;
class DaqMonitorBEInterface;

class VpspScanHistograms : public CommissioningHistograms {

 public:
  
  VpspScanHistograms( MonitorUserInterface* );
  VpspScanHistograms( DaqMonitorBEInterface* );
  virtual ~VpspScanHistograms();
  
  typedef SummaryHistogramFactory<VpspScanAnalysis> Factory;
  
  /** */
  void histoAnalysis( bool debug );

  /** */
  void createSummaryHisto( const sistrip::Monitorable&,
			   const sistrip::Presentation&,
			   const std::string& top_level_dir,
			   const sistrip::Granularity& );

 protected:

  std::map<uint32_t,VpspScanAnalysis> data_;

  std::auto_ptr<Factory> factory_;

};

#endif // DQM_SiStripCommissioningClients_VpspScanHistograms_H


