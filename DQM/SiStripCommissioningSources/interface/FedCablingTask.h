#ifndef DQM_SiStripCommissioningSources_FedCablingTask_h
#define DQM_SiStripCommissioningSources_FedCablingTask_h

#include "DQM/SiStripCommissioningSources/interface/CommissioningTask.h"
#include <vector>

/**
   @class FedCablingTask

   This object is stored in the TaskMap using FecKey as the key,
   rather than FedKey as for the other commissioning tasks.
*/
class FedCablingTask : public CommissioningTask {

 public:
  
  FedCablingTask( DaqMonitorBEInterface*, const FedChannelConnection& );
  virtual ~FedCablingTask();
  
 private:
  
  virtual void book();
  virtual void fill( const SiStripEventSummary&, 
		     const uint16_t& fed_id,
		     const std::map<uint16_t,float>& fed_ch );
  virtual void update();
  
  /** HistoSet for FED cabling. First element contains histo info for
      FED id, second element contains histo info for FED channel. */
  std::vector<HistoSet> cabling_;
  
};

#endif // DQM_SiStripCommissioningSources_FedCablingTask_h

