#include "DQM/SiStripCommissioningSources/interface/LatencyTask.h"
#include "DataFormats/SiStripCommon/interface/SiStripHistoTitle.h"
#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#define NBINS (100)
#define LOWBIN (-2500)
#define HIGHBIN (0)

// -----------------------------------------------------------------------------
//
std::map<std::string, CommissioningTask::HistoSet> LatencyTask::timingMap_;

// -----------------------------------------------------------------------------
//
LatencyTask::LatencyTask( DaqMonitorBEInterface* dqm,
			      const FedChannelConnection& conn ) :
  CommissioningTask( dqm, conn, "LatencyTask" ),timing_(0)
{
  LogDebug("Commissioning") << "[LatencyTask::LatencyTask] Constructing object...";
}

// -----------------------------------------------------------------------------
//
LatencyTask::~LatencyTask() {
  LogDebug("Commissioning") << "[LatencyTask::LatencyTask] Destructing object...";
}

// -----------------------------------------------------------------------------
//
void LatencyTask::book() {
  LogDebug("Commissioning") << "[LatencyTask::book]";

  // construct the histo title
  // by setting the granularity to sistrip::TRACKER, the title will be identical for all detkeys.
  // therefore, only one histo will be booked/analyzed
  std::string title = SiStripHistoTitle( sistrip::EXPERT_HISTO, 
					 sistrip::APV_LATENCY, 
  					 sistrip::DET_KEY, 
					 0,//connection().detId(),
					 sistrip::TRACKER, 
					 0 ).title(); 
  // look if such an histogram is already booked
  if(timingMap_.find(title)!=timingMap_.end()) {
    // if already booked, use it
    LogDebug("Commissioning") << "[LatencyTask::book] using existing histogram.";
  } else {
    // if not, book it
    timingMap_[title] = HistoSet();
    int nBins = NBINS;
    LogDebug("Commissioning") << "[LatencyTask::book] booking a new histogram.";
    timingMap_[title].histo_ = dqm()->bookProfile( title, title,    // name and title
  				         nBins, LOWBIN, HIGHBIN,   // binning + range
				         100, 0., -1. );  // Y range : automatic
  
    timingMap_[title].vNumOfEntries_.resize(nBins,0);
    timingMap_[title].vSumOfContents_.resize(nBins,0);
    timingMap_[title].vSumOfSquares_.resize(nBins,0);
  }
  timing_ = &(timingMap_[title]);
  LogDebug("Commissioning") << "Binning is " << timing_->vNumOfEntries_.size();
  LogDebug("Commissioning") << "[LatencyTask::book] done";
}

// -----------------------------------------------------------------------------
//
void LatencyTask::fill( const SiStripEventSummary& summary,
			  const edm::DetSet<SiStripRawDigi>& digis ) {
  LogDebug("Commissioning") << "[LatencyTask::fill]";
  // retrieve the delay from the EventSummary
  uint32_t delay = const_cast<SiStripEventSummary&>(summary).latency();
  float correctedDelay = 0.;
  LogDebug("Commissioning") << "[LatencyTask::fill]; the delay is " << delay;
  // loop on the strips to find the (maybe) non-zero digi
  for(unsigned int strip=0;strip<digis.data.size();strip++) {
    if(digis.data[strip].adc()!=0) {
      // apply the TOF correction
      float tof = (digis.data[strip].adc()>>8)/10.;
      correctedDelay = delay*(-25.) - tof;
      if((digis.data[strip].adc()>>8)==255) continue; // skip hit if TOF is in overflow
      // compute the bin
      int bin = int((correctedDelay-LOWBIN)/((HIGHBIN-LOWBIN)/NBINS));
      LogDebug("Commissioning") << "[LatencyTask::fill]; using a hit with value " << ( digis.data[strip].adc()&0xff )
                                << " at corrected delay of " << correctedDelay
				<< " in bin " << bin << "  (tof is " << tof << "( since adc = " << digis.data[strip].adc() << "))";
      updateHistoSet( *timing_,bin,digis.data[strip].adc()&0xff);
      //break;
    }
  }
}

// -----------------------------------------------------------------------------
//
void LatencyTask::update() {
  LogDebug("Commissioning") << "[LatencyTask::update]";
  updateHistoSet( *timing_ );
}

