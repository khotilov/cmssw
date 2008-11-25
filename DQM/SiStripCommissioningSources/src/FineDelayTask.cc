#include "DQM/SiStripCommissioningSources/interface/FineDelayTask.h"
#include "DataFormats/SiStripCommon/interface/SiStripHistoTitle.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <DQMServices/Core/interface/MonitorElement.h>

#define NBINS (500)
#define LOWBIN (-125)
#define HIGHBIN (125)


// -----------------------------------------------------------------------------
//
std::map<std::string, CommissioningTask::HistoSet> FineDelayTask::timingMap_;
MonitorElement* FineDelayTask::mode_ = NULL;

// -----------------------------------------------------------------------------
//
FineDelayTask::FineDelayTask( DQMStore* dqm,
			      const FedChannelConnection& conn ) :
  CommissioningTask( dqm, conn, "FineDelayTask" ), timing_(0)
{
  LogDebug("Commissioning") << "[FineDelayTask::FineDelayTask] Constructing object...";
  // compute the fiber length correction
  float length=conn.fiberLength();
  // convert cm to ns
  float c=30; //speed of light in cm/ns
  float refractionIndex = 1.4; // refraction index of the optical fibers
  fiberLengthCorrection_ =  length/c*refractionIndex;
}

// -----------------------------------------------------------------------------
//
FineDelayTask::~FineDelayTask() {
  LogDebug("Commissioning") << "[FineDelayTask::FineDelayTask] Destructing object...";
}

// -----------------------------------------------------------------------------
//
void FineDelayTask::book() {
  LogDebug("Commissioning") << "[FineDelayTask::book]";

  // construct the histo title
  // by setting the granularity to sistrip::TRACKER, the title will be identical for all detkeys.
  // therefore, only one histo will be booked/analyzed
  std::string title = SiStripHistoTitle( sistrip::EXPERT_HISTO, 
					 sistrip::FINE_DELAY, 
  					 sistrip::DET_KEY, 
					 0,
					 sistrip::TRACKER, 
					 0,
					 sistrip::extrainfo::clusterCharge_ ).title(); 
  // look if such an histogram is already booked
  if(timingMap_.find(title)!=timingMap_.end()) {
    // if already booked, use it
    LogDebug("Commissioning") << "[FineDelayTask::book] using existing histogram.";
  } else {
    // if not, book it
    int nBins = NBINS;
    LogDebug("Commissioning") << "[LatencyTask::book] booking a new histogram.";
    timingMap_[title].histo( dqm()->bookProfile( title, title,    // name and title
						 nBins, LOWBIN, HIGHBIN,   // binning + range
						 100, 0., -1., "s" ) );  // Y range : automatic

    timingMap_[title].vNumOfEntries_.resize(nBins,0);
    timingMap_[title].vSumOfContents_.resize(nBins,0);
    timingMap_[title].vSumOfSquares_.resize(nBins,0);
  }
  timing_ = &(timingMap_[title]);
  LogDebug("Commissioning") << "Binning is " << timing_->vNumOfEntries_.size();
  LogDebug("Commissioning") << "[FineDelayTask::book] done";
  if(!mode_) {
    std::string pwd = dqm()->pwd();
    std::string rootDir = pwd.substr(0,pwd.find(sistrip::root_ + "/")+sistrip::root_.size());
    dqm()->setCurrentFolder( rootDir );
    mode_ = dqm()->bookInt("latencyCode");
  }

}

// -----------------------------------------------------------------------------
//
void FineDelayTask::fill( const SiStripEventSummary& summary,
			  const edm::DetSet<SiStripRawDigi>& digis ) {
  LogDebug("Commissioning") << "[FineDelayTask::fill]";
  // retrieve the delay from the EventSummary
  float delay = const_cast<SiStripEventSummary&>(summary).ttcrx();
  uint32_t latencyCode = (const_cast<SiStripEventSummary&>(summary).layerScanned()>>24)&0xff;
  LogDebug("Commissioning") << "[FineDelayTask::fill]: layerScanned() is " << const_cast<SiStripEventSummary&>(summary).layerScanned();
  int latencyShift = latencyCode & 0x3f;             // number of bunch crossings between current value and start of scan... must be positive
  if((latencyCode>>6)==2) latencyShift -= 3;         // layer in deconv, rest in peak
  if((latencyCode>>6)==1) latencyShift += 3;         // layer in peak, rest in deconv
  float correctedDelay = delay - (latencyShift*25.); // shifts the delay so that 0 corresponds to the current settings.

  LogDebug("Commissioning") << "[FineDelayTask::fill]; the delay is " << delay;
  // loop on the strips to find the (maybe) non-zero digi
  for(unsigned int strip=0;strip<digis.data.size();strip++) {
    if(digis.data[strip].adc()!=0) {
      // apply the TOF correction
      float tof = (digis.data[strip].adc()>>8)/10.;
      correctedDelay = delay - (latencyShift*25.) - tof;
      if((digis.data[strip].adc()>>8)==255) continue; // skip hit if TOF is in overflow
      // compute the bin
      float nbins = NBINS;
      float lowbin = LOWBIN;
      float highbin = HIGHBIN;
      int bin = int((correctedDelay-lowbin)/((highbin-lowbin)/nbins));
      LogDebug("Commissioning") << "[FineDelayTask::fill]; using a hit with value " << ( digis.data[strip].adc()&0xff )
                                << " at corrected delay of " << correctedDelay
				<< " in bin " << bin << "  (tof is " << tof << "( since adc = " << digis.data[strip].adc() << "))";
      updateHistoSet( *timing_,bin,digis.data[strip].adc()&0xff);
      if(mode_) mode_->Fill(latencyCode);
    }
  }
}

// -----------------------------------------------------------------------------
//
void FineDelayTask::update() {
  LogDebug("Commissioning") << "[FineDelayTask::update]";
  updateHistoSet( *timing_ );
}

