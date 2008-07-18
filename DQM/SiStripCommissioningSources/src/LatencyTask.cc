#include "DQM/SiStripCommissioningSources/interface/LatencyTask.h"
#include "DataFormats/SiStripCommon/interface/SiStripHistoTitle.h"
#include <DataFormats/SiStripDetId/interface/SiStripDetId.h>
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#define NBINS (192)
#define LOWBIN (-4800)
#define HIGHBIN (0)

// -----------------------------------------------------------------------------
//
std::map<std::string, CommissioningTask::HistoSet> LatencyTask::timingMap_;
std::map<std::string, CommissioningTask::HistoSet> LatencyTask::clusterMap_;

// -----------------------------------------------------------------------------
//
LatencyTask::LatencyTask( DQMStore* dqm,
			      const FedChannelConnection& conn ) :
  CommissioningTask( dqm, conn, "LatencyTask" ),timing_(0),cluster_(0),firstReading_(-1)
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
					 0,
					 sistrip::TRACKER, 
					 0,
					 sistrip::extrainfo::clusterCharge_).title(); 
  // look if such an histogram is already booked
  if(timingMap_.find(title)!=timingMap_.end()) {
    // if already booked, use it
    LogDebug("Commissioning") << "[LatencyTask::book] using existing histogram.";
  } else {
    // if not, book it
    timingMap_[title] = HistoSet();
    int nBins = NBINS;
    std::string pwd = dqm()->pwd();
    std::string rootDir = pwd.substr(0,pwd.find(sistrip::root_ + "/")+sistrip::root_.size());
    rootDir += "/"; rootDir += sistrip::controlView_;
    dqm()->setCurrentFolder( rootDir );
    LogDebug("Commissioning") << "[LatencyTask::book] booking a new histogram in " << dqm()->pwd();
    timingMap_[title].histo( dqm()->bookProfile( title, title,    // name and title
						 nBins, LOWBIN, HIGHBIN,   // binning + range
						 100, 0., -1.,"" ) );  // Y range : automatic
    dqm()->setCurrentFolder( pwd );
  
    timingMap_[title].vNumOfEntries_.resize(nBins,0);
    timingMap_[title].vSumOfContents_.resize(nBins,0);
    timingMap_[title].vSumOfSquares_.resize(nBins,0);
  }
  timing_ = &(timingMap_[title]);
  // same histo at the partition level
  title = SiStripHistoTitle( sistrip::EXPERT_HISTO, 
					 sistrip::APV_LATENCY, 
  					 sistrip::DET_KEY, 
					 int(SiStripDetId(connection().detId()).subDetector()),
					 sistrip::PARTITION, 
					 0,
					 sistrip::extrainfo::clusterCharge_).title(); 
  // look if such an histogram is already booked
  if(timingMap_.find(title)!=timingMap_.end()) {
    // if already booked, use it
    LogDebug("Commissioning") << "[LatencyTask::book] using existing histogram.";
  } else {
    // if not, book it
    timingMap_[title] = HistoSet();
    int nBins = NBINS;
    LogDebug("Commissioning") << "[LatencyTask::book] booking a new histogram in " << dqm()->pwd();
    timingMap_[title].histo( dqm()->bookProfile( title, title,    // name and title
						 nBins, LOWBIN, HIGHBIN,   // binning + range
						 100, 0., -1.,"" ) );  // Y range : automatic
  
    timingMap_[title].vNumOfEntries_.resize(nBins,0);
    timingMap_[title].vSumOfContents_.resize(nBins,0);
    timingMap_[title].vSumOfSquares_.resize(nBins,0);
  }
  timingPartition_ = &(timingMap_[title]);
  // construct the histo title
  // by setting the granularity to sistrip::TRACKER, the title will be identical for all detkeys.
  // therefore, only one histo will be booked/analyzed
  title = SiStripHistoTitle( sistrip::EXPERT_HISTO, 
			     sistrip::APV_LATENCY, 
                             sistrip::DET_KEY, 
                             0,
                             sistrip::TRACKER, 
                             0,
                             sistrip::extrainfo::occupancy_).title(); 
  // look if such an histogram is already booked
  if(clusterMap_.find(title)!=clusterMap_.end()) {
    // if already booked, use it
    LogDebug("Commissioning") << "[LatencyTask::book] using existing histogram.";
  } else {
    // if not, book it
    clusterMap_[title] = HistoSet();
    int nBins = NBINS;
    std::string pwd = dqm()->pwd();
    std::string rootDir = pwd.substr(0,pwd.find(sistrip::root_ + "/")+sistrip::root_.size());
    rootDir += "/"; rootDir += sistrip::controlView_;
    dqm()->setCurrentFolder( rootDir );
    LogDebug("Commissioning") << "[LatencyTask::book] booking a new histogram in " << dqm()->pwd();
    clusterMap_[title].histo( dqm()->book1D( title, title,    // name and title
                                             nBins, LOWBIN, HIGHBIN ));  // binning + range
    dqm()->setCurrentFolder( pwd );
  
    clusterMap_[title].isProfile_=false;
    clusterMap_[title].vNumOfEntries_.resize(nBins,0);
    clusterMap_[title].vSumOfContents_.resize(nBins,0);
    clusterMap_[title].vSumOfSquares_.resize(nBins,0);
  }
  cluster_ = &(clusterMap_[title]);
  // same histo at the partition level
  title = SiStripHistoTitle( sistrip::EXPERT_HISTO, 
			     sistrip::APV_LATENCY, 
                             sistrip::DET_KEY, 
                             int(SiStripDetId(connection().detId()).subDetector()),
                             sistrip::PARTITION,
                             0,
                             sistrip::extrainfo::occupancy_).title(); 
  // look if such an histogram is already booked
  if(clusterMap_.find(title)!=clusterMap_.end()) {
    // if already booked, use it
    LogDebug("Commissioning") << "[LatencyTask::book] using existing histogram.";
  } else {
    // if not, book it
    clusterMap_[title] = HistoSet();
    int nBins = NBINS;
    LogDebug("Commissioning") << "[LatencyTask::book] booking a new histogram in " << dqm()->pwd();
    clusterMap_[title].histo( dqm()->book1D( title, title,    // name and title
                                             nBins, LOWBIN, HIGHBIN ));  // binning + range
  
    clusterMap_[title].isProfile_=false;
    clusterMap_[title].vNumOfEntries_.resize(nBins,0);
    clusterMap_[title].vSumOfContents_.resize(nBins,0);
    clusterMap_[title].vSumOfSquares_.resize(nBins,0);
  }
  clusterPartition_ = &(clusterMap_[title]);

  LogDebug("Commissioning") << "[LatencyTask::book] done";
}

// -----------------------------------------------------------------------------
//
void LatencyTask::fill( const SiStripEventSummary& summary,
			const edm::DetSet<SiStripRawDigi>& digis ) {
  LogDebug("Commissioning") << "[LatencyTask::fill]";
  // retrieve the delay from the EventSummary
  int32_t delay = static_cast<int32_t>( const_cast<SiStripEventSummary&>(summary).latency() );
  if(firstReading_==-1) firstReading_ = delay;
  float correctedDelay = 0.;
  LogDebug("Commissioning") << "[LatencyTask::fill]; the delay is " << delay;
  // loop on the strips to find the (maybe) non-zero digi
  unsigned int nclusters = 0;
  for(unsigned int strip=0;strip<digis.data.size();strip++) {
    if(digis.data[strip].adc()!=0) {
      // count the "cluster"
      ++nclusters;
      // no TOF correction is applied.
      // 2 reasons: the effect is a priori to thin to be seen with 25ns steps
      // and it biases the result by one clock due to the 25bins in the HistoSet
      correctedDelay = delay*(-25.); // no TOF correction is applied. 
      // compute the bin
      int bin = int((correctedDelay-LOWBIN)/((HIGHBIN-LOWBIN)/NBINS));
      LogDebug("Commissioning") << "[LatencyTask::fill]; using a hit with value " << ( digis.data[strip].adc()&0xff )
                                << " at corrected delay of " << correctedDelay
				<< " in bin " << bin ;
      updateHistoSet( *timing_,bin,digis.data[strip].adc()&0xff);
      updateHistoSet( *timingPartition_,bin,digis.data[strip].adc()&0xff);
    }
  }
  // set the occupancy
  int bin = int((delay*(-25.)-LOWBIN)/((HIGHBIN-LOWBIN)/NBINS));
  LogDebug("Commissioning") << "[LatencyTask::fill]; occupancy is " << nclusters;
  updateHistoSet( *cluster_,bin,nclusters );
  updateHistoSet( *clusterPartition_,bin,nclusters );
}

// -----------------------------------------------------------------------------
//
void LatencyTask::update() {
  LogDebug("Commissioning") << "[LatencyTask::update]";
  updateHistoSet( *timing_ );
  updateHistoSet( *timingPartition_ );
  updateHistoSet( *cluster_ );
  updateHistoSet( *clusterPartition_ );
}

