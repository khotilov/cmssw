#include "DQM/SiStripCommissioningClients/interface/ApvTimingHistograms.h"
#include "DQM/SiStripCommissioningSummary/interface/SummaryGenerator.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>
#include <sstream>
#include <iomanip>
 
using namespace std;
using namespace sistrip;

// -----------------------------------------------------------------------------
/** */
ApvTimingHistograms::ApvTimingHistograms( MonitorUserInterface* mui ) 
  : CommissioningHistograms( mui, sistrip::APV_TIMING ),
    factory_( new Factory )
{
  cout << endl // LogTrace(mlDqmClient_) 
       << "[ApvTimingHistograms::" << __func__ << "]"
       << " Constructing object...";
}

// -----------------------------------------------------------------------------
/** */
ApvTimingHistograms::~ApvTimingHistograms() {
  cout << endl // LogTrace(mlDqmClient_) 
       << "[ApvTimingHistograms::" << __func__ << "]"
       << " Destructing object...";
}

// -----------------------------------------------------------------------------	 
/** */	 
void ApvTimingHistograms::histoAnalysis( bool debug ) {
  
  std::map<uint32_t,ApvTimingAnalysis*>::iterator ianal;

  // Clear std::map holding analysis objects
  for ( ianal = data_.begin(); ianal != data_.end(); ianal++ ) { 
    if ( ianal->second ) { delete ianal->second; }
  } 
  
  // Reset minimum / maximum delays
  float time_min =  1. * sistrip::invalid_;
  float time_max = -1. * sistrip::invalid_;
  uint32_t device_min = sistrip::invalid_;
  uint32_t device_max = sistrip::invalid_;
  
  // Iterate through std::map containing std::vectors of profile histograms
  CollationsMap::const_iterator iter = collations().begin();
  for ( ; iter != collations().end(); iter++ ) {
    
    // Check std::vector of histos is not empty (should be 1 histo)
    if ( iter->second.empty() ) {
      cerr << endl // edm::LogWarning(mlDqmClient_) 
	   << "[ApvTimingHistograms::" << __func__ << "]"
	   << " Zero collation histograms found!";
      continue;
    }
    
    // Retrieve pointers to profile histos for this FED channel 
    std::vector<TH1*> profs;
    Collations::const_iterator ihis = iter->second.begin(); 
    for ( ; ihis != iter->second.end(); ihis++ ) {
      TProfile* prof = ExtractTObject<TProfile>().extract( ihis->second->getMonitorElement() );
      if ( prof ) { profs.push_back(prof); }
    } 
    
    // Perform histo analysis
    ApvTimingAnalysis* anal = new ApvTimingAnalysis( iter->first );
    anal->analysis( profs );
    data_[iter->first] = anal; 
    
    // Check tick height is valid
    if ( anal->height() < 100. ) { 
      cerr << endl // edm::LogWarning(mlDqmClient_) 
	   << "[ApvTimingHistograms::" << __func__ << "]"
	   << " Tick mark height too small: " << anal->height();
      continue; 
    }

    // Check time of rising edge
    if ( anal->time() > sistrip::maximum_ ) { continue; }
    
    // Find maximum time
    if ( anal->time() > time_max ) { 
      time_max = anal->time(); 
      device_max = iter->first;
    }
    
    // Find minimum time
    if ( anal->time() < time_min ) { 
      time_min = anal->time(); 
      device_min = iter->first;
    }
    
  }
  
  cout << endl // LogTrace(mlDqmClient_) 
       << "[ApvTimingHistograms::" << __func__ << "]"
       << " Analyzed histograms for " 
       << collations().size() 
       << " FED channels";

  // Check max time
  if ( time_max > sistrip::maximum_ ||
       time_max < -1.*sistrip::maximum_ ) { 
    cerr << endl // edm::LogWarning(mlDqmClient_) 
	 << "[ApvTimingHistograms::" << __func__ << "]"
	 << " Unable to set maximum time! Found unexpected value: "
	 << time_max;
    return; 
  }
  
  SiStripFecKey max( device_max );
  cout << endl // LogTrace(mlDqmClient_) 
       << "[ApvTimingHistograms::" << __func__ << "]"
       << " Device (FEC/slot/ring/CCU/module/channel) " 
       << max.fecCrate_ << "/" 
       << max.fecSlot_ << "/" 
       << max.fecRing_ << "/" 
       << max.ccuAddr_ << "/"
       << max.ccuChan_ << "/"
       << max.channel() 
       << " has maximum delay (rising edge) [ns]:" << time_max;
  
  SiStripFecKey min( device_min );
  cout << endl // LogTrace(mlDqmClient_) 
       << "[ApvTimingHistograms::" << __func__ << "]"
       << " Device (FEC/slot/ring/CCU/module/channel) " 
       << min.fecCrate_ << "/" 
       << min.fecSlot_ << "/" 
       << min.fecRing_ << "/" 
       << min.ccuAddr_ << "/"
       << min.ccuChan_ << "/"
       << max.channel() 
       << " has minimum delay (rising edge) [ns]:" << time_min;
  
  // Set maximum time for all analysis objects
  for ( ianal = data_.begin(); ianal != data_.end(); ianal++ ) { 
    ianal->second->maxTime( time_max ); 
    static uint16_t cntr = 0;
    if ( debug ) {
      std::stringstream ss;
      ianal->second->print( ss ); 
      cout << endl // LogTrace(mlDqmClient_) 
	   << ss.str();
      cntr++;
    }
  }
  
}

// -----------------------------------------------------------------------------
/** */
void ApvTimingHistograms::createSummaryHisto( const sistrip::Monitorable& mon, 
					      const sistrip::Presentation& pres, 
					      const std::string& dir,
					      const sistrip::Granularity& gran ) {
  cout << endl // LogTrace(mlDqmClient_)
       << "[ApvTimingHistograms::" << __func__ << "]";
  
  // Check view 
  sistrip::View view = SiStripEnumsAndStrings::view(dir);
  if ( view == sistrip::UNKNOWN_VIEW ) { return; }
  
  // Analyze histograms
  histoAnalysis( false );
  
  // Extract data to be histogrammed
  uint32_t bins = factory_->init( mon, pres, view, dir, gran, data_ );
  
  // Create summary histogram (if it doesn't already exist)
  TH1* summary = histogram( mon, pres, view, dir, bins );
  
  // Fill histogram with data
  factory_->fill( *summary );
  
}
