#ifndef DQM_SiStripCommissioningAnalysis_OptoScanAlgorithm_H
#define DQM_SiStripCommissioningAnalysis_OptoScanAlgorithm_H

#include "DQM/SiStripCommissioningAnalysis/interface/CommissioningAlgorithm.h"
#include <boost/cstdint.hpp>
#include <vector>

class OptoScanAnalysis;
class TProfile;
class TH1;

/** 
   @class OptoScanAnalysis
   @author M. Wingham, R.Bainbridge
   @brief Histogram-based analysis for opto bias/gain scan.
*/
class OptoScanAlgorithm : public CommissioningAlgorithm {
  
 public:
  
  OptoScanAlgorithm( OptoScanAnalysis* const );
  
  virtual ~OptoScanAlgorithm() {;}

  /** Histogram pointer and title. */
  Histo histo( const uint16_t& gain, 
	       const uint16_t& digital_level ) const;
  
 private:

  OptoScanAlgorithm() {;}
  
  /** Extracts and organises histograms. */
  void extract( const std::vector<TH1*>& );

  /** Performs histogram anaysis. */
  void analyse();

 private:
  
  /** Pointers and titles for histograms. */
  std::vector< std::vector<Histo> > histos_;
  
};

#endif // DQM_SiStripCommissioningAnalysis_OptoScanAlgorithm_H

