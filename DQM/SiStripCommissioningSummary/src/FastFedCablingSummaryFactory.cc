#include "DQM/SiStripCommissioningSummary/interface/FastFedCablingSummaryFactory.h"
#include "CondFormats/SiStripObjects/interface/CommissioningAnalysis.h"
#include "DataFormats/SiStripCommon/interface/SiStripEnumsAndStrings.h"
#include "DQM/SiStripCommissioningSummary/interface/SummaryGenerator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>
#include <sstream>

using namespace sistrip;

// -----------------------------------------------------------------------------
//
uint32_t FastFedCablingSummaryFactory::init( const sistrip::Monitorable& mon, 
					     const sistrip::Presentation& pres,
					     const sistrip::View& view, 
					     const std::string& level, 
					     const sistrip::Granularity& gran,
					     const std::map<uint32_t,CommissioningAnalysis*>& data ) {
  
// // -----------------------------------------------------------------------------
// //
// uint32_t SummaryPlotFactory<FastFedCablingAnalysis*>::init( const sistrip::Monitorable& mon, 
// 							    const sistrip::Presentation& pres,
// 							    const sistrip::View& view, 
// 							    const std::string& level, 
// 							    const sistrip::Granularity& gran,
// 							    const std::map<uint32_t,FastFedCablingAnalysis*>& data ) {
  
//   // Some initialisation
//   SummaryPlotFactoryBase::init( mon, pres, view, level, gran );
  
//   // Check if generator object exists
//   if ( !SummaryPlotFactoryBase::generator_ ) { return 0; }
  
//   // Extract monitorable
//   std::map<uint32_t,FastFedCablingAnalysis*>::const_iterator iter = data.begin();
//   for ( ; iter != data.end(); iter++ ) {
//     if ( !iter->second ) { continue; }
//     float value = 1. * sistrip::invalid_;
//     float error = 1. * sistrip::invalid_;
//     if ( SummaryPlotFactoryBase::mon_ == sistrip::FAST_CABLING_HIGH_LEVEL ) { 
//       value = iter->second->highLevel(); 
//       error = iter->second->highRms(); 
//     } else if ( SummaryPlotFactoryBase::mon_ == sistrip::FAST_CABLING_LOW_LEVEL ) { 
//       value = iter->second->lowLevel(); 
//       error = iter->second->lowRms(); 
// //     } else if ( SummaryPlotFactoryBase::mon_ == sistrip::FAST_CABLING_ADC_RANGE ) { 
// //       value = iter->second->highLevel() - iter->second->lowLevel(); 
//     } else if ( SummaryPlotFactoryBase::mon_ == sistrip::FAST_CABLING_MAX ) { 
//       value = iter->second->max(); 
//     } else if ( SummaryPlotFactoryBase::mon_ == sistrip::FAST_CABLING_MIN ) { 
//       value = iter->second->min(); 
//     } else if ( SummaryPlotFactoryBase::mon_ == sistrip::FAST_CABLING_CONNS_PER_FED ) { 
//       value = 1. * static_cast<uint16_t>( iter->second->isValid() ); 
//     } else { 
//       edm::LogWarning(mlSummaryPlots_)
// 	<< "[SummaryPlotFactory::" << __func__ << "]" 
// 	<< " Unexpected monitorable: "
// 	<< SiStripEnumsAndStrings::monitorable( SummaryPlotFactoryBase::mon_ );
//       continue; 
//     }
    
//     SummaryPlotFactoryBase::generator_->fillMap( SummaryPlotFactoryBase::level_, 
// 						 SummaryPlotFactoryBase::gran_, 
// 						 iter->first, 
// 						 value,
// 						 error );

//   }
  
  return SummaryPlotFactoryBase::generator_->nBins();
  
}

// -----------------------------------------------------------------------------
//
void FastFedCablingSummaryFactory::fill( TH1& summary_histo ) {

// //------------------------------------------------------------------------------
// //
// void SummaryPlotFactory<FastFedCablingAnalysis*>::fill( TH1& summary_histo ) {
  
//   // Histogram filling and formating
//   SummaryPlotFactoryBase::fill( summary_histo );
  
//   if ( !SummaryPlotFactoryBase::generator_ ) { return; }
  
//   // Histogram formatting
//   if ( SummaryPlotFactoryBase::mon_ == sistrip::FAST_CABLING_HIGH_LEVEL ) {
//     SummaryPlotFactoryBase::generator_->axisLabel( "\"High\" light level [ADC]" );
//   } else if ( SummaryPlotFactoryBase::mon_ == sistrip::FAST_CABLING_LOW_LEVEL ) {
//     SummaryPlotFactoryBase::generator_->axisLabel( "\"Low\" light level [ADC]" );
//   } else if ( SummaryPlotFactoryBase::mon_ == sistrip::FAST_CABLING_MAX ) {
//     SummaryPlotFactoryBase::generator_->axisLabel( "Maximum light level [ADC]" );
//   } else if ( SummaryPlotFactoryBase::mon_ == sistrip::FAST_CABLING_MIN ) {
//     SummaryPlotFactoryBase::generator_->axisLabel( "Minumum light level [ADC]" );
//   } else if ( SummaryPlotFactoryBase::mon_ == sistrip::FAST_CABLING_CONNS_PER_FED ) { 
//     SummaryPlotFactoryBase::generator_->axisLabel( "Connected channels per FED" );
//   } else { 
//     edm::LogWarning(mlSummaryPlots_)
//       << "[SummaryPlotFactory::" << __func__ << "]" 
//       << " Unexpected SummaryHisto value:"
//       << SiStripEnumsAndStrings::monitorable( SummaryPlotFactoryBase::mon_ );
//   } 
  
}

// -----------------------------------------------------------------------------
//
//template class FastFedCablingSummaryFactory;//<CommissiniongAnalysis*>;

