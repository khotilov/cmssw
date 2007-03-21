// Last commit: $Id: $

#ifndef DataFormats_SiStripCommon_ConstantsForHardwareSystems_H
#define DataFormats_SiStripCommon_ConstantsForHardwareSystems_H

#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/SiStripCommon/interface/Constants.h"
#include "boost/cstdint.hpp"

/**
   @file ConstantsForHardwareSystems.h
   @brief Constants and enumerated types for FED/FEC systems
*/

namespace sistrip { 

  // -------------------- FED ids --------------------
  
  static const uint16_t FED_ID_MIN     = static_cast<uint16_t>( FEDNumbering::getSiStripFEDIds().first );
  static const uint16_t FED_ID_MAX     = static_cast<uint16_t>( FEDNumbering::getSiStripFEDIds().second );
  static const uint16_t CMS_FED_ID_MAX = static_cast<uint16_t>( FEDNumbering::lastFEDId() );
  static const uint16_t NUMBER_OF_FEDS = static_cast<uint16_t>( FED_ID_MAX - FED_ID_MIN );
  
  // -------------------- FEDs to channels --------------------

  static const uint16_t FEDCH_PER_FEUNIT = 12;
  static const uint16_t FEUNITS_PER_FED  = 8;
  static const uint16_t FEDCH_PER_FED    = FEDCH_PER_FEUNIT * FEUNITS_PER_FED; // 96

  // -------------------- Front-end devices --------------------

  static const uint16_t APVS_PER_FEDCH   = 2;
  static const uint16_t APVS_PER_FEUNIT  = APVS_PER_FEDCH * FEDCH_PER_FEUNIT; // 24
  static const uint16_t APVS_PER_FED     = APVS_PER_FEUNIT * FEUNITS_PER_FED; // 194

  static const uint16_t APVS_PER_CHAN    = 2;
  static const uint16_t CHANS_PER_LLD    = 3;
  
  // -------------------- Detector strips -------------------- 

  static const uint16_t STRIPS_PER_APV    = 128;
  static const uint16_t STRIPS_PER_FEDCH  = STRIPS_PER_APV * APVS_PER_FEDCH;
  static const uint16_t STRIPS_PER_FEUNIT = STRIPS_PER_FEDCH * FEDCH_PER_FEUNIT; // 3072
  static const uint16_t STRIPS_PER_FED    = STRIPS_PER_FEUNIT * FEUNITS_PER_FED; // 24576

  // -------------------- FED buffers --------------------

  static const uint16_t DAQ_HDR_SIZE        = 8;
  static const uint16_t TRK_HDR_SIZE        = 8;
  static const uint16_t FE_HDR_SIZE         = 16;
  static const uint16_t APV_ERROR_HDR_SIZE  = 24;
  static const uint16_t FULL_DEBUG_HDR_SIZE = 8 * FE_HDR_SIZE;

  // -------------------- Control system --------------------

  static const uint16_t FEC_RING_MIN    =   1;
  static const uint16_t FEC_RING_MAX    =   8;

  static const uint16_t CCU_ADDR_MIN    =   1;
  static const uint16_t CCU_ADDR_MAX    = 127;

  static const uint16_t CCU_CHAN_MIN    =  16;
  static const uint16_t CCU_CHAN_MAX    =  31;

  static const uint16_t LLD_CHAN_MIN    =   1;
  static const uint16_t LLD_CHAN_MAX    =   3;

  static const uint16_t APV_I2C_MIN     =  32;
  static const uint16_t APV_I2C_MAX     =  37;

  // -------------------- VME crates --------------------

  static const uint16_t SLOTS_PER_CRATE    =  20;

  static const uint16_t CRATE_SLOT_MIN     =   2; // slot 1 is reserved for VME controller
  static const uint16_t CRATE_SLOT_MAX     =  21;

  static const uint16_t MAX_FEDS_PER_CRATE =  16;
  static const uint16_t MAX_FECS_PER_CRATE =  20;

  static const uint16_t FED_CRATE_MIN      =   1;
  static const uint16_t FED_CRATE_MAX      =  60;
  
  static const uint16_t FEC_CRATE_MIN      =   1;
  static const uint16_t FEC_CRATE_MAX      =   4;
  


  // -------------------- String constants -------------------- 

  static const std::string unknownApvReadoutMode_   = "UnknownApvReadoutMode";
  static const std::string undefinedApvReadoutMode_ = "UndefinedApvReadoutMode";

  static const std::string apvPeakMode_ = "ApvPeakMode";
  static const std::string apvDeconMode_ = "ApvDeconMode";
  static const std::string apvMultiMode_ = "ApvMultiMode";

  static const std::string unknownFedReadoutMode_   = "UnknownFedReadoutMode";
  static const std::string undefinedFedReadoutMode_ = "UndefinedFedReadoutMode";

  static const std::string fedScopeMode_     = "FedScopeMode";
  static const std::string fedVirginRaw_     = "FedVirginRaw";
  static const std::string fedProcRaw_       = "FedProcessedRaw";
  static const std::string fedZeroSuppr_     = "FedZeroSuppressed";
  static const std::string fedZeroSupprLite_ = "FedZeroSupprressedLite";
  
  // -------------------- Enumerators --------------------
  
  enum ApvReadoutMode { UNKNOWN_APV_READOUT_MODE = sistrip::unknown_,
			UNDEFINED_APV_READOUT_MODE = sistrip::invalid_,
			APV_PEAK_MODE = 1, 
			APV_DECON_MODE = 2, 
			APV_MULTI_MODE = 3
  };
  
  enum FedReadoutMode { UNKNOWN_FED_READOUT_MODE = sistrip::unknown_,
			UNDEFINED_FED_READOUT_MODE = sistrip::invalid_,
			FED_SCOPE_MODE = 1, 
			FED_VIRGIN_RAW = 2, 
			FED_PROC_RAW = 6, 
			FED_ZERO_SUPPR = 10,
			FED_ZERO_SUPPR_LITE = 12
  };

  enum FedReadoutPath { UNKNOWN_FED_READOUT_PATH = sistrip::unknown_,
			UNDEFINED_FED_READOUT_PATH = sistrip::invalid_,
			VME_READOUT = 1, 
			SLINK_READOUT = 2 
  };
  
  enum FedBufferFormat { UNKNOWN_FED_BUFFER_FORMAT = sistrip::unknown_,
			 UNDEFINED_FED_BUFFER_FORMAT = sistrip::invalid_,
			 FULL_DEBUG_FORMAT = 1, 
			 APV_ERROR_FORMAT = 2 
  };
  
  enum FedSuperMode { UNKNOWN_FED_SUPER_MODE = sistrip::unknown_,
		      UNDEFINED_FED_SUPER_MODE = sistrip::invalid_,
		      REAL = 0, 
		      FAKE = 1 
  };

}
  
#endif // DataFormats_SiStripCommon_ConstantsForHardwareSystems_H


