#ifndef DataFormats_SiStripDigi_Classes_H
#define DataFormats_SiStripDigi_Classes_H

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include <vector>
#include <string>

#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
namespace {
  namespace {
    edm::Wrapper<SiStripDigi > zs0;
    edm::Wrapper<std::vector<SiStripDigi> > zs1;
    edm::Wrapper<edm::DetSet<SiStripDigi> > zs2;
    edm::Wrapper<std::vector<edm::DetSet<SiStripDigi> > > zs3;
    edm::Wrapper<edm::DetSetVector<SiStripDigi> > zs4;
  }
}

#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
namespace {
  namespace {
    edm::Wrapper<SiStripRawDigi > raw0;
    edm::Wrapper<std::vector<SiStripRawDigi> > raw1;
    edm::Wrapper<edm::DetSet<SiStripRawDigi> > raw2;
    edm::Wrapper<std::vector<edm::DetSet<SiStripRawDigi> > > raw3;
    edm::Wrapper<edm::DetSetVector<SiStripRawDigi> > raw4;
  }
}

#include "DataFormats/SiStripDigi/interface/SiStripDigiCollection.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include "boost/cstdint.hpp"
namespace {
  namespace {

    edm::Wrapper<SiStripDigiCollection> digis;

    edm::Wrapper<std::vector<bool> > list_of_fed_id;

    edm::Wrapper<uint8_t*> ptr_to_buffer;
    edm::Wrapper<std::vector<uint8_t*> > ptrs_to_buffers;

    edm::Wrapper<uint32_t> size_of_buffer;
    edm::Wrapper<std::vector<uint32_t> > size_of_buffers;

    edm::Wrapper<sistrip::FedReadoutMode> fed_readout_mode;
    edm::Wrapper<std::vector<sistrip::FedReadoutMode> > fed_readout_modes;

    edm::Wrapper<sistrip::FedReadoutPath> fed_readout_path;
    edm::Wrapper<std::vector<sistrip::FedReadoutPath> > fed_readout_paths;

    edm::Wrapper<uint16_t> payload_position;
    edm::Wrapper<std::vector<uint16_t> > payload_positions;
    edm::Wrapper<std::vector<std::vector<uint16_t> > > more_payload_positions;

  }
}

#include "DataFormats/SiStripDigi/interface/SiStripEventSummary.h"
namespace {
  namespace {
    edm::Wrapper< sistrip::RunType > run_type; 
    edm::Wrapper< sistrip::ApvReadoutMode > apv_mode; 
    edm::Wrapper< sistrip::FedReadoutMode > fed_mode; 
    edm::Wrapper< SiStripEventSummary > summary;
  }
}

#endif // DataFormats_SiStripDigi_Classes_H


 
