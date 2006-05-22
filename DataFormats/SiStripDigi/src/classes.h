#ifndef DataFormats_SiStripDigi_Classes_H
#define DataFormats_SiStripDigi_Classes_H

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include <vector>

#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
namespace {
  namespace {
    edm::Wrapper<SiStripDigi> zs0;
    edm::Wrapper< std::vector<SiStripDigi>  > zs1;
    edm::Wrapper< edm::DetSet<SiStripDigi> > zs2;
    edm::Wrapper< std::vector<edm::DetSet<SiStripDigi> > > zs3;
    edm::Wrapper< edm::DetSetVector<SiStripDigi> > zs4;
  }
}

#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
namespace {
  namespace {
    edm::Wrapper<SiStripRawDigi> raw0;
    edm::Wrapper< std::vector<SiStripRawDigi>  > raw1;
    edm::Wrapper< edm::DetSet<SiStripRawDigi> > raw2;
    edm::Wrapper< std::vector<edm::DetSet<SiStripRawDigi> > > raw3;
    edm::Wrapper< edm::DetSetVector<SiStripRawDigi> > raw4;
  }
}
    
#include "DataFormats/SiStripDigi/interface/SiStripEventSummary.h"
#include "DataFormats/SiStripCommon/interface/SiStripEnumeratedTypes.h"
namespace {
  namespace {
    edm::Wrapper<sistrip::Task> task;
    edm::Wrapper<sistrip::FedReadoutMode> fed_mode;
    edm::Wrapper<SiStripEventSummary> summary;

  }
}

#include "DataFormats/SiStripDigi/interface/Profile.h"
#include "DataFormats/SiStripDigi/interface/Histo.h"
namespace {
  namespace {
    edm::Wrapper<Profile> profile;
    edm::Wrapper< edm::DetSetVector<Profile> > profiles;
    edm::Wrapper<Histo> histo;
    edm::Wrapper< edm::DetSetVector<Histo> > histos;

  }
}

#endif // DataFormats_SiStripDigi_Classes_H


 
