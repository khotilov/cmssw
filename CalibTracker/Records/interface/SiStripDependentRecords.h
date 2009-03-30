#ifndef CALIBTRACKER_RECORDS_SISTRIPDEPENDENTRECORDS_H
#define CALIBTRACKER_RECORDS_SISTRIPDEPENDENTRECORDS_H

#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "boost/mpl/vector.hpp"

#include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

class SiStripFecCablingRcd : public edm::eventsetup::DependentRecordImplementation<SiStripFecCablingRcd,
  boost::mpl::vector<SiStripFedCablingRcd> > {};

class SiStripDetCablingRcd : public edm::eventsetup::DependentRecordImplementation<SiStripDetCablingRcd,
  boost::mpl::vector<SiStripFedCablingRcd> > {};

class SiStripRegionCablingRcd : public edm::eventsetup::DependentRecordImplementation<SiStripRegionCablingRcd,
  boost::mpl::vector<SiStripDetCablingRcd,TrackerDigiGeometryRecord> > {};

class SiStripGainRcd : public edm::eventsetup::DependentRecordImplementation<SiStripGainRcd, boost::mpl::vector<SiStripApvGainRcd> > {};
class SiStripGainSimRcd : public edm::eventsetup::DependentRecordImplementation<SiStripGainSimRcd, boost::mpl::vector<SiStripApvGainSimRcd> > {};

class SiStripQualityRcd : public edm::eventsetup::DependentRecordImplementation<SiStripQualityRcd, boost::mpl::vector<SiStripBadModuleRcd, SiStripBadFiberRcd, SiStripBadChannelRcd, SiStripBadStripRcd, SiStripDetCablingRcd, SiStripDCSStatusRcd, SiStripDetVOffRcd> > {};


class SiStripHashedDetIdRcd : public edm::eventsetup::DependentRecordImplementation<SiStripHashedDetIdRcd, boost::mpl::vector<TrackerDigiGeometryRecord> > {};

#endif 

