// Last commit: $Id: SiStripGainBuilderFromDb.cc,v 1.11 2008/07/17 10:27:59 giordano Exp $
// Latest tag:  $Name: V03-00-03 $
// Location:    $Source: /cvs_server/repositories/CMSSW/CMSSW/OnlineDB/SiStripESSources/src/SiStripGainBuilderFromDb.cc,v $

#include "OnlineDB/SiStripESSources/interface/SiStripGainBuilderFromDb.h"
#include "CondFormats/SiStripObjects/interface/SiStripApvGain.h"
#include "OnlineDB/SiStripESSources/interface/SiStripCondObjBuilderFromDb.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <vector>
#include <string>

using namespace std;
using namespace sistrip;

// -----------------------------------------------------------------------------
/** */
SiStripGainBuilderFromDb::SiStripGainBuilderFromDb( const edm::ParameterSet& pset ) 
  : SiStripGainESSource( pset )
{
  LogTrace(mlESSources_) 
    << "[SiStripGainBuilderFromDb::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
/** */
SiStripGainBuilderFromDb::~SiStripGainBuilderFromDb() {
  LogTrace(mlESSources_)
    << "[SiStripGainBuilderFromDb::" << __func__ << "]"
    << " Destructing object...";
}

// -----------------------------------------------------------------------------
/** */
SiStripApvGain* SiStripGainBuilderFromDb::makeGain() {
  LogTrace(mlESSources_) 
    << "[SiStripGainBuilderFromDb::" << __func__ << "]"
    << " Constructing Gain object...";

  /** Service to access onlineDB and extract pedestal/gain */
  edm::Service<SiStripCondObjBuilderFromDb> condObjBuilder_;
  
  // Create Gain object 
  SiStripApvGain* gain;
  condObjBuilder_->getValue(gain);
  return gain;
  
}

