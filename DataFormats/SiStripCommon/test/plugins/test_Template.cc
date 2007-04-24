// Last commit: $Id: test_Template.cc,v 1.1 2007/04/04 07:35:50 bainbrid Exp $

#include "DataFormats/SiStripCommon/test/plugins/test_Template.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <boost/cstdint.hpp>
#include <iostream>
#include <sstream>

using namespace sistrip;

// -----------------------------------------------------------------------------
// 
test_Template::test_Template( const edm::ParameterSet& pset ) 
{
  LogTrace(mlTest_)
    << "[test_Template::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
// 
test_Template::~test_Template() {
  LogTrace(mlTest_)
    << "[test_Template::" << __func__ << "]"
    << " Destructing object...";
}

// -----------------------------------------------------------------------------
// 
void test_Template::beginJob( const edm::EventSetup& setup ) {
  
  std::stringstream ss;
  ss << "[test_Template::" << __func__ << "]"
     << " Initializing...";
  LogTrace(mlTest_) << ss.str();

}

// -----------------------------------------------------------------------------
// 
void test_Template::analyze( const edm::Event& event, 
			     const edm::EventSetup& setup ) {
  LogTrace(mlTest_) 
    << "[test_Template::" << __func__ << "]"
    << " Analyzing run/event "
    << event.id().run() << "/"
    << event.id().event();
}


