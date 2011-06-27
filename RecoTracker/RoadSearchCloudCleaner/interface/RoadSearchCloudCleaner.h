#ifndef RoadSearchCloudCleaner_h
#define RoadSearchCloudCleaner_h

//
// Package:         RecoTracker/RoadSearchCloudCleaner
// Class:           RoadSearchCloudCleaner
// 
// Description:     Calls RoadSeachCloudCleanerAlgorithm
//                  to find RoadSearchClouds.
//
// Original Author: Steve Wagner, stevew@pizero.colorado.edu
// Created:         Sat Feb 19 22:00:00 UTC 2006
//
// $Author: gutsche $
// $Date: 2006/03/28 23:12:10 $
// $Revision: 1.1 $
//

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/EDProduct.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoTracker/RoadSearchCloudCleaner/interface/RoadSearchCloudCleanerAlgorithm.h"

namespace cms
{
  class RoadSearchCloudCleaner : public edm::EDProducer
  {
  public:

    explicit RoadSearchCloudCleaner(const edm::ParameterSet& conf);

    virtual ~RoadSearchCloudCleaner();

    virtual void produce(edm::Event& e, const edm::EventSetup& c);

  private:
    RoadSearchCloudCleanerAlgorithm roadSearchCloudCleanerAlgorithm_;
    edm::ParameterSet conf_;

  };
}


#endif
