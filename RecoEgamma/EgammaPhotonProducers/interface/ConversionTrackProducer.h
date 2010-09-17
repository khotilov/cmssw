#ifndef ConversionTrackProducer_h
#define ConversionTrackProducer_h

//
// Package:         RecoTracker/FinalTrackSelectors
// Class:           ConversionTrackProducer
// 
// Description:     Hit Dumper
//
// Original Author: Steve Wagner, stevew@pizero.colorado.edu
// Created:         Sat Jan 14 22:00:00 UTC 2006
//
// $Author: vlimant $
// $Date: 2009/05/14 15:06:08 $
// $Revision: 1.6 $
//

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/EDProduct.h"

#include "DataFormats/EgammaTrackReco/interface/ConversionTrack.h"
#include "DataFormats/EgammaTrackReco/interface/ConversionTrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

  class ConversionTrackProducer : public edm::EDProducer
  {
  public:

    explicit ConversionTrackProducer(const edm::ParameterSet& conf);

    virtual ~ConversionTrackProducer();

    virtual void produce(edm::Event& e, const edm::EventSetup& c);

  private:
    edm::ParameterSet conf_;

    std::auto_ptr<reco::ConversionTrackCollection> outputTrks;
  };
#endif
