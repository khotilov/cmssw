#ifndef TrackHitPositions_H
#define TrackHitPositions_H

// system include files
//#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "CommonTools/TrackerMap/interface/TrackerMap.h"
#include "CalibTracker/SiStripCommon/interface/SiStripDetInfoFileReader.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DQM/SiStripCommon/interface/TkHistoMap.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"

#include <sstream>

class TrackHitPositions : public edm::EDAnalyzer {

 public:
  explicit TrackHitPositions( const edm::ParameterSet& );
  ~TrackHitPositions(){}

  void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob();

 private:

  void initializeMap(const edm::ESHandle<SiStripQuality> & SiStripQuality_, TrackerMap * tkMap);

  unsigned long long m_cacheID_;
  std::string dataLabel_;
  std::string TkMapFileName_;
  edm::FileInPath fp_;
  bool saveTkHistoMap_;
  //Global Info
  int NTkBadComponent[4]; //k: 0=BadModule, 1=BadFiber, 2=BadApv, 3=BadStrips
  int NBadComponent[4][19][4];
  //legend: NBadComponent[i][j][k]= SubSystem i, layer/disk/wheel j, BadModule/Fiber/Apv k
  //     i: 0=TIB, 1=TID, 2=TOB, 3=TEC
  //     k: 0=BadModule, 1=BadFiber, 2=BadApv, 3=BadStrips
  std::stringstream ssV[4][19];

  TrackerMap * tkMap, *tkMapFullIOVs;
  SiStripDetInfoFileReader* reader;
  TkHistoMap* tkhisto;
  double ptCut_;
};
#endif
