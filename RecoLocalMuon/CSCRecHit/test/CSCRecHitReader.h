#ifndef RecoLocalMuon_CSCRecHitReader_H
#define RecoLocalMuon_CSCRecHitReader_H

/** \class CSCRecHitReader
 *  Basic analyzer class which accesses 2D CSCRecHits
 *  and plot resolution comparing them with muon simhits
 *
 *  Author: D. Fortin  - UC Riverside
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Handle.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include <Geometry/CSCGeometry/interface/CSCChamber.h>
#include <Geometry/CSCGeometry/interface/CSCLayer.h>

#include "CSCRecHitHistograms.h"

#include <vector>
#include <map>
#include <string>

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

class PSimHit;
class TFile;
class CSCLayer;
class CSCDetId;

class CSCRecHitReader : public edm::EDAnalyzer {
public:
  /// Constructor
  CSCRecHitReader(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~CSCRecHitReader();

  // Operations

  /// Perform the real analysis
  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);


protected:

private: 

  // Histograms
  H2DRecHit *hRHPAll;
  H2DRecHit *hRHPME1a;
  H2DRecHit *hRHPME1b;
  H2DRecHit *hRHPME12;
  H2DRecHit *hRHPME13;
  H2DRecHit *hRHPME21;
  H2DRecHit *hRHPME22;
  H2DRecHit *hRHPME31;
  H2DRecHit *hRHPME32;
  H2DRecHit *hRHPME4;


  // The file which will store the histos
  TFile *theFile;
  // Switch for debug output
  bool debug;
  // Root file name
  std::string rootFileName;
  std::string simHitLabel;
  std::string recHitLabel;
  int minRechitChamber;
  int maxRechitChamber;
  int WhichEndCap;

};


#endif




