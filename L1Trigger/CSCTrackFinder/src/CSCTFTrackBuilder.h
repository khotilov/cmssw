#ifndef CSCTrackFinder_CSCTFTrackBuilder_h
#define CSCTrackFinder_CSCTFTrackBuilder_h

#include <vector>
#include <DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h>
#include <DataFormats/L1CSCTrackFinder/interface/CSCTriggerContainer.h>
#include <DataFormats/L1CSCTrackFinder/interface/TrackStub.h>
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h>
#include <DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h>
#include <L1Trigger/CSCTrackFinder/src/CSCTFDTReceiver.h>
///KK
#include <string.h>
#include <FWCore/Framework/interface/EventSetup.h>
///

#include <FWCore/ParameterSet/interface/ParameterSet.h>

class CSCMuonPortCard;
class CSCTFSectorProcessor;

class CSCTFTrackBuilder
{
 public:
///KK
  void initialize(const edm::EventSetup& c);
///

  enum { nEndcaps = 2, nSectors = 6};

  CSCTFTrackBuilder(const edm::ParameterSet& pset);
  ~CSCTFTrackBuilder();

  void buildTracks(const CSCCorrelatedLCTDigiCollection*, const L1MuDTChambPhContainer*,
		   L1CSCTrackCollection*, CSCTriggerContainer<csctf::TrackStub>*);

 private:

  CSCTFDTReceiver* my_dtrc;
  CSCTFSectorProcessor* my_SPs[nEndcaps][nSectors];
};

#endif
