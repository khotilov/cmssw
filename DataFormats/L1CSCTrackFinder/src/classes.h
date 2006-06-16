#include <vector>
#include <DataFormats/L1CSCTrackFinder/interface/L1TrackId.h>
#include <DataFormats/L1CSCTrackFinder/interface/L1Track.h>
#include <DataFormats/Common/interface/Wrapper.h>

namespace
{
  namespace
    {
      csc::L1Track cL1TRK;
      csc::L1TrackId cL1TRKID;

      std::pair<csc::L1Track,MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi> > vL1Trk;

      std::vector<csc::L1Track> vL1TRK;
      std::vector<csc::L1TrackId> vL1TRKID;
      std::vector<std::pair<csc::L1Track,MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi> > > vL1TrkColl;
      
      edm::Wrapper<std::vector<csc::L1Track> > wL1TRK;
      edm::Wrapper<std::vector<std::pair<csc::L1Track,MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi> > > > wL1TrkColl;
    }
}
