// Read in strip digi collection and apply calibrations to ADC counts


#include <RecoLocalMuon/CSCRecHitB/src/CSCStripGain.h>

#include <FWCore/Utilities/interface/Exception.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>

#include <DataFormats/MuonDetId/interface/CSCDetId.h>

#include <CondFormats/CSCObjects/interface/CSCGains.h>
#include <CondFormats/DataRecord/interface/CSCGainsRcd.h>
#include <CondFormats/CSCObjects/interface/CSCReadoutMappingFromFile.h>
#include <CondFormats/CSCObjects/interface/CSCReadoutMappingForSliceTest.h>

#include <map>
#include <vector>

CSCStripGain::CSCStripGain( const edm::ParameterSet & ps ) {

  debug                  = ps.getUntrackedParameter<bool>("CSCDebug");  
  theCSCMap              = CSCReadoutMappingFromFile( ps );                                              
}

CSCStripGain::~CSCStripGain() {

}


/* getStripWeights
 *
 */
void CSCStripGain::getStripGain( const CSCDetId& id, float* weights ) {

  // Compute channel id used for retrieving gains from database
  bool isME1a = false;
  int ec = id.endcap();
  int st = id.station();
  int rg = id.ring();
  int ch = id.chamber();
  int la = id.layer();

  // Note that ME-1a constants are stored in ME-11 (ME-1b)
  if (id.station() == 1 && id.ring() == 4 ) {
    rg = 1;
    isME1a = true;
  }

  int chId=220000000 + ec*100000 + st*10000 + rg*1000 + ch*10 + la;

  int nStrips = 80;
  if ( st == 1 && rg == 1) nStrips = 64;
  if ( st == 1 && rg == 3) nStrips = 64;

  // Note that ME1/a constants are stored in ME1/1 (ME1/b) starting at entry 64
  if ( st == 1  && rg == 4 ) {
    rg = 1;
    isME1a = true;
  }
  
  int me1a_id = 0;
  int strip1 = 0;
  
  if (isME1a) strip1 = 64;

  float w = 1.0;

//  CSCGains *Gains_n = const_cast<CSCGains*> (Gains);

  for ( int sid = strip1; sid < nStrips; sid++ ) {
    if (Gains_->gains.find(chId) != Gains_->gains.end( ) ) {
      w = globalGainAvg/Gains_->gains[chId][sid].gain_slope;
    } else {
      w = 1.0;
    }
    if (w > 2.0) w = 2.0;
    if (w < 0.5) w = 0.5;
 
    if ( !isME1a ) {
      weights[sid] = w;
    } else {
      if (sid > 63 ) {
        weights[me1a_id]    = w;
        weights[me1a_id+16] = w;
        weights[me1a_id+32] = w;
        me1a_id++;
      }
    }
  }
}
