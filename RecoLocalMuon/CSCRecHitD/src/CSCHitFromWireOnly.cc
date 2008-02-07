/* This is CSCHitFromWireOnly
 *
 * Finds wiregroup with hits, and fill in CSCWireHitCollection
 * which includes only DetId and wiregroup #
 *
 * \author  Dominique Fortin - UCR
 *
 */
//---- Some changes from Stoyan Stoynev - NU 

#include <RecoLocalMuon/CSCRecHitD/interface/CSCHitFromWireOnly.h>
#include <RecoLocalMuon/CSCRecHitD/interface/CSCWireHit.h>

#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>

#include <DataFormats/CSCDigi/interface/CSCWireDigi.h>

#include <FWCore/MessageLogger/interface/MessageLogger.h>

#include <iostream>


CSCHitFromWireOnly::CSCHitFromWireOnly( const edm::ParameterSet& ps ) {
  
  debug                  = ps.getUntrackedParameter<bool>("CSCDebug");
  deltaT                 = ps.getUntrackedParameter<int>("CSCWireClusterDeltaT");
  //clusterSize            = ps.getUntrackedParameter<int>("CSCWireClusterMaxSize");
}


CSCHitFromWireOnly::~CSCHitFromWireOnly(){}


/* runWire
 *
 */
std::vector<CSCWireHit> CSCHitFromWireOnly::runWire( const CSCDetId& id, const CSCLayer* layer, const CSCWireDigiCollection::Range& rwired ) {
  
  std::vector<CSCWireHit> hitsInLayer;
    
  layer_ = layer;
  layergeom_ = layer->geometry();
  bool any_digis = true;
  int n_wgroup = 0;


  // Loop over wire digi collection
  for ( CSCWireDigiCollection::const_iterator it = rwired.first; it != rwired.second; ++it ) {
    
    const CSCWireDigi wdigi = *it;


    if ( any_digis ) {
      any_digis = false;
      makeWireCluster( wdigi );
      n_wgroup = 1;
    } else {
      if ( !addToCluster( wdigi ) ) {
	// Make Wire Hit from cluster, delete old cluster and start new one
	float whit_pos = findWireHitPosition();
	CSCWireHit whit(id, whit_pos, wire_in_cluster, theTime);
	hitsInLayer.push_back( whit );	
	makeWireCluster( wdigi );
	n_wgroup = 1;
      } else {
	n_wgroup++;
      }
    }
    // Don't forget to fill last wire hit !!!
    if ( rwired.second - it == 1) {           
      float whit_pos = findWireHitPosition();
      CSCWireHit whit(id, whit_pos, wire_in_cluster, theTime);
      hitsInLayer.push_back( whit );
      n_wgroup++;
    }
  }

  return hitsInLayer;
}


/* makeWireCluster
 *
 */
void CSCHitFromWireOnly::makeWireCluster(const CSCWireDigi & digi) {
  wire_cluster.clear();
  wire_in_cluster.clear();
  theLastChannel  = digi.getWireGroup();
  theTime         = digi.getTimeBin();
  wire_cluster.push_back( digi );
}


/* addToCluster
 *
 */
bool CSCHitFromWireOnly::addToCluster(const CSCWireDigi & digi) {

  
  int iwg = digi.getWireGroup();
  
  if ( iwg == theLastChannel ){
    return true;  // Same wire group but different tbin -> ignore
  }
  else{
    if ( (iwg == theLastChannel+1) && (abs(digi.getTimeBin()-theTime)<= deltaT) ) {
      theLastChannel = iwg;
      wire_cluster.push_back( digi );
      return true;
    }
  }
  
  return false;
}


/* findWireHitPosition
 *
 * This position is expressed in terms of wire #... is a float since it may be a fraction.
 */
float CSCHitFromWireOnly::findWireHitPosition() {
  
  // Again use center of mass to determine position of wire hit
  // To do so, need to know wire spacing and # of wires
  
  float y = 0.0;
  
  for ( unsigned i = 0; i < wire_cluster.size(); ++i ) {
    CSCWireDigi wdigi = wire_cluster[i];
    int wgroup = wdigi.getWireGroup();
    wire_in_cluster.push_back( wgroup );
    y += float( wgroup );
  }       

  float wiregpos = y /wire_cluster.size() ;

  return wiregpos;

}

