// This is  CSCHitFromStripOnly.cc

#include <RecoLocalMuon/CSCRecHitB/src/CSCHitFromStripOnly.h>
#include <RecoLocalMuon/CSCRecHitB/src/CSCPeakBinOfStripPulse.h>
#include <RecoLocalMuon/CSCRecHitB/src/CSCStripGain.h>
#include <RecoLocalMuon/CSCRecHitB/src/CSCStripData.h>
#include <RecoLocalMuon/CSCRecHitB/src/CSCStripHitData.h>
#include <RecoLocalMuon/CSCRecHitB/interface/CSCStripHit.h>
//#include <RecoLocalMuon/CSCRecHitB/interface/CSCStripHitCollection.h>

#include <DataFormats/MuonDetId/interface/CSCDetId.h>
//#include <DataFormats/CSCDigi/interface/CSCStripDigi.h>

#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCChamberSpecs.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>

#include <CondFormats/CSCObjects/interface/CSCGains.h>
#include <CondFormats/DataRecord/interface/CSCGainsRcd.h>

#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/Exception.h>

#include <cmath>
#include <string>
#include <vector>
#include <iostream>


CSCHitFromStripOnly::CSCHitFromStripOnly( const edm::ParameterSet& ps ) {
  
  debug                      = ps.getUntrackedParameter<bool>("CSCDebug");
  useCalib                   = ps.getUntrackedParameter<bool>("CSCUseCalibrations");
  theClusterSize             = ps.getUntrackedParameter<int>("CSCStripClusterSize");
  theThresholdForAPeak       = ps.getUntrackedParameter<double>("CSCStripPeakThreshold");
  theThresholdForCluster     = ps.getUntrackedParameter<double>("CSCStripClusterChargeCut");
  useCleanStripCollection    = ps.getUntrackedParameter<bool>("CSCuseCleanStripCollection");

  pulseheightOnStripFinder_  = new CSCPeakBinOfStripPulse( ps );
  stripGain_                 = new CSCStripGain( ps );
  
}


CSCHitFromStripOnly::~CSCHitFromStripOnly() {
  delete pulseheightOnStripFinder_;
  delete stripGain_;
}


/* runStrip
 *
 * Search for strip with ADC output exceeding theThresholdForAPeak.  For each of these strips,
 * build a cluster of strip of size theClusterSize (typically 5 strips).  Finally, make
 * a Strip Hit out of these clusters by finding the center-of-mass position of the hit
 */
std::vector<CSCStripHit> CSCHitFromStripOnly::runStrip( const CSCDetId& id, const CSCLayer* layer,
                                                        const CSCStripDigiCollection::Range& rstripd, 
                                                        const CSCCLCTDigiCollection* clcts ) {	

  std::vector<CSCStripHit> hitsInLayer;
  
  // cache layer info for ease of access
  id_        = id;
  layer_     = layer;
  layergeom_ = layer_->geometry();
  specs_     = layer->chamber()->specs();  
  Nstrips    = specs_->nStrips();

  // Get strips which are part of the CLCT
  std::vector<int> clctStrips;
  if ( useCleanStripCollection ) {
    clctStrips = getCLCTStrips( id, clcts );
    if ( clctStrips.size() < 1) return hitsInLayer;
  }

  TmaxOfCluster = 5;

  // Get gains for cluster and compute correction weights:

  // Fill all the gain weights at once, so save calling this function many times.
  // N.B. in database, strip_id starts at 0, whereas it starts at 1 in detId
  // Initialize weights to 1. and crosstalk to 0., in case database isn't populated or using MC

  if ( useCalib ) {
    stripGain_->setCalibration( globalGainAvg, gains_ );
    stripGain_->getStripGain( id_, gainWeight );
  }
  
  // Fill adc map and find maxima (potential hits)
  fillPulseHeights( rstripd, clctStrips );
  findMaxima();      

  // Make a Strip Hit out of each strip local maximum
  for ( unsigned imax = 0; imax < theMaxima.size(); ++imax ) {

    // Initialize parameters entering the CSCStripHit
    ClusterSize = theClusterSize;
    theStrips.clear();
    strips_adc.clear();

    // This is where centroid position is determined
    // The strips_adc vector is also filled here
    // Remember, the array starts at 0, but the stripId starts at 1...
    float strippos = makeCluster( theMaxima[imax]+1 );    
    
    if ( strippos < 0 || TmaxOfCluster < 3 ) continue;
    
    CSCStripHit striphit( id, strippos, TmaxOfCluster, theStrips, strips_adc );
    hitsInLayer.push_back( striphit ); 
  }

  return hitsInLayer;
}


/* makeCluster
 *
 */
float CSCHitFromStripOnly::makeCluster( int centerStrip ) {
  
  float strippos = -1.;
  ClusterSize = theClusterSize;
  std::vector<CSCStripHitData> stripDataV;
 
  // We only want to use strip position in terms of strip # for the strip hit.
    
  // If the cluster size is such that you go beyond the edge of detector, shrink cluster appropriatly
  for ( int i = 1; i < theClusterSize/2 + 1; ++i) {
 
    if ( centerStrip - i < 1 || centerStrip + i > specs_->nStrips() ) {

      // Shrink cluster size, but keep it an odd number of strips.
      ClusterSize = 2*i - 1;  
    }
  }

  for ( int i = -ClusterSize/2; i <= ClusterSize/2; ++i ) {
    CSCStripHitData data = makeStripData(centerStrip, i);
    stripDataV.push_back( data );
    theStrips.push_back( centerStrip + i );
  }
  
  strippos = findHitOnStripPosition( stripDataV, centerStrip );
  
  return strippos;
}


/* makeStripData
 *
 * Gain corrections are applied here
 */
CSCStripHitData CSCHitFromStripOnly::makeStripData(int centerStrip, int offset) {
  
  CSCStripHitData prelimData(-1.,0.,0.,0.,0.,0);
  int thisStrip = centerStrip+offset;

  int tmax      = thePulseHeightMap[centerStrip-1].t();
  TmaxOfCluster = tmax;

  float adc[8];

  if ( tmax == 3 ) {
    adc[0] = thePulseHeightMap[thisStrip-1].y2();
    adc[1] = thePulseHeightMap[thisStrip-1].y3();
    adc[2] = thePulseHeightMap[thisStrip-1].y4();
    adc[3] = thePulseHeightMap[thisStrip-1].y5();
  } else if ( tmax == 4 ) {
    adc[0] = thePulseHeightMap[thisStrip-1].y3();
    adc[1] = thePulseHeightMap[thisStrip-1].y4();
    adc[2] = thePulseHeightMap[thisStrip-1].y5();
    adc[3] = thePulseHeightMap[thisStrip-1].y6();
  } else if ( tmax == 5 ) {
    adc[0] = thePulseHeightMap[thisStrip-1].y4();
    adc[1] = thePulseHeightMap[thisStrip-1].y5();
    adc[2] = thePulseHeightMap[thisStrip-1].y6();
    adc[3] = thePulseHeightMap[thisStrip-1].y7();
  } else if ( tmax == 6 ) {
    adc[0] = thePulseHeightMap[thisStrip-1].y5();
    adc[1] = thePulseHeightMap[thisStrip-1].y6();
    adc[2] = thePulseHeightMap[thisStrip-1].y7();
    adc[3] = 0.1;
  } else {
    adc[0] = 0.1;
    adc[1] = 0.1;
    adc[2] = 0.1;
    adc[3] = 0.1;
    std::cout << "[CSCHitFromStripOnly::makeStripData] Tmax out of range: contact CSC expert!!!" << std::endl;
  }
  
  if ( offset == 0 ) {
    prelimData = CSCStripHitData(thisStrip, adc[0], adc[1], adc[2], adc[3], TmaxOfCluster);
  } else {
    int sign = offset>0 ? 1 : -1;
    
    // If there's another maximum that would like to use part of this cluster, 
    // it gets shared in proportion to the height of the maxima
    for ( int i = 1; i <= ClusterSize/2; ++i ) {

      // Find the direction of the offset
      int testStrip = thisStrip + sign*i;

      std::vector<int>::iterator otherMax = find(theMaxima.begin(), theMaxima.end(), testStrip);

      // No other maxima found, so just store
      if ( otherMax == theMaxima.end() ) {
        prelimData = CSCStripHitData(thisStrip, adc[0], adc[1], adc[2], adc[3], TmaxOfCluster);
      } else {
        if ( tmax == 3 ) {
          adc[4] = thePulseHeightMap[testStrip-1].y2();
          adc[5] = thePulseHeightMap[testStrip-1].y3();
          adc[6] = thePulseHeightMap[testStrip-1].y4();
          adc[7] = thePulseHeightMap[testStrip-1].y5();
        } else if ( tmax == 4 ) {
          adc[4] = thePulseHeightMap[testStrip-1].y3();
          adc[5] = thePulseHeightMap[testStrip-1].y4();
          adc[6] = thePulseHeightMap[testStrip-1].y5();
          adc[7] = thePulseHeightMap[testStrip-1].y6();
         } else if ( tmax == 5 ) {
          adc[4] = thePulseHeightMap[testStrip-1].y4();
          adc[5] = thePulseHeightMap[testStrip-1].y5();
          adc[6] = thePulseHeightMap[testStrip-1].y6();
          adc[7] = thePulseHeightMap[testStrip-1].y7();
         } else if ( tmax == 6 ) {
          adc[4] = thePulseHeightMap[testStrip-1].y5();
          adc[5] = thePulseHeightMap[testStrip-1].y6();
          adc[6] = thePulseHeightMap[testStrip-1].y7();
          adc[7] = 0.1;
        } else {
          adc[4] = 0.1;
          adc[5] = 0.1;
          adc[6] = 0.1;
          adc[7] = 0.1;
        }
        for (int k = 0; k < 4; ++k) adc[k]   = adc[k] * adc[k] / (adc[k]+adc[k+4]);
        prelimData = CSCStripHitData(thisStrip, adc[0], adc[1], adc[2], adc[3], TmaxOfCluster);
      }
    }
  }
  return prelimData;
}
 

/* fillPulseHeights
 *
 */
void CSCHitFromStripOnly::fillPulseHeights( const CSCStripDigiCollection::Range& rstripd, std::vector<int> clctStrips ) {
  
  // Loop over strip digis and fill the pulseheight map
  
  thePulseHeightMap.clear();
  thePulseHeightMap.resize(100);

  // std::cout << "MaxStripADC " ;
  float maxADC = 0;  
  float maxADCCluster = 0;  
  float adc_cluster[3];  
  adc_cluster[0] = 0.;  
  adc_cluster[1] = 0.;  
  adc_cluster[2] = 0.;  

  for ( CSCStripDigiCollection::const_iterator it = rstripd.first; it != rstripd.second; ++it ) {
    bool fill            = true;
    int  thisChannel     = (*it).getStrip(); 
    std::vector<int> sca = (*it).getADCCounts();
    float height[6];
    float hmax = 0.;
    int   tmax = -1;

    // If using CLCT, test if strip is near CLCT strips (within +/-2 strips)
    if ( useCleanStripCollection ) {
      if ( !foundCLCTMatch( thisChannel, clctStrips) ) {
        // Don't forget that the ME_11/a strips are ganged !!!
        // Have to loop 2 more times to populate strips 17-48.
        if ( id_.station() == 1 && id_.ring() == 4 ) {
          for ( int j = 0; j < 3; ++j ) {
            thePulseHeightMap[thisChannel+16*j-1] = CSCStripData( float(thisChannel+16*j), hmax, tmax, 0., 0., 0., 0., 0., 0.);
          }
        } else {
          thePulseHeightMap[thisChannel-1] = CSCStripData( float(thisChannel), hmax, tmax, 0., 0., 0., 0., 0., 0.);
        }
        continue;
      }
    }

    
    fill = pulseheightOnStripFinder_->peakAboveBaseline( (*it), hmax, tmax, height );
    // std::cout << " " << hmax ;
    float maxCluster = 0;

    adc_cluster[2] = adc_cluster[1];
    adc_cluster[1] = adc_cluster[0];
    adc_cluster[0] = hmax;  

    for (int j = 0; j < 3; ++j ) maxCluster += adc_cluster[j];   
 
    if (hmax > maxADC) maxADC = hmax;
    if (maxCluster > maxADCCluster ) maxADCCluster = maxCluster;

    // Don't forget that the ME_11/a strips are ganged !!!
    // Have to loop 2 more times to populate strips 17-48.
    
    if ( id_.station() == 1 && id_.ring() == 4 ) {
      for ( int j = 0; j < 3; ++j ) {
        thePulseHeightMap[thisChannel+16*j-1] = CSCStripData( float(thisChannel+16*j), hmax, tmax, height[0], height[1], height[2], height[3], height[4], height[5]);
        if ( useCalib ) thePulseHeightMap[thisChannel+16*j-1] *= gainWeight[thisChannel-1];
      }
    } else {
      thePulseHeightMap[thisChannel-1] = CSCStripData( float(thisChannel), hmax, tmax, height[0], height[1], height[2], height[3], height[4], height[5]);
      if ( useCalib ) thePulseHeightMap[thisChannel-1] *= gainWeight[thisChannel-1];
    }
  }
  // std::cout << std::endl;
  // std::cout << "MAXADC " << maxADC << std::endl;
  // std::cout << "MAXADCCluster " << maxADCCluster << std::endl;
}


/* findMaxima
 *
 * fills vector theMaxima with the local maxima in the pulseheight distribution
 * of the strips. The threshold defining a maximum is a configurable parameter.
 * A typical value is 30.
 */
void CSCHitFromStripOnly::findMaxima() {
  
  theMaxima.clear();
  for ( size_t i = 0; i < thePulseHeightMap.size(); ++i ) {

    float heightPeak = thePulseHeightMap[i].ymax();

    // sum 3 strips so that hits between strips are not suppressed
    float heightCluster;
    if ( i == 0 ) {
      heightCluster = thePulseHeightMap[i].ymax()+thePulseHeightMap[i+1].ymax();
    } else if ( i == thePulseHeightMap.size()-1) {  
      heightCluster = thePulseHeightMap[i-1].ymax()+thePulseHeightMap[i].ymax();
    } else {
      heightCluster = thePulseHeightMap[i-1].ymax()+thePulseHeightMap[i].ymax()+thePulseHeightMap[i+1].ymax();
    }
   
    // Have found a strip Hit if...
    if (( heightPeak                   > theThresholdForAPeak         ) &&  
        ( heightCluster                > theThresholdForCluster       ) && 
        ( thePulseHeightMap[i].ymax()  > thePulseHeightMap[i-1].ymax()) &&
	( thePulseHeightMap[i].ymax() >= thePulseHeightMap[i+1].ymax()) &&
        ( thePulseHeightMap[i].t()     > 2                            ) &&
        ( thePulseHeightMap[i].t()     < 7                            )) {
      theMaxima.push_back(i);
    }
  }
}



/* findHitOnStripPosition
 *
 */
float CSCHitFromStripOnly::findHitOnStripPosition( const std::vector<CSCStripHitData>& data, const int& centerStrip ) {
  
  float strippos = -1.;
  
  int nstrips = data.size();
  if ( nstrips == 0 ) return strippos;
  
  // biggestStrip is strip with largest pulse height 
  // Use pointer subtraction

  int biggestStrip = max_element(data.begin(), data.end()) - data.begin();
  strippos = data[biggestStrip].x() * 1.;
  
  // If more than one strip:  use centroid to find center of cluster
  // but only use time bin == tmax (otherwise, bias centroid).
  float sum  = 0.;
  float sum_w= 0.;

   
  // CFEB trigger problem:
  // Make sure there are non-zero entries on either side of the central strip
  // otherwise, centroid will be biased.
  float stripLeft = 0.;
  float stripRight = 0.;

  
  for ( unsigned i = 0; i != data.size(); ++i ) {
    float w0 = data[i].y0();
    float w1 = data[i].y1();
    float w2 = data[i].y2();
    float w3 = data[i].y3();

    if (w0 < 0.) w0 = 0.001;
    if (w1 < 0.) w1 = 0.001;
    if (w2 < 0.) w2 = 0.001;
    if (w3 < 0.) w3 = 0.001;
    

    if (i == data.size()/2 -1) stripLeft = w1;
    if (i == data.size()/2 +1) stripRight = w1;


     // Fill the adcs to the strip hit --> needed for Gatti fitter
    strips_adc.push_back( w0 );
    strips_adc.push_back( w1 );
    strips_adc.push_back( w2 );
    strips_adc.push_back( w3 );
 
    sum_w += w1;
    sum   += w1 * data[i].x();
  }

  // CFEB trigger problem:
  // Test that we have readout entries on either side of the central strip
  // What's the minimum ADC count ???  Try 2 for now..
  if (int(strippos)%16 < 2) {
    if ( stripLeft > 2. && stripRight > 2. && sum_w > 0.) strippos = sum / sum_w;
  } else if ( sum_w > 0.) {
    strippos = sum / sum_w;    
  } 

  return strippos;
}


/*
 * Get strips which are part of CLCT(s) for this DetId (layer)
 */
std::vector<int> CSCHitFromStripOnly::getCLCTStrips( const CSCDetId& id, const CSCCLCTDigiCollection* clcts ) {
  std::vector<int> strips;

  for ( CSCCLCTDigiCollection::DigiRangeIterator it = clcts->begin(); it != clcts->end(); ++it ) {
  
    // Test if within same chamber
    // note that CLCT wiregroup is average over all 6 layers !!!
    const CSCDetId& clctId = (*it).first;
    if ( (clctId.chamber() == id.chamber()) &&
         (clctId.station() == id.station()) &&
         (clctId.ring()    == id.ring())    &&
         (clctId.endcap()  == id.endcap())  ) {
      for ( std::vector<CSCCLCTDigi>::const_iterator itr = (*it).second.first; itr != (*it).second.second; ++itr) {
        strips.push_back( int(itr->getKeyStrip()) );      
      }
    }  
  }    
  return strips;
}


/* Look at strips which matches CLCT(s)
 *
 */
bool CSCHitFromStripOnly::foundCLCTMatch( int stripId, std::vector<int> clctStrips ) {

  bool foundMatch = false;

  // Note: CSCT Key strip is 1/2 strip # and start at 0 whereas stripId start at 1 and is full strip Id
  for (unsigned i = 0; i < clctStrips.size(); ++i )
    if (abs(stripId - clctStrips[i]/2) < 6 ) return true;

  return foundMatch;
}

