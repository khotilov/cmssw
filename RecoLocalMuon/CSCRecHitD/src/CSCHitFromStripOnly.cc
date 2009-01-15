// This is  CSCHitFromStripOnly.cc

#include "RecoLocalMuon/CSCRecHitD/src/CSCHitFromStripOnly.h"
#include "RecoLocalMuon/CSCRecHitD/src/CSCStripData.h"
#include "RecoLocalMuon/CSCRecHitD/src/CSCStripHitData.h"
#include "RecoLocalMuon/CSCRecHitD/src/CSCStripHit.h"
#include "RecoLocalMuon/CSCRecHitD/src/CSCPedestalChoice.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCChamberSpecs.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <algorithm>
#include <string>
#include <vector>


CSCHitFromStripOnly::CSCHitFromStripOnly( const edm::ParameterSet& ps ) : recoConditions_(0), calcped_(0) {
  
  useCalib                   = ps.getUntrackedParameter<bool>("CSCUseCalibrations");
  bool useStaticPedestals    = ps.getParameter<bool>("CSCUseStaticPedestals");
  int noOfTimeBinsForDynamicPed = ps.getParameter<int>("CSCNoOfTimeBinsForDynamicPedestal");
  //theClusterSize             = ps.getUntrackedParameter<int>("CSCStripClusterSize");
  theThresholdForAPeak       = ps.getUntrackedParameter<double>("CSCStripPeakThreshold");
  theThresholdForCluster     = ps.getUntrackedParameter<double>("CSCStripClusterChargeCut");

  LogTrace("CSCRecHit")  << "[CSCHitFromStripOnly] CSCUseStaticPedestals = " << useStaticPedestals;
  if ( !useStaticPedestals ) LogTrace("CSCRecHit")  << "[CSCHitFromStripOnly] CSCNoOfTimeBinsForDynamicPedestal = " 
						    << noOfTimeBinsForDynamicPed;

  if ( useStaticPedestals ) {
    calcped_ = new CSCStaticPedestal();
  }
  else {
    if ( noOfTimeBinsForDynamicPed == 1 ) {
      calcped_ = new CSCDynamicPedestal1();
    }
    else {
      calcped_ = new CSCDynamicPedestal2(); // NORMAL DEFAULT!
    }
  }

}


CSCHitFromStripOnly::~CSCHitFromStripOnly() {
  delete calcped_;
}


/* runStrip
 *
 * Search for strip with ADC output exceeding theThresholdForAPeak.  For each of these strips,
 * build a cluster of strip of size theClusterSize (typically 5 strips).  Finally, make
 * a Strip Hit out of these clusters by finding the center-of-mass position of the hit
 */
std::vector<CSCStripHit> CSCHitFromStripOnly::runStrip( const CSCDetId& id, const CSCLayer* layer,
                                               const CSCStripDigiCollection::Range& rstripd ) 
{	

  std::vector<CSCStripHit> hitsInLayer;
  
  // cache layer info for ease of access
  id_        = id;
  layer_     = layer;
  nstrips_   = layer->chamber()->specs()->nStrips();

  tmax_cluster = 5;

  // Get gain correction weights for all strips in layer, and cache in gainWeight.
  // They're used in fillPulseHeights below.
  if ( useCalib ) {
    recoConditions_->stripWeights( id, gainWeight );
  }
  
  // Store pulseheights from SCA and find maxima (potential hits)
  fillPulseHeights( rstripd );
  findMaxima();      

  // Make a Strip Hit out of each strip local maximum
  for ( size_t imax = 0; imax != theMaxima.size(); ++imax ) {

    // Initialize parameters entering the CSCStripHit
    clusterSize = theClusterSize;
    theStrips.clear();
    strips_adc.clear();
    strips_adcRaw.clear();

    // This is where centroid position is determined
    // The strips_adc vector is also filled here
    // Remember, the array starts at 0, but the stripId starts at 1...
    float strippos = makeCluster( theMaxima[imax]+1 );    
    
    if ( strippos < 0 || tmax_cluster < 3 ) continue;

    //---- If two maxima are too close the error assigned will be width/sqrt(12) - see CSCXonStrip_MatchGatti.cc
    int maximum_to_left = 99; //---- If there is one maximum - the distance is set to 99 (strips)
    int maximum_to_right = 99;
    if(imax<theMaxima.size()-1){
      maximum_to_right = theMaxima.at(imax+1) - theMaxima.at(imax);
    }
    if(imax>0 ){
      maximum_to_left =  theMaxima.at(imax-1) - theMaxima.at(imax);
    }
    if(fabs(maximum_to_right) < fabs(maximum_to_left)){
      theClosestMaximum.push_back(maximum_to_right);
    }
    else{
      theClosestMaximum.push_back(maximum_to_left);
    }
    
    //---- Check if a neighbouring strip is a dead strip
    bool deadStrip = isNearDeadStrip(id, theMaxima.at(imax)); 
    
    CSCStripHit striphit( id, strippos, tmax_cluster, theStrips, strips_adc, strips_adcRaw,
			  theConsecutiveStrips.at(imax), theClosestMaximum.at(imax), deadStrip);
    hitsInLayer.push_back( striphit ); 
  }

  return hitsInLayer;
}


/* makeCluster
 *
 */
float CSCHitFromStripOnly::makeCluster( int centerStrip ) {
  
  float strippos = -1.;
  clusterSize = theClusterSize;
  std::vector<CSCStripHitData> stripDataV;
 
  // We only want to use strip position in terms of strip # for the strip hit. //@@ What other choice is there?
    
  // If the cluster size is such that you go beyond the edge of detector, shrink cluster appropriately
  for ( int i = 1; i < theClusterSize/2 + 1; ++i ) {
 
    if ( centerStrip - i < 1 || centerStrip + i > int(nstrips_) ) {

      // Shrink cluster size, but keep it an odd number of strips.
      clusterSize = 2*i - 1;  
    }
  }

  for ( int i = -clusterSize/2; i <= clusterSize/2; ++i ) {
    CSCStripHitData data = makeStripData(centerStrip, i);
    stripDataV.push_back( data );
    theStrips.push_back( centerStrip + i );
  }
  
  strippos = findHitOnStripPosition( stripDataV, centerStrip );
  
  return strippos;
}


/** makeStripData
 *
 */
CSCStripHitData CSCHitFromStripOnly::makeStripData(int centerStrip, int offset) {
  
  CSCStripHitData prelimData;
  int thisStrip = centerStrip+offset;

  int tmax      = thePulseHeightMap[centerStrip-1].tmax();
  tmax_cluster  = tmax;

  std::vector<float> adc(4);
  std::vector<float> adcRaw(4);

  // Fill adc & adcRaw
	
  int istart = tmax-1;
  int istop  = std::min( tmax+2, 7 ) ; // there are only time bins 0-7
  adc[3] = 0.1; // in case it isn't filled
	
  if ( tmax > 2 && tmax < 7 ) { // for time bins 3-6
    int ibin = thisStrip-1;
		
    std::copy( thePulseHeightMap[ibin].ph().begin()+istart, 
	 thePulseHeightMap[ibin].ph().begin()+istop+1, adc.begin() );
			
    std::copy( thePulseHeightMap[ibin].phRaw().begin()+istart, 
	 thePulseHeightMap[ibin].phRaw().begin()+istop+1, adcRaw.begin() );
  } 
  else {
    adc[0] = 0.1;
    adc[1] = 0.1;
    adc[2] = 0.1;
    adc[3] = 0.1;
    adcRaw = adc;
    LogTrace("CSCRecHit")  << "[CSCHitFromStripOnly::makeStripData] Tmax out of range: contact CSC expert!";
  }
  
  if ( offset == 0 ) {
    prelimData = CSCStripHitData(thisStrip, tmax_cluster, adcRaw, adc);
  } else {
    int sign = offset>0 ? 1 : -1;
    // If there's another maximum that would like to use part of this cluster, 
    // it gets shared in proportion to the height of the maxima
    for ( int i = 1; i <= clusterSize/2; ++i ) {

      // Find the direction of the offset
      int testStrip = thisStrip + sign*i;
      std::vector<int>::iterator otherMax = find(theMaxima.begin(), theMaxima.end(), testStrip-1);

      // No other maxima found, so just store
      if ( otherMax == theMaxima.end() ) {
        prelimData = CSCStripHitData(thisStrip, tmax_cluster, adcRaw, adc);      } 
      else {
	
     	// Another maximum found - share       
	std::vector<float> adc1(4);
	std::vector<float> adcRaw1(4);
	std::vector<float> adc2(4);
	std::vector<float> adcRaw2(4);
	// In case we only copy (below) into 3 of the 4 bins i.e. when istart=5, istop=7
	adc1[3]    = 0.1; 
        adc2[3]    = 0.1; 
        adcRaw1[3] = 0.1; 
        adcRaw2[3] = 0.1; 

        // Fill adcN with content of time bins tmax-1 to tmax+2 (if it exists!)
        if ( tmax > 2 && tmax < 7 ) { // for time bin tmax from 3-6
	  int ibin = testStrip-1;
	  int jbin = centerStrip-1;
	  std::copy(thePulseHeightMap[ibin].ph().begin()+istart, 
	       thePulseHeightMap[ibin].ph().begin()+istop+1, adc1.begin());
				
	  std::copy(thePulseHeightMap[ibin].phRaw().begin()+istart, 
	       thePulseHeightMap[ibin].phRaw().begin()+istop+1, adcRaw1.begin());
										
	  std::copy(thePulseHeightMap[jbin].ph().begin()+istart, 
	       thePulseHeightMap[jbin].ph().begin()+istop+1, adc2.begin());  
						
	  std::copy(thePulseHeightMap[jbin].phRaw().begin()+istart, 
	       thePulseHeightMap[jbin].phRaw().begin()+istop+1, adcRaw2.begin());
	} 
	else {
	  adc1.assign(4, 0.1);
	  adcRaw1 = adc1;
	  adc2.assign(4, 0.1);
	  adcRaw2 = adc2;
	}
				
        // Scale shared strip B ('adc') by ratio of peak of ADC counts from central strip A ('adc2')
        // to sum of A and neighbouring maxima C ('adc1')

	for (size_t k = 0; k < 4; ++k){
	  if(adc1[k]>0    && adc2[k]>0)    adc[k]    = adc[k] * adc2[k] / ( adc1[k]+adc2[k] );
	  if(adcRaw1[k]>0 && adcRaw2[k]>0) adcRaw[k] = adcRaw[k] * adcRaw2[k] / ( adcRaw1[k]+adcRaw2[k] );     
	}
        prelimData = CSCStripHitData(thisStrip, tmax_cluster, adcRaw, adc);
      }
    }
  }
  return prelimData;
}
 

/* fillPulseHeights
 *
 */
void CSCHitFromStripOnly::fillPulseHeights( const CSCStripDigiCollection::Range& rstripd ) {
  
  // Loop over strip digis in one CSCLayer and fill PulseHeightMap with pedestal-subtracted
  // SCA pulse heights.
  
  thePulseHeightMap.clear();
  thePulseHeightMap.resize(100); //@@ WHY NOT JUST 80?

  // for storing sca pulseheights once they may no longer be integer (e.g. after ped subtraction)
  std::vector<float> sca;
  sca.reserve(8);
	
  for ( CSCStripDigiCollection::const_iterator it = rstripd.first; it != rstripd.second; ++it ) {
    int  thisChannel        = (*it).getStrip(); 
    std::vector<int> scaRaw = (*it).getADCCounts();
    sca.clear();
    // Fill sca from scaRaw, implicitly converting to float
    std::copy( scaRaw.begin(), scaRaw.end(), std::back_inserter( sca ));

    //@@ Find bin with largest pulseheight (_before_ ped subtraction - shouldn't matter, right?)
    int tmax =  std::max_element( sca.begin(), sca.end() ) - sca.begin(); // counts from 0

    // get pedestal - calculated as appropriate - for this sca pulse
    float ped = calcped_->pedestal(sca, recoConditions_, id_, thisChannel );

    // subtract the pedestal (from BOTH sets of sca pulseheights)
    std::for_each( sca.begin(), sca.end(), CSCSubtractPedestal( ped ) );
    std::for_each( scaRaw.begin(), scaRaw.end(), CSCSubtractPedestal( ped ) );

    //@@ Max in first 3 or last time bins is unacceptable, if so set to zero (why?)
    float phmax = 0.;
    if ( tmax > 2 && tmax < 7 ) {
      phmax = sca[tmax];
    }
		
    // Fill the map, possibly apply gains from cond data, and unfold ME1A channels
    // (To apply gains use CSCStripData::op*= which scales only the non-raw sca ph's;
    // but note that both sca & scaRaw are pedestal-subtracted.)

    if ( id_.ring() != 4 ) { // non-ME1a
      thePulseHeightMap[thisChannel-1] = CSCStripData( thisChannel, phmax, tmax, scaRaw, sca );
      if ( useCalib ) thePulseHeightMap[thisChannel-1] *= gainWeight[thisChannel-1];
    } 
    else { // ME1a, so unfold its 16 channels to its 48 strips
      for ( int j = 0; j < 3; ++j ) {
	thePulseHeightMap[thisChannel-1+16*j] = CSCStripData( thisChannel+16*j, phmax, tmax, scaRaw, sca );
	if ( useCalib ) thePulseHeightMap[thisChannel-1+16*j] *= gainWeight[thisChannel-1];
      }
    }

  }
}


/* findMaxima
 *
 * fills vector 'theMaxima' with the local maxima in the pulseheight distribution
 * of the strips. The threshold defining a maximum is a configurable parameter.
 * A typical value is 30.
 */
void CSCHitFromStripOnly::findMaxima() {
  
  theMaxima.clear();
  theConsecutiveStrips.clear();
  theClosestMaximum.clear();
  for ( size_t i=0; i!=thePulseHeightMap.size(); ++i ) {

    float heightPeak = thePulseHeightMap[i].phmax();

    // sum 3 strips so that hits between strips are not suppressed
    float heightCluster;

    bool maximumFound = false;
    // Left edge of chamber
    if ( i == 0 ) {
      heightCluster = thePulseHeightMap[i].phmax()+thePulseHeightMap[i+1].phmax();
      // Have found a strip Hit if...
      if (( heightPeak                   > theThresholdForAPeak           ) &&  
          ( heightCluster                > theThresholdForCluster         ) && 
	        ( thePulseHeightMap[i].phmax() >= thePulseHeightMap[i+1].phmax()) &&
          ( thePulseHeightMap[i].tmax()  > 2                              ) &&
          ( thePulseHeightMap[i].tmax()  < 7                              )) {
        theMaxima.push_back(i);
	      maximumFound = true;
      }
    // Right edge of chamber
    } else if ( i == thePulseHeightMap.size()-1) {  
      heightCluster = thePulseHeightMap[i-1].phmax()+thePulseHeightMap[i].phmax();
      // Have found a strip Hit if...
      if (( heightPeak                   > theThresholdForAPeak          ) &&  
          ( heightCluster                > theThresholdForCluster        ) && 
          ( thePulseHeightMap[i].phmax() > thePulseHeightMap[i-1].phmax()) &&
          ( thePulseHeightMap[i].tmax()  > 2                             ) &&
          ( thePulseHeightMap[i].tmax()  < 7                             )) {
        theMaxima.push_back(i);
	      maximumFound = true;
      }
    // Any other strips
    } else {
      heightCluster = thePulseHeightMap[i-1].phmax()+thePulseHeightMap[i].phmax()+thePulseHeightMap[i+1].phmax();
      // Have found a strip Hit if...
      if (( heightPeak                   > theThresholdForAPeak           ) &&  
          ( heightCluster                > theThresholdForCluster         ) && 
          ( thePulseHeightMap[i].phmax() > thePulseHeightMap[i-1].phmax() ) &&
	        ( thePulseHeightMap[i].phmax() >= thePulseHeightMap[i+1].phmax()) &&
          ( thePulseHeightMap[i].tmax()  > 2                              ) &&
          ( thePulseHeightMap[i].tmax()  < 7                             )) {
        theMaxima.push_back(i);
      	maximumFound = true;
      }
    }
    //---- Consecutive strips with charge (real cluster); if too wide - measurement is not accurate 
    if(maximumFound){
      int numberOfConsecutiveStrips = 1;
      float testThreshold = 10.;//---- ADC counts; 
                                //---- this is not XTalk corrected so it is correct in first approximation only
      int j = 0;
      for(int l = 0; l<8; ++l){
        if(j<0) edm::LogError("CSCRecHit") << "This should never occur!!! Contact CSC expert!";
        ++j;
        bool signalPresent = false;
        for(int k = 0; k<2; ++k){
          j*= -1;//---- check from left and right
          int anotherConsecutiveStrip = i+j;
          if(anotherConsecutiveStrip>=0 && anotherConsecutiveStrip<int( thePulseHeightMap.size() )){
            if(thePulseHeightMap[anotherConsecutiveStrip].phmax()>testThreshold){
              ++numberOfConsecutiveStrips;
              signalPresent = true;
            }
          }
        }
        if(!signalPresent){
          break;
        }
      }
      theConsecutiveStrips.push_back(numberOfConsecutiveStrips);
    }
  }
}



/* findHitOnStripPosition
 *
 */
float CSCHitFromStripOnly::findHitOnStripPosition( const std::vector<CSCStripHitData>& data, const int& centerStrip ) {
  
  float strippos = -1.;
  
  if ( data.size() < 1 ) return strippos;
  
  // biggestStrip is strip with largest pulse height
  // Use pointer subtraction

  int biggestStrip = max_element(data.begin(), data.end()) - data.begin();
  strippos = data[biggestStrip].strip() * 1.;
  
  // If more than one strip:  use centroid to find center of cluster
  // but only use time bin == tmax (otherwise, bias centroid).
  float sum  = 0.;
  float sum_w= 0.;

  float stripLeft = 0.;
  float stripRight = 0.;

  std::vector<float> w(4);
  std::vector<float> wRaw(4);
  
  for ( size_t i = 0; i != data.size(); ++i ) {
    w = data[i].ph();
    wRaw = data[i].phRaw();

    // Require ADC to be > 0.
    for ( size_t j = 0; j != w.size(); ++j ) {
      if ( w[j] < 0. ) w[j] = 0.001;
    }

    if (i == data.size()/2 -1) stripLeft = w[1];
    if (i == data.size()/2 +1) stripRight = w[1];


    // Fill the data members 
    std::copy( w.begin(), w.end(), std::back_inserter(strips_adc));
    std::copy( wRaw.begin(), wRaw.end(), std::back_inserter(strips_adcRaw));

    if ( data[i].strip() < 1 ){
      LogTrace("CSCRecHit") << "problem in indexing of strip, strip id is: " << data[i].strip();
    } 
    sum_w += w[1];
    sum   += w[1] * data[i].strip();
  }

  if ( sum_w > 0.) strippos = sum / sum_w;    

  return strippos;
}

bool CSCHitFromStripOnly::isNearDeadStrip(const CSCDetId& id, int centralStrip){

  //@@ Tim says: not sure I understand this properly... but just moved code to CSCRecoConditions
  // where it can handle the conversion from strip to channel etc.
  return recoConditions_->nearBadStrip( id, centralStrip );

} 

// Define space for static
const int CSCHitFromStripOnly::theClusterSize;
