// This is CSCMake2DRecHit

#include <RecoLocalMuon/CSCRecHitB/src/CSCMake2DRecHit.h>
#include <RecoLocalMuon/CSCRecHitB/src/CSCStripCrosstalk.h>
#include <RecoLocalMuon/CSCRecHitB/src/CSCStripNoiseMatrix.h>
#include <RecoLocalMuon/CSCRecHitB/src/CSCFitXonStripWithGatti.h>
#include <RecoLocalMuon/CSCRecHitB/src/probab.h>
#include <RecoLocalMuon/CSCRecHitB/interface/CSCStripHit.h>
#include <RecoLocalMuon/CSCRecHitB/interface/CSCWireHit.h>

#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>

#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCChamberSpecs.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>

#include <CondFormats/CSCObjects/interface/CSCGains.h>
#include <CondFormats/DataRecord/interface/CSCGainsRcd.h>
#include <CondFormats/CSCObjects/interface/CSCcrosstalk.h>
#include <CondFormats/DataRecord/interface/CSCcrosstalkRcd.h>
#include <CondFormats/CSCObjects/interface/CSCNoiseMatrix.h>
#include <CondFormats/DataRecord/interface/CSCNoiseMatrixRcd.h>

#include <Geometry/Vector/interface/LocalPoint.h>

#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/Exception.h>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>


/* Constructor
 *
 */
CSCMake2DRecHit::CSCMake2DRecHit(const edm::ParameterSet& ps){
    
  debug                      = ps.getUntrackedParameter<bool>("CSCDebug");
  useCalib                   = ps.getUntrackedParameter<bool>("CSCUseCalibrations");
  stripWireDeltaTime         = ps.getUntrackedParameter<int>("CSCstripWireDeltaTime");
  useGatti                   = ps.getUntrackedParameter<bool>("CSCUseGattiFit");
  maxGattiChi2               = ps.getUntrackedParameter<double>("CSCMaxGattiChi2");

  stripCrosstalk_            = new CSCStripCrosstalk( ps );
  stripNoiseMatrix_          = new CSCStripNoiseMatrix( ps );
  xFitWithGatti_             = new CSCFitXonStripWithGatti( ps );
}   


/* Destructor
 *
 */
CSCMake2DRecHit::~CSCMake2DRecHit() {
  delete stripCrosstalk_;
  delete stripNoiseMatrix_;
  delete xFitWithGatti_;
}


/* hitFromStripAndWire
 *
 */
CSCRecHit2D CSCMake2DRecHit::hitFromStripAndWire(const CSCDetId& id, const CSCLayer* layer,
                                                 const CSCWireHit& wHit, const CSCStripHit& sHit){

  if (debug) std::cout <<"[CSCMake2DRecHit::hitFromStripAndWire] creating 2-D hit" << std::endl;  

  // Cache layer info for ease of access
  layer_        = layer;
  layergeom_    = layer_->geometry();
  specs_        = layer->chamber()->specs();
  id_           = id;

  double sigma, chisq, prob;
  sigma = chisq = 0.00;
  float dx2, dy2, dxy, tpeak;
  dx2 = dy2 = dxy = tpeak = 0.;
  prob  = 1.00;
    
  CSCRecHit2D::ChannelContainer channels;
  CSCRecHit2D::ADCContainer adcMap;
  CSCRecHit2D::ChannelContainer wgroups;

  // In case hit falls outside fiducial volume of chamber:
  bool keepHit = true;
  LocalPoint lpFailed(-999., -999.);
  LocalError localerrFailed(dx2, dxy, dy2);  
  CSCRecHit2D failedHit( id, lpFailed, localerrFailed, channels, adcMap, wgroups, tpeak, chisq, prob );

  float slopeRight[100];
  float slopeLeft[100];
  float interRight[100];
  float interLeft[100];
  std::vector<float> nMatrix;

  // Fill x-talk and noise matrix at once:
  if ( useCalib ) {
    stripCrosstalk_->setCrossTalk( xtalk_ );
    stripCrosstalk_->getCrossTalk( id, slopeLeft, interLeft, slopeRight, interRight );
    stripNoiseMatrix_->setNoiseMatrix( gains_, noise_ );
    stripNoiseMatrix_->getNoiseMatrix( id, nMatrix ); 
  } else {
    // FIXME:  NOT SUPPOSED TO HARDWIRE ANYTHING !!!
    nMatrix.push_back(10.06);  // 3,3
    nMatrix.push_back(3.85);   // 3,4
    nMatrix.push_back(2.96);   // 3,5
    nMatrix.push_back(8.53);   // 4,4
    nMatrix.push_back(3.12);   // 4,5
    nMatrix.push_back(2.05);   // 4,6
    nMatrix.push_back(8.72);   // 5,5
    nMatrix.push_back(3.29);   // 5,6
    nMatrix.push_back(1.83);   // 5,7
    nMatrix.push_back(7.65);   // 6,6
    nMatrix.push_back(2.98);   // 6,7
    nMatrix.push_back(0.);     // 6,8  <- does not exist
    nMatrix.push_back(7.95);   // 7,7
    nMatrix.push_back(0.);     // 7,8  <- does not exist
    nMatrix.push_back(0.);     // 7,9  <- does not exist
  }

  
  // Find wire hit position and wire properties
  float wireg = wHit.wHitPos();                            // Position in terms of wire group #
  wgroups = wHit.wgroups();
  int wireg1 = int( wireg );
  int wireg2 = wireg1;
  if ( (wireg - wireg1) != 0. ) {
    wireg2++;
  }  
  float wgoffset   = wireg - wireg1;
  float the_wire1  = layergeom_->middleWireOfGroup(wireg1);
  float the_wire2  = layergeom_->middleWireOfGroup(wireg2);
  int n_wires      = layergeom_->numberOfWiresPerGroup(wireg1); 
  double wSpace    = layergeom_->wirePitch();
  float wAngle     = layergeom_->wireAngle();
  float dy         = wSpace * (n_wires) / sqrt(12.);       // Error on y'  (on y if wAngle = 0)


  // Find strip position and properties
  CSCRecHit2D::ChannelContainer strips = sHit.strips();
  tpeak = sHit.tmax() * 50.;
  int stripClusterSize = strips.size();
  float strip_pos = sHit.sHitPos();
  int ch = int(strip_pos);
  float strip_offset = float(strip_pos) - ch;
  if ( strip_offset >= 0.5 ) ch++;
  int centerStrip = ch;
  double sAngle   = layergeom_->stripAngle(ch);

  std::vector<float> adcs = sHit.s_adc();
  std::vector<float> adc2;
  for ( int t = 0; t < stripClusterSize*4; t++ ) {
    if ( t%4 == 0 ) adc2.clear();
    adc2.push_back(adcs[t]);
    if ( (t+1)%4 == 0) adcMap.put( strips[t/4], adc2.begin(), adc2.end() ); 
  }


  // CFEB readout problem:
  // If the strip positon and the channel # are the same, it means I'm either at the edge of the
  // chamber, or that one of the strip next to the central strip wasn't readout.  Flag these hits as
  // we don't want them to go through the Gatti fit procedure.
  bool hasCFEBProblem = false;
  if (strip_offset == 0. && ch%16 < 2) { 
    hasCFEBProblem = true;
    if (debug &&  ch > 1 && ch < specs_->nStrips() ) 
      std::cout << "Have found potential CFEB problem for strip: " << ch << std::endl; 
  }
  

  // CFEB trigger problem:
  // Also use the following if the hit is flagged as deficient


  // If at the edge, then used 1 strip cluster only :
  if ( ch <= 1 || ch >= specs_->nStrips() || hasCFEBProblem) {

    LocalPoint lp1 = layergeom_->stripWireIntersection( centerStrip, the_wire1);
    LocalPoint lp2 = layergeom_->stripWireIntersection( centerStrip, the_wire2);
    float x1 = lp1.x();
    float x2 = lp2.x();
    float y1 = lp1.y();
    float y2 = lp2.y();

    float x = ( x1 + wgoffset * (x2 - x1) );
    float y = ( y1 + wgoffset * (y2 - y1) );

    LocalPoint lp3(x, y);

    // Ensure that y position is within active area (ME_11 chambers):
    LocalPoint lp0;
    keepHit = keepHitInFiducial( lp3, lp0 );

    if ( !keepHit ) { 
      if (debug) std::cout <<"[CSCMake2DRecHit::hitFromStripAndWire] failedHit" << std::endl;
      return failedHit;
    }
    sigma =  layergeom_->stripPitch(lp0)/sqrt(12);  
    double dx = sigma/sin(sAngle);
    
    float sinangdif = sin(sAngle - wAngle);
    float sin2angdif = sinangdif * sinangdif;   
    
    float wsins = dy * sin(sAngle);
    float wcoss = dy * cos(sAngle);
    float ssinw = dx * sin(wAngle);
    float scosw = dx * cos(wAngle);
    
    dx2 = (scosw*scosw + wcoss*wcoss)/sin2angdif;
    dy2 = (ssinw*ssinw + wsins*wsins)/sin2angdif;
    dxy = (scosw*ssinw + wcoss*wsins)/sin2angdif;
    
    LocalError localerr(dx2, dxy, dy2);

    CSCRecHit2D rechit( id, lp0, localerr, strips, adcMap, wgroups, tpeak, chisq, prob );

    if (debug) std::cout << "Found rechit in chamber ME" << id.station() << "/" << id.ring() << "-" << id.chamber() << " and layer " << id.layer() << " with local position:  x = " << x << "  y = " << y << std::endl; 

    nMatrix.clear();

    return rechit;  
  } 


  // If not at the edge, used cluster of size ClusterSize:

  int ch0 = int(strip_pos);
  LocalPoint lp11  = layergeom_->stripWireIntersection( ch0,   the_wire1);
  LocalPoint lp111 = layergeom_->stripWireIntersection( ch0,   the_wire2);
  LocalPoint lp22  = layergeom_->stripWireIntersection( ch0+1, the_wire1);
  LocalPoint lp222 = layergeom_->stripWireIntersection( ch0+1, the_wire2);

  float x11  = lp11.x();
  float x111 = lp111.x();
  float x1   = (x11 + wgoffset * (x111 - x11) );
  float y11  = lp11.y();
  float y111 = lp111.y();
  float y1   = (y11 + wgoffset * (y111 - y11) );
  LocalPoint lp1(x1, y1);

  float x22  = lp22.x();
  float x222 = lp222.x();
  float x2   = (x22 + wgoffset * (x222 - x22) );
  float y22  = lp22.y();
  float y222 = lp222.y();
  float y2   = (y22 + wgoffset * (y222 - y22) );
  LocalPoint lp2(x2, y2);
  

  // Ensure that y position is within active area (ME_11 chambers):
  LocalPoint lp3;
  keepHit = keepHitInFiducial( lp1, lp3 );
  if ( !keepHit ) { 
    if (debug) std::cout <<"[CSCMake2DRecHit::hitFromStripAndWire] failedHit" << std::endl;
    return failedHit;
  }
  LocalPoint lp4;
  keepHit = keepHitInFiducial( lp2, lp4 );
  if ( !keepHit ) { 
    if (debug) std::cout <<"[CSCMake2DRecHit::hitFromStripAndWire] failedHit" << std::endl;
    return failedHit;
  }

  // Use center of gravity to determine local x and y position:
  float x3 = lp3.x();
  float x4 = lp4.x();
  float x  = (1. - strip_offset) * x3 + strip_offset * x4;
  float y3 = lp3.y();
  float y4 = lp4.y();
  float y  = (1. - strip_offset) * y3 + strip_offset * y4;

  // Build local point
  float stripWidth = fabs(x2 - x1);  
  sigma =  stripWidth/sqrt(12);              

  // Here try to improve the strip position by applying Gatti fitter

  if ( useGatti ) {   

    LocalPoint lp55  = layergeom_->stripWireIntersection( centerStrip, the_wire1);
    LocalPoint lp555 = layergeom_->stripWireIntersection( centerStrip, the_wire2);
    x1 = lp55.x();
    y1 = lp55.y();
    x2 = lp555.x();
    y2 = lp555.y();

    x1   = (x1 + wgoffset * (x2 - x1) );
    y1   = (y1 + wgoffset * (y2 - y1) );
  
    LocalPoint lp5(x1, y1);
                                   
    LocalPoint lp6;
    keepHit = keepHitInFiducial( lp5, lp6 );
    if ( !keepHit ) { 
      if (debug) std::cout <<"[CSCMake2DRecHit::hitFromStripAndWire] failedHit" << std::endl;
      nMatrix.clear();
      return failedHit;
    }	

    // Local position at center of strip
    float x_to_gatti = lp6.x();   
    float x_fit;
    double sigma_fit, chisq_fit;

    std::vector<float> xtalks, nmatrix;
    for ( int j = centerStrip -2; j < centerStrip + 1; j++) {  // Index of strip is starts at 1 whereas array 
      xtalks.push_back(slopeLeft[j]);                          // starts at 0 --> make use of this here
      xtalks.push_back(interLeft[j]);
      xtalks.push_back(slopeRight[j]);
      xtalks.push_back(interRight[j]);
      for ( int k = 0; k < 15; k++ ) nmatrix.push_back(nMatrix[j*15 + k]);
    }   

    xFitWithGatti_->findXOnStrip( layer_, sHit, x_to_gatti, stripWidth, xtalks, nmatrix, x_fit, tpeak, sigma_fit, chisq_fit );

    if ( chisq_fit < maxGattiChi2 ) {
      x     = x_fit;
      sigma = sigma_fit;
      chisq = chisq_fit;
      int ndof = 5;
      prob  = probab(ndof, chisq_fit);
    }
  }
  LocalPoint lp0(x, y);


  // Now compute the errors properly on local x and y

  double dx = sigma / sin(sAngle);         
  
  float sinangdif = sin(sAngle - wAngle);
  float sin2angdif = sinangdif * sinangdif;
  
  float wsins = dy * sin(sAngle); 
  float wcoss = dy * cos(sAngle);
  float ssinw = dx * sin(wAngle);
  float scosw = dx * cos(wAngle);
  
  dx2 = (scosw*scosw + wcoss*wcoss)/sin2angdif;
  dy2 = (ssinw*ssinw + wsins*wsins)/sin2angdif;
  dxy = (scosw*ssinw + wcoss*wsins)/sin2angdif;
  
  LocalError localerr(dx2, dxy, dy2);

  // store rechit    
  CSCRecHit2D rechit( id, lp0, localerr, strips, adcMap, wgroups, tpeak, chisq, prob );

  if (debug) std::cout << "Found rechit in chamber ME" << id.station() << "/" << id.ring() << "-" << id.chamber() << " and layer " << id.layer() << " with local position:  x = " << x << "  y = " << y << std::endl; 

  nMatrix.clear();

  return rechit;
}



/* hitFromWireOnly
 *
 */
CSCRecHit2D CSCMake2DRecHit::hitFromWireOnly(const CSCDetId& id, const CSCLayer* layer, const CSCWireHit& wHit){

  // Cache layer info for ease of access
  layer_        = layer;
  layergeom_    = layer_->geometry();
  specs_        = layer->chamber()->specs();
  id_           = id;
  float wAngle  = layergeom_->wireAngle();
  int this_layer= id.layer();

  double sigma, chisq, prob;
  sigma = chisq = 0.00;
  prob  = 1.00;
  float dx2, dy2, dxy;
  dx2 = dy2 = dxy = 0.;

  CSCRecHit2D::ChannelContainer strips;
  CSCRecHit2D::ChannelContainer wiregroups;
  CSCRecHit2D::ADCContainer adcMap;

  // Find wire hit position and wire properties
  float tpeak = wHit.tmax() * 50.;
  float wireg = wHit.wHitPos();                            // Position in terms of wire group #
  CSCRecHit2D::ChannelContainer wgroups = wHit.wgroups();
  int wireg1 = int( wireg );
  int wireg2 = wireg1;
  if ( (wireg - wireg1) != 0. ) {
    wireg2++;
  }
  float wgoffset   = wireg - wireg1;
    
  LocalPoint lp1 = layergeom_->localCenterOfWireGroup( wireg1 ); 
  LocalPoint lp2 = layergeom_->localCenterOfWireGroup( wireg2 );

  float x = lp1.x() + wgoffset * ( lp2.x() - lp1.x() ); 
  float y = lp1.y() + wgoffset * ( lp2.y() - lp1.y() ); 

  LocalPoint lp0(x,y);

  // Now compute errors properly:

  float width1 = layergeom_->lengthOfWireGroup( wireg1 );
  float width2 = layergeom_->lengthOfWireGroup( wireg1 );
  float width  = width1 + wgoffset * (width2 - width1 );

  sigma = width * cos(wAngle) / sqrt(12.);

  // What are the proper error on these ?????

  dx2 = sigma * sigma;
  dy2 = ( lp1.y()-lp2.y() )*( lp1.y()-lp2.y() ) /12.; 
  dxy = 0.;

  LocalError localerr(dx2, dxy, dy2);

  CSCRecHit2D rechit( id, lp0, localerr, strips, adcMap, wgroups, tpeak, chisq, prob );

  if (debug) std::cout << "Found wire hit in layer " << this_layer << " with local position:  x = " << x << "  y = " << y << std::endl; 

  return rechit;
}



/* hitFromStripOnly
 *
 * Don't bother applying Gatti fit here since we don't know the "y", i.e. the width of the strip.
 * Centroid will do.  Use y = 0 to find the "x", then assume the strip width at top/sqrt(12) is error (maximum deviation)
 */
CSCRecHit2D CSCMake2DRecHit::hitFromStripOnly(const CSCDetId& id, const CSCLayer* layer, const CSCStripHit& sHit){

  // Cache layer info for ease of access
  layer_        = layer;
  layergeom_    = layer_->geometry();
  specs_        = layer->chamber()->specs();
  id_           = id;

  int this_layer= id.layer();

  double sigma, chisq, prob;
  chisq = 0.00;
  prob  = 1.00;

  CSCRecHit2D::ChannelContainer wgroups;
  CSCRecHit2D::ADCContainer adcMap;

  // Strip position and properties
  CSCRecHit2D::ChannelContainer strips = sHit.strips();
  int stripClusterSize = strips.size();
  float tpeak = sHit.tmax()*50.;
  float strip_pos = sHit.sHitPos();
  int ch = int(strip_pos);

  float chamberHalfLength = layergeom_->length()/2.;
  float stripWidth = layergeom_->stripPitch(); 

  std::vector<float> adcs = sHit.s_adc();
  std::vector<float> adc2;
  for ( int t = 0; t < stripClusterSize*4; t++ ) {
    if ( t%4 == 0 ) adc2.clear();
    adc2.push_back(adcs[t]);
    if ( (t+1)%4 == 0) adcMap.put( strips[t/4], adc2.begin(), adc2.end() );
  }

  float x, xAtTop, xAtBottom;

  // Use central strip for this case (not worth using centroid)
  x         = layergeom_->xOfStrip( ch, 0.);
  xAtTop    = layergeom_->xOfStrip( ch, chamberHalfLength);
  xAtBottom = layergeom_->xOfStrip( ch, -chamberHalfLength);

  float y=0.;

  LocalPoint lp0(x, y);

  sigma = (0.5*(xAtTop+xAtBottom)+stripWidth)/sqrt(12.);     // Be on the conservative side...

  // Now compute the errors, but again, what are they ???

  double dx2 = sigma * sigma;
  double dy2 = chamberHalfLength*chamberHalfLength/12.;      // i.e. sigma_y = 0.5 * chamberLength/sqrt(12);
  double dxy = 0.;

  LocalError localerr(dx2, dxy, dy2);

  CSCRecHit2D rechit( id, lp0, localerr, strips, adcMap, wgroups, tpeak, chisq, prob );

  if (debug) std::cout << "Found strip hit in layer " << this_layer << " with local position:  x = " << x << "  y = " << y << std::endl; 

  return rechit;
}



/* keepHitInFiducial
 *
 * Only useful for ME11 chambers.
 */
bool CSCMake2DRecHit::keepHitInFiducial( LocalPoint& lp1, LocalPoint& lp0 ) {

  bool isInFiducial = true;

  // Allow extra margin for future tuning etc. For now 0.1 cm.
  float marginAtEdge = 0.1; 

  // Initialize parameters needed for correction:
  float x1 = lp1.x();
  float y1 = lp1.y();
  int centralStrip = layergeom_->nearestStrip( lp0 );
  float wAngle     = layergeom_->wireAngle();      
  float apothem    = layergeom_->length()/2.;

  // Test what type of chambers we are dealing with

  if (id_.station() == 1 && id_.ring() == 4 ) {
    // Within ME1a chambers, allow hit to be reconstructed below bottom, 
    // but not above top of chamber
    if (y1 > apothem + marginAtEdge) {
      lp0 = lp1;
      isInFiducial = false;
      return isInFiducial;
    } 
  } else if (id_.station() == 1 && id_.ring() == 1 ) {
    // Within ME1b chambers, allow hit to be reconstructed above top, 
    // but not below bottom of chamber
    if (y1 < - (apothem + marginAtEdge)) {
      lp0 = lp1;
      isInFiducial = false;
      return isInFiducial;
    }

  } else {
   // Within all other chambers, wires are parallel to top/bottom of chamber 
   // so hit should always be within fiducial volume
   lp0 = lp1;
   return isInFiducial;
  }

  // Now deal with ME11 chambers:

  // Test if beyond edge of chamber
  if ( fabs(y1) > apothem )  {                                

    // beyond bottom of chamber
    if ( y1 < 0. ) {       

      float yNearWire     = apothem + marginAtEdge/cos( wAngle );           
      float xNearWire     = x1 - yNearWire * sin( wAngle );
      LocalPoint lpGetWire( xNearWire, yNearWire );
      float closestWire   = layergeom_->nearestWire( lpGetWire );
      lp0                 = layergeom_->stripWireIntersection( centralStrip, closestWire );

   // beyond top of chamber      
    } else {            

      float yNearWire     = apothem - marginAtEdge/cos( wAngle );           
      float xNearWire     = x1 + yNearWire * sin( wAngle );
      LocalPoint lpGetwire( xNearWire, yNearWire );
      float closestWire   = layergeom_->nearestWire( lpGetwire );
      lp0                 = layergeom_->stripWireIntersection( centralStrip, closestWire );
    }

  // Within chamber, so keep local point as is
  } else {              
    lp0 = lp1;
  } 

  return isInFiducial;
}



/* isHitInFiducial
 *
 * Only useful for ME11 chambers.
 */
bool CSCMake2DRecHit::isHitInFiducial( const CSCLayer* layer, const CSCRecHit2D& rh ) {

  bool isInFiducial = true;

  // Allow extra margin for future tuning etc.
  float marginAtEdge = 0.1; 

  const CSCLayerGeometry* layergeom = layer->geometry();
  
  float y = rh.localPosition().y();
  float apothem = layergeom->length()/2.;

  if ( fabs(y) > (apothem+marginAtEdge) ) isInFiducial=false;

  return isInFiducial;
}
 
 

