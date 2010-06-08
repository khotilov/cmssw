// This is CSCMake2DRecHit
 
#include <RecoLocalMuon/CSCRecHitD/src/CSCMake2DRecHit.h>
#include <RecoLocalMuon/CSCRecHitD/src/CSCXonStrip_MatchGatti.h>
#include <RecoLocalMuon/CSCRecHitD/src/CSCStripHit.h>
#include <RecoLocalMuon/CSCRecHitD/src/CSCWireHit.h>
#include <RecoLocalMuon/CSCRecHitD/src/CSCRecoConditions.h>
 
#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>

#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCChamberSpecs.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>

#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/Exception.h>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>


CSCMake2DRecHit::CSCMake2DRecHit(const edm::ParameterSet& ps){
    
  useCalib            = ps.getParameter<bool>("CSCUseCalibrations");
  stripWireDeltaTime  = ps.getParameter<int>("CSCstripWireDeltaTime"); //@@ Non-standard  CSC*s*trip...
  useTimingCorrections= ps.getParameter<bool>("CSCUseTimingCorrections");

  xMatchGatti_        = new CSCXonStrip_MatchGatti( ps );

}   


CSCMake2DRecHit::~CSCMake2DRecHit() {
  delete xMatchGatti_;
}


CSCRecHit2D CSCMake2DRecHit::hitFromStripAndWire(const CSCDetId& id, const CSCLayer* layer,
                                                 const CSCWireHit& wHit, const CSCStripHit& sHit)
{
  // Cache layer info for ease of access
  layer_        = layer;
  layergeom_    = layer_->geometry();
  specs_        = layer->chamber()->specs();
  id_           = id;
  
  const float sqrt_12 = 3.4641;
  
  float tpeak = -99.;
  
  CSCRecHit2D::ADCContainer adcMap;
  CSCRecHit2D::ChannelContainer wgroups;
  
  // Find wire hit position and wire properties
  wgroups = wHit.wgroups();

  int wg_left = wgroups[0];;
  int wg_right = wgroups[wgroups.size()-1];
  
  int Nwires1 = layergeom_->numberOfWiresPerGroup( wg_left );
  int Nwires2 = layergeom_->numberOfWiresPerGroup( wg_right );

  float Mwire1 = layergeom_->middleWireOfGroup( wg_left );
  float Mwire2 = layergeom_->middleWireOfGroup( wg_right );
  
  int centerWire_left = (int) (Mwire1 - Nwires1 / 2. + 0.5);
  int centerWire_right = (int) (Mwire2 + Nwires2 / 2.);
  
  float centerWire = (centerWire_left + centerWire_right) / 2.;

  //---- WGs around dead HV segment regions may need special treatment...
  //---- This is not addressed here.
    
  float sigmaWire = 0.;
  if(true == wHit.isNearDeadWG() || wgroups.size()>2){
    //---- worst possible case; take most conservative approach
    for(unsigned int iWG=0;iWG<wgroups.size();iWG++){
      sigmaWire+=layergeom_->yResolution( wgroups[iWG] );
    }
  }
  else if(2==wgroups.size()){
    //---- 2 WGs - get the larger error (overestimation if a single track is passing
    //---- between the WGs; underestimation if there are two separate signal sources)
    if(layergeom_->yResolution( wgroups[0] ) > layergeom_->yResolution( wgroups[1] )){
      sigmaWire  = layergeom_->yResolution( wgroups[0]);
    }    
    else{
      sigmaWire  = layergeom_->yResolution( wgroups[1]);
    }
  }  
  else if(1==wgroups.size()){
    //---- simple - just 1 WG
    sigmaWire  = layergeom_->yResolution( wgroups[0]);
  }
  
  // Find strips position and properties
  
  CSCRecHit2D::ChannelContainer strips = sHit.strips();
  int tmax = sHit.tmax();
  int nStrip = strips.size();
  int idCenterStrip = nStrip/2;
  int centerStrip = strips[idCenterStrip];
  
  // Retrieve strip pulseheights from the CSCStripHit
  const std::vector<float>& adc    = sHit.s_adc();
  const std::vector<float>& adcRaw = sHit.s_adcRaw();

  std::vector<float> adc2;
  std::vector<float> adc2Raw;

  LogTrace("CSCRecHit") << "CSCMake2DRecHit: dump of adc values to be added to rechit follows...";

  for ( int iStrip = 0; iStrip < nStrip; ++iStrip) {

    adc2.clear();
    adc2Raw.clear();
    for ( int t = 0; t < 4; ++t ){
      adc2.push_back(adc[t+iStrip*4]);
      adc2Raw.push_back(adcRaw[t+iStrip*4]);
    }
    // OLD: Rechit takes _calibrated_ adc values
    // adcMap.put( strips[iStrip], adc2.begin(), adc2.end() ); 
    // NEW: Rechit takes _raw_ adc values
    adcMap.put( strips[iStrip], adc2Raw.begin(), adc2Raw.end() ); 

    LogTrace("CSCRecHit") << "CSCMake2DRecHit: strip = " << strips[iStrip] << 
      " adcs= " << adc2Raw[0] << " " << adc2Raw[1] << " " << adc2Raw[2] << " " << adc2Raw[3];

    if (iStrip == nStrip/2 ) 
      tpeak = 50. * ( adc2[0]*(tmax-1) + adc2[1]*tmax + adc2[2]*(tmax+1) ) / (adc2[0]+adc2[1]+adc2[2]);
  }

  float positionWithinTheStrip= -99.;
  float sigmaWithinTheStrip = -99.;
  int quality = -1;
  LocalPoint lp0(0., 0.);
  
  float stripWidth = -99.;
  // If at the edge, then used 1 strip cluster only
  if ( centerStrip == 1 || centerStrip == specs_->nStrips() || nStrip < 2 ) {
    lp0 = layergeom_->stripWireIntersection( centerStrip, centerWire);
    positionWithinTheStrip = 0.;
    stripWidth = layergeom_->stripPitch(lp0);
    sigmaWithinTheStrip = stripWidth / sqrt_12;
    quality = 2;
  }
  else {
    // If not at the edge, used cluster of size ClusterSize:
    LocalPoint lp11  = layergeom_->stripWireIntersection( centerStrip, centerWire);
    stripWidth = layergeom_->stripPitch( lp11 );
    
    //---- Calculate local position within the strip
    float xWithinChamber = lp11.x();
    quality = 0;
    if(layergeom_->inside(lp11 )){// save time; this hit is to be discarded anyway - see isHitInFiducial(...)

      xMatchGatti_->findXOnStrip( id, layer_, sHit, centerStrip, 
         			xWithinChamber,
				stripWidth, tpeak, positionWithinTheStrip, 
				sigmaWithinTheStrip, quality);
    }				
    lp0 = LocalPoint( xWithinChamber, layergeom_->yOfWire(centerWire, xWithinChamber) );
  }
  
  // compute the errors in local x and y
  LocalError localerr = layergeom_->localError( centerStrip, 
						sigmaWithinTheStrip, sigmaWire );

  // Before storing the recHit, take the opportunity to change its time
  if (useTimingCorrections){
    float chipCorrection = 170. - recoConditions_->chipCorrection(id,centerStrip);
    float phaseCorrection = (sHit.stripsl1a()[0]>> (15-0) & 0x1)*25.;
    float chamberCorrection = recoConditions_->chamberTimingCorrection(id);

    GlobalPoint gp0 = layer_->surface().toGlobal(lp0);
    float tofCorrection = gp0.mag()/29.9792458;

    // printf("RecHit in e:%d s:%d r:%d c:%d l:%d strip:%d \n",id.endcap(),id.station(), id.ring(),id.chamber(),id.layer(),centerStrip);
    // printf("\t tpeak before = %5.2f \t chipCorr %5.2f phaseCorr %5.2f chamberCorr %5.2f tofCorr %5.2f\n",
    // 	   tpeak,chipCorrection, phaseCorrection,chamberCorrection,tofCorrection);
    tpeak = tpeak + chipCorrection + phaseCorrection + chamberCorrection-tofCorrection;
    // to more or less zero out the chambers (TOF will still be visible)
    if(id.station() ==1 && (id.ring() ==1 || id.ring() == 4))
      tpeak = tpeak - 205;
    else
      tpeak = tpeak - 216;
    // printf("\t tpeak after = %5.2f\n",tpeak);

  }

  // store rechit

   /// Retrive the L1APhase+strips combination
   CSCRecHit2D::ChannelContainer L1A_and_strips = sHit.stripsTotal();        /// L1A
  // (sigmaWithinTheStrip/stripWidth) is in strip widths just like positionWithinTheStrip is!
     CSCRecHit2D rechit( id, lp0, localerr, L1A_and_strips,                  /// L1A
		      adcMap, wgroups, tpeak, positionWithinTheStrip, 
		      sigmaWithinTheStrip/stripWidth, quality);

  /// To see RecHit content (L1A feature included) (to be commented out)
  // rechit.print();

  LogTrace("CSCRecHit") << "CSCMake2DRecHit: rechit created in layer " << id << "... \n" << rechit << "\n";

  return rechit;

}


bool CSCMake2DRecHit::isHitInFiducial( const CSCLayer* layer, const CSCRecHit2D& rh ) {

  bool isInFiducial = true;
  const CSCLayerGeometry* layergeom_ = layer->geometry();
  LocalPoint rhPosition =  rh.localPosition();
  // is the rechit within the chamber? 
  //(the problem occurs in ME11a/b otherwise it is OK)
  // we could use also 
  //bool inside( const Local3DPoint&, const LocalError&, float scale=1.f ) const;
  if(!layergeom_->inside(rhPosition)){
    isInFiducial = false;
  }
  
  return isInFiducial;
}
 

void CSCMake2DRecHit::setConditions( const CSCRecoConditions* reco ) {
  xMatchGatti_->setConditions( reco );
  // And cache for use here
  recoConditions_ = reco;
} 

