// This is CSCLayerGeometry.cc

#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCChamberSpecs.h>
#include <Geometry/CSCGeometry/interface/CSCWireGeometry.h>

#include <Geometry/CSCGeometry/src/CSCUngangedStripTopology.h>
#include <Geometry/CSCGeometry/src/CSCGangedStripTopology.h>
#include <Geometry/CSCGeometry/src/CSCWireGroupPackage.h>

#include <DataFormats/GeometryVector/interface/LocalPoint.h>

#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/Exception.h>

#include <CLHEP/Units/SystemOfUnits.h>

#include <algorithm>
#include <iostream>
#include <cmath>

CSCLayerGeometry::CSCLayerGeometry( int iChamberType,
         const TrapezoidalPlaneBounds& bounds,
         int nstrips, float stripOffset, float stripPhiPitch,
	 float whereStripsMeet, float extentOfStripPlane, float yCentreOfStripPlane,
	 const CSCWireGroupPackage& wg, float wireAngleInDegrees, double yOfFirstWire )
  :   TrapezoidalPlaneBounds( bounds ), theWireTopology( 0 ),
      theStripTopology( 0 ), myName( "CSCLayerGeometry" ), 
      chamberType( iChamberType ) {

  LogTrace("CSC") << myName <<": being constructed, this=" << this;

  // For simplicity cache the base class TPB dims even though they could be
  // retrieved indirectly through the TPB public interface
  apothem     = bounds.length() / 2.;  
  hTopEdge    = bounds.width() / 2.;
  hBottomEdge = bounds.widthAtHalfLength() - hTopEdge; // t+b=2w

  // Ganged strips in ME1A?
  bool gangedME1A = ( iChamberType == 1 && CSCChamberSpecs::gangedStrips() );

  CSCStripTopology* aStripTopology = 
        new CSCUngangedStripTopology(nstrips, stripPhiPitch,
	    extentOfStripPlane, whereStripsMeet, stripOffset, yCentreOfStripPlane );

  if ( gangedME1A ) {
    theStripTopology = new CSCGangedStripTopology( *aStripTopology, 16 );
    delete aStripTopology;
  }
  else {
    theStripTopology = aStripTopology;
  }

  if ( ! CSCChamberSpecs::realWireGeometry() ) {
    // Approximate ORCA_8_8_0 and earlier calculated geometry...
    float wangler = wireAngleInDegrees*degree; // convert angle to radians
    float wireCos = cos(wangler);
    float wireSin = sin(wangler);
    float y2 = apothem * wireCos + hBottomEdge * fabs(wireSin);
    float wireSpacing = wg.wireSpacing/10.; // in cm
    float wireOffset = -y2 + wireSpacing/2.;
    yOfFirstWire = wireOffset/wireCos;
  }

  theWireTopology = new CSCWireTopology( wg, yOfFirstWire, wireAngleInDegrees );

} 

CSCLayerGeometry::CSCLayerGeometry(const CSCLayerGeometry& melg) :
  TrapezoidalPlaneBounds(melg.hBottomEdge, melg.hTopEdge, melg.apothem,
			 0.5 * melg.thickness() ),
  theWireTopology(0), theStripTopology(0), 
  hBottomEdge(melg.hBottomEdge), hTopEdge(melg.hTopEdge),
  apothem(melg.apothem)
{
  // CSCStripTopology is abstract, so need clone()
  if (melg.theStripTopology) theStripTopology = melg.theStripTopology->clone();
  // CSCWireTopology is concrete, so direct copy
  if (melg.theWireTopology) theWireTopology = new CSCWireTopology(*(melg.theWireTopology));
}

CSCLayerGeometry& CSCLayerGeometry::operator=(const CSCLayerGeometry& melg)
{
  if ( &melg != this ) {
    delete theStripTopology;
    if ( melg.theStripTopology )
      theStripTopology=melg.theStripTopology->clone();
    else
      theStripTopology=0;

    delete theWireTopology;
    if ( melg.theWireTopology )
      theWireTopology=new CSCWireTopology(*(melg.theWireTopology));
    else
      theWireTopology=0;

    hBottomEdge     = melg.hBottomEdge;
    hTopEdge        = melg.hTopEdge;
    apothem         = melg.apothem;
  }
  return *this;
}

CSCLayerGeometry::~CSCLayerGeometry()
{
  LogTrace("CSC") << myName << ": being destroyed, this=" << this << 
    "\nDeleting theStripTopology=" << theStripTopology << 
    " and theWireTopology=" << theWireTopology;
  delete theStripTopology;
  delete theWireTopology;
}


LocalPoint 
CSCLayerGeometry::stripWireIntersection( int strip, float wire ) const
{
  // This allows _float_ wire no. so that we can calculate the
  // intersection of a strip with the mid point of a wire group 
  // containing an even no. of wires (which is not an actual wire),
  // as well as for a group containing an odd no. of wires.

  // Equation of wire and strip as straight lines in local xy
  // y = mx + c where m = tan(angle w.r.t. x axis)
  // At the intersection x = -(cs-cw)/(ms-mw)
  // At y=0, 0 = ms * xOfStrip(strip) + cs => cs = -ms*xOfStrip
  // At x=0, yOfWire(wire) = 0 + cw => cw = yOfWire

  float ms = tan( stripAngle(strip) );
  float mw = tan( wireAngle() );
  float xs = xOfStrip(strip);
  float xi = ( ms * xs + yOfWire(wire) ) / ( ms - mw );
  float yi = ms * (xi - xs );

  return LocalPoint(xi, yi);
}

LocalPoint CSCLayerGeometry::stripWireGroupIntersection( int strip, int wireGroup) const
{
  // middleWire is only an actual wire for a group with an odd no. of wires
  float middleWire = middleWireOfGroup( wireGroup );
  return stripWireIntersection(strip, middleWire);
}

float CSCLayerGeometry::stripAngle(int strip) const 
{
  // Cleverly subtly change meaning of stripAngle once more.
  // In TrapezoidalStripTopology it is angle measured
  // counter-clockwise from y axis. 
  // In APTST and RST it is angle measured 
  // clockwise from y axis.
  // Output of this function is measured counter-clockwise 
  // from x-axis, so it is a conventional 2-dim azimuthal angle
  // in the (x,y) local coordinates

  // We want angle at centre of strip (strip N covers
  // *float* range N-1 to N-epsilon)

  return M_PI_2 - theStripTopology->stripAngle(strip-0.5);
}


LocalPoint CSCLayerGeometry::localCenterOfWireGroup( int wireGroup ) const {

  // It can use CSCWireTopology::yOfWireGroup for y,
  // But x requires mixing with 'extent' of wire plane

  // If the wires are NOT tilted, default to simple calculation...
  if ( fabs(wireAngle() ) < 1.E-6 )  {
    float y = yOfWireGroup( wireGroup );
    return LocalPoint( 0., y );
  }
  else {
    // w is "wire" at the center of the wire group
    float w = middleWireOfGroup( wireGroup );
    std::vector<float> store = theWireTopology->wireValues( w );
    return LocalPoint( store[0], store[1] );
  }
}

float CSCLayerGeometry::lengthOfWireGroup( int wireGroup ) const {
  // Return length of 'wire' in the middle of the wire group
   float w = middleWireOfGroup( wireGroup );
   std::vector<float> store = theWireTopology->wireValues( w );
   return store[2];
}
    
void CSCLayerGeometry::setTopology( CSCStripTopology * newTopology ) {
   delete theStripTopology;
   theStripTopology = newTopology;
}

std::ostream & operator<<(std::ostream & stream, const CSCLayerGeometry & lg) {
  stream << "LayerGeometry " << std::endl
         << "------------- " << std::endl
         << "numberOfStrips        " << lg.numberOfStrips() << std::endl
         << "numberOfWires         " << lg.numberOfWires() << std::endl
         << "numberOfWireGroups    " << lg.numberOfWireGroups() << std::endl
         << "wireAngle  (rad)      " << lg.wireAngle() << std::endl
    //         << "wireAngle  (deg)      " << lg.theWireAngle << std::endl
    //         << "sin(wireAngle)        " << lg.theWireSin << std::endl
    //         << "cos(wireAngle)        " << lg.theWireCos << std::endl
         << "wirePitch             " << lg.wirePitch() << std::endl
         << "stripPitch            " << lg.stripPitch() << std::endl
    //         << "numberOfWiresPerGroup " << lg.theNumberOfWiresPerGroup << std::endl
    //         << "numberOfWiresInLastGroup " << lg.theNumberOfWiresInLastGroup << std::endl
    //         << "wireOffset            " << lg.theWireOffset << std::endl
    //         << "whereStripsMeet       " << lg.whereStripsMeet << std::endl;
         << "hBottomEdge           " << lg.hBottomEdge << std::endl
         << "hTopEdge              " << lg.hTopEdge << std::endl
         << "apothem               " << lg.apothem << std::endl;
    return stream;
}

