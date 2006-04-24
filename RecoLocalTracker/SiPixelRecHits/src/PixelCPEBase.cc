//#include "Utilities/Configuration/interface/Architecture.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
//#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"

//#include "CommonDet/BasicDet/interface/Topology.h"
//#include "CommonDet/BasicDet/interface/Det.h"
//#include "CommonDet/BasicDet/interface/DetUnit.h"
//#include "CommonDet/BasicDet/interface/DetType.h"
//#include "Tracker/SiPixelDet/interface/PixelDetType.h"
//#include "Tracker/SiPixelDet/interface/PixelDigi.h"
//#include "Tracker/SiPixelDet/interface/PixelTopology.h"
//  #include "CommonDet/DetGeometry/interface/ActiveMediaShape.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEBase.h"

//#define DEBUG

// MessageLogger
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Magnetic field
#include "MagneticField/Engine/interface/MagneticField.h"

#include <iostream>
using namespace std;

const float PI = 3.141593;
const float degsPerRad = 57.29578;


//-----------------------------------------------------------------------------
//  A fairly boring constructor.  All quantities are DetUnit-dependent, and
//  will be initialized in setTheDet().
//-----------------------------------------------------------------------------
PixelCPEBase::PixelCPEBase(edm::ParameterSet const & conf, const MagneticField *mag) 
{
  //--- Lorentz angle tangent per Tesla
  theTanLorentzAnglePerTesla =
    conf.getParameter<double>("TanLorentzAnglePerTesla");

  //--- Algorithm's verbosity
  theVerboseLevel = 
    conf.getUntrackedParameter<int>("VerboseLevel",20);

  //-- Magnetic Field
  magfield_ = mag;
}


//-----------------------------------------------------------------------------
//  One function to cache the variables common for one DetUnit.
//-----------------------------------------------------------------------------
void
PixelCPEBase::setTheDet( const GeomDetUnit & det )const 
{
  if ( theDet == &det )
    return;       // we have already seen this det unit

  //--- This is a new det unit, so cache it
  theDet = dynamic_cast<const PixelGeomDetUnit*>( &det );
  if (! theDet) {
    // &&& Fatal error!  TO DO: throw an exception!
    assert(0);
  }

  //--- theDet->type() returns a GeomDetType, which implements subDetector()
  thePart = theDet->type().subDetector();
  switch (thePart) {
  case GeomDetType::PixelBarrel:
    // A barrel!  A barrel!
    break;
  case GeomDetType::PixelEndcap:
    // A forward!  A forward!
    break;
  default:
    LogDebug("PixelCPEBase") 
      << "PixelCPEBase:: a non-pixel detector type in here? Yuck!" ;
    //  &&& Should throw an exception here!
    assert(0);
  }
       
  //--- The location in of this DetUnit in a cyllindrical coord system (R,Z)
  //--- The call goes via BoundSurface, returned by theDet->surface(), but
  //--- position() is implemented in GloballyPositioned<> template
  //--- ( BoundSurface : Surface : GloballyPositioned<float> )
  theDetR = theDet->surface().position().perp();
  theDetZ = theDet->surface().position().z();


  //--- Define parameters for chargewidth calculation

  //--- bounds() is implemented in BoundSurface itself.
  theThickness = theDet->surface().bounds().thickness();

  //--- Cache the topology.
  theTopol
    = dynamic_cast<const RectangularPixelTopology*>( & (theDet->specificTopology()) );

  //---- The geometrical description of one module/plaquette
  theNumOfRow = theTopol->nrows();      // rows in x
  theNumOfCol = theTopol->ncolumns();   // cols in y
  std::pair<float,float> pitchxy = theTopol->pitch();
  thePitchX = pitchxy.first;            // pitch along x
  thePitchY = pitchxy.second;           // pitch along y

  //--- Find the offset
  MeasurementPoint  offset = 
    theTopol->measurementPosition( LocalPoint(0., 0.) );  
  theOffsetX = offset.x();
  theOffsetY = offset.y();

  //--- Find if the E field is flipped: i.e. whether it points
  //--- from the beam, or towards the beam.  (The voltages are
  //--- applied backwards on every other module in barrel and
  //--- blade in forward.)
  theSign = isFlipped() ? -1 : 1;

  //--- The Lorentz shift.
  theLShift = lorentzShift();

  if (theVerboseLevel > 5) {
    LogDebug("PixelCPEBase") << "***** PIXEL LAYOUT *****" << " theThickness = " << theThickness
    << " thePitchX  = " << thePitchX 
    << " thePitchY  = " << thePitchY 
    << " theOffsetX = " << theOffsetX 
    << " theOffsetY = " << theOffsetY 
    << " theLShift  = " << theLShift;
  }
}



//-----------------------------------------------------------------------------
//  Compute alpha_ and beta_ from the position of this DetUnit.
//  &&& Needs to be consolidated with estimatedAlphaForBarrel() below.
//  &&& Does not work for Forward!
//-----------------------------------------------------------------------------
void PixelCPEBase::
computeAnglesFromDetPosition(const SiPixelCluster & cl, 
			     const GeomDetUnit    & det ) const
{
  //calculate center
  float xmin = float(cl.minPixelRow()) + 0.5;
  float xmax = float(cl.maxPixelRow()) + 0.5;
  float xcenter = 0.5*( xmin + xmax );
  
  alpha_ = estimatedAlphaForBarrel(xcenter);
  beta_  = 0.0;                             // &&& ????
}


//-----------------------------------------------------------------------------
//  Compute alpha_ and beta_ from the LocalTrajectoryParameters.
//  Note: should become const after both localParameters() become const.
//-----------------------------------------------------------------------------
void PixelCPEBase::
computeAnglesFromTrajectory(const SiPixelCluster & cl,
			    const GeomDetUnit    & det, 
			    const LocalTrajectoryParameters & ltp)
{
  LocalVector localDir = ltp.momentum()/ltp.momentum().mag();

  // &&& Or, maybe we need to move to the local frame ???
  //  LocalVector localDir( theDet->toLocal(theState.globalDirection()));
  //thePart = theDet->type().part();

  float locx = localDir.x();
  float locy = localDir.y();
  float locz = localDir.z();

  alpha_ = acos(locx/sqrt(locx*locx+locz*locz));
  if ( isFlipped() )                    // &&& check for FPIX !!!
    alpha_ = PI - alpha_ ;

  beta_ = acos(locy/sqrt(locy*locy+locz*locz));
}



//-----------------------------------------------------------------------------
//  Estimate theAlpha for barrel, based on the det position.
//  &&& Needs to be consolidated from the above.
//-----------------------------------------------------------------------------
float 
PixelCPEBase::estimatedAlphaForBarrel(float centerx) const
{
  float tanalpha = theSign * (centerx-theOffsetX) * thePitchX / theDetR;
  return PI/2.0 - atan(tanalpha);
}





//-----------------------------------------------------------------------------
//  The local position.
//-----------------------------------------------------------------------------
LocalPoint
PixelCPEBase::localPosition(const SiPixelCluster& cluster, const GeomDetUnit & det) const
{
  setTheDet( det );
  //return theTopol->localPosition(measurementPosition(cluster, det)); 
  MeasurementPoint mp = measurementPosition(cluster, det);
  LocalPoint lp = theTopol->localPosition(mp);
  return lp;
}




//-----------------------------------------------------------------------------
//  Once we have the position, feed it to the topology to give us
//  the error.  
//  &&& APPARENTLY THIS METHOD IS NOT BEING USED ??? (Delete it?)
//-----------------------------------------------------------------------------
MeasurementError  
PixelCPEBase::measurementError( const SiPixelCluster& cluster, const GeomDetUnit & det) const 
{
  LocalPoint lp( localPosition(cluster, det) );
  LocalError le( localError(   cluster, det) );
  return theTopol->measurementError( lp, le );
}




//-----------------------------------------------------------------------------
//  Takes the cluster, calculates xpos() and ypos(), applies the Lorentz
//  shift, and then makes a MeasurementPoint.  This could really be
//  folded back into the localPosition().
//-----------------------------------------------------------------------------
MeasurementPoint 
PixelCPEBase::measurementPosition( const SiPixelCluster& cluster, const GeomDetUnit & det) const 
{
  if (theVerboseLevel > 15) {
    LogDebug("PixelCPEBase") <<
      "X-pos = " << xpos(cluster) << 
      " Y-pos = " << ypos(cluster) << 
      " Lshf = " << theLShift ;
  }
  return MeasurementPoint( xpos(cluster)-theLShift, 
  			   ypos(cluster));
}








//-----------------------------------------------------------------------------
// The isFlipped() is a silly way to determine which detectors are inverted.
// In the barrel for every 2nd ladder the E field direction is in the
// global r direction (points outside from the z axis), every other
// ladder has the E field inside. Something similar is in the 
// forward disks (2 sides of the blade). This has to be recognised
// because the charge sharing effect is different.
//
// The isFliped does it by looking and the relation of the local (z always
// in the E direction) to global coordinates. There is probably a much 
// better way.
//
// Plan: ignore it for the moment
//-----------------------------------------------------------------------------
bool 
PixelCPEBase::isFlipped() const
{
// &&& Not sure what the need is -- ask Danek.
  //  float tmp1 = theDet->toGlobal( Local3DPoint(0.,0.,0.) ).perp();
  //   float tmp2 = theDet->toGlobal( Local3DPoint(0.,0.,1.) ).perp();
  //   if ( tmp2<tmp1 ) return true;
  // else return false;
  float tmp1 = theDet->surface().toGlobal(Local3DPoint(0.,0.,0.)).perp();
  float tmp2 = theDet->surface().toGlobal(Local3DPoint(0.,0.,1.)).perp();
  //cout << " 1: " << tmp1 << " 2: " << tmp2 << endl;
  if ( tmp2<tmp1 ) return true;
  else return false;    
}




//-----------------------------------------------------------------------------
// From Danek: "geomCorrection() is sort of second order effect, ignore it for 
// the moment. I have to to derive it again and understand better what it means."
//-----------------------------------------------------------------------------
float 
PixelCPEBase::geomCorrection() const
{ 
  //@@ the geometrical correction are calculated only
  //@@ for the barrel part (am I right?)  &&& ??????????????????
  if (thePart == GeomDetType::PixelEndcap) return 0;
  else return theThickness / theDetR;
}



//-----------------------------------------------------------------------------
//  Lorentz shift, but valid only for the barrel.
//-----------------------------------------------------------------------------
float 
PixelCPEBase::lorentzShift() const
{
  // Implement only the x direction shift for now (OK for barrel)
  // &&& TEMPORARY 
  //LocalVector dir = theDet->driftDirection( LocalPoint(0,0));
  //LocalVector dir(0,0,1);  // &&& TEMPORARY HACK
  LocalVector dir = driftDirection(magfield_->inTesla(theDet->surface().position()) );

  // max shift in cm 
  float xdrift = dir.x()/dir.z() * theThickness;  
  // express the shift in units of pitch, 
  // divide by 2 to get the average correction
  float lshift = xdrift / thePitchX / 2.; 

  //cout << "Lorentz Drift = " << lshift << endl;
  //cout << "X Drift = " << dir.x() << endl;
  //cout << "Z Drift = " << dir.z() << endl;
 
  return lshift;  
}






//-----------------------------------------------------------------------------
//  Sum the pixels in the first and the last row, and the total.  Returns
//  a vector of three elements with q_first, q_last and q_total.
//  &&& Really need a simpler & cleaner way, this is very confusing...
//-----------------------------------------------------------------------------
vector<float> 
PixelCPEBase::xCharge(const vector<SiPixelCluster::Pixel>& pixelsVec, 
				    const float& xmin, const float& xmax)const
{
  vector<float> charge; 

  //calculate charge in the first and last pixel in y
  // and the total cluster charge
  float q1 = 0, q2 = 0, qm=0;
  int isize = pixelsVec.size();
  for (int i = 0;  i < isize; ++i) {
    if (pixelsVec[i].x == xmin) q1 += pixelsVec[i].adc;
    else if (pixelsVec[i].x == xmax) q2 += pixelsVec[i].adc;
    else qm += pixelsVec[i].adc;
  }
  charge.clear();
  charge.push_back(q1); 
  charge.push_back(q2); 
  charge.push_back(qm);
  return charge;
} 



//-----------------------------------------------------------------------------
//  Sum the pixels in the first and the last column, and the total.  Returns
//  a vector of three elements with q_first, q_last and q_total.
//  &&& Really need a simpler & cleaner way, this is very confusing...
//-----------------------------------------------------------------------------
vector<float> 
PixelCPEBase::yCharge(const vector<SiPixelCluster::Pixel>& pixelsVec,
				    const float& ymin, const float& ymax)const
{
  vector<float> charge; 

  //calculate charge in the first and last pixel in y
  // and the inner cluster charge
  float q1 = 0, q2 = 0, qm=0;
  int isize = pixelsVec.size();
  for (int i = 0;  i < isize; ++i) {
    if (pixelsVec[i].y == ymin) q1 += pixelsVec[i].adc;
    else if (pixelsVec[i].y == ymax) q2 += pixelsVec[i].adc;
    else if (pixelsVec[i].y < ymax && 
	     pixelsVec[i].y > ymin ) qm += pixelsVec[i].adc;
  }
  charge.clear();
  charge.push_back(q1); 
  charge.push_back(q2); 
  charge.push_back(qm);

  return charge;
} 




//-----------------------------------------------------------------------------
//  Drift direction.
//  NB: it's correct only for the barrel!  &&& Need to fix it for the forward.
//  Assumption: setTheDet() has been called already.
//-----------------------------------------------------------------------------
LocalVector 
PixelCPEBase::driftDirection( GlobalVector bfield )const 
{
  Frame detFrame(theDet->surface().position(), theDet->surface().rotation());
  LocalVector Bfield = detFrame.toLocal(bfield);


  //  if    (DetId(detID).subdetId()==  PixelSubdetector::PixelBarrel){
  float dir_x =  theTanLorentzAnglePerTesla * Bfield.y();
  float dir_y = -theTanLorentzAnglePerTesla * Bfield.x();
  float dir_z = 1.; // E field always in z direction
  LocalVector theDriftDirection = LocalVector(dir_x,dir_y,dir_z);

  if (theVerboseLevel > 0) {
    LogDebug("PixelCPEBase") << " The drift direction in local coordinate is " 
  	 << theDriftDirection    ;
  }

  return theDriftDirection;
  //  }
}




