#ifndef RecoLocalTracker_SiPixelRecHits_PixelCPEGeneric_H
#define RecoLocalTracker_SiPixelRecHits_PixelCPEGeneric_H

// \class PixelCPEGeneric  -- a generalized CPE reco for the idealized detector
//
// The basic idea of this class is to use generic formulae in order
// to achieve clean and minimal code.  It should work for
// - both normal and big pixels
// - both barrel and forward
// - both "FromDetPosition" and "FromTrackAngles" (i.e. by the track fit)
//
// This is possible since, in its nature, the original "ORCA" algorithm by 
// Danek and Susana is the same in both X and Y directions, provided that
// one correctly computes angles alpha_ and beta_ up front.  Thus, all
// geometrical and special corrections are dropped, since the presumption
// is that alpha_ and beta_ are determined as best as possible.  That means
// that they either come from the track, or, if they come from the 
// position of the DetUnit, they include all geometrical information 
// possible for this DetUnit:
// - for both the barrel and the forward, we use the cluster position 
//   instead of the center of the module/plaquette
// - for the forward, the tilt of the blades is included too
//
// In addtion, anything which is special for the computation of the lorentz
// angle is done in setTheDet() method.  So the algorithm per se does not
// need to worry about it.  This includes extra E*B term (a.k.a. "alpha2Order")
// and extra tilt in the forward.
//
// Thus, the formula for the computation of the hit position is very
// simple, and is described in Morris's note (IN ???) on the generalizaton
// of the pixel algorithm.

#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEBase.h"
#include "CondTools/SiPixel/interface/SiPixelDBErrorParametrization.h"

// Already defined in the base class
//#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
//#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
//#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
//#include "Geometry/CommonDetAlgo/interface/MeasurementPoint.h"
//#include "Geometry/CommonDetAlgo/interface/MeasurementError.h"
//#include "Geometry/Surface/interface/GloballyPositioned.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <utility>
#include <vector>


#if 0
/** \class PixelCPEGeneric
 * Perform the position and error evaluation of pixel hits using 
 * the Det angle to estimate the track impact angle 
*/
#endif

class MagneticField;
class SiPixelCPEParmErrors;
class PixelCPEGeneric : public PixelCPEBase
{
 public:
  // PixelCPEGeneric( const DetUnit& det );
  PixelCPEGeneric(edm::ParameterSet const& conf, const MagneticField*, const SiPixelCPEParmErrors*);
  ~PixelCPEGeneric() {;}

  LocalPoint localPosition (const SiPixelCluster& cluster, const GeomDetUnit & det) const; 
  
  // However, we do need to implement localError().
  LocalError localError   (const SiPixelCluster& cl, const GeomDetUnit & det) const;
  
  MeasurementPoint measurementPosition ( const SiPixelCluster&, 
					 const GeomDetUnit & det) const;
/*   MeasurementError measurementError    ( const SiPixelCluster&,  */
/* 					  const GeomDetUnit & det) const; */


 private:
  //--------------------------------------------------------------------
  //  Methods.
  //------------------------------------------------------------------
  double
    generic_position_formula( int size,                //!< Size of this projection.
			      double Q_f,              //!< Charge in the first pixel.
			      double Q_l,              //!< Charge in the last pixel.
			      double upper_edge_first_pix, //!< As the name says.
			      double lower_edge_last_pix,  //!< As the name says.
			      double half_lorentz_shift,   //!< L-shift at half thickness
			      double cot_angle,            //!< cot of alpha_ or beta_
			      double pitch,            //!< thePitchX or thePitchY
			      bool first_is_big,       //!< true if the first is big
			      bool last_is_big,        //!< true if the last is big
			      double eff_charge_cut_low, //!< Use edge if > W_eff (in pix) &&&
			      double eff_charge_cut_high,//!< Use edge if < W_eff (in pix) &&&
			      double size_cut,           //!< Use edge when size == cuts
			      float & cot_angle_from_length  //!< Aux output: angle from len
			      ) const;

  void
    collect_edge_charges(const SiPixelCluster& cluster,  //!< input, the cluster
			 float & Q_f_X,              //!< output, Q first  in X 
			 float & Q_l_X,              //!< output, Q last   in X
			 float & Q_m_X,              //!< output, Q middle in X
			 float & Q_f_Y,              //!< output, Q first  in Y 
			 float & Q_l_Y,              //!< output, Q last   in Y
			 float & Q_m_Y               //!< output, Q middle in Y
			 ) const;
  
  
  //--- Errors squared in x and y.  &&& Need to be revisited.
  float err2X(bool&, int&) const;
  float err2Y(bool&, int&) const;

  //--- Cuts made externally settable
  double the_eff_charge_cut_lowX;
  double the_eff_charge_cut_lowY;
  double the_eff_charge_cut_highX;
  double the_eff_charge_cut_highY;
  double the_size_cutX;
  double the_size_cutY;

	SiPixelDBErrorParametrization * dbErrors_;

 protected:
  //--- These functions are no longer needed, yet they are declared 
  //--- pure virtual in the base class.
  float xpos( const SiPixelCluster& ) const { return -999000.0; }  // &&& should abort
  float ypos( const SiPixelCluster& ) const { return -999000.0; }  // &&& should abort

};

#endif




