#ifndef Math_Point3D_h
#define Math_Point3D_h
// $Id: Point3D.h,v 1.6 2006/06/12 10:40:31 llista Exp $
#include <Rtypes.h>
#include <Math/Cartesian3D.h>
#include <Math/Point3D.h>
#include <Math/GenVector/CoordinateSystemTags.h>

namespace math {
  /// point in space with cartesian internal representation
  typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > XYZPointD;
  /// point in space with cartesian internal representation
  typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > XYZPointF;
  /// point in space with cartesian internal representation
  typedef XYZPointD XYZPoint;

  /// point in local coordinate system
  typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>, ROOT::Math::LocalCoordinateSystemTag> LocalPoint;
  /// point in global coordinate system
  typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>, ROOT::Math::GlobalCoordinateSystemTag> GlobalPoint;
}

#endif
