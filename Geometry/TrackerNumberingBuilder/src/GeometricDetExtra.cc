//#include "CondFormats/GeometryObjects/interface/PGeometricDetExtra.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDetExtra.h"
//#include "DetectorDescription/Core/interface/DDExpandedView.h"
/* #include "DetectorDescription/Base/interface/DDRotationMatrix.h" */
/* #include "DetectorDescription/Base/interface/DDTranslation.h" */
/* #include "DetectorDescription/Core/interface/DDSolidShapes.h" */
/* #include "DataFormats/GeometrySurface/interface/Surface.h" */
/* #include "DataFormats/GeometrySurface/interface/Bounds.h" */
//#include "DataFormats/DetId/interface/DetId.h"

//#include <vector>
//#include "FWCore/ParameterSet/interface/types.h"

/**
 * Constructors to be used when looping over DDD
 */
GeometricDetExtra::GeometricDetExtra( GeometricDet const * gd, DetId id, GeoHistory& gh,  double vol, double dens, double wgt, double cpy, std::string& mat, bool dd )
  : _mygd(gd), _geographicalId(id), _parents(gh), _volume(vol), _density(dens), _weight(wgt), _copy(cpy), _material(mat), _fromDD(dd) 
{ 
  //  std::cout << " made _mygd = " << _mygd << " to the geographicalId() " << _geographicalId << std::endl;

}

GeometricDetExtra::~GeometricDetExtra()
{ }

