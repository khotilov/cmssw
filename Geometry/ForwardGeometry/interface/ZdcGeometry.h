#ifndef Geometry_ForwardGeometry_ZdcGeometry_h
#define Geometry_ForwardGeometry_ZDcGeometry_h

#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/ForwardGeometry/interface/ZdcTopology.h"

class ZdcGeometry : public CaloSubdetectorGeometry {
public:

  explicit ZdcGeometry(const ZdcTopology * topology);
  virtual ~ZdcGeometry();
  
  virtual const std::vector<DetId>& getValidDetIds( DetId::Detector det    = DetId::Detector ( 0 ) ,
						    int             subdet = 0   ) const;

  virtual DetId getClosestCell(const GlobalPoint& r) const ;

private:
  const ZdcTopology * theTopology;
  mutable DetId::Detector lastReqDet_;
  mutable int lastReqSubdet_;
      mutable std::vector<DetId> m_validIds;

};


#endif

