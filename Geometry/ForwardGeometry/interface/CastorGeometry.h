#ifndef Geometry_ForwardGeometry_CastorGeometry_h
#define Geometry_ForwardGeometry_CastorGeometry_h 1

#include "CondFormats/AlignmentRecord/interface/CastorAlignmentRcd.h"
#include "DataFormats/HcalDetId/interface/HcalCastorDetId.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/ForwardGeometry/interface/CastorTopology.h"
#include "Geometry/Records/interface/CastorGeometryRecord.h"
#include "Geometry/Records/interface/PCastorRcd.h"

#include <vector>

class CastorGeometry : public CaloSubdetectorGeometry 
{
   public:

      typedef CastorAlignmentRcd   AlignmentRecord ;
      typedef CastorGeometryRecord AlignedRecord   ;
      typedef PCastorRcd           PGeometryRecord ;
      typedef HcalCastorDetId      DetIdType       ;

      enum { k_NumberOfCellsForCorners = HcalCastorDetId::kSizeForDenseIndexing } ;

      enum { k_NumberOfShapes = 4 } ;

      enum { k_NumberOfParametersPerShape = 6 } ;

      static std::string dbString() { return "PCastorRcd" ; }

      virtual unsigned int numberOfTransformParms() const { return 3 ; }

      virtual unsigned int numberOfShapes() const { return k_NumberOfShapes ; }
      virtual unsigned int numberOfParametersPerShape() const { return k_NumberOfParametersPerShape ; }

      CastorGeometry() ;

      explicit CastorGeometry(const CastorTopology * topology);
      virtual ~CastorGeometry();

      virtual const std::vector<DetId>& getValidDetIds(
	 DetId::Detector det    = DetId::Detector(0) ,
	 int             subdet = 0 ) const ;

      virtual DetId getClosestCell(const GlobalPoint& r) const ;

      static std::string producerTag() { return "CASTOR" ; }

      static unsigned int numberOfAlignments() { return 0 ; }

      static unsigned int alignmentTransformIndexLocal( const DetId& id ) ;

      static unsigned int alignmentTransformIndexGlobal( const DetId& id ) ;

      static std::vector<HepPoint3D> localCorners( const double* pv, 
						   unsigned int  i,
						   HepPoint3D&   ref ) ;

      static CaloCellGeometry* newCell( const GlobalPoint& f1 ,
					const GlobalPoint& f2 ,
					const GlobalPoint& f3 ,
					CaloCellGeometry::CornersMgr* mgr,
					const double*      parm,
					const DetId&       detId     ) ;

private:

      const CastorTopology * theTopology;
      mutable DetId::Detector lastReqDet_;
      mutable int lastReqSubdet_;
      mutable std::vector<DetId> m_validIds;
      bool m_ownsTopology ;
};


#endif
