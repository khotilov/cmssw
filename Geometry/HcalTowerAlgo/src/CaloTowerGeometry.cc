#include "Geometry/CaloGeometry/interface/IdealObliquePrism.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGenericDetId.h"
#include "Geometry/HcalTowerAlgo/interface/CaloTowerGeometry.h"

CaloTowerGeometry::CaloTowerGeometry() {
}
  

CaloTowerGeometry::~CaloTowerGeometry() {}


unsigned int
CaloTowerGeometry::alignmentTransformIndexLocal( const DetId& id )
{
   const CaloGenericDetId gid ( id ) ;

   assert( gid.isCaloTower() ) ;

   const CaloTowerDetId cid ( id ) ;

   const unsigned int iea ( cid.ietaAbs() ) ;

   const unsigned int ip ( ( cid.iphi() - 1 )/4 ) ;

   const int izoff ( ( cid.zside() + 1 )/2 ) ;

   const unsigned int offset ( izoff*3*18) ;

   return ( offset + ip + ( CaloTowerDetId::kEndIEta < iea ? 36 :
			    ( CaloTowerDetId::kBarIEta < iea ? 18 : 0 ) ) ) ;
}

unsigned int
CaloTowerGeometry::alignmentTransformIndexGlobal( const DetId& id )
{
   return (unsigned int) DetId::Calo - 1 ;
}

std::vector<HepGeom::Point3D<double> > 
CaloTowerGeometry::localCorners( const double* pv,
				 unsigned int  i,
				 HepGeom::Point3D<double> &   ref )
{
   return ( calogeom::IdealObliquePrism::localCorners( pv, ref ) ) ;
}

CaloCellGeometry* 
CaloTowerGeometry::newCell( const GlobalPoint& f1 ,
			    const GlobalPoint& f2 ,
			    const GlobalPoint& f3 ,
			    CaloCellGeometry::CornersMgr* mgr,
			    const double*      parm ,
			    const DetId&       detId   ) 
{
   const CaloGenericDetId cgid ( detId ) ;

   assert( cgid.isCaloTower() ) ;

   return ( new calogeom::IdealObliquePrism( f1, mgr, parm ) ) ;
}
