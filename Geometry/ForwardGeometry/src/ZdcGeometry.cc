#include "Geometry/CaloGeometry/interface/CaloGenericDetId.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/ForwardGeometry/interface/ZdcGeometry.h"
#include "Geometry/ForwardGeometry/interface/IdealZDCTrapezoid.h"
#include "ZdcHardcodeGeometryData.h"

ZdcGeometry::ZdcGeometry() :
   theTopology( new ZdcTopology ),
   lastReqDet_(DetId::Detector(0)), 
   lastReqSubdet_(0),
   m_ownsTopology ( true )
{
}

ZdcGeometry::ZdcGeometry( const ZdcTopology* topology) :
   theTopology(topology),
   lastReqDet_(DetId::Detector(0)), 
   lastReqSubdet_(0),
   m_ownsTopology ( false )
{
}

ZdcGeometry::~ZdcGeometry() 
{
   if( m_ownsTopology ) delete theTopology ;
}

const std::vector<DetId>& 
ZdcGeometry::getValidDetIds( DetId::Detector det,
			     int             subdet ) const 
{
   const std::vector<DetId>& baseIds ( CaloSubdetectorGeometry::getValidDetIds() ) ;
   if( det    == DetId::Detector( 0 ) &&
       subdet == 0                        )
   {
      return baseIds ;
   }
   
   if( lastReqDet_    != det    ||
       lastReqSubdet_ != subdet    ) 
   {
      lastReqDet_     = det    ;
      lastReqSubdet_  = subdet ;
      m_validIds.clear();
      m_validIds.reserve( baseIds.size() ) ;
   }

   if( m_validIds.empty() ) 
   {
      for( int i ( 0 ) ; i != baseIds.size() ; ++i ) 
      {
	 const DetId id ( baseIds[i] );
	 if( id.det()      == det    &&
	     id.subdetId() == subdet    )
	 { 
	    m_validIds.push_back( id ) ;
	 }
      }
      std::sort(m_validIds.begin(),m_validIds.end());
   }
   return m_validIds;
}

DetId ZdcGeometry::getClosestCell(const GlobalPoint& r) const
{
  // first find the side
  double z = r.z();
  double x = r.x();
  double y = r.y();
  double dz = 0.;
  double zt = 0.;

  int zside = 0;
  if(z >= 0)
    zside = 1;
  else
    zside =-1;

  bool isPositive = false;
  if(z>0)isPositive = true;
  z = fabs(z);
  
  // figure out if is closer to EM, HAD or LUM section
  HcalZDCDetId::Section section = HcalZDCDetId::Unknown;
  if(z<= theZSectionBoundaries[1])section = HcalZDCDetId::EM;
  if(theZSectionBoundaries[1]<z<= theZSectionBoundaries[2])section = HcalZDCDetId::LUM;
  if(z>theZSectionBoundaries[2])section = HcalZDCDetId::HAD;

  // figure out channel
  int channel = -1;
  if(section ==HcalZDCDetId::EM){
    if(x < theXChannelBoundaries[1]) channel = 1;
    if(theXChannelBoundaries[1]<= x <theXChannelBoundaries[2])channel = 2;
    if(theXChannelBoundaries[2]<= x <theXChannelBoundaries[3])channel = 3;
    if(theXChannelBoundaries[3]<= x <theXChannelBoundaries[4])channel = 4;
    if(x > theXChannelBoundaries[4])channel = 5;
  }
  
  if(section == HcalZDCDetId::LUM){
    if(z <= theZLUMChannelBoundaries[1])channel = 1;
    if(z > theZLUMChannelBoundaries[1])channel = 2;
  } 
  if(section == HcalZDCDetId::HAD){
    if(fabs(y) > dYPlate*sin(tiltangle))
      dz = (y > 0.) ?  dYPlate*cos(tiltangle) : -  dYPlate*sin(tiltangle);
    else
      dz = (y > 0.) ?  y/tan(tiltangle) : -y/tan(tiltangle); 
    zt = z - dz;
    if(zt< theZHadChannelBoundaries[1]) channel = 1;
    if(theZHadChannelBoundaries[1]<=  zt <theZHadChannelBoundaries[2])channel = 2;
    if(theZHadChannelBoundaries[2]<=  zt <theZHadChannelBoundaries[3])channel = 3;
    if(zt > theZHadChannelBoundaries[4])channel = 4;
  }
  
  HcalZDCDetId bestId  = HcalZDCDetId(section,isPositive,channel);
  return bestId;
}

unsigned int
ZdcGeometry::alignmentTransformIndexLocal( const DetId& id )
{
   const CaloGenericDetId gid ( id ) ;

   assert( gid.isZDC() ) ;

   unsigned int index ( 0 ) ;// to be implemented

   return index ;
}

unsigned int
ZdcGeometry::alignmentTransformIndexGlobal( const DetId& id )
{
   return (unsigned int)DetId::Calo ;
}

std::vector<HepPoint3D> 
ZdcGeometry::localCorners( const double* pv,
			   unsigned int  i,
			   HepPoint3D&   ref )
{
   return ( calogeom::IdealZDCTrapezoid::localCorners( pv, ref ) ) ;
}

CaloCellGeometry* 
ZdcGeometry::newCell( const GlobalPoint& f1 ,
		      const GlobalPoint& f2 ,
		      const GlobalPoint& f3 ,
		      CaloCellGeometry::CornersMgr* mgr,
		      const double*      parm ,
		      const DetId&       detId   ) 
{
   const CaloGenericDetId cgid ( detId ) ;

   assert( cgid.isZDC() ) ;

   return ( new calogeom::IdealZDCTrapezoid( f1, mgr, parm ) ) ;
}

