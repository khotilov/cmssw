#ifndef GEOMETRY_CALOGEOMETRY_CALOGEOMETRYDBEP_H
#define GEOMETRY_CALOGEOMETRY_CALOGEOMETRYDBEP_H 1

// system include files
#include <memory>
#include "boost/shared_ptr.hpp"

// user include files
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "CondFormats/AlignmentRecord/interface/GlobalPositionRcd.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/CaloGeometry/interface/PreshowerStrip.h"
#include "Geometry/CaloGeometry/interface/CaloGenericDetId.h"
#include "Geometry/CaloEventSetup/interface/CaloGeometryLoader.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"

#include "CondFormats/GeometryObjects/interface/PCaloGeometry.h"


#include "CondFormats/Alignment/interface/Alignments.h"

#include <Math/Transform3D.h>
#include <Math/EulerAngles.h>

//Forward declaration

//
// class declaration
//

template <class T, class U>
class CaloGeometryDBEP : public edm::ESProducer 
{
   public:

      typedef boost::shared_ptr< CaloSubdetectorGeometry > PtrType ;
      typedef CaloSubdetectorGeometry::TrVec  TrVec      ;
      typedef CaloSubdetectorGeometry::DimVec DimVec     ;
      typedef CaloSubdetectorGeometry::IVec   IVec       ;
      typedef HepGeom::Point3D<double> HepPoint3D;
      typedef HepGeom::Transform3D  HepTransform3D;

      CaloGeometryDBEP<T,U>( const edm::ParameterSet& ps ) :
	 m_applyAlignment ( ps.getUntrackedParameter<bool>("applyAlignment", false) )
      {
	 setWhatProduced( this,
			  &CaloGeometryDBEP<T,U>::produceAligned,
			  edm::es::Label( T::producerName() ) ) ;//+std::string("TEST") ) ) ;
      }

      virtual ~CaloGeometryDBEP<T,U>() {}
      PtrType produceAligned( const typename T::AlignedRecord& iRecord ) 
      {
	 const Alignments* alignPtr  ( 0 ) ;
	 const Alignments* globalPtr ( 0 ) ;
	 if( m_applyAlignment ) // get ptr if necessary
	 {
	    edm::ESHandle< Alignments >                                      alignments ;
	    iRecord.template getRecord< typename T::AlignmentRecord >().get( alignments ) ;

	    assert( alignments.isValid() && // require valid alignments and expected size
		    ( alignments->m_align.size() == T::numberOfAlignments() ) ) ;
	    alignPtr = alignments.product() ;

	    edm::ESHandle< Alignments >                          globals   ;
	    iRecord.template getRecord<GlobalPositionRcd>().get( globals ) ;

	    assert( globals.isValid() ) ;
	    globalPtr = globals.product() ;
	 }

	 TrVec  tvec ;
	 DimVec dvec ;
	 IVec   ivec ;

	 if( U::writeFlag() )
	 {
	    edm::ESHandle<CaloSubdetectorGeometry> pG ;
	    iRecord.get( T::producerName() + std::string("_master"), pG ) ; 

	    const CaloSubdetectorGeometry* pGptr ( pG.product() ) ;

	    pGptr->getSummary( tvec, ivec, dvec ) ;

	    U::write( tvec, dvec, ivec, T::dbString() ) ;
	 }
	 else
	 {
	    edm::ESHandle<PCaloGeometry> pG ;
	    iRecord.template getRecord<typename T::PGeometryRecord >().get( pG ) ; 

	    tvec = pG->getTranslation() ;
	    dvec = pG->getDimension() ;
	    ivec = pG->getIndexes() ;
	 }	 
/*
	 if( toDB )
	 {
	    edm::ESHandle<CaloSubdetectorGeometry> pG ;
	    iRecord.get( T::producerName() + std::string("xml"), pG ) ; 

	    const CaloSubdetectorGeometry* pGptr ( pG.product() ) ;

	    pGptr->getSummary( tvec, ivec, dvec ) ;
	 }
	 else
	 {
	    // from DB

	    edm::ESHandle<PCaloGeometry> peG;   
	    iRecord.template getRecord<typename T::PGeometryRecord >().get( peG ) ; 
         
	    tvec = peG->getTranslation() ;
	    dvec = peG->getDimension() ;
	    ivec = peG->getIndexes() ;
	    }*/
//*********************************************************************************************

	 const unsigned int nTrParm ( tvec.size()/T::k_NumberOfCellsForCorners ) ;

	 assert( dvec.size() == T::k_NumberOfShapes * T::k_NumberOfParametersPerShape ) ;

	 PtrType ptr ( new T ) ;

	 ptr->fillDefaultNamedParameters() ;

	 ptr->allocateCorners( T::k_NumberOfCellsForCorners ) ;

	 ptr->allocatePar(    dvec.size() ,
			      T::k_NumberOfParametersPerShape ) ;

	 for( unsigned int i ( 0 ) ; i != T::k_NumberOfCellsForCorners ; ++i )
	 {
	    const unsigned int nPerShape ( T::k_NumberOfParametersPerShape ) ;
	    DimVec dims ;
	    dims.reserve( nPerShape ) ;

	    const unsigned int indx ( ivec.size()==1 ? 0 : i ) ;

	    DimVec::const_iterator dsrc ( dvec.begin() + ivec[indx]*nPerShape ) ;

	    for( unsigned int j ( 0 ) ; j != nPerShape ; ++j )
	    {
	       dims.push_back( *dsrc ) ;
	       ++dsrc ;
	    }

	    const double* myParm ( CaloCellGeometry::getParmPtr( dims, 
								 ptr->parMgr(), 
								 ptr->parVecVec() ) ) ;


	    const DetId id ( T::DetIdType::detIdFromDenseIndex( i ) ) ;
    
	    const unsigned int iGlob ( T::alignmentTransformIndexGlobal( id ) ) ;

	    const AlignTransform* gt ( 0 == globalPtr ? 0 :
				       ( iGlob < globalPtr->m_align.size() ?
					 &globalPtr->m_align[ iGlob ] : 0 ) ) ;

	    const unsigned int iLoc ( T::alignmentTransformIndexLocal( id ) ) ;

	    assert( 0 == alignPtr ||
		    iLoc < alignPtr->m_align.size() ) ;

	    const AlignTransform* at ( 0 == alignPtr ? 0 :
				       &alignPtr->m_align[ iLoc ] ) ;

	    assert( 0 == at || ( T::alignmentTransformIndexLocal( DetId( at->rawId() ) ) == iLoc ) ) ;

	    const CaloGenericDetId gId ( id ) ;

	    HepPoint3D lRef ;
	    const std::vector<HepPoint3D> lc ( T::localCorners( &dims.front(), i, lRef ) ) ;

	    const HepPoint3D lBck ( 0.25*(lc[4]+lc[5]+lc[6]+lc[7] ) ) ; // ctr rear  face in local
	    const HepPoint3D lCor ( lc[0] ) ;

	    const unsigned int jj ( i*nTrParm ) ;
	    HepTransform3D tr ;
	    const ROOT::Math::Translation3D tl ( tvec[jj], tvec[jj+1], tvec[jj+2] ) ;
	    const ROOT::Math::EulerAngles ea (
	       6==nTrParm ?
	       ROOT::Math::EulerAngles( tvec[jj+3], tvec[jj+4], tvec[jj+5] ) :
	       ROOT::Math::EulerAngles() ) ;
	    const ROOT::Math::Transform3D rt ( ea, tl ) ;
	    double xx,xy,xz,dx,yx,yy,yz,dy,zx,zy,zz,dz;
	    rt.GetComponents(xx,xy,xz,dx,yx,yy,yz,dy,zx,zy,zz,dz) ;
	    tr = HepTransform3D( CLHEP::HepRep3x3( xx, xy, xz,
						   yx, yy, yz,
						   zx, zy, zz ), 
				 CLHEP::Hep3Vector(dx,dy,dz)     );

	    const HepTransform3D atr ( 0 == at ? tr :
				       ( 0 == gt ? at->transform()*tr :
					 gt->transform()*at->transform()*tr ) ) ;

	    const HepPoint3D  gRef ( atr*lRef ) ;
	    const GlobalPoint fCtr ( gRef.x(), gRef.y(), gRef.z() ) ;
	    const HepPoint3D  gBck ( atr*lBck ) ;
	    const GlobalPoint fBck ( gBck.x(), gBck.y(), gBck.z() ) ;
	    const HepPoint3D  gCor ( atr*lCor ) ;
	    const GlobalPoint fCor ( gCor.x(), gCor.y(), gCor.z() ) ;

	    CaloCellGeometry* cell ( T::newCell(  fCtr, fBck, fCor,
						  ptr->cornersMgr() , myParm ) ) ;

	    ptr->addCell( id, cell ) ;    
	 }

	 ptr->initializeParms() ; // initializations; must happen after cells filled

	 return ptr ; 
      }

   private:

      bool        m_applyAlignment ;
};

#endif
