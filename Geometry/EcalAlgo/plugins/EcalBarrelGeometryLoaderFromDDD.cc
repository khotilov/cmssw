#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"

#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"

#include "Geometry/CaloEventSetup/interface/CaloGeometryLoader.h"
#include "Geometry/CaloEventSetup/interface/CaloGeometryLoader.icc"

template class CaloGeometryLoader< EcalBarrelGeometry > ;

//#include "DetectorDescription/Core/interface/DDInit.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"


using namespace std;

typedef CaloGeometryLoader< EcalBarrelGeometry > EcalBGL ;

template <>
void 
EcalBGL::fillGeom( EcalBarrelGeometry*     geom ,
		   const EcalBGL::ParmVec& vv ,
		   const HepGeom::Transform3D&   tr ,
		   const DetId&            id    )
{
   std::vector<double> pv ;
   pv.reserve( vv.size() ) ;
   for( unsigned int i ( 0 ) ; i != vv.size() ; ++i )
   {
      const double factor ( 1==i || 2==i || 6==i || 10==i ? 1 : k_ScaleFromDDDtoGeant ) ;
      pv.push_back( factor*vv[i] ) ;
   }

   CaloCellGeometry::CornersVec corners ( geom->cornersMgr() ) ;
   corners.resize() ;

   TruncatedPyramid::createCorners( pv, tr, corners ) ;

   const double* parmPtr ( CaloCellGeometry::getParmPtr( pv, 
							 geom->parMgr(), 
							 geom->parVecVec() ) ) ;

   TruncatedPyramid* cell ( new TruncatedPyramid( corners , parmPtr ) ) ;

   geom->addCell( id, cell );
}

template <>
void 
EcalBGL::fillNamedParams( DDFilteredView      fv,
			  EcalBarrelGeometry* geom )
{
   bool doSubDets = fv.firstChild();

   while( doSubDets )
   {
      DDsvalues_type sv(fv.mergedSpecifics());
        
      //nxtalPhi
      DDValue valnPhi("nxtalPhi");
      if( DDfetch( &sv, valnPhi ) )
      {
	 const vector<double>& fvec = valnPhi.doubles();

	 // this parameter can only appear once
	 assert(fvec.size() == 1);
	 geom->setNumXtalsPhiDirection((int)fvec[0]);
      }
      else
	 continue;
          

      DDValue valnEta("nxtalEta");
      if( DDfetch( &sv, valnEta ) ) 
      {
	 const vector<double>& fmvec = valnEta.doubles();

	 // there can only be one such value
	 assert(fmvec.size() == 1);
            
	 geom->setNumXtalsEtaDirection((int)fmvec[0]);
      }
      else
	 // once we find nxtalPhi, the rest must also be defined
	 assert( 1 == 0 );

      //EtaBaskets
      DDValue valEtaB("EtaBaskets");
      if( DDfetch( &sv, valEtaB ) ) 
      {
	 const vector<double>& ebvec = valEtaB.doubles();
	 assert(ebvec.size() > 0);
	 vector<int> EtaBaskets;
	 EtaBaskets.resize(ebvec.size());
	 for (unsigned i = 0; i<ebvec.size(); ++i) 
	    EtaBaskets[i] = (int) ebvec[i];
	 geom->setEtaBaskets(EtaBaskets);
      }
      else
	 // once we find nxtalPhi, the rest must also be defined
	 assert( 1 == 0 );

      //PhiBaskets
      DDValue valPhi("PhiBaskets");
      if (DDfetch(&sv,valPhi)) 
      {
	 const vector<double> & pvec = valPhi.doubles();
	 assert(pvec.size() > 0);
	 geom->setBasketSizeInPhi((int)pvec[0]);
      }
      else
	 // once we find nxtalPhi, the rest must also be defined
	 assert( 1 == 0 );

      break;
    
      doSubDets = fv.nextSibling(); // go to next layer
   }
}

