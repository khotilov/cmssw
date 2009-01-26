#ifndef GEOMETRY_CALOEVENTSETUP_CALOGEOMETRYDBWRITER_H
#define GEOMETRY_CALOEVENTSETUP_CALOGEOMETRYDBWRITER_H 1

#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Framework/interface/Event.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
//#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/GeometryObjects/interface/PCaloGeometry.h"
#include "Geometry/Records/interface/PEcalBarrelRcd.h"
#include "Geometry/Records/interface/PEcalEndcapRcd.h"
#include "Geometry/Records/interface/PEcalPreshowerRcd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

class CaloGeometryDBWriter
{
   public:

      typedef CaloSubdetectorGeometry::TrVec  TrVec      ;
      typedef CaloSubdetectorGeometry::DimVec DimVec     ;
      typedef CaloSubdetectorGeometry::IVec   IVec       ;

      static bool writeFlag() { return true ; }

      static void write( const TrVec&  tvec, 
			 const DimVec& dvec, 
			 const IVec&   ivec,
			 std::string   tag   )
      {
	 PCaloGeometry* peg = new PCaloGeometry( tvec , dvec, ivec );  
  
	 edm::Service<cond::service::PoolDBOutputService> mydbservice;
	 if( !mydbservice.isAvailable() )
	 {
	    edm::LogError("PCaloDBGeometryBuilder")<<"PoolDBOutputService unavailable";
	 }
	 else
	 {
	    if ( mydbservice->isNewTagRequest( tag ) ) 
	    {
	       mydbservice->createNewIOV<PCaloGeometry>( 
		  peg,
		  mydbservice->beginOfTime(),
		  mydbservice->endOfTime(),
		  tag ) ;
	    }
	    else 
	    {
	       mydbservice->appendSinceTime<PCaloGeometry>(
		  peg,
		  mydbservice->currentTime(),
		  tag ) ;
	    }
	 }
      }

      CaloGeometryDBWriter() {}
      virtual ~CaloGeometryDBWriter() {}
};

#endif
