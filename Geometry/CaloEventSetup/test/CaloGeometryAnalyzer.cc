// -*- C++ -*-
//
// Package:    CaloGeometryAnalyzer
// Class:      CaloGeometryAnalyzer
// 
/**\class CaloGeometryAnalyzer CaloGeometryAnalyzer.cc test/CaloGeometryAnalyzer/src/CaloGeometryAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/CaloGeometry/interface/PreshowerStrip.h"
#include "Geometry/EcalPreshowerAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include <fstream>
#include <iomanip>

//
// class decleration
//

class CaloGeometryAnalyzer : public edm::EDAnalyzer {
   public:
      explicit CaloGeometryAnalyzer( const edm::ParameterSet& );
      ~CaloGeometryAnalyzer();


      virtual void analyze( const edm::Event&, const edm::EventSetup& );
   private:
      // ----------member data ---------------------------
      void build( const CaloGeometry& cg , 
		  DetId::Detector     det, 
		  int                 subdetn, 
		  const char*         name);

      void ctrcor( const DetId::Detector   det     , 
		   const int               subdetn , 
		   const DetId&            did     ,
		   const CaloCellGeometry& cell ,
		   std::fstream&           fCtr    ,
		   std::fstream&           fCor       );
  int pass_;
  //  bool fullEcalDump_;

      EEDetId gid( unsigned int ix, 
		   unsigned int iy,
		   unsigned int iz,
		   const EEDetId& did ) const ;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CaloGeometryAnalyzer::CaloGeometryAnalyzer( const edm::ParameterSet& iConfig )
{
   //now do what ever initialization is needed
  pass_=0;
  //  fullEcalDump_=iConfig.getUntrackedParameter<bool>("fullEcalDump",false);
}


CaloGeometryAnalyzer::~CaloGeometryAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

EEDetId
CaloGeometryAnalyzer::gid( unsigned int ix, 
			   unsigned int iy,
			   unsigned int iz,
			   const EEDetId& did ) const
{
   return ( EEDetId::validDetId( ix, iy, iz ) ? 
	    EEDetId( ix, iy, iz ) : 
	    ( EEDetId::validDetId( ix+1, iy  , iz ) ? 
	      EEDetId( ix+1, iy  , iz ) :
	      ( EEDetId::validDetId( ix-1 , iy, iz ) ? 
		EEDetId( ix-1  , iy, iz ) :
		( EEDetId::validDetId( ix , iy+1, iz ) ? 
		  EEDetId( ix  , iy+1, iz ) :
		  ( EEDetId::validDetId( ix , iy-1, iz ) ? 
		    EEDetId( ix, iy-1, iz ) :
		    ( EEDetId::validDetId( ix+1, iy+1, iz ) ? 
		      EEDetId( ix+1, iy+1, iz ) :
		      ( EEDetId::validDetId( ix+1 , iy-1, iz ) ? 
			EEDetId( ix+1, iy-1, iz ) :
			( EEDetId::validDetId( ix-1, iy+1, iz ) ? 
			  EEDetId( ix-1, iy+1, iz ) :
			  ( EEDetId::validDetId( ix-1, iy-1, iz ) ? 
			    EEDetId( ix-1, iy-1, iz ) : did )  ) ) ) ) ) ) ) ) ;
}

void 
CaloGeometryAnalyzer::ctrcor( const DetId::Detector   det     , 
			      const int               subdetn , 
			      const DetId&            did     ,
			      const CaloCellGeometry& cell    ,
			      std::fstream&           fCtr    ,
			      std::fstream&           fCor       )
{
   if( det     == DetId::Ecal &&
       subdetn == EcalBarrel     )
   {
      const EBDetId ebid ( did ) ;
      const int ie  ( ebid.ieta() ) ;
      const int ip  ( ebid.iphi() ) ;
      fCtr << std::setw(4) << ie
	   << std::setw(4) << ip ;
      fCor << std::setw(4) << ie
	   << std::setw(4) << ip ;
   }
   if( det     == DetId::Ecal &&
       subdetn == EcalEndcap     )
   {
      const EEDetId eeid ( did ) ;
      const int ix ( eeid.ix() ) ;
      const int iy ( eeid.iy() ) ;
//      const int iz ( eeid.zside() ) ;
      fCtr << std::setw(4) << ix
	   << std::setw(4) << iy ;
      fCor << std::setw(4) << ix
	   << std::setw(4) << iy ;
   }
   if( det     == DetId::Ecal &&
       subdetn == EcalPreshower     )
   {
      const ESDetId esid ( did ) ;
      const int pl ( esid.plane() ) ;
      const int ix ( esid.six() ) ;
      const int iy ( esid.siy() ) ;
//      const int iz ( esid.zside() ) ;
      const int st ( esid.strip() ) ;
      fCtr << std::setw(4) << pl
	   << std::setw(4) << ix
	   << std::setw(4) << iy
	   << std::setw(4) << st ;
      fCor << std::setw(4) << pl
	   << std::setw(4) << ix
	   << std::setw(4) << iy
	   << std::setw(4) << st ;
   }
   if( det     == DetId::Hcal    ) 
   {
      const HcalDetId hcid ( did ) ;
      const int ie ( hcid.ieta() ) ;
      const int ip ( hcid.iphi() ) ;
      const int de ( hcid.depth() ) ;
//      const int iz ( hcid.zside() ) ;
      fCtr << std::setw(4) << ie
	   << std::setw(4) << ip
	   << std::setw(4) << de ;
   }

   fCtr << std::fixed << std::setw(12) << std::setprecision(4)
	<< cell.getPosition().x()
	<< std::fixed << std::setw(12) << std::setprecision(4)
	<< cell.getPosition().y()
	<< std::fixed << std::setw(12) << std::setprecision(4)
	<< cell.getPosition().z()
	<< std::endl ;

   const CaloCellGeometry::CornersVec& co ( cell.getCorners() ) ;

   for( unsigned int j ( 0 ) ; j < co.size() ; ++(++(++j)) )
   {
      fCor << std::fixed << std::setw(10) << std::setprecision(4)
	   << co[j].x()
	   << std::fixed << std::setw(10) << std::setprecision(4)
	   << co[j].y()
	   << std::fixed << std::setw(10) << std::setprecision(4)
	   << co[j].z() ;
   }
   fCor << std::endl ;
}

void 
CaloGeometryAnalyzer::build( const CaloGeometry& cg, 
			     DetId::Detector det, 
			     int subdetn, 
			     const char* name) 
{
   std::cout<<"Name now is "<<name<<std::endl ;

   const std::string fnameCtr  ( std::string( name ) + ".ctr" ) ;
   const std::string fnameCor  ( std::string( name ) + ".cor" ) ;
   const std::string fnameRoot ( std::string( name ) + ".C" ) ;
   std::fstream fCtr(fnameCtr.c_str() ,std::ios_base::out);
   std::fstream fCor(fnameCor.c_str() ,std::ios_base::out);
   std::fstream f   (fnameRoot.c_str(),std::ios_base::out);
   const CaloSubdetectorGeometry* geom=cg.getSubdetectorGeometry(det,subdetn);

   f << "{" << std::endl;
   f << "  TGeoManager* geoManager = new TGeoManager(\"ROOT\", \"" << name << "\");" << std::endl;
   f << "  TGeoMaterial* dummyMaterial = new TGeoMaterial(\"Vacuum\", 0,0,0); " << std::endl;
   f << "  TGeoMedium* dummyMedium =  new TGeoMedium(\"Vacuum\",1,dummyMaterial);" << std::endl;
   f << "  TGeoVolume* world=geoManager->MakeBox(\"world\",dummyMedium, 8000.0, 8000.0, 14000.0); " << std::endl;
   f << "  geoManager->SetTopVolume(world); " << std::endl;
   f << "  TGeoVolume* box; " << std::endl;
   int n=0;
   std::vector<DetId> ids=geom->getValidDetIds(det,subdetn);
   for (std::vector<DetId>::iterator i=ids.begin(); i!=ids.end(); i++) {
      n++;
      const CaloCellGeometry* cell=geom->getGeometry(*i);

      ctrcor( det,
	      subdetn,
	      *i,
	      *cell,
	      fCtr,
	      fCor  ) ;

      if (det == DetId::Ecal)
      {
	 if (subdetn == EcalBarrel)
	 {
	    f << "  // " << EBDetId(*i) << std::endl;
	    
	    f << "  // Checking getClosestCell for position " << dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(0.) << std::endl;
	    EBDetId closestCell = EBDetId(geom->getClosestCell(dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(0.))) ;
	    f << "  // Return position is " << closestCell << std::endl;
	    assert (closestCell == EBDetId(*i) );
	 }
	 if (subdetn == EcalEndcap)
	 {
	    const EEDetId did ( *i ) ;
	    const int ix ( did.ix() ) ;
	    const int iy ( did.iy() ) ;
	    const int iz ( did.zside() ) ;
	    const TruncatedPyramid* tp ( dynamic_cast<const TruncatedPyramid*>(cell) ) ;
	    f << "  // Checking getClosestCell for position " << tp->getPosition(0.) << std::endl;

	    const EEDetId closestCell ( geom->getClosestCell( tp->getPosition(0.) ) ) ;
	    f << "  // Return position is " << closestCell << std::endl;
	    assert( closestCell == did ) ;

	    const GlobalVector xx ( 2.5,   0, 0 ) ;
	    const GlobalVector yy (   0, 2.5, 0 ) ;
	    const GlobalVector zz (   0,   0, 1 ) ;
	    const GlobalPoint pointIn ( tp->getPosition(  1.) ) ; 
	    const GlobalPoint pointFr ( tp->getPosition( -1.) ) ; 
	    const GlobalPoint pointBk ( tp->getPosition( 25.) ) ; 
	    const GlobalPoint pointXP ( tp->getPosition(1.) + xx ) ; 
	    const GlobalPoint pointXM ( tp->getPosition(1.) - xx ) ; 
	    const GlobalPoint pointYP ( tp->getPosition(1.) + yy ) ; 
	    const GlobalPoint pointYM ( tp->getPosition(1.) - yy ) ; 
	    const GlobalPoint pointPP ( tp->getPosition(1.) + xx + yy ) ; 
	    const GlobalPoint pointPM ( tp->getPosition(1.) + xx - yy ) ; 
	    const GlobalPoint pointMP ( tp->getPosition(1.) - xx + yy ) ; 
	    const GlobalPoint pointMM ( tp->getPosition(1.) - xx - yy ) ; 
	    const EEDetId didXP ( gid( ix+1, iy  , iz, did ) ) ;
	    const EEDetId didXM ( gid( ix-1, iy  , iz, did ) ) ;
	    const EEDetId didYP ( gid( ix  , iy+1, iz, did ) ) ;
	    const EEDetId didYM ( gid( ix  , iy-1, iz, did ) ) ;
	    const EEDetId didPP ( gid( ix+1, iy+1, iz, did ) ) ;
	    const EEDetId didPM ( gid( ix+1, iy-1, iz, did ) ) ;
	    const EEDetId didMP ( gid( ix-1, iy+1, iz, did ) ) ;
	    const EEDetId didMM ( gid( ix-1, iy-1, iz, did ) ) ;

	    assert(  cell->inside( pointIn ) ) ;
	    assert( !cell->inside( pointFr ) ) ;
	    assert( !cell->inside( pointBk ) ) ;
	    assert( !cell->inside( pointXP ) ) ;
	    assert( !cell->inside( pointXM ) ) ;
	    assert( !cell->inside( pointYP ) ) ;
	    assert( !cell->inside( pointYM ) ) ;
	    assert( !cell->inside( pointPP ) ) ;
	    assert( !cell->inside( pointPM ) ) ;
	    assert( !cell->inside( pointMP ) ) ;
	    assert( !cell->inside( pointMM ) ) ;

	    const EEDetId ccBk ( geom->getClosestCell( pointBk ) ) ;
	    const EEDetId ccIn ( geom->getClosestCell( pointIn ) ) ;
	    const EEDetId ccFr ( geom->getClosestCell( pointFr ) ) ;
	    const EEDetId ccXP ( geom->getClosestCell( pointXP ) ) ;
	    const EEDetId ccXM ( geom->getClosestCell( pointXM ) ) ;
	    const EEDetId ccYP ( geom->getClosestCell( pointYP ) ) ;
	    const EEDetId ccYM ( geom->getClosestCell( pointYM ) ) ;
	    const EEDetId ccPP ( geom->getClosestCell( pointPP ) ) ;
	    const EEDetId ccPM ( geom->getClosestCell( pointPM ) ) ;
	    const EEDetId ccMP ( geom->getClosestCell( pointMP ) ) ;
	    const EEDetId ccMM ( geom->getClosestCell( pointMM ) ) ;

	    assert( ccIn == did ) ;
	    assert( ccFr == did ) ;
	    assert( ccBk == did ) ;
	    assert( ccXP == didXP ||
		    !geom->getGeometry(didXP)->inside( pointXP ) ) ;
	    assert( ccXM == didXM ||
		    !geom->getGeometry(didXM)->inside( pointXM ) ) ;
	    assert( ccYP == didYP  ||
		    !geom->getGeometry(didYP)->inside( pointYP ) ) ;
	    assert( ccYM == didYM  ||
		    !geom->getGeometry(didYM)->inside( pointYM )) ;
	    assert( ccPP == didPP  ||
		    !geom->getGeometry(didPP)->inside( pointPP ) ) ;
	    assert( ccPM == didPM ||
		    !geom->getGeometry(didPM)->inside( pointPM ) ) ;
	    assert( ccMP == didMP ||
		    !geom->getGeometry(didMP)->inside( pointMP ) ) ;
	    assert( ccMM == didMM ||
		    !geom->getGeometry(didMM)->inside( pointMM ) ) ;
	 }
	 if (subdetn == EcalPreshower) 
	 {
	    f << "  // " << ESDetId(*i) << std::endl;
	    f << "  // Checking getClosestCell for position " << cell->getPosition() << " in plane " << ESDetId(*i).plane() << std::endl;
	    ESDetId closestCell=ESDetId((dynamic_cast<const EcalPreshowerGeometry*>(geom))->getClosestCellInPlane(cell->getPosition(),ESDetId(*i).plane()));
	    f << "  // Return position is " << closestCell << std::endl;
	    //sanity checks
            int o_zside = ESDetId(*i).zside();
            //int o_plane = ESDetId(*i).plane();
            int o_six   = ESDetId(*i).six();
            int o_siy   = ESDetId(*i).siy();
            //int o_strip = ESDetId(*i).strip();

            assert ((o_six <= 20 && cell->getPosition().x() < 0.) || (o_six > 20 && cell->getPosition().x() > 0.));
            assert ((o_siy <= 20 && cell->getPosition().y() < 0.) || (o_siy > 20 && cell->getPosition().y() > 0.));
            assert ((o_zside < 0 && cell->getPosition().z() < 0.) || (o_zside > 0 && cell->getPosition().z() > 0.));
	    assert (closestCell == ESDetId(*i) );
	 }
      }
      else if (det == DetId::Hcal)
      {
	 f << "  // " << HcalDetId(*i) << std::endl;
      }
    
      if (det == DetId::Hcal && subdetn==HcalForward) 
	 f << "  box=geoManager->MakeBox(\"point\",dummyMedium,1.0,1.0,1.0);" << std::endl;
      else
	 f << "  box=geoManager->MakeBox(\"point\",dummyMedium,3.0,3.0,3.0);" << std::endl;
      f << "  world->AddNode(box,"<< n << ",new TGeoHMatrix(TGeoTranslation(" << 
	 cell->getPosition().x() << "," << cell->getPosition().y() << "," << cell->getPosition().z() << ")));" << std::endl;
      //   f << (HcalDetId)(*i) << " " << cell->getPosition() << std::endl;
   }
   f << "  geoManager->CloseGeometry();" << std::endl;
   f << "world->Voxelize(\"\"); // now the new geometry is valid for tracking, so you can do \n // even raytracing \n //  if (!canvas) { \n    TCanvas* canvas=new TCanvas(\"EvtDisp\",\"EvtDisp\",500,500); \n //  } \n  canvas->Modified(); \n  canvas->Update();      \n  world->Draw(); \n";
   f << "}" << std::endl;
   f.close();
   fCtr.close();
   fCor.close();
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
CaloGeometryAnalyzer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
   using namespace edm;

   edm::ESHandle<CaloGeometry> pG;
   iSetup.get<IdealGeometryRecord>().get(pG);     
   //
   // get the ecal & hcal geometry
   //
   if (pass_==0) {
     build(*pG,DetId::Ecal,EcalBarrel,"eb");
     build(*pG,DetId::Ecal,EcalEndcap,"ee");
     //Test eeGetClosestCell in Florian Point
     std::cout << "Checking getClosestCell for position" << GlobalPoint(-38.9692,-27.5548,-317) << std::endl;
     std::cout << "Position of Closest Cell in EE " << dynamic_cast<const TruncatedPyramid*>(pG->getGeometry(EEDetId((*pG).getSubdetectorGeometry(DetId::Ecal,EcalEndcap)->getClosestCell(GlobalPoint(-38.9692,-27.5548,-317)))))->getPosition(0.) << std::endl;
     build(*pG,DetId::Ecal,EcalPreshower,"es");
     build(*pG,DetId::Hcal,HcalBarrel,"hb");
     build(*pG,DetId::Hcal,HcalEndcap,"he");
     build(*pG,DetId::Hcal,HcalOuter,"ho");
     build(*pG,DetId::Hcal,HcalForward,"hf");
     
   }

   pass_++;
      
}

//define this as a plug-in
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(CaloGeometryAnalyzer);
