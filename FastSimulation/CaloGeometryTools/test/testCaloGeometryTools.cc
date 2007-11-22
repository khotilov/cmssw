// -*- C++ -*-
//
// Package:    testCaloCaloGeometryTools
// Class:      testCaloCaloGeometryTools
// 
/**\class testCaloCaloGeometryTools testCaloGeometryTools.cc test/testCaloGeometryTools/src/testCaloGeometryTools.cc

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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/EcalBarrelAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalEndcapAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"

#include "FastSimulation/CaloGeometryTools/interface/CaloGeometryHelper.h"
#include "FastSimulation/CaloHitMakers/interface/EcalHitMaker.h"
#include "FastSimulation/CaloGeometryTools/interface/Crystal.h"
#include "FastSimulation/Utilities/interface/Histos.h"
#include "FastSimulation/Utilities/interface/RandomEngine.h"

#include <TCanvas.h>
#include <TVirtualPad.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TText.h>
#include <TFile.h>
#include <TArrow.h>
#include <TBox.h>
#include <TPolyLine3D.h>
#include <TMarker.h>
#include <iostream>
#include <sstream>
//
// class decleration
//

typedef math::XYZVector XYZVector;
typedef math::XYZVector XYZPoint;

class testCaloGeometryTools : public edm::EDAnalyzer {
public:
  explicit testCaloGeometryTools( const edm::ParameterSet& );
  ~testCaloGeometryTools();
  
  
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
private:
  // ----------member data ---------------------------
  void testpoint(const XYZPoint& , std::string name, bool barrel);
  void checkSM();
  void checkSC();
  void testBorderCrossing();
  int pass_;

  Histos * myHistos;
  CaloGeometryHelper myGeometry;

  const RandomEngine* random;

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
testCaloGeometryTools::testCaloGeometryTools( const edm::ParameterSet& iConfig )
{
  myHistos = Histos::instance();
  myHistos->book("h100",150,0.,1.5,100,0.,35.);
}


testCaloGeometryTools::~testCaloGeometryTools()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  std::cout << " Writing the histo " << std::endl;
  myHistos->put("Grid.root");
  std::cout << " done " << std::endl;
}

// ------------ method called to produce the data  ------------
void
testCaloGeometryTools::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
   using namespace edm;
   
  // Initialize the random number generator service
  edm::Service<edm::RandomNumberGenerator> rng;
  if ( ! rng.isAvailable() ) {
    throw cms::Exception("Configuration")
      << "prod requires the RandomGeneratorService\n"
         "which is not present in the configuration file.\n"
         "You must add the service in the configuration file\n"
         "or remove the module that requires it";
  }
  random = new RandomEngine(&(*rng));

   edm::ESHandle<CaloTopology> theCaloTopology;
   iSetup.get<CaloTopologyRecord>().get(theCaloTopology);     

   edm::ESHandle<CaloGeometry> pG;
   iSetup.get<IdealGeometryRecord>().get(pG);     



   // Setup the tools
   double bField000 = 4.;
   myGeometry.setupGeometry(*pG);
   myGeometry.setupTopology(*theCaloTopology);
   myGeometry.initialize(bField000);
   
   // Take a point in the barrel
   XYZPoint p1(129,0.,-50);
   testpoint(p1,"barrel",true);
   XYZPoint p2(60,60,-317);
   testpoint(p1,"endcap",false);

   checkSM();
   checkSC();
   testBorderCrossing();
}

void testCaloGeometryTools::checkSM()
{
  std::vector<DetId> vec(myGeometry.getEcalBarrelGeometry()->getValidDetIds(DetId::Ecal,EcalBarrel));
  unsigned size=vec.size();
  for(unsigned ic=0;ic<size;++ic)
    {
      const CaloCellGeometry * geom=myGeometry.getEcalBarrelGeometry()->getGeometry(vec[ic]);
      GlobalPoint p=geom->getPosition();
      XYZPoint pp(p.x(),p.y(),p.z());
       // Build the name of the object
      std::ostringstream oss,oss2;
      oss << "iM"<< ic;
      oss2 <<"iSM" << ic;
      TMarker * myMarker= new TMarker(pp.eta(),pp.phi(),1);
      myMarker->SetMarkerColor(EBDetId(vec[ic]).im());
      TMarker * myMarker2= new TMarker(pp.eta(),pp.phi(),1);
      myMarker2->SetMarkerColor(EBDetId(vec[ic]).ism());
      myHistos->addObject(oss.str(),myMarker);
      myHistos->addObject(oss2.str(),myMarker2);
    }
}

void testCaloGeometryTools::checkSC()
{
  std::vector<DetId> vec(myGeometry.getEcalEndcapGeometry()->getValidDetIds(DetId::Ecal,EcalEndcap));
  unsigned size=vec.size();
  for(unsigned ic=0;ic<size;++ic)
    {
      const CaloCellGeometry * geom=myGeometry.getEcalEndcapGeometry()->getGeometry(vec[ic]);
      GlobalPoint p=geom->getPosition();
      XYZPoint pp(p.x(),p.y(),p.z());
       // Build the name of the object
      std::ostringstream oss,oss2;
      if(p.z()>0)
	oss << "iSCP"<< ic;
      else
	oss << "iSCN" << ic ;
      TMarker * myMarker= new TMarker(pp.x(),pp.y(),1);
      if(pp.perp2()<100.)
	std::cout << EEDetId(vec[ic]) << " " << pp.x() << " " << pp.y() << std::endl;
      myMarker->SetMarkerColor(EEDetId(vec[ic]).isc()%100);
      myMarker->SetMarkerStyle(22);
      myHistos->addObject(oss.str(),myMarker);
    }
}



void testCaloGeometryTools::testpoint(const XYZPoint& point, std::string name, bool barrel)
{
   DetId myCell = myGeometry.getClosestCell(point,true,barrel);
   EcalHitMaker myGrid(&myGeometry,point,myCell,1,7,0,random);
   
   std::vector<Crystal> myCrystals=myGrid.getCrystals();
   
   float xmin,ymin,zmin,xmax,ymax,zmax;
   xmin=ymin=zmin=99999;
   xmax=ymax=zmax=-9999;
   unsigned nxtal = myCrystals.size();

   std::vector<float> xp,yp,zp;
   
   for(unsigned ic=0;ic<nxtal;++ic)
     {
       XYZPoint p= myCrystals[ic].getCenter();

       myCrystals[ic].getDrawingCoordinates(xp,yp,zp);
       TPolyLine3D * myxtal= new TPolyLine3D(xp.size(),&xp[0],&yp[0],&zp[0]);
       

       // Build the name of the object
       std::ostringstream oss;
       oss << name << ic;

       myHistos->addObject(oss.str(),myxtal);

       if(xmin > p.x()) xmin=p.x();
       if(ymin > p.y()) ymin=p.y();
       if(zmin > p.z()) zmin=p.z();
       if(xmax < p.x()) xmax=p.x();
       if(ymax < p.y()) ymax=p.y();
       if(zmax < p.z()) zmax=p.z();
     }
   TH3F * frame = new TH3F(std::string(name+"frame").c_str(),"",100,xmin*0.9,xmax*1.1,100,ymin*0.9,ymax*1.1,100,zmin*0.9,zmax*1.1);
   myHistos->addObject("frame"+name,frame);
  
}

void testCaloGeometryTools::testBorderCrossing()
{
  // Barrel 
  std::vector<DetId> vec(myGeometry.getEcalBarrelGeometry()->getValidDetIds(DetId::Ecal,EcalBarrel));
  unsigned size=vec.size();
  unsigned counter=0;
  for(unsigned ic=0;ic<size;++ic)
    {
      std::vector<DetId> neighbours=myGeometry.getNeighbours(vec[ic]);
      for(unsigned in=0;in<8;++in)
	{
	  if(neighbours[in].null()) continue;
	  if(myGeometry.borderCrossing(vec[ic],neighbours[in]))
	    {
	      const CaloCellGeometry * geom=myGeometry.getEcalBarrelGeometry()->getGeometry(vec[ic]);
	      GlobalPoint p1=geom->getPosition();
	      XYZPoint pp1(p1.x(),p1.y(),p1.z());
	      geom=myGeometry.getEcalBarrelGeometry()->getGeometry(neighbours[in]);
	      GlobalPoint p2=geom->getPosition();
	      XYZPoint pp2(p2.x(),p2.y(),p2.z());
	      TMarker * myMarker= new TMarker((pp1+pp2).eta()*0.5,((pp1+pp2)*0.5).phi(),22);
	      std::ostringstream oss;
	      oss << "iBCB"<< counter;
	      myHistos->addObject(oss.str(),myMarker);
	      ++counter;
	    }
	}
    }
  
  // Endcap 
  vec=myGeometry.getEcalEndcapGeometry()->getValidDetIds(DetId::Ecal,EcalEndcap);
  size=vec.size();
  counter=0;
  for(unsigned ic=0;ic<size;++ic)
    {
      std::vector<DetId> neighbours=myGeometry.getNeighbours(vec[ic]);
      for(unsigned in=0;in<8;++in)
	{
	  if(neighbours[in].null()) continue;
	  if(myGeometry.borderCrossing(vec[ic],neighbours[in]))
	    {
	      const CaloCellGeometry * geom=myGeometry.getEcalEndcapGeometry()->getGeometry(vec[ic]);
	      GlobalPoint p1=geom->getPosition();
	      XYZPoint pp1(p1.x(),p1.y(),p1.z());
	      geom=myGeometry.getEcalEndcapGeometry()->getGeometry(neighbours[in]);
	      GlobalPoint p2=geom->getPosition();
	      XYZPoint pp2(p2.x(),p2.y(),p2.z());
	      TMarker * myMarker= new TMarker((pp1+pp2).x()*0.5,(pp1+pp2).y()*0.5,22);
	      std::ostringstream oss;
	      if(p1.z()>0)
		oss << "iBCEN"<< counter;
	      else
		oss << "iBCEP" << counter; 
	      
	      myHistos->addObject(oss.str(),myMarker);
	      ++counter;
	    }
	}
    }

}


//define this as a plug-in
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(testCaloGeometryTools);
