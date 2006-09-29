// -*- C++ -*-
//
// Package:    TestTranslation
// Class:      TestTranslation
// 
//
// Description: Module to test the Alignment software
//
//


// system include files
#include <string>
#include <TTree.h>
#include <TFile.h>
#include <TRotMatrix.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Alignment/MuonAlignment/interface/AlignableMuon.h"
#include "Alignment/MuonAlignment/interface/AlignableDTChamber.h"
#include "Alignment/MuonAlignment/interface/AlignableCSCChamber.h"

#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTChamber.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"

#include "Geometry/Surface/interface/Surface.h"

#include <vector>

//
//
// class declaration
//

class TestTranslation : public edm::EDAnalyzer {
public:
  explicit TestTranslation( const edm::ParameterSet& );
  ~TestTranslation();
  
  
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
private:
  // ----------member data ---------------------------
  TTree* theTree;
  TFile* theFile;
  float x,y,z,phi,theta,length,thick,width;
  TRotMatrix* dir;

  typedef Surface::RotationType    RotationType;
  typedef Surface::PositionType    PositionType;


  void apply( Alignable* );

};

//
// constructors and destructor
//
TestTranslation::TestTranslation( const edm::ParameterSet& iConfig ) 
{ 

  // Open root file and define tree
  std::string fileName = iConfig.getUntrackedParameter<std::string>("fileName","test.root");
  theFile = new TFile( fileName.c_str(), "RECREATE" );
  theTree = new TTree( "theTree", "Detector units positions" );
  
  theTree->Branch("x",      &x,      "x/F"      );
  theTree->Branch("y",      &y,      "y/F"      );
  theTree->Branch("z",      &z,      "z/F"      );
  theTree->Branch("phi",    &phi,    "phi/F"    );
  theTree->Branch("theta",  &theta,  "theta/F"  );
  theTree->Branch("length", &length, "length/F" );
  theTree->Branch("width",  &width,  "width/F"  );
  theTree->Branch("thick",  &thick,  "thick/F"  );
  dir = 0;
  theTree->Branch("dir",    "TRotMatrix", &dir  );

}


TestTranslation::~TestTranslation()
{ 
  
  theTree->Write();
  theFile->Close();
  
}


void
TestTranslation::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

   
  edm::LogInfo("MuonAlignment") << "Starting!";

/*
  //
  // Build alignable muon geometry from event setup
  //
  edm::ESHandle<DDCompactView> cpv;
  iRecord.getRecord<IdealGeometryRecord>().get( cpv );

  DTGeometryBuilderFromDDD  DTGeometryBuilder;
  CSCGeometryBuilderFromDDD CSCGeometryBuilder;

  theDTGeometry   = boost::shared_ptr<DTGeometry>( DTGeometryBuilder.build( &(*cpv) ) );
  theCSCGeometry  = boost::shared_ptr<CSCGeometry>( CSCGeometryBuilder.build( &(*cpv) ) );

  AlignableMuon* theAlignableMuon = new AlignableMuon( &(*theDTGeometry) , &(*theCSCGeometry) );
*/

  //
  // Retrieve ideal geometry from event setup
  //
  edm::ESHandle<DTGeometry> dtGeometry;
  edm::ESHandle<CSCGeometry> cscGeometry;
  iSetup.get<MuonGeometryRecord>().get( dtGeometry );
  iSetup.get<MuonGeometryRecord>().get( cscGeometry );

  AlignableMuon* theAlignableMuon = new AlignableMuon( &(*dtGeometry), &(*cscGeometry) );



  // Apply  alignment

  std::vector<Alignable*> theDTAlignables = theAlignableMuon->DTChambers();
  std::vector<Alignable*> theCSCAlignables = theAlignableMuon->CSCEndcaps();


  for ( std::vector<Alignable*>::iterator iter = theDTAlignables.begin();
		                          iter != theDTAlignables.end(); iter++ ){ 
    apply( *iter ); 
  }

  theDTAlignables.clear();

  for ( std::vector<Alignable*>::iterator iter = theCSCAlignables.begin();
		                          iter != theCSCAlignables.end(); iter++ ){ 
    apply( *iter ); 
  }

  theCSCAlignables.clear();


edm::LogInfo("MuonAlignment") << "Done!";

}



void TestTranslation::apply( Alignable* it )
{
	  std::cout << "------------------------" << std::endl
		    << " BEFORE TRANSLATION " << std::endl;

	  GlobalPoint  pos_i  = (it)->globalPosition() ;
          RotationType dir_i  = (it)->globalRotation();

	  std::cout << "x=" << pos_i.x() << ",  y=" << pos_i.y() << ",  z=" << pos_i.z() << std::endl; 

	  float dx = 1.0;
          float dy = 2.0;
          float dz = 3.0;
          GlobalVector dr( dx, dy, dz );
	  it->move( dr );

	  std::cout << "------------------------" << std::endl
		    << " AFTER TRANSLATION " << std::endl;

          GlobalPoint  pos_f  = (it)->globalPosition() ;
          RotationType dir_f = (it)->globalRotation();

          std::cout << "x=" << pos_f.x() << ",  y=" << pos_f.y() << ",  z=" << pos_f.z()  << std::endl ;

	  std::cout << "------------------------" << std::endl;

}
  

//  delete theAlignableMuon ;


//define this as a plug-in
DEFINE_FWK_MODULE(TestTranslation)
