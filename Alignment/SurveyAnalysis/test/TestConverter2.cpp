// -*- C++ -*-
//
// Package:    TestConverter2
// Class:      TestConverter2
// 
//
// Description: Module to test the SurveyConverter software
//
//
// Original Author:  Roberto Covarelli
//         Created:  March 16, 2006
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

#include "Alignment/TrackerAlignment/interface/AlignableTracker.h"
#include "Alignment/TrackerAlignment/interface/AlignableTrackerBarrelLayer.h"
#include "Alignment/TrackerAlignment/interface/AlignableTrackerRod.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/GeometrySurface/interface/TkRotation.h"
#include "CondFormats/Alignment/interface/Alignments.h"
#include "CondFormats/Alignment/interface/AlignTransform.h"
#include "CondFormats/AlignmentRecord/interface/TrackerAlignmentRcd.h"
#include "CondFormats/Alignment/interface/AlignmentErrors.h"
#include "CondFormats/Alignment/interface/AlignTransformError.h"
#include "CondFormats/AlignmentRecord/interface/TrackerAlignmentErrorRcd.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"

#include "Alignment/SurveyAnalysis/interface/SurveyDataReader.h"
//
//
// class declaration
//
using namespace std;

class TestConverter2 : public edm::EDAnalyzer {

public:
  explicit TestConverter2( const edm::ParameterSet& );
  ~TestConverter2();
  
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
private:
  // ----------member data ---------------------------
  TTree* theTree;
  TFile* theFile;
  edm::ParameterSet theParameterSet;
  float dx_,dy_,dz_;
  float dtx_,dty_,dtz_;
  float dkx_,dky_,dkz_;
  float dnx_,dny_,dnz_;
  float errx_,erry_,errz_; 
  int subdid_, fwbw_, layerdisk_, frontback_, stringrod_, petal_, module_;
  // TRotMatrix* rot_;

};

//
// constructors and destructor
//
TestConverter2::TestConverter2( const edm::ParameterSet& iConfig ) :
  theParameterSet( iConfig )
{ 
  
  // Open root file and define tree
  std::string fileName = theParameterSet.getUntrackedParameter<std::string>("fileName","test.root");
  theFile = new TFile( fileName.c_str(), "RECREATE" );
  theTree = new TTree( "theTree", "Detector units positions" );
  
  theTree->Branch("subDetId",     &subdid_,     "subDetId/I"     );
  theTree->Branch("Fw_Bw",        &fwbw_,       "Fw_Bw/I"        );
  theTree->Branch("FrBa_InEx",    &frontback_,  "FrBa_inEx/I"    );
  theTree->Branch("LayerDisk",    &layerdisk_,  "LayerDisk/I"    );
  theTree->Branch("Petal",        &petal_,      "Petal/I"        );
  theTree->Branch("StRingRod",    &stringrod_,  "StRingRod/I"    );
  theTree->Branch("Module",       &module_,     "Module/I"       );
  theTree->Branch("dx",     &dx_,     "dx/F"     );
  theTree->Branch("dy",     &dy_,     "dy/F"     );
  theTree->Branch("dz",     &dz_,     "dz/F"     );
  theTree->Branch("dtx",    &dtx_,    "dtx/F"    );
  theTree->Branch("dty",    &dty_,    "dty/F"    );
  theTree->Branch("dtz",    &dtz_,    "dtz/F"    );
  theTree->Branch("dkx",    &dkx_,    "dkx/F"    );
  theTree->Branch("dky",    &dky_,    "dky/F"    );
  theTree->Branch("dkz",    &dkz_,    "dkz/F"    );
  theTree->Branch("dnx",    &dnx_,    "dnx/F"    );
  theTree->Branch("dny",    &dny_,    "dny/F"    );
  theTree->Branch("dnz",    &dnz_,    "dnz/F"    ); 
  theTree->Branch("errx",   &errx_,   "errx/F"   );
  theTree->Branch("erry",   &erry_,   "erry/F"   );
  theTree->Branch("errz",   &errz_,   "errz/F"   );
}


TestConverter2::~TestConverter2()
{ 
  
  theTree->Write();
  theFile->Close();
  
}


void
TestConverter2::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
   
  edm::LogInfo("TrackerAlignment") << "Starting!";

  //
  // Retrieve tracker geometry from event setup
  edm::ESHandle<TrackerGeometry> trackerGeometry;
  iSetup.get<TrackerDigiGeometryRecord>().get( trackerGeometry );

  // Retrieve alignment[Error]s from DBase
  edm::ESHandle<Alignments> alignments;
  iSetup.get<TrackerAlignmentRcd>().get( alignments );
  edm::ESHandle<AlignmentErrors> alignmentErrors;
  iSetup.get<TrackerAlignmentErrorRcd>().get( alignmentErrors );
  
  std::vector<AlignTransformError> alignErrors = alignmentErrors->m_alignError;

  // Now loop on detector units, and store difference position and orientation w.r.t. survey
  
   for ( std::vector<AlignTransform>::const_iterator iGeomDet = alignments->m_align.begin();
		iGeomDet != alignments->m_align.end(); iGeomDet++ )
	{
	  
	  for ( std::vector<GeomDet*>::const_iterator iDet = trackerGeometry->dets().begin();
		iDet != trackerGeometry->dets().end(); iDet++ )
	    {
              
	      // std::cout << (*iDet)->geographicalId().rawId() << " " << (*iGeomDet).rawId() << std::endl;
              if ((*iDet)->geographicalId().rawId() == (*iGeomDet).rawId()) {
                
		for ( std::vector<AlignTransformError>::const_iterator it = alignErrors.begin();
		      it != alignErrors.end(); it++ ) {
		  
		  if ((*it).rawId() == (*iGeomDet).rawId()) {
		  
		    DetId * thisId = new DetId( (*iGeomDet).rawId() );

		    if (thisId->subdetId() == int(StripSubdetector::TIB) ) {
		      TIBDetId * thisTIBid = new TIBDetId( *thisId );
		      subdid_ = 3;
		      layerdisk_ = thisTIBid->layer(); 
		      std::vector<unsigned int> theString = thisTIBid->string();
		      fwbw_ = theString[0];
		      frontback_ = theString[1];
		      stringrod_ = theString[2];
		      petal_ = 0;
		      module_ = thisTIBid->module();

		    } else if (thisId->subdetId() == int(StripSubdetector::TID) ) {
		      TIDDetId * thisTIDid = new TIDDetId( *thisId );
		      subdid_ = 4;
		      layerdisk_ = thisTIDid->wheel(); 
		      std::vector<unsigned int> theModule = thisTIDid->module();
		      frontback_ = theModule[0];
		      module_ = theModule[1];
		      petal_ = 0;
		      fwbw_ = thisTIDid->side();
		      stringrod_ = thisTIDid->ring();
		     
		    } else if (thisId->subdetId() == int(StripSubdetector::TOB) ) {
		      TOBDetId * thisTOBid = new TOBDetId( *thisId );
		      subdid_ = 5;
		      layerdisk_ = thisTOBid->layer(); 
		      std::vector<unsigned int> theRod = thisTOBid->rod();
		      fwbw_ = theRod[0];
		      stringrod_ = theRod[1];
		      petal_ = 0;
		      frontback_ = 0;
		      module_ = thisTOBid->module();
		      
		    } else if (thisId->subdetId() == int(StripSubdetector::TEC) ) {
		      TECDetId * thisTECid = new TECDetId( *thisId );
		      subdid_ = 6;
		      layerdisk_ = thisTECid->wheel(); 
		      std::vector<unsigned int> thePetal = thisTECid->petal();
		      frontback_ = thePetal[0];
		      petal_ = thePetal[1];
		      fwbw_ = thisTECid->side();
		      stringrod_ = thisTECid->ring();
		      module_ = thisTECid->module();
		         
		    } else {
		      std::cout << "WARNING!!! this DetId (" << thisId->rawId() << ") does not belong to SiStrip tracker" << std::endl;
		      break;
		    }

		    // HepRotation fromAngles( (*iGeomDet).eulerAngles()  );
                    HepRotation fromAngles( (*iGeomDet).rotation()  );
		    Surface::RotationType rotation( fromAngles.xx(), fromAngles.xy(), fromAngles.xz(),
						    fromAngles.yx(), fromAngles.yy(), fromAngles.yz(),
						    fromAngles.zx(), fromAngles.zy(), fromAngles.zz() );
		    
		    dx_      = (*iGeomDet).translation().x() - (*iDet)->position().x(); 
		    dy_      = (*iGeomDet).translation().y() - (*iDet)->position().y();
		    dz_      = (*iGeomDet).translation().z() - (*iDet)->position().z();
		    dtx_     = rotation.xx() - (*iDet)->rotation().xx();
		    dty_     = rotation.xy() - (*iDet)->rotation().xy();
		    dtz_     = rotation.xz() - (*iDet)->rotation().xz();
		    dkx_     = rotation.yx() - (*iDet)->rotation().yx();
		    dky_     = rotation.yy() - (*iDet)->rotation().yy();
		    dkz_     = rotation.yz() - (*iDet)->rotation().yz();
		    dnx_     = rotation.zx() - (*iDet)->rotation().zx();
		    dny_     = rotation.zy() - (*iDet)->rotation().zy();
		    dnz_     = rotation.zz() - (*iDet)->rotation().zz();
                    HepSymMatrix errMat = (*it).matrix();
                    errx_    = sqrt(errMat[0][0]); 
		    erry_    = sqrt(errMat[1][1]);
		    errz_    = sqrt(errMat[2][2]); 
		    
		    theTree->Fill();
		    
		    cout << "DetId = " << (*iGeomDet).rawId() << " " << endl;
		    cout << "DetId decodified = " << subdid_ << " " << fwbw_ << " " << frontback_ << " " << layerdisk_ << " " << stringrod_ << " " << petal_ << " " << module_ << endl;
		    cout << "X pos : TRACKER_MOVED = " << std::fixed << std::setprecision(2) << (*iGeomDet).translation().x() << " / TRACKER_ORIGINAL = " << (*iDet)->position().x() << endl;
		    cout << "Y pos : TRACKER_MOVED = " << std::fixed << std::setprecision(2) << (*iGeomDet).translation().y() << " / TRACKER_ORIGINAL = " << (*iDet)->position().y() << endl;
		    cout << "Z pos : TRACKER_MOVED = " << std::fixed << std::setprecision(2) << (*iGeomDet).translation().z() << " / TRACKER_ORIGINAL = " << (*iDet)->position().z() << endl;
		    cout << "SPATIAL DISTANCE = " << std::fixed << std::setprecision(3) << sqrt(pow(dx_,2)+pow(dy_,2)+pow(dz_,2)) << endl;
		    cout << "Trans vect X : TRACKER_MOVED = " << std::fixed << std::setprecision(5) << rotation.xx() << " / TRACKER_ORIGINAL = " << (*iDet)->rotation().xx() << endl;
		    cout << "Trans vect Y : TRACKER_MOVED = " << std::fixed << std::setprecision(5) << rotation.xy() << " / TRACKER_ORIGINAL = " << (*iDet)->rotation().xy() << endl;
		    cout << "Trans vect Z : TRACKER_MOVED = " << std::fixed << std::setprecision(5) << rotation.xz() << " / TRACKER_ORIGINAL = " << (*iDet)->rotation().xz() << endl; 
		    cout << "Long vect X : TRACKER_MOVED = " << std::fixed << std::setprecision(5) << rotation.yx() << " / TRACKER_ORIGINAL = " << (*iDet)->rotation().yx() << endl;	
		    cout << "Long vect Y : TRACKER_MOVED = " << std::fixed << std::setprecision(5) << rotation.yy() << " / TRACKER_ORIGINAL = " << (*iDet)->rotation().yy() << endl;  
		    cout << "Long vect Z : TRACKER_MOVED = " << std::fixed << std::setprecision(5) << rotation.yz() << " / TRACKER_ORIGINAL = " << (*iDet)->rotation().yz() << endl;
		    cout << "Norm vect X : TRACKER_MOVED = " << std::fixed << std::setprecision(5) << rotation.zx() << " / TRACKER_ORIGINAL = " << (*iDet)->rotation().zx() << endl;
		    cout << "Norm vect Y : TRACKER_MOVED = " << std::fixed << std::setprecision(5) << rotation.zy() << " / TRACKER_ORIGINAL = " << (*iDet)->rotation().zy() << endl; 
		    cout << "Norm vect Z : TRACKER_MOVED = " << std::fixed << std::setprecision(5) << rotation.zz() << " / TRACKER_ORIGINAL = " << (*iDet)->rotation().zz() << endl; 
		  }
		}
	      }
	    }
	}	   
   
  edm::LogInfo("TrackerAlignment") << "Done!";

}

//define this as a plug-in
DEFINE_FWK_MODULE(TestConverter2);
