// -*- C++ -*-
//
// Package:    TestConverter
// Class:      TestConverter
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
#include "Geometry/Surface/interface/TkRotation.h"
#include "CondFormats/Alignment/interface/Alignments.h"
#include "CondFormats/Alignment/interface/AlignTransform.h"
#include "CondFormats/DataRecord/interface/TrackerAlignmentRcd.h"
#include "CondFormats/Alignment/interface/AlignmentErrors.h"
#include "CondFormats/Alignment/interface/AlignTransformError.h"
#include "CondFormats/DataRecord/interface/TrackerAlignmentErrorRcd.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"

#include "Alignment/SurveyAnalysis/interface/SurveyDataReader.h"
//
//
// class declaration
//
using namespace std;
static const int NFILES = 2;

class TestConverter : public edm::EDAnalyzer {

  typedef unsigned int DetIdType;
  typedef std::map< DetIdType, std::vector<float> > MapType;
  typedef std::pair< DetIdType, std::vector<float> > PairType;
  typedef std::map< std::vector<int>, std::vector<float> > MapTypeOr;
  typedef std::pair< std::vector<int>, std::vector<float> > PairTypeOr;

public:
  explicit TestConverter( const edm::ParameterSet& );
  ~TestConverter();
  
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
  int Id_;
  // TRotMatrix* rot_;

};

//
// constructors and destructor
//
TestConverter::TestConverter( const edm::ParameterSet& iConfig ) :
  theParameterSet( iConfig )
{ 
  
  // Open root file and define tree
  std::string fileName = theParameterSet.getUntrackedParameter<std::string>("fileName","test.root");
  theFile = new TFile( fileName.c_str(), "RECREATE" );
  theTree = new TTree( "theTree", "Detector units positions" );
  
  theTree->Branch("Id",     &Id_,     "Id/I"     );
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


TestConverter::~TestConverter()
{ 
  
  theTree->Write();
  theFile->Close();
  
}


void
TestConverter::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
   
  edm::LogInfo("TrackerAlignment") << "Starting!";

  //
  // Read in the survey information from the text files
  // 
  edm::ParameterSet textFiles = theParameterSet.getParameter<edm::ParameterSet>( "textFileNames" );
  std::string textFileNames[NFILES]; 
  std::string fileType[NFILES];    
  textFileNames[0] = textFiles.getUntrackedParameter<std::string>("forTIB","NONE");  
  fileType[0] = "TIB";
  textFileNames[1] = textFiles.getUntrackedParameter<std::string>("forTID","NONE");
  fileType[1] = "TID";

  SurveyDataReader dataReader;
  for (int ii=0 ; ii<NFILES ;ii++) {
    if ( textFileNames[ii] == "NONE" )
      throw cms::Exception("BadConfig") << fileType[ii] << " input file not found in configuration";
    dataReader.readFile( textFileNames[ii], fileType[ii] );
  } 

  edm::LogInfo("TrackerAlignment") << "Files read";

  const MapTypeOr theSurveyMap = dataReader.surveyMap();

  edm::LogInfo("TrackerAlignment") << "Map written";

  //
  // Retrieve tracker geometry from event setup
  //
  // edm::ESHandle<TrackerGeometry> trackerGeometry;
  // iSetup.get<TrackerDigiGeometryRecord>().get( trackerGeometry );

  // Retrieve alignment[Error]s from DBase
  edm::ESHandle<Alignments> alignments;
  iSetup.get<TrackerAlignmentRcd>().get( alignments );
  edm::ESHandle<AlignmentErrors> alignmentErrors;
  iSetup.get<TrackerAlignmentErrorRcd>().get( alignmentErrors );
  
  std::vector<AlignTransformError> alignErrors = alignmentErrors->m_alignError;
  int countDet = 0;

  // Now loop on detector units, and store difference position and orientation w.r.t. survey
  
   for ( std::vector<AlignTransform>::const_iterator iGeomDet = alignments->m_align.begin();
		iGeomDet != alignments->m_align.end(); iGeomDet++ )
	{
          
          if (countDet == 0) countDet = countDet + 3;
	  unsigned int comparisonVect[6] = {0,0,0,0,0,0};
	  
	  DetId * thisId = new DetId( (*iGeomDet).rawId() );
	  if (thisId->subdetId() == int(StripSubdetector::TIB)) {
	    
	    comparisonVect[0] = int(StripSubdetector::TIB);
	    TIBDetId * thisTIBid = new TIBDetId( *thisId );
	    comparisonVect[1] = thisTIBid->layer();  
            if (comparisonVect[1] < 3) {countDet--;} else {countDet = countDet - 3;}
	    std::vector<unsigned int> theString = thisTIBid->string();
	    comparisonVect[2] = theString[0];
	    comparisonVect[3] = theString[1];
	    comparisonVect[4] = theString[2];
	    comparisonVect[5] = thisTIBid->module();
	    
	  } else if (thisId->subdetId() == int(StripSubdetector::TID)) {
	    
	    comparisonVect[0] = int(StripSubdetector::TID);
	    TIDDetId * thisTIDid = new TIDDetId( *thisId );
	    comparisonVect[1] = thisTIDid->side();
	    comparisonVect[2] = thisTIDid->wheel();
	    comparisonVect[3] = thisTIDid->ring(); 
            if (comparisonVect[3] < 3) {countDet--;} else {countDet = countDet - 3;}
	    std::vector<unsigned int> theModule = thisTIDid->module();
	    comparisonVect[4] = theModule[0];
	    comparisonVect[5] = theModule[1];
	    
	  }
	  
	  if (countDet == 0) { // Store only r-phi for double-sided modules

	    for ( MapTypeOr::const_iterator it = theSurveyMap.begin(); it != theSurveyMap.end(); it++ ) {
	      std::vector<int> locPos = (it)->first;
	      std::vector<float> align_params = (it)->second;
	      
	      if (locPos[0] == int(comparisonVect[0]) &&
		  locPos[1] == int(comparisonVect[1]) &&
		  locPos[2] == int(comparisonVect[2]) &&
		  locPos[3] == int(comparisonVect[3]) &&
		  locPos[4] == int(comparisonVect[4]) &&
		  locPos[5] == int(comparisonVect[5]) ) {
		
                
		for ( std::vector<AlignTransformError>::const_iterator it = alignErrors.begin();
		      it != alignErrors.end(); it++ ) {
		  
		  if ((*it).rawId() == (*iGeomDet).rawId()) {
		  
		    const CLHEP::HepRotation& rot = (*iGeomDet).rotation();
		    align::RotationType rotation( rot.xx(), rot.xy(), rot.xz(),
						  rot.yx(), rot.yy(), rot.yz(),
						  rot.zx(), rot.zy(), rot.zz() );
		    
		    Id_     = (*iGeomDet).rawId();    
		    dx_      = (*iGeomDet).translation().x() - align_params[15]; 
		    dy_      = (*iGeomDet).translation().y() - align_params[16];
		    dz_      = (*iGeomDet).translation().z() - align_params[17];
		    dtx_     = rotation.xx() - align_params[21];
		    dty_     = rotation.xy() - align_params[22];
		    dtz_     = rotation.xz() - align_params[23];
		    dkx_     = rotation.yx() - align_params[24];
		    dky_     = rotation.yy() - align_params[25];
		    dkz_     = rotation.yz() - align_params[26];
		    dnx_     = rotation.zx() - align_params[18];
		    dny_     = rotation.zy() - align_params[19];
		    dnz_     = rotation.zz() - align_params[20];
                    HepSymMatrix errMat = (*it).matrix();
                    errx_    = sqrt(errMat[0][0]); 
		    erry_    = sqrt(errMat[1][1]);
		    errz_    = sqrt(errMat[2][2]); 
		    
		    theTree->Fill();
		    
		    if (dkx_ > 0.04) {
		      cout << "DetId = " << Id_ << " " << endl;
		      cout << "DetId decodified = " << comparisonVect[0] << " " << comparisonVect[1] << " " << comparisonVect[2] << " " << comparisonVect[3] << " " << comparisonVect[4] << " " << comparisonVect[5] << endl;
		      cout << "X pos : TRACKER_MOVED = " << std::fixed << std::setprecision(2) << (*iGeomDet).translation().x() << " / SURVEY RICCARDO = " << align_params[15] << endl;
		      cout << "Y pos : TRACKER_MOVED = " << std::fixed << std::setprecision(2) << (*iGeomDet).translation().y() << " / SURVEY RICCARDO = " << align_params[16] << endl;
		      cout << "Z pos : TRACKER_MOVED = " << std::fixed << std::setprecision(2) << (*iGeomDet).translation().z() << " / SURVEY RICCARDO = " << align_params[17] << endl;
		      cout << "SPATIAL DISTANCE = " << std::fixed << std::setprecision(3) << sqrt(pow(dx_,2)+pow(dy_,2)+pow(dz_,2)) << endl;
		      cout << "Trans vect X : TRACKER_MOVED = " << std::fixed << std::setprecision(3) << rotation.xx() << " / SURVEY RICCARDO = " << align_params[21] << endl;
		      cout << "Trans vect Y : TRACKER_MOVED = " << std::fixed << std::setprecision(3) << rotation.xy() << " / SURVEY RICCARDO = " << align_params[22] << endl;
		      cout << "Trans vect Z : TRACKER_MOVED = " << std::fixed << std::setprecision(3) << rotation.xz() << " / SURVEY RICCARDO = " << align_params[23] << endl; 
		      cout << "Long vect X : TRACKER_MOVED = " << std::fixed << std::setprecision(3) << rotation.yx() << " / SURVEY RICCARDO = " << align_params[24] << endl;	
		      cout << "Long vect Y : TRACKER_MOVED = " << std::fixed << std::setprecision(3) << rotation.yy() << " / SURVEY RICCARDO = " << align_params[25] << endl;  
		      cout << "Long vect Z : TRACKER_MOVED = " << std::fixed << std::setprecision(3) << rotation.yz() << " / SURVEY RICCARDO = " << align_params[26] << endl;
		      cout << "Norm vect X : TRACKER_MOVED = " << std::fixed << std::setprecision(3) << rotation.zx() << " / SURVEY RICCARDO = " << align_params[18] << endl;
		      cout << "Norm vect Y : TRACKER_MOVED = " << std::fixed << std::setprecision(3) << rotation.zy() << " / SURVEY RICCARDO = " << align_params[19] << endl; 
		      cout << "Norm vect Z : TRACKER_MOVED = " << std::fixed << std::setprecision(3) << rotation.zz() << " / SURVEY RICCARDO = " << align_params[20] << endl; 
		    }
		  }
		}
	      }
	    }	  
	  } 
	}
  edm::LogInfo("TrackerAlignment") << "Done!";

}

//define this as a plug-in
DEFINE_FWK_MODULE(TestConverter);
