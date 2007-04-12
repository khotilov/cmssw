/** \file LaserAlignment.cc
 *  LAS reconstruction module
 *
 *  $Date: 2007/04/05 13:20:11 $
 *  $Revision: 1.8 $
 *  \author Maarten Thomas
 */

#include "Alignment/LaserAlignment/plugins/LaserAlignment.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "DataFormats/GeometrySurface/interface/BoundSurface.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeomBuilderFromGeometricDet.h"
#include "Geometry/TrackingGeometryAligner/interface/GeometryAligner.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"

#include "DataFormats/LaserAlignment/interface/LASBeamProfileFit.h"
#include "DataFormats/LaserAlignment/interface/LASBeamProfileFitCollection.h"
#include "DataFormats/LaserAlignment/interface/LASAlignmentParameter.h"
#include "DataFormats/LaserAlignment/interface/LASAlignmentParameterCollection.h"

// Conditions database
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"


LaserAlignment::LaserAlignment(edm::ParameterSet const& theConf) 
  : theEvents(0), 
    theStoreToDB(theConf.getUntrackedParameter<bool>("saveToDbase", false)),
    theSaveHistograms(theConf.getUntrackedParameter<bool>("saveHistograms",false)),
    theDebugLevel(theConf.getUntrackedParameter<int>("DebugLevel",0)),
    theNEventsPerLaserIntensity(theConf.getUntrackedParameter<int>("NumberOfEventsPerLaserIntensity",100)),
    theNEventsForAllIntensities(theConf.getUntrackedParameter<int>("NumberOfEventsForAllIntensities",100)),
    theDoAlignmentAfterNEvents(theConf.getUntrackedParameter<int>("DoAlignmentAfterNEvents",1000)),
    theAlignPosTEC(theConf.getUntrackedParameter<bool>("AlignPosTEC",false)),
    theAlignNegTEC(theConf.getUntrackedParameter<bool>("AlignNegTEC",false)), 
    theAlignTEC2TEC(theConf.getUntrackedParameter<bool>("AlignTECTIBTOBTEC",false)),
    theIsGoodFit(false),
    theSearchPhiTIB(theConf.getUntrackedParameter<double>("SearchWindowPhiTIB",0.05)),
    theSearchPhiTOB(theConf.getUntrackedParameter<double>("SearchWindowPhiTOB",0.05)),
    theSearchPhiTEC(theConf.getUntrackedParameter<double>("SearchWindowPhiTEC",0.05)),
    theSearchZTIB(theConf.getUntrackedParameter<double>("SearchWindowZTIB",1.0)),
    theSearchZTOB(theConf.getUntrackedParameter<double>("SearchWindowZTOB",1.0)),
    thePhiErrorScalingFactor(theConf.getUntrackedParameter<double>("PhiErrorScalingFactor",1.0)),
    theDigiProducersList(theConf.getParameter<Parameters>("DigiProducersList")),
    theFile(),
    theCompression(theConf.getUntrackedParameter<int>("ROOTFileCompression",1)),
    theFileName(theConf.getUntrackedParameter<std::string>("ROOTFileName","test.root")),
    theBeamFitPS(theConf.getParameter<edm::ParameterSet>("BeamProfileFitter")),
    theAlignmentAlgorithmPS(theConf.getParameter<edm::ParameterSet>("AlignmentAlgorithm")),
    theMinAdcCounts(theConf.getUntrackedParameter<int>("MinAdcCounts",0)),
    theHistogramNames(), theHistograms(),
    theLaserPhi(),
    theLaserPhiError(),
		thePosTECR4BeamPositions(),
		thePosTECR6BeamPositions(),
		theNegTECR4BeamPositions(),
		theNegTECR6BeamPositions(),
		theTEC2TECBeamPositions(),
		thePosTECR4BeamPositionErrors(),
		thePosTECR6BeamPositionErrors(),
		theNegTECR4BeamPositionErrors(),
		theNegTECR6BeamPositionErrors(),
		theTEC2TECBeamPositionErrors(),
    theNumberOfIterations(0), theNumberOfAlignmentIterations(0),
    theBeamFitter(),
    theLASAlignPosTEC(),
    theLASAlignNegTEC(),
    theLASAlignTEC2TEC(),
		theAlignmentAlgorithmBW(),
		theUseBSFrame(theConf.getUntrackedParameter<bool>("UseBeamSplitterFrame", true)),
    theDigiStore(),
    theBeamProfileFitStore(),
    theDigiVector(),
    theAlignableTracker(),
	  theAlignRecordName( "TrackerAlignmentRcd" ),
	  theErrorRecordName( "TrackerAlignmentErrorRcd" )
{
  // load the configuration from the ParameterSet  
  edm::LogInfo("LaserAlignment") << "==========================================================="
				  << "\n===                Start configuration                  ==="
				  << "\n    theDebugLevel               = " << theDebugLevel
				  << "\n    theAlignPosTEC              = " << theAlignPosTEC
				  << "\n    theAlignNegTEC              = " << theAlignNegTEC
				  << "\n    theAlignTEC2TEC             = " << theAlignTEC2TEC
				  << "\n    theSearchPhiTIB             = " << theSearchPhiTIB
				  << "\n    theSearchPhiTOB             = " << theSearchPhiTOB
				  << "\n    theSearchPhiTEC             = " << theSearchPhiTEC 
				  << "\n    theSearchZTIB               = " << theSearchZTIB
				  << "\n    theSearchZTOB               = " << theSearchZTOB
				  << "\n    theMinAdcCounts             = " << theMinAdcCounts
				  << "\n    theNEventsPerLaserIntensity = " << theNEventsPerLaserIntensity
				  << "\n    theNEventsForAllIntensiteis = " << theNEventsForAllIntensities
				  << "\n    theDoAlignmentAfterNEvents  = " << theDoAlignmentAfterNEvents
				  << "\n    ROOT filename               = " << theFileName
				  << "\n    compression                 = " << theCompression
				  << "\n===========================================================";

  // alias for the Branches in the root files
  std::string alias ( theConf.getParameter<std::string>("@module_label") );  

  // declare the product to produce
  produces<edm::DetSetVector<SiStripDigi> >().setBranchAlias( alias + "siStripDigis" );
  produces<LASBeamProfileFitCollection>().setBranchAlias( alias + "LASBeamProfileFits" );
	produces<LASAlignmentParameterCollection>().setBranchAlias( alias + "LASAlignmentParameters" );

  // the beam profile fitter
  theBeamFitter = new BeamProfileFitter(theBeamFitPS);

  // the alignable tracker parts
  theLASAlignPosTEC = new LaserAlignmentPosTEC;

  theLASAlignNegTEC = new LaserAlignmentNegTEC;

  theLASAlignTEC2TEC = new LaserAlignmentTEC2TEC;

	// the alignment algorithm from Bruno
	theAlignmentAlgorithmBW = new AlignmentAlgorithmBW;
  
  // counter for the number of iterations, i.e. the number of BeamProfile fits and
  // local Millepede fits
  theNumberOfIterations = 0;
}

LaserAlignment::~LaserAlignment()
{
  if (theSaveHistograms)
    {
      // close the rootfile
      closeRootFile();
    }

  if (theFile != 0) { delete theFile; }

  if (theBeamFitter != 0) { delete theBeamFitter; }

  if (theLASAlignPosTEC != 0) { delete theLASAlignPosTEC; }
  if (theLASAlignNegTEC != 0) { delete theLASAlignNegTEC; }
  if (theLASAlignTEC2TEC != 0) { delete theLASAlignTEC2TEC; }
  if (theAlignableTracker != 0) { delete theAlignableTracker; }
	if (theAlignmentAlgorithmBW != 0) { delete theAlignmentAlgorithmBW; }
}

double LaserAlignment::angle(double theAngle)
{
  return (theAngle >= 0.0) ? theAngle : theAngle + 2.0*M_PI;
}

void LaserAlignment::beginJob(const edm::EventSetup& theSetup)
{
  // creating a new file
  theFile = new TFile(theFileName.c_str(),"RECREATE","CMS ROOT file");
  theFile->SetCompressionLevel(theCompression);
      
  // initialize the histograms
  if (theFile) 
    {
      this->initHistograms();
    }
  else 
    {
      throw cms::Exception("LaserAlignment") << "<LaserAlignment::beginJob()>: ERROR!!! something wrong with the RootFile" << std::endl;
    } 

  LogDebug("LaserAlignment:beginJob()") << " access the Tracker Geometry ";
  // access the tracker
  theSetup.get<TrackerDigiGeometryRecord>().get( theTrackerGeometry );
  theSetup.get<IdealGeometryRecord>().get( gD );

  // Create the alignable hierarchy
  LogDebug("LaserAlignment:beginJob()") << " create the alignable hierarchy ";
  theAlignableTracker = new AlignableTracker( &(*gD),
					      &(*theTrackerGeometry) );

}

std::vector<int> LaserAlignment::checkBeam(std::vector<std::string>::const_iterator iHistName, std::map<std::string, std::pair<DetId, TH1D*> >::iterator iHist)
{
	std::vector<int> result;
	std::string stringDisc;
  std::string stringRing;
  std::string stringBeam;
  bool isTEC2TEC = false;
  int theDisc = 0;
  int theRing = 0;
  int theBeam = 0;
  int theTECSide = 0;
  
  // check if we are in the Endcap
  switch (((iHist->second).first).subdetId())
    {
    case StripSubdetector::TIB:
      {
	break;
      }
    case StripSubdetector::TOB:
      {
	break;
      }
    case StripSubdetector::TEC:
      {
	TECDetId theTECDetId(((iHist->second).first).rawId());

	theTECSide = theTECDetId.side(); // 1 for TEC-, 2 for TEC+

	stringBeam = (*iHistName).at(4);
	stringRing = (*iHistName).at(9);
	stringDisc = (*iHistName).at(14);
	isTEC2TEC = ( (*iHistName).size() > 21 ) ? true : false;
	break;
      }
    }

  if ( stringRing == "4" ) { theRing = 4; }
  else if ( stringRing == "6" ) { theRing = 6; }

  if ( stringDisc == "1" ) { theDisc = 0; }
  else if ( stringDisc == "2" ) { theDisc = 1; }
  else if ( stringDisc == "3" ) { theDisc = 2; } 
  else if ( stringDisc == "4" ) { theDisc = 3; } 
  else if ( stringDisc == "5" ) { theDisc = 4; }
  else if ( stringDisc == "6" ) { theDisc = 5; } 
  else if ( stringDisc == "7" ) { theDisc = 6; } 
  else if ( stringDisc == "8" ) { theDisc = 7; } 
  else if ( stringDisc == "9" ) { theDisc = 8; } 

  if ( theRing == 4 )
    {
      if ( stringBeam == "0" ) { theBeam = 0; } 
      else if ( stringBeam == "1" ) { theBeam = 1; } 
      else if ( stringBeam == "2" ) { theBeam = 2; }
      else if ( stringBeam == "3" ) { theBeam = 3; } 
      else if ( stringBeam == "4" ) { theBeam = 4; }
      else if ( stringBeam == "5" ) { theBeam = 5; } 
      else if ( stringBeam == "6" ) { theBeam = 6; } 
      else if ( stringBeam == "7" ) { theBeam = 7; } 
    }
  else if ( theRing == 6 )
    {
      if ( stringBeam == "0" ) { theBeam = 0; } 
      else if ( stringBeam == "1" ) { theBeam = 1; } 
      else if ( stringBeam == "2" ) { theBeam = 2; }
      else if ( stringBeam == "3" ) { theBeam = 3; } 
      else if ( stringBeam == "4" ) { theBeam = 4; }
      else if ( stringBeam == "5" ) { theBeam = 5; } 
      else if ( stringBeam == "6" ) { theBeam = 6; } 
      else if ( stringBeam == "7" ) { theBeam = 7; } 
    }
		result.push_back(theTECSide);
		result.push_back(theRing);
		result.push_back(theBeam);
		result.push_back(theDisc);
		
		return result;
}

void LaserAlignment::produce(edm::Event& theEvent, edm::EventSetup const& theSetup) 
{
  LogDebug("LaserAlignment") << "==========================================================="
			      << "\n   Private analysis of event #"<< theEvent.id().event() 
			      << " in run #" << theEvent.id().run();
  
  // ----- check if the actual event can be used -----
  /* here we can later on add some criteria for good alignment events!? */
  theEvents++;
  
  LogDebug("LaserAlignment") << "     Total Event number = " << theEvents;
  
  
  // create an empty output collection
  std::auto_ptr<LASBeamProfileFitCollection> theFitOutput(new LASBeamProfileFitCollection);
  std::auto_ptr<LASAlignmentParameterCollection> theAlignmentParameterOutput(new LASAlignmentParameterCollection);

  theDigiVector.reserve(10000);
  theDigiVector.clear();
  
  // do the Tracker Statistics
  trackerStatistics(theEvent, theSetup);
  
  LogDebug("LaserAlignment") << "===========================================================";
  
  // the fit of the beam profiles
  if (theEvents % theNEventsPerLaserIntensity == 0)
    {
      // do the beam profile fit
      fit(theSetup);
    }

  if (theEvents % theNEventsForAllIntensities == 0)
    {
      // increase the counter for the iterations
      theNumberOfIterations++;
      
      // put the digis of the beams into the StripDigiCollection
      for (std::map<DetId, std::vector<SiStripDigi> >::const_iterator p = theDigiStore.begin();
	   p != theDigiStore.end(); ++p)
	{
	  edm::DetSet<SiStripDigi> collector((p->first).rawId());
	  
	  if ((p->second).size()>0)
	    {
	      collector.data = (p->second);
	      
	      theDigiVector.push_back(collector);
	    }
	}
      // clear the map to fill new digis for the next theNEventsForAllIntensities number of events
      theDigiStore.clear();

      // put the LASBeamProfileFits into the LASBeamProfileFitCollection
      // loop over the map with the histograms and lookup the LASBeamFits
      for (std::vector<std::string>::const_iterator iHistName = theHistogramNames.begin(); iHistName != theHistogramNames.end(); ++iHistName)
	{
	  std::map<std::string, std::pair<DetId, TH1D*> >::iterator iHist = theHistograms.find(*iHistName);
	  if ( iHist != theHistograms.end() )
	    {
	      std::map<std::string, std::vector<LASBeamProfileFit> >::iterator iBeamFit = theBeamProfileFitStore.find(*iHistName);
	      if ( iBeamFit != theBeamProfileFitStore.end() )
		{
		  // get the DetId from the map with histograms
		  unsigned int theDetId = ((iHist->second).first).rawId();
		  // get the information for the LASBeamProfileFitCollection
		  LASBeamProfileFitCollection::Range inputRange;
		  inputRange.first = (iBeamFit->second).begin();
		  inputRange.second = (iBeamFit->second).end();
	      
		  theFitOutput->put(inputRange,theDetId);

		  // now fill the fitted phi position and error into the appropriate vectors
		  // which will be used by the Alignment Algorithm
		  LASBeamProfileFit theFit = (iBeamFit->second).at(0);
		  theLaserPhi.push_back(theFit.phi());
		  theLaserPhiError.push_back(thePhiErrorScalingFactor * theFit.phiError());
		
			// fill also the LASvec2D for Bruno's algorithm
			std::vector<int> sectorLocation = checkBeam(iHistName, iHist);
			
			if (sectorLocation.at(0) == 1) // TEC- beams
			{
				if (sectorLocation.at(1) == 4) // Ring 4
				{
					theNegTECR4BeamPositions[sectorLocation.at(2)][sectorLocation.at(3)] = theFit.pitch() * theFit.mean();
					theNegTECR4BeamPositionErrors[sectorLocation.at(2)][sectorLocation.at(3)] = theFit.pitch() * theFit.meanError();
				}
				else if (sectorLocation.at(1) == 6) // Ring 6
				{
					theNegTECR6BeamPositions[sectorLocation.at(2)][sectorLocation.at(3)] = theFit.pitch() * theFit.mean();
					theNegTECR6BeamPositionErrors[sectorLocation.at(2)][sectorLocation.at(4)] = theFit.pitch() * theFit.meanError();
				}
			}
			else if (sectorLocation.at(0) == 2) // TEC+ beams
			{
				if (sectorLocation.at(1) == 4) // Ring 4
				{
					thePosTECR4BeamPositions[sectorLocation.at(2)][sectorLocation.at(3)] = theFit.pitch() * theFit.mean();
					thePosTECR6BeamPositionErrors[sectorLocation.at(2)][sectorLocation.at(3)] = theFit.pitch() * theFit.meanError();
				}
				else if (sectorLocation.at(1) == 6) // Ring 6
				{
					thePosTECR6BeamPositions[sectorLocation.at(2)][sectorLocation.at(3)] = theFit.pitch() * theFit.mean();
					thePosTECR6BeamPositionErrors[sectorLocation.at(2)][sectorLocation.at(3)] = theFit.pitch() * theFit.meanError();
 				}				
			}
			// implementation of TEC2TEC beams has to be done!!!!!!
		}
	      else
		{
		  // no BeamFit found for this layer, use the nominal phi position of the Det for the Alignment Algorithm

// 		  // access the tracker
// 		  edm::ESHandle<TrackerGeometry> theTrackerGeometry;
// 		  theSetup.get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);

		  double thePhi = angle(theTrackerGeometry->idToDet((iHist->second).first)->surface().position().phi());
		  double thePhiError = 0.1;

		  LogDebug("LaserAlignment") << " no LASBeamProfileFit found for " << (*iHistName) << "! Use nominal phi position (" 
					      << thePhi << ") for alignment ";
		  theLaserPhi.push_back(thePhi);
		  theLaserPhiError.push_back(thePhiErrorScalingFactor * thePhiError);
		
			/// for Bruno's algorithm get the pitch at strip 255.5 and multiply it with 255.5 to get the position in cm at the middle of the module
			const StripGeomDetUnit* const theStripDet = dynamic_cast<const StripGeomDetUnit*>(theTrackerGeometry->idToDet((iHist->second).first));
			
			double thePosition = 255.5 * theStripDet->specificTopology().localPitch(theStripDet->specificTopology().localPosition(255.5));
			double thePositionError = 0.05; // set the error to half a milllimeter in this case

			// fill also the LASvec2D for Bruno's algorithm
			std::vector<int> sectorLocation = checkBeam(iHistName, iHist);
			
			if (sectorLocation.at(0) == 1) // TEC- beams
			{
				if (sectorLocation.at(1) == 4) // Ring 4
				{
					theNegTECR4BeamPositions[sectorLocation.at(2)][sectorLocation.at(3)] = thePosition;
					theNegTECR4BeamPositionErrors[sectorLocation.at(2)][sectorLocation.at(3)] = thePositionError;
				}
				else if (sectorLocation.at(1) == 6) // Ring 6
				{
					theNegTECR6BeamPositions[sectorLocation.at(2)][sectorLocation.at(3)] = thePosition;
					theNegTECR6BeamPositionErrors[sectorLocation.at(2)][sectorLocation.at(4)] = thePositionError;
				}
			}
			else if (sectorLocation.at(0) == 2) // TEC+ beams
			{
				if (sectorLocation.at(1) == 4) // Ring 4
				{
					thePosTECR4BeamPositions[sectorLocation.at(2)][sectorLocation.at(3)] = thePosition;
					thePosTECR6BeamPositionErrors[sectorLocation.at(2)][sectorLocation.at(3)] = thePositionError;
				}
				else if (sectorLocation.at(1) == 6) // Ring 6
				{
					thePosTECR6BeamPositions[sectorLocation.at(2)][sectorLocation.at(3)] = thePosition;
					thePosTECR6BeamPositionErrors[sectorLocation.at(2)][sectorLocation.at(3)] = thePositionError;
 				}				
			}
			// implementation of TEC2TEC beams has to be done!!!!!
		}
	    }
	  else 
	    { 
	      throw cms::Exception("LaserAlignment") << " You are in serious trouble!!! no entry for " << (*iHistName) << " found. "
						      << " To avoid the calculation of wrong alignment corrections the program will abort! "; 
	    }
	}
      // clear the map to fill new LASBeamProfileFits for the next theNEventsForAllIntensities number of events
      theBeamProfileFitStore.clear();
    }
  
  // do the alignment of the tracker
  if (theEvents % theDoAlignmentAfterNEvents == 0)
    {

//       // Create the alignable hierarchy
//       theAlignableTracker = new AlignableTracker( &(*gD),
// 						  &(*theTrackerGeometry) );

      // do the Alignment of the Tracker (with Millepede) ...
      alignmentAlgorithm(theAlignmentAlgorithmPS, theAlignableTracker);

			// do the Alignment with Bruno's algorithm
			std::vector<LASAlignmentParameter> theAPPosTECR4;
			std::vector<LASAlignmentParameter> theAPPosTECR6;
			std::vector<LASAlignmentParameter> theAPNegTECR4;
			std::vector<LASAlignmentParameter> theAPNegTECR6;
			// which return type do we need to calculate the final correction ((corr_R4 + corr_R6)/2)??
			// get them from the returned std::vector<LASAlignmentParameter>!
			theAPPosTECR4 = theAlignmentAlgorithmBW->run("TEC+Ring4",thePosTECR4BeamPositions, 
																										thePosTECR4BeamPositionErrors, theUseBSFrame, 4);
			theAPPosTECR6 = theAlignmentAlgorithmBW->run("TEC+Ring6",thePosTECR6BeamPositions, 
																										thePosTECR6BeamPositionErrors, theUseBSFrame, 6);
			theAPNegTECR4 = theAlignmentAlgorithmBW->run("TEC-Ring4",theNegTECR4BeamPositions, 
																										theNegTECR4BeamPositionErrors, theUseBSFrame, 4);
			theAPNegTECR6 = theAlignmentAlgorithmBW->run("TEC-Ring6",theNegTECR6BeamPositions, 
																										theNegTECR6BeamPositionErrors, theUseBSFrame, 6);
			// implementation of TEC2TEC has to be done!!!!!
			// theAlignmentAlgorithmBW->run(theTEC2TECBeamPositions, theTEC2TECBeamPositionErrors);

			// put the alignment parameters in the output collection
			// first TEC+ Ring 4
		  LASAlignmentParameterCollection::Range inputRange;
		  inputRange.first = theAPPosTECR4.begin();
		  inputRange.second = theAPPosTECR4.end();
		  theAlignmentParameterOutput->put(inputRange,"TEC+Ring4");

			// now TEC+ Ring 6
			inputRange.first = theAPPosTECR6.begin();
			inputRange.second = theAPPosTECR6.end();
			theAlignmentParameterOutput->put(inputRange,"TEC+Ring6");
			
			// now TEC- Ring 4
			inputRange.first = theAPNegTECR4.begin();
			inputRange.second = theAPNegTECR6.end();
			theAlignmentParameterOutput->put(inputRange,"TEC-Ring4");
			
			// finally TEC- Ring 6
			inputRange.first = theAPNegTECR6.begin();
			inputRange.second = theAPNegTECR6.end();
			theAlignmentParameterOutput->put(inputRange,"TEC-Ring6");
			
		  // set the number of iterations to zero for the next alignment round *** NOT NEEDED FOR BRUNO'S ALGO!
      theNumberOfIterations = 0;

			// apply the calculated corrections to the alignable Tracker. This will propagate the corrections from the
			// higher level hierarchy (e.g. the TEC discs) to the lowest level in the Tracker Hierarchy: the Dets
			/**
			 * TO BE DONE!!!!!!!!!!
			**/

			// store the estimated alignment parameters into the DB
      // first get them
      Alignments* alignments =  theAlignableTracker->alignments();
      AlignmentErrors* alignmentErrors = theAlignableTracker->alignmentErrors();
  
      // Write alignments to DB: have to sort beforhand!
      if ( theStoreToDB )
        {
	  LogDebug("LaserAlignment") << " storing the calculated alignment parameters to the DataBase";
	      // Call service
	      edm::Service<cond::service::PoolDBOutputService> poolDbService;
	      if( !poolDbService.isAvailable() ) // Die if not available
	        throw cms::Exception("NotAvailable") << "PoolDBOutputService not available";

		  // Store
	      if ( poolDbService->isNewTagRequest(theAlignRecordName) )
	        poolDbService->createNewIOV<Alignments>( alignments, poolDbService->endOfTime(), 
	                                                 theAlignRecordName );
	      else
	        poolDbService->appendSinceTime<Alignments>( alignments, poolDbService->currentTime(), 
	                                                   theAlignRecordName );
	      if ( poolDbService->isNewTagRequest(theErrorRecordName) )
	        poolDbService->createNewIOV<AlignmentErrors>( alignmentErrors,
	                                                      poolDbService->endOfTime(), 
	                                                      theErrorRecordName );
	      else
	        poolDbService->appendSinceTime<AlignmentErrors>( alignmentErrors,
	                                                         poolDbService->currentTime(), 
	                                                         theErrorRecordName );
	}

//       // Store result to EventSetup
//       GeometryAligner aligner;
//       aligner.applyAlignments<TrackerGeometry>( &(*theTrackerGeometry), &(*alignments), &(*alignmentErrors) );
    }
  
  // create the output collection for the DetSetVector
  std::auto_ptr<edm::DetSetVector<SiStripDigi> > theDigiOutput(new edm::DetSetVector<SiStripDigi>(theDigiVector));

  // write output to file
  theEvent.put(theDigiOutput);
  theEvent.put(theFitOutput);
  theEvent.put(theAlignmentParameterOutput);

//   // clear the vector with pairs of DetIds and Histograms
//   theHistograms.clear();
}


void LaserAlignment::closeRootFile()
{
  theFile->Write();
}

void LaserAlignment::fillAdcCounts(TH1D * theHistogram, DetId theDetId,
				    edm::DetSet<SiStripDigi>::const_iterator digiRangeIterator,
				    edm::DetSet<SiStripDigi>::const_iterator digiRangeIteratorEnd)
{
  if (theDebugLevel > 4) std::cout << "<LaserAlignment::fillAdcCounts()>: DetUnit: " << theDetId.rawId() << std::endl;

  // loop over all the digis in this det
  for (; digiRangeIterator != digiRangeIteratorEnd; ++digiRangeIterator) 
    {
      const SiStripDigi *digi = &*digiRangeIterator;

      // store the digis from the laser beams. They are later on used to create
      // clusters and RecHits. In this way some sort of "Laser Tracks" can be
      // reconstructed, which are useable for Track Based Alignment      
      theDigiStore[theDetId].push_back((*digi));
      
      if ( theDebugLevel > 5 ) 
	{ std::cout << " Channel " << digi->channel() << " has " << digi->adc() << " adc counts " << std::endl; }

      // fill the number of adc counts in the histogram
      if (digi->channel() < 512)
	{
	  Double_t theBinContent = theHistogram->GetBinContent(digi->channel()) + digi->adc();
	  theHistogram->SetBinContent(digi->channel(), theBinContent);
	}
    }
}
// define the SEAL module
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(LaserAlignment);
