#ifndef LaserAlignment_LaserAlignment_H
#define LaserAlignment_LaserAlignment_H

/** \class LaserAlignment
 *  Main reconstruction module for the Laser Alignment System
 *
 *  $Date: 2007/12/04 23:51:42 $
 *  $Revision: 1.9 $
 *  \author Maarten Thomas
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Alignment/LaserAlignment/interface/BeamProfileFitter.h"

#include "Alignment/LaserAlignment/interface/LaserAlignmentPosTEC.h"
#include "Alignment/LaserAlignment/interface/LaserAlignmentNegTEC.h"
#include "Alignment/LaserAlignment/interface/LaserAlignmentTEC2TEC.h"
#include "Alignment/LaserAlignment/interface/AlignmentAlgorithmBW.h"

#include "Alignment/LaserAlignment/interface/LASvector.h"
#include "Alignment/LaserAlignment/interface/LASvector2D.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"

// Alignable Tracker needed to propagate the alignment corrections calculated 
// for the disks down to the lowest levels
#include "Alignment/TrackerAlignment/interface/AlignableTracker.h"


// ROOT
#include "TH1.h"
#include "TObject.h"
class TFile;
class TH1D;

#include <iostream>

class LaserAlignment : public edm::EDProducer, public TObject 
{
 public:
  typedef std::vector<edm::ParameterSet> Parameters;
	/// define vector and matrix formats for easier calculation of the alignment corrections
  typedef LASvector<double> LASvec;
  typedef LASvector2D<double> LASvec2D;

	/// constructor
  explicit LaserAlignment(edm::ParameterSet const& theConf);
	/// destructor
  ~LaserAlignment();
  
	/// begin job
  virtual void beginJob(const edm::EventSetup&);
	/// produce LAS products
  virtual void produce(edm::Event& theEvent, edm::EventSetup const& theSetup);

 private:
	/// return angle in radian between 0 and 2*pi
  double angle(double theAngle);
	/// check in which subdetector/sector/disc we currently are
	std::vector<int> checkBeam(std::vector<std::string>::const_iterator iHistName, std::map<std::string, std::pair<DetId, TH1D*> >::iterator iHist);
	/// write the ROOT file with histograms
  void closeRootFile();
  /// fill adc counts from the laser profiles into a histogram
  void fillAdcCounts(TH1D * theHistogram, DetId theDetId,
		     edm::DetSet<SiStripDigi>::const_iterator digiRangeIterator,
		     edm::DetSet<SiStripDigi>::const_iterator digiRangeIteratorEnd);
	/// initialize the histograms
  void initHistograms();
	/// search for dets which are hit by a laser beam and fill the profiles into a histogram
  void trackerStatistics(edm::Event const& theEvent, edm::EventSetup const& theSetup);
  /// do the beam profile fit
  void fit(edm::EventSetup const& theSetup);
  /// calculate alignment corrections
  void alignmentAlgorithm(edm::ParameterSet const& theAlgorithmConf, 
			  AlignableTracker * theAlignableTracker);
    
 private:
  int theEvents;
  bool theStoreToDB;
  bool theSaveHistograms;
  int theDebugLevel;
  int theNEventsPerLaserIntensity;
  int theNEventsForAllIntensities;
  int theDoAlignmentAfterNEvents;
  bool theAlignPosTEC;
  bool theAlignNegTEC;
  bool theAlignTEC2TEC;
  bool theUseBrunosAlgorithm;
  bool theIsGoodFit;
  double theSearchPhiTIB;
  double theSearchPhiTOB;
  double theSearchPhiTEC;
  double theSearchZTIB;
  double theSearchZTOB;

  double thePhiErrorScalingFactor;

  // digi producer
  Parameters theDigiProducersList;

  // Tree stuff
  TFile * theFile;
  int theCompression;
  std::string theFileName;

  // parameter set for BeamProfileFitter
  edm::ParameterSet theBeamFitPS;

  // parameter set for the AlignmentAlgorithm
  edm::ParameterSet theAlignmentAlgorithmPS;

  // minimum number of AdcCounts to fill
  int theMinAdcCounts;

  // vector with HistogramNames and the map with a pair<DetId, TH1D*> for each beam
  std::vector<std::string> theHistogramNames;
  std::map<std::string, std::pair<DetId, TH1D*> > theHistograms;

  // vectors to store the Phi positions of the
  // fitted laser beams plus the error on Phi
  std::vector<double> theLaserPhi;
  std::vector<double> theLaserPhiError;

  // counter for the iterations
  int theNumberOfIterations;
  // counter for the number of Alignment Iterations
  int theNumberOfAlignmentIterations;

  // the Beam Profile Fitter
  BeamProfileFitter * theBeamFitter;

  // the alignable Tracker parts
  LaserAlignmentPosTEC * theLASAlignPosTEC;
  LaserAlignmentNegTEC * theLASAlignNegTEC;
  LaserAlignmentTEC2TEC * theLASAlignTEC2TEC;

	/// Bruno's alignment algorithm
	AlignmentAlgorithmBW * theAlignmentAlgorithmBW;
	/// use the BS frame in the alignment algorithm (i.e. BS at z = 0)
	bool theUseBSFrame;

  // the map to store digis for cluster creation
  std::map<DetId, std::vector<SiStripDigi> > theDigiStore;
  // map to store temporary the LASBeamProfileFits
  std::map<std::string, std::vector<LASBeamProfileFit> > theBeamProfileFitStore;

  // the vector which contains the Digis
  std::vector<edm::DetSet<SiStripDigi> > theDigiVector;

  // tracker geometry;
  edm::ESHandle<GeometricDet> gD;
  edm::ESHandle<TrackerGeometry> theTrackerGeometry;


/*  edm::ESHandle<DDCompactView> cpv;
  edm::ESHandle<GeometricDet> gD;
  
  boost::shared_ptr<TrackerGeometry> theTrackerGeometry;
*/
  // alignable tracker
  AlignableTracker * theAlignableTracker;

  std::string theAlignRecordName, theErrorRecordName;

  // monitor histograms

  /* Laser Beams in TEC+ */
  // Adc counts for Beam 0 in Ring 4
  TH1D * theBeam0Ring4Disc1PosAdcCounts;
  TH1D * theBeam0Ring4Disc2PosAdcCounts;
  TH1D * theBeam0Ring4Disc3PosAdcCounts;
  TH1D * theBeam0Ring4Disc4PosAdcCounts;
  TH1D * theBeam0Ring4Disc5PosAdcCounts;
  TH1D * theBeam0Ring4Disc6PosAdcCounts;
  TH1D * theBeam0Ring4Disc7PosAdcCounts;
  TH1D * theBeam0Ring4Disc8PosAdcCounts;
  TH1D * theBeam0Ring4Disc9PosAdcCounts;

  // Adc counts for Beam 1 in Ring 4
  TH1D * theBeam1Ring4Disc1PosAdcCounts;
  TH1D * theBeam1Ring4Disc2PosAdcCounts;
  TH1D * theBeam1Ring4Disc3PosAdcCounts;
  TH1D * theBeam1Ring4Disc4PosAdcCounts;
  TH1D * theBeam1Ring4Disc5PosAdcCounts;
  TH1D * theBeam1Ring4Disc6PosAdcCounts;
  TH1D * theBeam1Ring4Disc7PosAdcCounts;
  TH1D * theBeam1Ring4Disc8PosAdcCounts;
  TH1D * theBeam1Ring4Disc9PosAdcCounts;

  // plots for TEC2TEC
  TH1D * theBeam1Ring4Disc1PosTEC2TECAdcCounts;
  TH1D * theBeam1Ring4Disc2PosTEC2TECAdcCounts;
  TH1D * theBeam1Ring4Disc3PosTEC2TECAdcCounts;
  TH1D * theBeam1Ring4Disc4PosTEC2TECAdcCounts;
  TH1D * theBeam1Ring4Disc5PosTEC2TECAdcCounts;

  // Adc counts for Beam 2 in Ring 4
  TH1D * theBeam2Ring4Disc1PosAdcCounts;
  TH1D * theBeam2Ring4Disc2PosAdcCounts;
  TH1D * theBeam2Ring4Disc3PosAdcCounts;
  TH1D * theBeam2Ring4Disc4PosAdcCounts;
  TH1D * theBeam2Ring4Disc5PosAdcCounts;
  TH1D * theBeam2Ring4Disc6PosAdcCounts;
  TH1D * theBeam2Ring4Disc7PosAdcCounts;
  TH1D * theBeam2Ring4Disc8PosAdcCounts;
  TH1D * theBeam2Ring4Disc9PosAdcCounts;

  // plots for TEC2TEC
  TH1D * theBeam2Ring4Disc1PosTEC2TECAdcCounts;
  TH1D * theBeam2Ring4Disc2PosTEC2TECAdcCounts;
  TH1D * theBeam2Ring4Disc3PosTEC2TECAdcCounts;
  TH1D * theBeam2Ring4Disc4PosTEC2TECAdcCounts;
  TH1D * theBeam2Ring4Disc5PosTEC2TECAdcCounts;

  // Adc counts for Beam 3 in Ring 4
  TH1D * theBeam3Ring4Disc1PosAdcCounts;
  TH1D * theBeam3Ring4Disc2PosAdcCounts;
  TH1D * theBeam3Ring4Disc3PosAdcCounts;
  TH1D * theBeam3Ring4Disc4PosAdcCounts;
  TH1D * theBeam3Ring4Disc5PosAdcCounts;
  TH1D * theBeam3Ring4Disc6PosAdcCounts;
  TH1D * theBeam3Ring4Disc7PosAdcCounts;
  TH1D * theBeam3Ring4Disc8PosAdcCounts;
  TH1D * theBeam3Ring4Disc9PosAdcCounts;

  // Adc counts for Beam 4 in Ring 4
  TH1D * theBeam4Ring4Disc1PosAdcCounts;
  TH1D * theBeam4Ring4Disc2PosAdcCounts;
  TH1D * theBeam4Ring4Disc3PosAdcCounts;
  TH1D * theBeam4Ring4Disc4PosAdcCounts;
  TH1D * theBeam4Ring4Disc5PosAdcCounts;
  TH1D * theBeam4Ring4Disc6PosAdcCounts;
  TH1D * theBeam4Ring4Disc7PosAdcCounts;
  TH1D * theBeam4Ring4Disc8PosAdcCounts;
  TH1D * theBeam4Ring4Disc9PosAdcCounts;

  // plots for TEC2TEC
  TH1D * theBeam4Ring4Disc1PosTEC2TECAdcCounts;
  TH1D * theBeam4Ring4Disc2PosTEC2TECAdcCounts;
  TH1D * theBeam4Ring4Disc3PosTEC2TECAdcCounts;
  TH1D * theBeam4Ring4Disc4PosTEC2TECAdcCounts;
  TH1D * theBeam4Ring4Disc5PosTEC2TECAdcCounts;

  // Adc counts for Beam 5 in Ring 4
  TH1D * theBeam5Ring4Disc1PosAdcCounts;
  TH1D * theBeam5Ring4Disc2PosAdcCounts;
  TH1D * theBeam5Ring4Disc3PosAdcCounts;
  TH1D * theBeam5Ring4Disc4PosAdcCounts;
  TH1D * theBeam5Ring4Disc5PosAdcCounts;
  TH1D * theBeam5Ring4Disc6PosAdcCounts;
  TH1D * theBeam5Ring4Disc7PosAdcCounts;
  TH1D * theBeam5Ring4Disc8PosAdcCounts;
  TH1D * theBeam5Ring4Disc9PosAdcCounts;

  // Adc counts for Beam 6 in Ring 4
  TH1D * theBeam6Ring4Disc1PosAdcCounts;
  TH1D * theBeam6Ring4Disc2PosAdcCounts;
  TH1D * theBeam6Ring4Disc3PosAdcCounts;
  TH1D * theBeam6Ring4Disc4PosAdcCounts;
  TH1D * theBeam6Ring4Disc5PosAdcCounts;
  TH1D * theBeam6Ring4Disc6PosAdcCounts;
  TH1D * theBeam6Ring4Disc7PosAdcCounts;
  TH1D * theBeam6Ring4Disc8PosAdcCounts;
  TH1D * theBeam6Ring4Disc9PosAdcCounts;

  // plots for TEC2TEC
  TH1D * theBeam6Ring4Disc1PosTEC2TECAdcCounts;
  TH1D * theBeam6Ring4Disc2PosTEC2TECAdcCounts;
  TH1D * theBeam6Ring4Disc3PosTEC2TECAdcCounts;
  TH1D * theBeam6Ring4Disc4PosTEC2TECAdcCounts;
  TH1D * theBeam6Ring4Disc5PosTEC2TECAdcCounts;

  // Adc counts for Beam 7 in Ring 4
  TH1D * theBeam7Ring4Disc1PosAdcCounts;
  TH1D * theBeam7Ring4Disc2PosAdcCounts;
  TH1D * theBeam7Ring4Disc3PosAdcCounts;
  TH1D * theBeam7Ring4Disc4PosAdcCounts;
  TH1D * theBeam7Ring4Disc5PosAdcCounts;
  TH1D * theBeam7Ring4Disc6PosAdcCounts;
  TH1D * theBeam7Ring4Disc7PosAdcCounts;
  TH1D * theBeam7Ring4Disc8PosAdcCounts;
  TH1D * theBeam7Ring4Disc9PosAdcCounts;

  // plots for TEC2TEC
  TH1D * theBeam7Ring4Disc1PosTEC2TECAdcCounts;
  TH1D * theBeam7Ring4Disc2PosTEC2TECAdcCounts;
  TH1D * theBeam7Ring4Disc3PosTEC2TECAdcCounts;
  TH1D * theBeam7Ring4Disc4PosTEC2TECAdcCounts;
  TH1D * theBeam7Ring4Disc5PosTEC2TECAdcCounts;

  // Adc counts for Beam 0 in Ring 6
  TH1D * theBeam0Ring6Disc1PosAdcCounts;
  TH1D * theBeam0Ring6Disc2PosAdcCounts;
  TH1D * theBeam0Ring6Disc3PosAdcCounts;
  TH1D * theBeam0Ring6Disc4PosAdcCounts;
  TH1D * theBeam0Ring6Disc5PosAdcCounts;
  TH1D * theBeam0Ring6Disc6PosAdcCounts;
  TH1D * theBeam0Ring6Disc7PosAdcCounts;
  TH1D * theBeam0Ring6Disc8PosAdcCounts;
  TH1D * theBeam0Ring6Disc9PosAdcCounts;

  // Adc counts for Beam 1 in Ring 6
  TH1D * theBeam1Ring6Disc1PosAdcCounts;
  TH1D * theBeam1Ring6Disc2PosAdcCounts;
  TH1D * theBeam1Ring6Disc3PosAdcCounts;
  TH1D * theBeam1Ring6Disc4PosAdcCounts;
  TH1D * theBeam1Ring6Disc5PosAdcCounts;
  TH1D * theBeam1Ring6Disc6PosAdcCounts;
  TH1D * theBeam1Ring6Disc7PosAdcCounts;
  TH1D * theBeam1Ring6Disc8PosAdcCounts;
  TH1D * theBeam1Ring6Disc9PosAdcCounts;

  // Adc counts for Beam 2 in Ring 6
  TH1D * theBeam2Ring6Disc1PosAdcCounts;
  TH1D * theBeam2Ring6Disc2PosAdcCounts;
  TH1D * theBeam2Ring6Disc3PosAdcCounts;
  TH1D * theBeam2Ring6Disc4PosAdcCounts;
  TH1D * theBeam2Ring6Disc5PosAdcCounts;
  TH1D * theBeam2Ring6Disc6PosAdcCounts;
  TH1D * theBeam2Ring6Disc7PosAdcCounts;
  TH1D * theBeam2Ring6Disc8PosAdcCounts;
  TH1D * theBeam2Ring6Disc9PosAdcCounts;

  // Adc counts for Beam 3 in Ring 6
  TH1D * theBeam3Ring6Disc1PosAdcCounts;
  TH1D * theBeam3Ring6Disc2PosAdcCounts;
  TH1D * theBeam3Ring6Disc3PosAdcCounts;
  TH1D * theBeam3Ring6Disc4PosAdcCounts;
  TH1D * theBeam3Ring6Disc5PosAdcCounts;
  TH1D * theBeam3Ring6Disc6PosAdcCounts;
  TH1D * theBeam3Ring6Disc7PosAdcCounts;
  TH1D * theBeam3Ring6Disc8PosAdcCounts;
  TH1D * theBeam3Ring6Disc9PosAdcCounts;

  // Adc counts for Beam 4 in Ring 6
  TH1D * theBeam4Ring6Disc1PosAdcCounts;
  TH1D * theBeam4Ring6Disc2PosAdcCounts;
  TH1D * theBeam4Ring6Disc3PosAdcCounts;
  TH1D * theBeam4Ring6Disc4PosAdcCounts;
  TH1D * theBeam4Ring6Disc5PosAdcCounts;
  TH1D * theBeam4Ring6Disc6PosAdcCounts;
  TH1D * theBeam4Ring6Disc7PosAdcCounts;
  TH1D * theBeam4Ring6Disc8PosAdcCounts;
  TH1D * theBeam4Ring6Disc9PosAdcCounts;

  // Adc counts for Beam 5 in Ring 6
  TH1D * theBeam5Ring6Disc1PosAdcCounts;
  TH1D * theBeam5Ring6Disc2PosAdcCounts;
  TH1D * theBeam5Ring6Disc3PosAdcCounts;
  TH1D * theBeam5Ring6Disc4PosAdcCounts;
  TH1D * theBeam5Ring6Disc5PosAdcCounts;
  TH1D * theBeam5Ring6Disc6PosAdcCounts;
  TH1D * theBeam5Ring6Disc7PosAdcCounts;
  TH1D * theBeam5Ring6Disc8PosAdcCounts;
  TH1D * theBeam5Ring6Disc9PosAdcCounts;

  // Adc counts for Beam 6 in Ring 6
  TH1D * theBeam6Ring6Disc1PosAdcCounts;
  TH1D * theBeam6Ring6Disc2PosAdcCounts;
  TH1D * theBeam6Ring6Disc3PosAdcCounts;
  TH1D * theBeam6Ring6Disc4PosAdcCounts;
  TH1D * theBeam6Ring6Disc5PosAdcCounts;
  TH1D * theBeam6Ring6Disc6PosAdcCounts;
  TH1D * theBeam6Ring6Disc7PosAdcCounts;
  TH1D * theBeam6Ring6Disc8PosAdcCounts;
  TH1D * theBeam6Ring6Disc9PosAdcCounts;

  // Adc counts for Beam 7 in Ring 6
  TH1D * theBeam7Ring6Disc1PosAdcCounts;
  TH1D * theBeam7Ring6Disc2PosAdcCounts;
  TH1D * theBeam7Ring6Disc3PosAdcCounts;
  TH1D * theBeam7Ring6Disc4PosAdcCounts;
  TH1D * theBeam7Ring6Disc5PosAdcCounts;
  TH1D * theBeam7Ring6Disc6PosAdcCounts;
  TH1D * theBeam7Ring6Disc7PosAdcCounts;
  TH1D * theBeam7Ring6Disc8PosAdcCounts;
  TH1D * theBeam7Ring6Disc9PosAdcCounts;

  /* Laser Beams in TEC- */
  // Adc counts for Beam 0 in Ring 4
  TH1D * theBeam0Ring4Disc1NegAdcCounts;
  TH1D * theBeam0Ring4Disc2NegAdcCounts;
  TH1D * theBeam0Ring4Disc3NegAdcCounts;
  TH1D * theBeam0Ring4Disc4NegAdcCounts;
  TH1D * theBeam0Ring4Disc5NegAdcCounts;
  TH1D * theBeam0Ring4Disc6NegAdcCounts;
  TH1D * theBeam0Ring4Disc7NegAdcCounts;
  TH1D * theBeam0Ring4Disc8NegAdcCounts;
  TH1D * theBeam0Ring4Disc9NegAdcCounts;

  // Adc counts for Beam 1 in Ring 4
  TH1D * theBeam1Ring4Disc1NegAdcCounts;
  TH1D * theBeam1Ring4Disc2NegAdcCounts;
  TH1D * theBeam1Ring4Disc3NegAdcCounts;
  TH1D * theBeam1Ring4Disc4NegAdcCounts;
  TH1D * theBeam1Ring4Disc5NegAdcCounts;
  TH1D * theBeam1Ring4Disc6NegAdcCounts;
  TH1D * theBeam1Ring4Disc7NegAdcCounts;
  TH1D * theBeam1Ring4Disc8NegAdcCounts;
  TH1D * theBeam1Ring4Disc9NegAdcCounts;

  // plots for TEC2TEC
  TH1D * theBeam1Ring4Disc1NegTEC2TECAdcCounts;
  TH1D * theBeam1Ring4Disc2NegTEC2TECAdcCounts;
  TH1D * theBeam1Ring4Disc3NegTEC2TECAdcCounts;
  TH1D * theBeam1Ring4Disc4NegTEC2TECAdcCounts;
  TH1D * theBeam1Ring4Disc5NegTEC2TECAdcCounts;

  // Adc counts for Beam 2 in Ring 4
  TH1D * theBeam2Ring4Disc1NegAdcCounts;
  TH1D * theBeam2Ring4Disc2NegAdcCounts;
  TH1D * theBeam2Ring4Disc3NegAdcCounts;
  TH1D * theBeam2Ring4Disc4NegAdcCounts;
  TH1D * theBeam2Ring4Disc5NegAdcCounts;
  TH1D * theBeam2Ring4Disc6NegAdcCounts;
  TH1D * theBeam2Ring4Disc7NegAdcCounts;
  TH1D * theBeam2Ring4Disc8NegAdcCounts;
  TH1D * theBeam2Ring4Disc9NegAdcCounts;

  // plots for TEC2TEC
  TH1D * theBeam2Ring4Disc1NegTEC2TECAdcCounts;
  TH1D * theBeam2Ring4Disc2NegTEC2TECAdcCounts;
  TH1D * theBeam2Ring4Disc3NegTEC2TECAdcCounts;
  TH1D * theBeam2Ring4Disc4NegTEC2TECAdcCounts;
  TH1D * theBeam2Ring4Disc5NegTEC2TECAdcCounts;

  // Adc counts for Beam 3 in Ring 4
  TH1D * theBeam3Ring4Disc1NegAdcCounts;
  TH1D * theBeam3Ring4Disc2NegAdcCounts;
  TH1D * theBeam3Ring4Disc3NegAdcCounts;
  TH1D * theBeam3Ring4Disc4NegAdcCounts;
  TH1D * theBeam3Ring4Disc5NegAdcCounts;
  TH1D * theBeam3Ring4Disc6NegAdcCounts;
  TH1D * theBeam3Ring4Disc7NegAdcCounts;
  TH1D * theBeam3Ring4Disc8NegAdcCounts;
  TH1D * theBeam3Ring4Disc9NegAdcCounts;

  // Adc counts for Beam 4 in Ring 4
  TH1D * theBeam4Ring4Disc1NegAdcCounts;
  TH1D * theBeam4Ring4Disc2NegAdcCounts;
  TH1D * theBeam4Ring4Disc3NegAdcCounts;
  TH1D * theBeam4Ring4Disc4NegAdcCounts;
  TH1D * theBeam4Ring4Disc5NegAdcCounts;
  TH1D * theBeam4Ring4Disc6NegAdcCounts;
  TH1D * theBeam4Ring4Disc7NegAdcCounts;
  TH1D * theBeam4Ring4Disc8NegAdcCounts;
  TH1D * theBeam4Ring4Disc9NegAdcCounts;

  // plots for TEC2TEC
  TH1D * theBeam4Ring4Disc1NegTEC2TECAdcCounts;
  TH1D * theBeam4Ring4Disc2NegTEC2TECAdcCounts;
  TH1D * theBeam4Ring4Disc3NegTEC2TECAdcCounts;
  TH1D * theBeam4Ring4Disc4NegTEC2TECAdcCounts;
  TH1D * theBeam4Ring4Disc5NegTEC2TECAdcCounts;

  // Adc counts for Beam 5 in Ring 4
  TH1D * theBeam5Ring4Disc1NegAdcCounts;
  TH1D * theBeam5Ring4Disc2NegAdcCounts;
  TH1D * theBeam5Ring4Disc3NegAdcCounts;
  TH1D * theBeam5Ring4Disc4NegAdcCounts;
  TH1D * theBeam5Ring4Disc5NegAdcCounts;
  TH1D * theBeam5Ring4Disc6NegAdcCounts;
  TH1D * theBeam5Ring4Disc7NegAdcCounts;
  TH1D * theBeam5Ring4Disc8NegAdcCounts;
  TH1D * theBeam5Ring4Disc9NegAdcCounts;

  // Adc counts for Beam 6 in Ring 4
  TH1D * theBeam6Ring4Disc1NegAdcCounts;
  TH1D * theBeam6Ring4Disc2NegAdcCounts;
  TH1D * theBeam6Ring4Disc3NegAdcCounts;
  TH1D * theBeam6Ring4Disc4NegAdcCounts;
  TH1D * theBeam6Ring4Disc5NegAdcCounts;
  TH1D * theBeam6Ring4Disc6NegAdcCounts;
  TH1D * theBeam6Ring4Disc7NegAdcCounts;
  TH1D * theBeam6Ring4Disc8NegAdcCounts;
  TH1D * theBeam6Ring4Disc9NegAdcCounts;

  // plots for TEC2TEC
  TH1D * theBeam6Ring4Disc1NegTEC2TECAdcCounts;
  TH1D * theBeam6Ring4Disc2NegTEC2TECAdcCounts;
  TH1D * theBeam6Ring4Disc3NegTEC2TECAdcCounts;
  TH1D * theBeam6Ring4Disc4NegTEC2TECAdcCounts;
  TH1D * theBeam6Ring4Disc5NegTEC2TECAdcCounts;

  // Adc counts for Beam 7 in Ring 4
  TH1D * theBeam7Ring4Disc1NegAdcCounts;
  TH1D * theBeam7Ring4Disc2NegAdcCounts;
  TH1D * theBeam7Ring4Disc3NegAdcCounts;
  TH1D * theBeam7Ring4Disc4NegAdcCounts;
  TH1D * theBeam7Ring4Disc5NegAdcCounts;
  TH1D * theBeam7Ring4Disc6NegAdcCounts;
  TH1D * theBeam7Ring4Disc7NegAdcCounts;
  TH1D * theBeam7Ring4Disc8NegAdcCounts;
  TH1D * theBeam7Ring4Disc9NegAdcCounts;

  // plots for TEC2TEC
  TH1D * theBeam7Ring4Disc1NegTEC2TECAdcCounts;
  TH1D * theBeam7Ring4Disc2NegTEC2TECAdcCounts;
  TH1D * theBeam7Ring4Disc3NegTEC2TECAdcCounts;
  TH1D * theBeam7Ring4Disc4NegTEC2TECAdcCounts;
  TH1D * theBeam7Ring4Disc5NegTEC2TECAdcCounts;

  // Adc counts for Beam 0 in Ring 6
  TH1D * theBeam0Ring6Disc1NegAdcCounts;
  TH1D * theBeam0Ring6Disc2NegAdcCounts;
  TH1D * theBeam0Ring6Disc3NegAdcCounts;
  TH1D * theBeam0Ring6Disc4NegAdcCounts;
  TH1D * theBeam0Ring6Disc5NegAdcCounts;
  TH1D * theBeam0Ring6Disc6NegAdcCounts;
  TH1D * theBeam0Ring6Disc7NegAdcCounts;
  TH1D * theBeam0Ring6Disc8NegAdcCounts;
  TH1D * theBeam0Ring6Disc9NegAdcCounts;

  // Adc counts for Beam 1 in Ring 6
  TH1D * theBeam1Ring6Disc1NegAdcCounts;
  TH1D * theBeam1Ring6Disc2NegAdcCounts;
  TH1D * theBeam1Ring6Disc3NegAdcCounts;
  TH1D * theBeam1Ring6Disc4NegAdcCounts;
  TH1D * theBeam1Ring6Disc5NegAdcCounts;
  TH1D * theBeam1Ring6Disc6NegAdcCounts;
  TH1D * theBeam1Ring6Disc7NegAdcCounts;
  TH1D * theBeam1Ring6Disc8NegAdcCounts;
  TH1D * theBeam1Ring6Disc9NegAdcCounts;

  // Adc counts for Beam 2 in Ring 6
  TH1D * theBeam2Ring6Disc1NegAdcCounts;
  TH1D * theBeam2Ring6Disc2NegAdcCounts;
  TH1D * theBeam2Ring6Disc3NegAdcCounts;
  TH1D * theBeam2Ring6Disc4NegAdcCounts;
  TH1D * theBeam2Ring6Disc5NegAdcCounts;
  TH1D * theBeam2Ring6Disc6NegAdcCounts;
  TH1D * theBeam2Ring6Disc7NegAdcCounts;
  TH1D * theBeam2Ring6Disc8NegAdcCounts;
  TH1D * theBeam2Ring6Disc9NegAdcCounts;

  // Adc counts for Beam 3 in Ring 6
  TH1D * theBeam3Ring6Disc1NegAdcCounts;
  TH1D * theBeam3Ring6Disc2NegAdcCounts;
  TH1D * theBeam3Ring6Disc3NegAdcCounts;
  TH1D * theBeam3Ring6Disc4NegAdcCounts;
  TH1D * theBeam3Ring6Disc5NegAdcCounts;
  TH1D * theBeam3Ring6Disc6NegAdcCounts;
  TH1D * theBeam3Ring6Disc7NegAdcCounts;
  TH1D * theBeam3Ring6Disc8NegAdcCounts;
  TH1D * theBeam3Ring6Disc9NegAdcCounts;

  // Adc counts for Beam 4 in Ring 6
  TH1D * theBeam4Ring6Disc1NegAdcCounts;
  TH1D * theBeam4Ring6Disc2NegAdcCounts;
  TH1D * theBeam4Ring6Disc3NegAdcCounts;
  TH1D * theBeam4Ring6Disc4NegAdcCounts;
  TH1D * theBeam4Ring6Disc5NegAdcCounts;
  TH1D * theBeam4Ring6Disc6NegAdcCounts;
  TH1D * theBeam4Ring6Disc7NegAdcCounts;
  TH1D * theBeam4Ring6Disc8NegAdcCounts;
  TH1D * theBeam4Ring6Disc9NegAdcCounts;

  // Adc counts for Beam 5 in Ring 6
  TH1D * theBeam5Ring6Disc1NegAdcCounts;
  TH1D * theBeam5Ring6Disc2NegAdcCounts;
  TH1D * theBeam5Ring6Disc3NegAdcCounts;
  TH1D * theBeam5Ring6Disc4NegAdcCounts;
  TH1D * theBeam5Ring6Disc5NegAdcCounts;
  TH1D * theBeam5Ring6Disc6NegAdcCounts;
  TH1D * theBeam5Ring6Disc7NegAdcCounts;
  TH1D * theBeam5Ring6Disc8NegAdcCounts;
  TH1D * theBeam5Ring6Disc9NegAdcCounts;

  // Adc counts for Beam 6 in Ring 6
  TH1D * theBeam6Ring6Disc1NegAdcCounts;
  TH1D * theBeam6Ring6Disc2NegAdcCounts;
  TH1D * theBeam6Ring6Disc3NegAdcCounts;
  TH1D * theBeam6Ring6Disc4NegAdcCounts;
  TH1D * theBeam6Ring6Disc5NegAdcCounts;
  TH1D * theBeam6Ring6Disc6NegAdcCounts;
  TH1D * theBeam6Ring6Disc7NegAdcCounts;
  TH1D * theBeam6Ring6Disc8NegAdcCounts;
  TH1D * theBeam6Ring6Disc9NegAdcCounts;

  // Adc counts for Beam 7 in Ring 6
  TH1D * theBeam7Ring6Disc1NegAdcCounts;
  TH1D * theBeam7Ring6Disc2NegAdcCounts;
  TH1D * theBeam7Ring6Disc3NegAdcCounts;
  TH1D * theBeam7Ring6Disc4NegAdcCounts;
  TH1D * theBeam7Ring6Disc5NegAdcCounts;
  TH1D * theBeam7Ring6Disc6NegAdcCounts;
  TH1D * theBeam7Ring6Disc7NegAdcCounts;
  TH1D * theBeam7Ring6Disc8NegAdcCounts;
  TH1D * theBeam7Ring6Disc9NegAdcCounts;

  // TOB Beams
  // Adc counts for Beam 0
  TH1D * theBeam0TOBPosition1AdcCounts;
  TH1D * theBeam0TOBPosition2AdcCounts;
  TH1D * theBeam0TOBPosition3AdcCounts;
  TH1D * theBeam0TOBPosition4AdcCounts;
  TH1D * theBeam0TOBPosition5AdcCounts;
  TH1D * theBeam0TOBPosition6AdcCounts;

  // Adc counts for Beam 1
  TH1D * theBeam1TOBPosition1AdcCounts;
  TH1D * theBeam1TOBPosition2AdcCounts;
  TH1D * theBeam1TOBPosition3AdcCounts;
  TH1D * theBeam1TOBPosition4AdcCounts;
  TH1D * theBeam1TOBPosition5AdcCounts;
  TH1D * theBeam1TOBPosition6AdcCounts;

  // Adc counts for Beam 2
  TH1D * theBeam2TOBPosition1AdcCounts;
  TH1D * theBeam2TOBPosition2AdcCounts;
  TH1D * theBeam2TOBPosition3AdcCounts;
  TH1D * theBeam2TOBPosition4AdcCounts;
  TH1D * theBeam2TOBPosition5AdcCounts;
  TH1D * theBeam2TOBPosition6AdcCounts;

  // Adc counts for Beam 3
  TH1D * theBeam3TOBPosition1AdcCounts;
  TH1D * theBeam3TOBPosition2AdcCounts;
  TH1D * theBeam3TOBPosition3AdcCounts;
  TH1D * theBeam3TOBPosition4AdcCounts;
  TH1D * theBeam3TOBPosition5AdcCounts;
  TH1D * theBeam3TOBPosition6AdcCounts;

  // Adc counts for Beam 4
  TH1D * theBeam4TOBPosition1AdcCounts;
  TH1D * theBeam4TOBPosition2AdcCounts;
  TH1D * theBeam4TOBPosition3AdcCounts;
  TH1D * theBeam4TOBPosition4AdcCounts;
  TH1D * theBeam4TOBPosition5AdcCounts;
  TH1D * theBeam4TOBPosition6AdcCounts;

  // Adc counts for Beam 5
  TH1D * theBeam5TOBPosition1AdcCounts;
  TH1D * theBeam5TOBPosition2AdcCounts;
  TH1D * theBeam5TOBPosition3AdcCounts;
  TH1D * theBeam5TOBPosition4AdcCounts;
  TH1D * theBeam5TOBPosition5AdcCounts;
  TH1D * theBeam5TOBPosition6AdcCounts;

  // Adc counts for Beam 6
  TH1D * theBeam6TOBPosition1AdcCounts;
  TH1D * theBeam6TOBPosition2AdcCounts;
  TH1D * theBeam6TOBPosition3AdcCounts;
  TH1D * theBeam6TOBPosition4AdcCounts;
  TH1D * theBeam6TOBPosition5AdcCounts;
  TH1D * theBeam6TOBPosition6AdcCounts;

  // Adc counts for Beam 7
  TH1D * theBeam7TOBPosition1AdcCounts;
  TH1D * theBeam7TOBPosition2AdcCounts;
  TH1D * theBeam7TOBPosition3AdcCounts;
  TH1D * theBeam7TOBPosition4AdcCounts;
  TH1D * theBeam7TOBPosition5AdcCounts;
  TH1D * theBeam7TOBPosition6AdcCounts;

  // TIB Beams
  // Adc counts for Beam 0
  TH1D * theBeam0TIBPosition1AdcCounts;
  TH1D * theBeam0TIBPosition2AdcCounts;
  TH1D * theBeam0TIBPosition3AdcCounts;
  TH1D * theBeam0TIBPosition4AdcCounts;
  TH1D * theBeam0TIBPosition5AdcCounts;
  TH1D * theBeam0TIBPosition6AdcCounts;

  // Adc counts for Beam 1
  TH1D * theBeam1TIBPosition1AdcCounts;
  TH1D * theBeam1TIBPosition2AdcCounts;
  TH1D * theBeam1TIBPosition3AdcCounts;
  TH1D * theBeam1TIBPosition4AdcCounts;
  TH1D * theBeam1TIBPosition5AdcCounts;
  TH1D * theBeam1TIBPosition6AdcCounts;

  // Adc counts for Beam 2
  TH1D * theBeam2TIBPosition1AdcCounts;
  TH1D * theBeam2TIBPosition2AdcCounts;
  TH1D * theBeam2TIBPosition3AdcCounts;
  TH1D * theBeam2TIBPosition4AdcCounts;
  TH1D * theBeam2TIBPosition5AdcCounts;
  TH1D * theBeam2TIBPosition6AdcCounts;

  // Adc counts for Beam 3
  TH1D * theBeam3TIBPosition1AdcCounts;
  TH1D * theBeam3TIBPosition2AdcCounts;
  TH1D * theBeam3TIBPosition3AdcCounts;
  TH1D * theBeam3TIBPosition4AdcCounts;
  TH1D * theBeam3TIBPosition5AdcCounts;
  TH1D * theBeam3TIBPosition6AdcCounts;

  // Adc counts for Beam 4
  TH1D * theBeam4TIBPosition1AdcCounts;
  TH1D * theBeam4TIBPosition2AdcCounts;
  TH1D * theBeam4TIBPosition3AdcCounts;
  TH1D * theBeam4TIBPosition4AdcCounts;
  TH1D * theBeam4TIBPosition5AdcCounts;
  TH1D * theBeam4TIBPosition6AdcCounts;

  // Adc counts for Beam 5
  TH1D * theBeam5TIBPosition1AdcCounts;
  TH1D * theBeam5TIBPosition2AdcCounts;
  TH1D * theBeam5TIBPosition3AdcCounts;
  TH1D * theBeam5TIBPosition4AdcCounts;
  TH1D * theBeam5TIBPosition5AdcCounts;
  TH1D * theBeam5TIBPosition6AdcCounts;

  // Adc counts for Beam 6
  TH1D * theBeam6TIBPosition1AdcCounts;
  TH1D * theBeam6TIBPosition2AdcCounts;
  TH1D * theBeam6TIBPosition3AdcCounts;
  TH1D * theBeam6TIBPosition4AdcCounts;
  TH1D * theBeam6TIBPosition5AdcCounts;
  TH1D * theBeam6TIBPosition6AdcCounts;

  // Adc counts for Beam 7
  TH1D * theBeam7TIBPosition1AdcCounts;
  TH1D * theBeam7TIBPosition2AdcCounts;
  TH1D * theBeam7TIBPosition3AdcCounts;
  TH1D * theBeam7TIBPosition4AdcCounts;
  TH1D * theBeam7TIBPosition5AdcCounts;
  TH1D * theBeam7TIBPosition6AdcCounts;

};
#endif
