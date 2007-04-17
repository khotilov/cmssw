/** SiPixelRecHitConverter.cc
 * ------------------------------------------------------
 * Description:  see SiPixelRecHitConverter.h
 * Authors:  P. Maksimovic (JHU), V.Chiochia (Uni Zurich)
 * History: Feb 27, 2006 -  initial version
 *          May 30, 2006 -  edm::DetSetVector and edm::Ref
 * ------------------------------------------------------
 */

// Our own stuff
#include "RecoLocalTracker/SiPixelRecHits/interface/SiPixelRecHitConverter.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/CPEFromDetPosition.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEInitial.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEParmError.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPETemplateReco.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEGeneric.h"

// Geometry
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

// Data Formats
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/Ref.h"

// Framework Already defined in the *.h file
//#include "DataFormats/Common/interface/Handle.h"
//#include "FWCore/Framework/interface/ESHandle.h"

// STL
#include <vector>
#include <memory>
#include <string>
#include <iostream>

// MessageLogger
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Magnetic Field
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

using namespace std;

namespace cms
{
  //---------------------------------------------------------------------------
  //!  Constructor: set the ParameterSet and defer all thinking to setupCPE().
  //---------------------------------------------------------------------------
  SiPixelRecHitConverter::SiPixelRecHitConverter(edm::ParameterSet const& conf) 
    : 
    conf_(conf),
    cpeName_("None"),     // bogus
    cpe_(0),              // the default, in case we fail to make one
    ready_(false),        // since we obviously aren't
    src_( conf.getParameter<edm::InputTag>( "src" ) )
  {
    //--- Declare to the EDM what kind of collections we will be making.
    produces<SiPixelRecHitCollection>();
    //--- Algorithm's verbosity
    theVerboseLevel = 
      conf.getUntrackedParameter<int>("VerboseLevel",0);
  }
  
  // Destructor
  SiPixelRecHitConverter::~SiPixelRecHitConverter() 
  { 
    delete cpe_;
  }  
  
  //---------------------------------------------------------------------------
  // Begin job: get magnetic field
  //---------------------------------------------------------------------------
  void SiPixelRecHitConverter::beginJob(const edm::EventSetup& c) 
  {
    edm::ESHandle<MagneticField> magfield;
    c.get<IdealMagneticFieldRecord>().get(magfield);
    setupCPE(magfield.product());
  
  }

  //---------------------------------------------------------------------------
  //! The "Event" entrypoint: gets called by framework for every event
  //---------------------------------------------------------------------------
  void SiPixelRecHitConverter::produce(edm::Event& e, const edm::EventSetup& es)
  {
    // Step A.1: get input data
    edm::Handle< edm::DetSetVector<SiPixelCluster> > input;
    e.getByLabel( src_, input);
    
    // Step A.2: get event setup
    edm::ESHandle<TrackerGeometry> geom;
    es.get<TrackerDigiGeometryRecord>().get( geom );
    
    // Step B: create empty output collection
    std::auto_ptr<SiPixelRecHitCollection> output(new SiPixelRecHitCollection);

    // Step C: Iterate over DetIds and invoke the strip CPE algorithm
    // on each DetUnit
    run( input, *output, geom );
    
    // Step D: write output to file
    e.put(output);
  }

  //---------------------------------------------------------------------------
  //!  Set up the specific algorithm we are going to use.  
  //---------------------------------------------------------------------------
  void SiPixelRecHitConverter::setupCPE(const MagneticField* mag) 
  {
    cpeName_ = conf_.getParameter<std::string>("CPE");
    
    if ( cpeName_ == "FromDetPosition" ) 
      {
	cout << endl;
	cout << "--------------------------- Using CPE FromDetPosition ---------------------------------------------" << endl;
	cpe_ = new CPEFromDetPosition(conf_,mag);
	ready_ = true;
      } 
    else if ( cpeName_ == "Initial" ) 
      {
	cout << endl;
	cout << "--------------------------- Using CPE Initial ---------------------------------------------" << endl;
	cpe_ = new PixelCPEInitial(conf_,mag);
	ready_ = true;
      } 
    else if ( cpeName_ == "ParmError" ) 
      {
	cout << endl;
	cout << "--------------------------- Using CPE ParmError ---------------------------------------------" << endl;
	cpe_ = new PixelCPEParmError(conf_,mag);
	ready_ = true;
      } 
    else if ( cpeName_ == "TemplateReco" ) 
      {
	cout << endl;
	cout << "--------------------------- Using CPE TemplateReco ---------------------------------------------" << endl;
	cpe_ = new PixelCPETemplateReco(conf_,mag);
	ready_ = true;
      } 
    else if ( cpeName_ == "Generic" ) 
      {
	cout << endl;
	cout << "--------------------------- Using CPE Generic --------------------------------------------------" << endl;
	cpe_ = new PixelCPEGeneric(conf_,mag);
	ready_ = true;
      } 
    else 
      {
	cout << endl;
	cout << "--------------------------- No valid CPE !!!!!!!!!! ---------------------------------------------" << endl;

	edm::LogError("SiPixelRecHitConverter") 
	  <<" Cluster parameter estimator " << cpeName_ << " is invalid.\n"
	  << "Possible choices:\n" 
	  << "    - FromDetPosition";
	ready_ = false;
      }
  }
  //---------------------------------------------------------------------------
  //!  Iterate over DetUnits, then over Clusters and invoke the CPE on each,
  //!  and make a RecHit to store the result.
  //!  New interface reading DetSetVector by V.Chiochia (May 30th, 2006)
  //---------------------------------------------------------------------------
  void SiPixelRecHitConverter::run(edm::Handle<edm::DetSetVector<SiPixelCluster> >  inputhandle,
				   SiPixelRecHitCollection &output,
				   edm::ESHandle<TrackerGeometry> & geom) 
  {
    if ( ! ready_ ) 
      {
	edm::LogError("SiPixelRecHitConverter") << " at least one CPE is not ready -- can't run!";
	// TO DO: throw an exception here?  The user may want to know...
	assert(0);
	return;   // clusterizer is invalid, bail out
      }
    
    int numberOfDetUnits = 0;
    int numberOfClusters = 0;
    
    const edm::DetSetVector<SiPixelCluster>& input = *inputhandle;
    
    edm::DetSetVector<SiPixelCluster>::const_iterator DSViter=input.begin();
    
    for ( ; DSViter != input.end() ; DSViter++) 
      {
	numberOfDetUnits++;
	unsigned int detid = DSViter->id;
	DetId detIdObject( detid );  
	const GeomDetUnit * genericDet = geom->idToDetUnit( detIdObject );
	const PixelGeomDetUnit * pixDet = dynamic_cast<const PixelGeomDetUnit*>(genericDet);
	assert(pixDet); 
	edm::OwnVector<SiPixelRecHit> recHitsOnDetUnit;
	
	edm::DetSet<SiPixelCluster>::const_iterator clustIt = DSViter->data.begin();
	for ( ; clustIt != DSViter->data.end(); clustIt++) 
	  {
	    numberOfClusters++;
	    std::pair<LocalPoint, LocalError> lv = 
	      cpe_->localParameters( *clustIt, *genericDet );
	    LocalPoint lp( lv.first );
	    LocalError le( lv.second );
	    
	    // Create a persistent edm::Ref to the cluster
	    edm::Ref< edm::DetSetVector<SiPixelCluster>, SiPixelCluster > cluster =
	      edm::makeRefTo( inputhandle, detid, clustIt);
	    
	    // Make a RecHit and add it to the temporary OwnVector.
	    //recHitsOnDetUnit.push_back( new SiPixelRecHit( lp, le, detIdObject, &*clustIt) );
	    recHitsOnDetUnit.push_back( new SiPixelRecHit( lp, le, detIdObject, cluster) );
	  } //  <-- End loop on Clusters
	
	if ( recHitsOnDetUnit.size()>0 ) 
	  {
	    output.put(detIdObject, 
		       recHitsOnDetUnit.begin(), recHitsOnDetUnit.end());
	    if (theVerboseLevel > 2) 
	      LogDebug("SiPixelRecHitConverter") << " Found " 
						 << recHitsOnDetUnit.size() << " RecHits on " << detid;	
	  }
	
      } //    <-- End loop on DetUnits
    
    if ( theVerboseLevel > 2 ) LogDebug ("SiPixelRecHitConverter") 
      << cpeName_ << " converted " << numberOfClusters 
      << " SiPixelClusters into SiPixelRecHits, in " 
      << numberOfDetUnits << " DetUnits."; 	   
  }
}  // end of namespace cms
