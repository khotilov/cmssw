/** SiPixelRecHitConverter.cc
 * ------------------------------------------------------
 * Description:  see SiPixelRecHitConverter.h
 * Authors:  P. Maksimovic (JHU), V.Chiochia (Uni Zurich)
 * History: Feb 27, 2006 -  initial version
 *          May 30, 2006 -  edm::DetSetVector and edm::Ref
 *          Aug 30, 2007 -  edmNew::DetSetVector
*			Jan 31, 2008 -  change to use Lorentz angle from DB (Lotte Wilke)
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
#include "DataFormats/Common/interface/DetSet2RangeMap.h"

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

// Error Parametrization DB Record
#include "CondFormats/DataRecord/interface/SiPixelCPEParmErrorsRcd.h"

// Lorentz Angle Record
#include "CondFormats/DataRecord/interface/SiPixelLorentzAngleRcd.h"


#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

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
    src_( conf.getParameter<edm::InputTag>( "src" ) ),
    theVerboseLevel(conf.getUntrackedParameter<int>("VerboseLevel",0)),
    m_newCont(conf.getUntrackedParameter<bool>("newContainer",false)),
    errorsFromDB_(conf.getUntrackedParameter<bool>("errorsFromDB",false))
  {
    //--- Declare to the EDM what kind of collections we will be making.
    produces<SiPixelRecHitCollection>();
    if (m_newCont)
      produces<SiPixelRecHitCollectionNew>();
   
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

	edm::ESHandle<SiPixelLorentzAngle> lorentzAngle;
    c.get<SiPixelLorentzAngleRcd>().get(lorentzAngle);
	
    edm::ESHandle<SiPixelCPEParmErrors> parmErrors;

    if(errorsFromDB_) {
      c.get<SiPixelCPEParmErrorsRcd>().get(parmErrors);
    }
	
    setupCPE(magfield.product(),parmErrors.product(),lorentzAngle.product());

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
    std::auto_ptr<SiPixelRecHitCollectionNew> output(new SiPixelRecHitCollectionNew);

    // Step C: Iterate over DetIds and invoke the strip CPE algorithm
    // on each DetUnit

    run( input, *output, geom );


    // Step Z: temporary write also the old collection
    std::auto_ptr<SiPixelRecHitCollection> old(new SiPixelRecHitCollection);
    edmNew::copy(*output,*old);
    e.put(old);
    
    
    // Step D: write output to file
    if (m_newCont)
      e.put(output);

  }

  //---------------------------------------------------------------------------
  //!  Set up the specific algorithm we are going to use.  
  //---------------------------------------------------------------------------
  void SiPixelRecHitConverter::setupCPE(const MagneticField* mag, const SiPixelCPEParmErrors* parmErrors, const SiPixelLorentzAngle* lorentzAngle) 
  {
    cpeName_ = conf_.getParameter<std::string>("CPE");
    if ( cpeName_ == "FromDetPosition" ) 
      {
	cpe_ = new CPEFromDetPosition(conf_,mag);
	ready_ = true;
      } 
    else if ( cpeName_ == "Initial" ) 
      {
	cpe_ = new PixelCPEInitial(conf_,mag,lorentzAngle);
	ready_ = true;
      } 
    else if ( cpeName_ == "ParmError" ) 
      {
	cpe_ = new PixelCPEParmError(conf_,mag,lorentzAngle);
	ready_ = true;
      } 
    else if ( cpeName_ == "TemplateReco" ) 
      {
	cpe_ = new PixelCPETemplateReco(conf_,mag);
	ready_ = true;
      } 
    else if ( cpeName_ == "Generic" ) 
      {
	cpe_ = new PixelCPEGeneric(conf_,mag,parmErrors);
	ready_ = true;
      } 
    else 
      {
	edm::LogError("SiPixelRecHitConverter") 
	  <<" Cluster parameter estimator " << cpeName_ << " is invalid.\n"
	  << "Possible choices:\n" 
	  << "    - FromDetPosition  (straight from ORCA)\n"
	  << "    - Initial          (ORCA algorithm fixed for bix pixels\n"
	  << "    - ParmError        (copy of FromTrackAngles from ORCA)\n"
	  << "    - TemplateReco     (fits to templates based on PIXELAV)\n"
	  << "    - Generic          (Initial rewritten for clarity and speed)\n";
	ready_ = false;
      }
    // &&& We should really throw a fatal error here!
  }
  //---------------------------------------------------------------------------
  //!  Iterate over DetUnits, then over Clusters and invoke the CPE on each,
  //!  and make a RecHit to store the result.
  //!  New interface reading DetSetVector by V.Chiochia (May 30th, 2006)
  //---------------------------------------------------------------------------
	void SiPixelRecHitConverter::run(edm::Handle<edm::DetSetVector<SiPixelCluster> >  inputhandle,
					SiPixelRecHitCollectionNew &output,
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
			SiPixelRecHitCollectionNew::FastFiller recHitsOnDetUnit(output,detid);
			
			edm::DetSet<SiPixelCluster>::const_iterator clustIt = DSViter->data.begin();

			for ( ; clustIt != DSViter->data.end(); clustIt++) 
			{
				numberOfClusters++;
				std::pair<LocalPoint, LocalError> lv = cpe_->localParameters( *clustIt, *genericDet );
				LocalPoint lp( lv.first );
				LocalError le( lv.second );
				// Create a persistent edm::Ref to the cluster
				edm::Ref< edm::DetSetVector<SiPixelCluster>, SiPixelCluster > cluster = edm::makeRefTo( inputhandle, detid, clustIt);
				// Make a RecHit and add it to the DetSet
				// old : recHitsOnDetUnit.push_back( new SiPixelRecHit( lp, le, detIdObject, &*clustIt) );
				SiPixelRecHit hit( lp, le, detIdObject, cluster);
				#ifdef SIPIXELRECHIT_HAS_EXTRA_INFO
				// Copy the extra stuff; unfortunately, only the derivatives of PixelCPEBase
				// are capable of doing that.  So until we get rid of CPEFromDetPosition
				// we'll have to dynamic_cast :(
				// &&& This cast can be moved to the setupCPE, so that it's done once per job.
				PixelCPEBase * cpeBase = dynamic_cast< PixelCPEBase* >( cpe_ );
				if (cpeBase) {
					hit.setProbabilityX( cpeBase->probabilityX() );
					hit.setProbabilityY( cpeBase->probabilityY() );
					hit.setQBin( cpeBase->qBin() );
					hit.setCotAlphaFromCluster( cpeBase->cotAlphaFromCluster() );
					hit.setCotBetaFromCluster ( cpeBase->cotBetaFromCluster()  );
				}
				#endif
				// 
				// Now save it =================
				recHitsOnDetUnit.push_back(hit);
				// =============================
			} //  <-- End loop on Clusters
			
			if ( recHitsOnDetUnit.size()>0 ) 
			{
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
