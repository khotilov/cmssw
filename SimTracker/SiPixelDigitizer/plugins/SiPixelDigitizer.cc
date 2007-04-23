// -*- C++ -*-
//
// Package:    SiPixelDigitizer
// Class:      SiPixelDigitizer
// 
/**\class SiPixelDigitizer SiPixelDigitizer.cc SimTracker/SiPixelDigitizer/src/SiPixelDigitizer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Michele Pioppi-INFN perugia
//         Created:  Mon Sep 26 11:08:32 CEST 2005
// $Id: SiPixelDigitizer.cc,v 1.26 2007/03/09 08:12:40 dkotlins Exp $
//
//


// system include files
#include <memory>
// user include files
#include "SimTracker/SiPixelDigitizer/interface/SiPixelDigitizer.h"
#include "SimTracker/SiPixelDigitizer/interface/SiPixelDigitizerAlgorithm.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
//using namespace std;


namespace cms
{
  SiPixelDigitizer::SiPixelDigitizer(const edm::ParameterSet& iConfig):
    conf_(iConfig),
    _pixeldigialgo(iConfig) 
  {
    edm::LogInfo ("PixelDigitizer ") <<"Enter the Pixel Digitizer";
    
    std::string alias ( iConfig.getParameter<std::string>("@module_label") ); 

    produces<edm::DetSetVector<PixelDigi> >().setBranchAlias( alias );
    produces<edm::DetSetVector<PixelDigiSimLink> >().setBranchAlias ( alias + "siPixelDigiSimLink");
    trackerContainers.clear();
    trackerContainers = iConfig.getParameter<std::vector<std::string> >("ROUList");
  }

  
  SiPixelDigitizer::~SiPixelDigitizer()
  {  edm::LogInfo ("PixelDigitizer ") <<"Destruct the Pixel Digitizer";}


  //
  // member functions
  //
  
  // ------------ method called to produce the data  ------------
  void
  SiPixelDigitizer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {
    // Step A: Get Inputs
    edm::Handle<CrossingFrame> cf;
    iEvent.getByType(cf);

    std::auto_ptr<MixCollection<PSimHit> > allPixelTrackerHits(new MixCollection<PSimHit>(cf.product(),trackerContainers));

    edm::ESHandle<TrackerGeometry> pDD;
    
    iSetup.get<TrackerDigiGeometryRecord> ().get(pDD);
 
    edm::ESHandle<MagneticField> pSetup;
    iSetup.get<IdealMagneticFieldRecord>().get(pSetup);

    //Loop on PSimHit
    SimHitMap.clear();

    MixCollection<PSimHit>::iterator isim;
    for (isim=allPixelTrackerHits->begin(); isim!= allPixelTrackerHits->end();isim++) {
      DetId detid=DetId((*isim).detUnitId());
      unsigned int subid=detid.subdetId();
      if ((subid==  PixelSubdetector::PixelBarrel) || (subid== PixelSubdetector::PixelEndcap)) {
	SimHitMap[(*isim).detUnitId()].push_back((*isim));
      }
    }

    // Step B: LOOP on PixelGeomDetUnit //
    for(TrackingGeometry::DetUnitContainer::const_iterator iu = pDD->detUnits().begin(); iu != pDD->detUnits().end(); iu ++){
      DetId idet=DetId((*iu)->geographicalId().rawId());
       unsigned int isub=idet.subdetId();
       
       
      if  ((isub==  PixelSubdetector::PixelBarrel) || (isub== PixelSubdetector::PixelEndcap)) {  
  
 
	//access to magnetic field in global coordinates
	GlobalVector bfield=pSetup->inTesla((*iu)->surface().position());
	LogDebug ("PixelDigitizer ") << "B-field(T) at "<<(*iu)->surface().position()<<"(cm): " 
				     << pSetup->inTesla((*iu)->surface().position());
	//

	edm::DetSet<PixelDigi> collector((*iu)->geographicalId().rawId());
	edm::DetSet<PixelDigiSimLink> linkcollector((*iu)->geographicalId().rawId());
	
	
 	collector.data=
 	  _pixeldigialgo.run(SimHitMap[(*iu)->geographicalId().rawId()],
 			     dynamic_cast<PixelGeomDetUnit*>((*iu)),
 			     bfield);
	if (collector.data.size()>0){
	  theDigiVector.push_back(collector);

	  //digisimlink
	  if(SimHitMap[(*iu)->geographicalId().rawId()].size()>0){
	    linkcollector.data=_pixeldigialgo.make_link();
	    if (linkcollector.data.size()>0) theDigiLinkVector.push_back(linkcollector);
       	  }
	
	}
      }

    }
    // Step C: create collection with the cache vector of DetSet 
    std::auto_ptr<edm::DetSetVector<PixelDigi> > 
      output(new edm::DetSetVector<PixelDigi>(theDigiVector) );
    std::auto_ptr<edm::DetSetVector<PixelDigiSimLink> > 
      outputlink(new edm::DetSetVector<PixelDigiSimLink>(theDigiLinkVector) );

    // Step D: write output to file 
    iEvent.put(output);
    iEvent.put(outputlink);
  }
}


