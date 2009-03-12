// File: StdHitNtuplizer.cc
// Description: see StdHitNtuplizer.h
// Authors: H. Cheung
//--------------------------------------------------------------


#include "SLHCUpgradeSimulations/Geometry/test/StdHitNtuplizer.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// DataFormats
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

// Geometry
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

// For ROOT
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

// STD
#include <memory>
#include <string>
#include <iostream>

using namespace std;
using namespace edm;
using namespace reco;

StdHitNtuplizer::StdHitNtuplizer(edm::ParameterSet const& conf) : 
  conf_(conf), 
  src_( conf.getParameter<edm::InputTag>( "src" ) ),
  rphiRecHits_( conf.getParameter<edm::InputTag>("rphiRecHits") ),
  stereoRecHits_( conf.getParameter<edm::InputTag>("stereoRecHits") ),
  matchedRecHits_( conf.getParameter<edm::InputTag>("matchedRecHits") ),
  tfile_(0), 
  pixeltree_(0), 
  striptree_(0),
  pixeltree2_(0)
{
}


StdHitNtuplizer::~StdHitNtuplizer() { }  

void StdHitNtuplizer::endJob() 
{
  std::cout << " StdHitNtuplizer::endJob" << std::endl;
  tfile_->Write();
  tfile_->Close();
}



void StdHitNtuplizer::beginJob(const edm::EventSetup& es)
{
  std::cout << " StdHitNtuplizer::beginJob" << std::endl;
  std::string outputFile = conf_.getParameter<std::string>("OutputFile");
 
  tfile_ = new TFile ( outputFile.c_str() , "RECREATE" );
  pixeltree_ = new TTree("PixelNtuple","Pixel hit analyzer ntuple");
  striptree_ = new TTree("StripNtuple","Strip hit analyzer ntuple");
  pixeltree2_ = new TTree("Pixel2Ntuple","Track Pixel hit analyzer ntuple");

  int bufsize = 64000;

  //Common Branch
  pixeltree_->Branch("evt",    &evt_,      "run/I:evtnum/I", bufsize);
  pixeltree_->Branch("pixel_recHit", &recHit_, 
    "x/F:y:xx:xy:yy:row:col:gx:gy:gz:subid/I:layer:nsimhit:hx/F:hy:tx:ty:theta:phi", bufsize);
  pixeltree2_->Branch("evt",    &evt_,      "run/I:evtnum/I", bufsize);
  pixeltree2_->Branch("pixel_recHit", &recHit_, 
    "x/F:y:xx:xy:yy:row:col:gx:gy:gz:subid/I:layer:nsimhit:hx/F:hy:tx:ty:theta:phi", bufsize);
  
  // Strip Branches 
  striptree_->Branch("evt",    &evt_,      "run/I:evtnum/I", bufsize);
  striptree_->Branch("strip_recHit", &striprecHit_,
    "x/F:y:xx:xy:yy:row:col:gx:gy:gz:subid/I:layer:nsimhit:hx/F:hy:tx:ty:theta:phi", bufsize);

  // geometry setup
  edm::ESHandle<TrackerGeometry>        geometry;

  es.get<TrackerDigiGeometryRecord>().get(geometry);

  theGeometry = &(*geometry);
}

// Functions that gets called by framework every event
void StdHitNtuplizer::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  // fastsim rechits
  //edm::Handle<SiTrackerGSRecHit2DCollection> theGSRecHits;
  //edm::InputTag hitProducer;
  //hitProducer = conf_.getParameter<edm::InputTag>("HitProducer");
  //e.getByLabel(hitProducer, theGSRecHits);

  edm::Handle<SiPixelRecHitCollection> recHitColl;
  e.getByLabel( src_, recHitColl);

  // for finding matched simhit
  TrackerHitAssociator associate( e, conf_ );

  //std::cout << " Step A: Standard RecHits found " << recHitColl->size() << std::endl;
  if(recHitColl->size() > 0) {
    //Loop over all rechits in SiPixelRecHitCollection (can also loop only over DetId)
    SiPixelRecHitCollection::const_iterator theRecHitRangeIteratorBegin = recHitColl->begin();
    SiPixelRecHitCollection::const_iterator theRecHitRangeIteratorEnd   = recHitColl->end();
    SiPixelRecHitCollection::const_iterator iterRecHit;

    std::string detname ;
    std::vector<PSimHit> matched;
    std::vector<PSimHit>::const_iterator closest_simhit;

    for ( iterRecHit = theRecHitRangeIteratorBegin; 
          iterRecHit != theRecHitRangeIteratorEnd; ++iterRecHit) {

      const DetId& detId =  iterRecHit->geographicalId();
      const GeomDet* geomDet( theGeometry->idToDet(detId) );

      unsigned int subdetId = detId.subdetId();
      int layerNumber=0;
      int ringNumber = 0;
      int stereo = 0;
      if ( subdetId == StripSubdetector::TIB) {
        detname = "TIB";
	TIBDetId tibid(detId.rawId());
	layerNumber = tibid.layer();
	stereo = tibid.stereo();
      } else if ( subdetId ==  StripSubdetector::TOB ) {
        detname = "TOB";
	TOBDetId tobid(detId.rawId());
	layerNumber = tobid.layer();
	stereo = tobid.stereo();
      } else if ( subdetId ==  StripSubdetector::TID) {
        detname = "TID";
	TIDDetId tidid(detId.rawId());
	layerNumber = tidid.wheel();
	ringNumber = tidid.ring();
	stereo = tidid.stereo();
      } else if ( subdetId ==  StripSubdetector::TEC ) {
        detname = "TEC";
	TECDetId tecid(detId.rawId());
	layerNumber = tecid.wheel();
	ringNumber = tecid.ring();
	stereo = tecid.stereo();
      } else if ( subdetId ==  PixelSubdetector::PixelBarrel ) {
        detname = "PXB";
	PXBDetId pxbid(detId.rawId());
	layerNumber = pxbid.layer();
	stereo = 1;
      } else if ( subdetId ==  PixelSubdetector::PixelEndcap ) {
        detname = "PXF";
	PXFDetId pxfid(detId.rawId());
	layerNumber = pxfid.disk();
	stereo = 1;
      }
      // get matched simhit
      matched.clear();
      matched = associate.associateHit(*iterRecHit);
      if ( !matched.empty() ) {
        float closest = 9999.9;
        std::vector<PSimHit>::const_iterator closestit = matched.begin();
        LocalPoint lp = iterRecHit->localPosition();
        float rechit_x = lp.x();
        float rechit_y = lp.y();
        //loop over simhits and find closest
        for (std::vector<PSimHit>::const_iterator m = matched.begin(); m<matched.end(); m++) 
        {
          float sim_x1 = (*m).entryPoint().x();
          float sim_x2 = (*m).exitPoint().x();
          float sim_xpos = 0.5*(sim_x1+sim_x2);
          float sim_y1 = (*m).entryPoint().y();
          float sim_y2 = (*m).exitPoint().y();
          float sim_ypos = 0.5*(sim_y1+sim_y2);
            
          float x_res = fabs(sim_xpos - rechit_x);
          float y_res = fabs(sim_ypos - rechit_y);
          float dist = sqrt(x_res*x_res + y_res*y_res);
          if ( dist < closest ) {
                closest = dist;
                closestit = m;
          }
        } // end of simhit loop
        closest_simhit = closestit;
      } // end matched emtpy
/////comment out begin
//      std::cout << "Found SiPixelRecHit in " << detname << " from detid " << detId.rawId()
//		<< " subdet = " << subdetId
//		<< " layer = " << layerNumber
//		<< " Stereo = " << stereo
//		<< std::endl;
//      std::cout << "Rechit global x/y/z/r : "
//                 << geomDet->surface().toGlobal(iterRecHit->localPosition()).x() << " " 
//                 << geomDet->surface().toGlobal(iterRecHit->localPosition()).y() << " " 
//                 << geomDet->surface().toGlobal(iterRecHit->localPosition()).z() << " " 
//                 << geomDet->surface().toGlobal(iterRecHit->localPosition()).perp() << std::endl;
//comment out end
      unsigned int subid = detId.subdetId();
      int layer_num = 0;
      if ( (subid==1)||(subid==2) ) {
        // 1 = PXB, 2 = PXF
        if ( subid ==  PixelSubdetector::PixelBarrel ) {
	  PXBDetId pxbid(detId.rawId());
	  layer_num   = pxbid.layer();
        } else if ( subid ==  PixelSubdetector::PixelEndcap ) {
	  PXFDetId pxfid(detId.rawId());
	  layer_num   = pxfid.disk();
        }
        int num_simhit = matched.size();
        fillPRecHit(subid, layer_num, iterRecHit, num_simhit, closest_simhit, geomDet);
        fillEvt(e);
        pixeltree_->Fill();
        init();
/*
        LocalPoint lp = iterRecHit->localPosition();
        LocalError le = iterRecHit->localPositionError();
        std::cout << "Found SiPixelRecHit in " << detname << " from detid " << detId.rawId()
		<< " subdet = " << subdetId
		<< " layer = " << layerNumber
                << "global x/y/z/r = "
                 << geomDet->surface().toGlobal(lp).x() << " " 
                 << geomDet->surface().toGlobal(lp).y() << " " 
                 << geomDet->surface().toGlobal(lp).z() << " " 
                 << geomDet->surface().toGlobal(lp).perp() 
                << " err x/y = " << sqrt(le.xx()) << " " << sqrt(le.yy()) << std::endl;
*/
        //std::cout << "   lp x,y = " << lp.x() << " " << lp.y() << " lpe xx,xy,yy = "
        //          << le.xx() << " " << le.xy() << " " << le.yy() << std::endl;
      }
    } // end of rechit loop
  } // end of loop test on recHitColl size

// Now loop over recotracks

  edm::Handle<View<reco::Track> >  trackCollection;
  edm::InputTag trackProducer;
  trackProducer = conf_.getParameter<edm::InputTag>("trackProducer");
  e.getByLabel(trackProducer, trackCollection);

/*
  std::cout << " num of reco::Tracks with "
            << trackProducer.process()<<":"
            << trackProducer.label()<<":"
            << trackProducer.instance()
            << ": " << trackCollection->size() << "\n";
*/
  int rT = 0;
  for(View<reco::Track>::size_type i=0; i<trackCollection->size(); ++i){
      ++rT;
      RefToBase<reco::Track> track(trackCollection, i);
//      std::cout << " num of hits for track " << rT << " = " << track->recHitsSize() << std::endl;
      for(trackingRecHit_iterator ih=track->recHitsBegin(); ih != track->recHitsEnd(); ++ih) {
        TrackingRecHit * hit = (*ih)->clone();
        const DetId& detId =  hit->geographicalId();
        const GeomDet* geomDet( theGeometry->idToDet(detId) );

        unsigned int subdetId = detId.subdetId();
        int layerNumber=0;
        int ringNumber = 0;
        int stereo = 0;
        std::string detname;
        if ( subdetId == StripSubdetector::TIB) {
          detname = "TIB";
          TIBDetId tibid(detId.rawId());
          layerNumber = tibid.layer();
          stereo = tibid.stereo();
        } else if ( subdetId ==  StripSubdetector::TOB ) {
          detname = "TOB";
	  TOBDetId tobid(detId.rawId());
          layerNumber = tobid.layer();
          stereo = tobid.stereo();
        } else if ( subdetId ==  StripSubdetector::TID) {
          detname = "TID";
          TIDDetId tidid(detId.rawId());
          layerNumber = tidid.wheel();
          ringNumber = tidid.ring();
          stereo = tidid.stereo();
        } else if ( subdetId ==  StripSubdetector::TEC ) {
          detname = "TEC";
          TECDetId tecid(detId.rawId());
          layerNumber = tecid.wheel();
          ringNumber = tecid.ring();
          stereo = tecid.stereo();
        } else if ( subdetId ==  PixelSubdetector::PixelBarrel ) {
          detname = "PXB";
          PXBDetId pxbid(detId.rawId());
          layerNumber = pxbid.layer();
          stereo = 1;
        } else if ( subdetId ==  PixelSubdetector::PixelEndcap ) {
          detname = "PXF";
          PXFDetId pxfid(detId.rawId());
          layerNumber = pxfid.disk();
          stereo = 1;
        }
/////comment out begin
//        std::cout << "RecHit in " << detname << " from detid " << detId.rawId()
//                  << " subdet = " << subdetId
//                  << " layer = " << layerNumber
//                  << " Stereo = " << stereo
//                  << std::endl;
        if(hit->isValid()) {
          unsigned int subid = detId.subdetId();
          if ( (subid==1)||(subid==2) ) {
            // 1 = PXB, 2 = PXF
            fillPRecHit(subid, ih, geomDet);
            fillEvt(e);
            pixeltree2_->Fill();
            init();
/*
            TrackingRecHit * hit = (*ih)->clone();
            LocalPoint lp = hit->localPosition();
            LocalError le = hit->localPositionError();
//            std::cout << "   lp x,y = " << lp.x() << " " << lp.y() << " lpe xx,xy,yy = "
//                  << le.xx() << " " << le.xy() << " " << le.yy() << std::endl;
            std::cout << "Found RecHit in " << detname << " from detid " << detId.rawId()
		<< " subdet = " << subdetId
		<< " layer = " << layerNumber
                << "global x/y/z/r = "
                 << geomDet->surface().toGlobal(lp).x() << " " 
                 << geomDet->surface().toGlobal(lp).y() << " " 
                 << geomDet->surface().toGlobal(lp).z() << " " 
                 << geomDet->surface().toGlobal(lp).perp() 
                << " err x/y = " << sqrt(le.xx()) << " " << sqrt(le.yy()) << std::endl;
*/
          }
        }
      } //end of loop on tracking rechits
  } // end of loop on recotracks

  // now for strip rechits
  edm::Handle<SiStripRecHit2DCollection> rechitsrphi;
  edm::Handle<SiStripRecHit2DCollection> rechitsstereo;
  edm::Handle<SiStripMatchedRecHit2DCollection> rechitsmatched;
  e.getByLabel(rphiRecHits_, rechitsrphi);
  e.getByLabel(stereoRecHits_, rechitsstereo);
  e.getByLabel(matchedRecHits_, rechitsmatched);

  //std::cout << " Step A: Standard Strip RPHI RecHits found " << rechitsrphi->size() << std::endl;
  //std::cout << " Step A: Standard Strip Stereo RecHits found " << rechitsstereo->size() << std::endl;
  //std::cout << " Step A: Standard Strip Matched RecHits found " << rechitsmatched->size() << std::endl;
  if(rechitsrphi->size() > 0) {
    //Loop over all rechits in RPHI collection (can also loop only over DetId)
    SiStripRecHit2DCollection::const_iterator theRecHitRangeIteratorBegin = rechitsrphi->begin();
    SiStripRecHit2DCollection::const_iterator theRecHitRangeIteratorEnd   = rechitsrphi->end();
    SiStripRecHit2DCollection::const_iterator iterRecHit;

    std::string detname ;

    for ( iterRecHit = theRecHitRangeIteratorBegin; 
          iterRecHit != theRecHitRangeIteratorEnd; ++iterRecHit) {

      const DetId& detId =  iterRecHit->geographicalId();
      const GeomDet* geomDet( theGeometry->idToDet(detId) );

      unsigned int subdetId = detId.subdetId();
      int layerNumber=0;
      int ringNumber = 0;
      int stereo = 0;
      if ( subdetId == StripSubdetector::TIB) {
        detname = "TIB";
	TIBDetId tibid(detId.rawId());
	layerNumber = tibid.layer();
	stereo = tibid.stereo();
      } else if ( subdetId ==  StripSubdetector::TOB ) {
        detname = "TOB";
	TOBDetId tobid(detId.rawId());
	layerNumber = tobid.layer();
	stereo = tobid.stereo();
      } else if ( subdetId ==  StripSubdetector::TID) {
        detname = "TID";
	TIDDetId tidid(detId.rawId());
	layerNumber = tidid.wheel();
	ringNumber = tidid.ring();
	stereo = tidid.stereo();
      } else if ( subdetId ==  StripSubdetector::TEC ) {
        detname = "TEC";
	TECDetId tecid(detId.rawId());
	layerNumber = tecid.wheel();
	ringNumber = tecid.ring();
	stereo = tecid.stereo();
      } else if ( subdetId ==  PixelSubdetector::PixelBarrel ) {
        detname = "PXB";
	PXBDetId pxbid(detId.rawId());
	layerNumber = pxbid.layer();
	stereo = 1;
      } else if ( subdetId ==  PixelSubdetector::PixelEndcap ) {
        detname = "PXF";
	PXFDetId pxfid(detId.rawId());
	layerNumber = pxfid.disk();
	stereo = 1;
      }
/////comment out begin
//      std::cout << "Found SiStripRPhiRecHit in " << detname << " from detid " << detId.rawId()
//                << " subdet = " << subdetId
//                << " layer = " << layerNumber
//                << " Stereo = " << stereo
//                << std::endl;
//      std::cout << "Rechit global x/y/z/r : "
//                << geomDet->surface().toGlobal(iterRecHit->localPosition()).x() << " " 
//                << geomDet->surface().toGlobal(iterRecHit->localPosition()).y() << " " 
//                << geomDet->surface().toGlobal(iterRecHit->localPosition()).z() << " " 
//                << geomDet->surface().toGlobal(iterRecHit->localPosition()).perp() << std::endl;
//comment out end
      unsigned int subid = detId.subdetId();
      fillSRecHit(subid, iterRecHit, geomDet);
      fillEvt(e);
      striptree_->Fill();
      init();
    } // end of rechit loop
  } // end of loop test on rechit size

  // now stereo hits
  if(rechitsstereo->size() > 0) {
    //Loop over all rechits in RPHI collection (can also loop only over DetId)
    SiStripRecHit2DCollection::const_iterator theRecHitRangeIteratorBegin = rechitsstereo->begin();
    SiStripRecHit2DCollection::const_iterator theRecHitRangeIteratorEnd   = rechitsstereo->end();
    SiStripRecHit2DCollection::const_iterator iterRecHit;

    std::string detname ;

    for ( iterRecHit = theRecHitRangeIteratorBegin; 
          iterRecHit != theRecHitRangeIteratorEnd; ++iterRecHit) {

      const DetId& detId =  iterRecHit->geographicalId();
      const GeomDet* geomDet( theGeometry->idToDet(detId) );

      unsigned int subdetId = detId.subdetId();
      int layerNumber=0;
      int ringNumber = 0;
      int stereo = 0;
      if ( subdetId == StripSubdetector::TIB) {
        detname = "TIB";
	TIBDetId tibid(detId.rawId());
	layerNumber = tibid.layer();
	stereo = tibid.stereo();
      } else if ( subdetId ==  StripSubdetector::TOB ) {
        detname = "TOB";
	TOBDetId tobid(detId.rawId());
	layerNumber = tobid.layer();
	stereo = tobid.stereo();
      } else if ( subdetId ==  StripSubdetector::TID) {
        detname = "TID";
	TIDDetId tidid(detId.rawId());
	layerNumber = tidid.wheel();
	ringNumber = tidid.ring();
	stereo = tidid.stereo();
      } else if ( subdetId ==  StripSubdetector::TEC ) {
        detname = "TEC";
	TECDetId tecid(detId.rawId());
	layerNumber = tecid.wheel();
	ringNumber = tecid.ring();
	stereo = tecid.stereo();
      } else if ( subdetId ==  PixelSubdetector::PixelBarrel ) {
        detname = "PXB";
	PXBDetId pxbid(detId.rawId());
	layerNumber = pxbid.layer();
	stereo = 1;
      } else if ( subdetId ==  PixelSubdetector::PixelEndcap ) {
        detname = "PXF";
	PXFDetId pxfid(detId.rawId());
	layerNumber = pxfid.disk();
	stereo = 1;
      }
/////comment out begin
//      std::cout << "Found SiStripStereoRecHit in " << detname << " from detid " << detId.rawId()
//                << " subdet = " << subdetId
//                << " layer = " << layerNumber
//                << " Stereo = " << stereo
//                << std::endl;
//      std::cout << "Rechit global x/y/z/r : "
//                << geomDet->surface().toGlobal(iterRecHit->localPosition()).x() << " " 
//                << geomDet->surface().toGlobal(iterRecHit->localPosition()).y() << " " 
//                << geomDet->surface().toGlobal(iterRecHit->localPosition()).z() << " " 
//                << geomDet->surface().toGlobal(iterRecHit->localPosition()).perp() << std::endl;
//comment out end
      unsigned int subid = detId.subdetId();
      fillSRecHit(subid, iterRecHit, geomDet);
      fillEvt(e);
      striptree_->Fill();
      init();
    } // end of rechit loop
  } // end of loop test on rechit size
            
  // now matched hits
  if(rechitsmatched->size() > 0) {
    //Loop over all rechits in RPHI collection (can also loop only over DetId)
    SiStripMatchedRecHit2DCollection::const_iterator theRecHitRangeIteratorBegin = rechitsmatched->begin();
    SiStripMatchedRecHit2DCollection::const_iterator theRecHitRangeIteratorEnd   = rechitsmatched->end();
    SiStripMatchedRecHit2DCollection::const_iterator iterRecHit;

    std::string detname ;

    for ( iterRecHit = theRecHitRangeIteratorBegin; 
          iterRecHit != theRecHitRangeIteratorEnd; ++iterRecHit) {

      const DetId& detId =  iterRecHit->geographicalId();
      const GeomDet* geomDet( theGeometry->idToDet(detId) );

      unsigned int subdetId = detId.subdetId();
      int layerNumber=0;
      int ringNumber = 0;
      int stereo = 0;
      if ( subdetId == StripSubdetector::TIB) {
        detname = "TIB";
	TIBDetId tibid(detId.rawId());
	layerNumber = tibid.layer();
	stereo = tibid.stereo();
      } else if ( subdetId ==  StripSubdetector::TOB ) {
        detname = "TOB";
	TOBDetId tobid(detId.rawId());
	layerNumber = tobid.layer();
	stereo = tobid.stereo();
      } else if ( subdetId ==  StripSubdetector::TID) {
        detname = "TID";
	TIDDetId tidid(detId.rawId());
	layerNumber = tidid.wheel();
	ringNumber = tidid.ring();
	stereo = tidid.stereo();
      } else if ( subdetId ==  StripSubdetector::TEC ) {
        detname = "TEC";
	TECDetId tecid(detId.rawId());
	layerNumber = tecid.wheel();
	ringNumber = tecid.ring();
	stereo = tecid.stereo();
      } else if ( subdetId ==  PixelSubdetector::PixelBarrel ) {
        detname = "PXB";
	PXBDetId pxbid(detId.rawId());
	layerNumber = pxbid.layer();
	stereo = 1;
      } else if ( subdetId ==  PixelSubdetector::PixelEndcap ) {
        detname = "PXF";
	PXFDetId pxfid(detId.rawId());
	layerNumber = pxfid.disk();
	stereo = 1;
      }
/////comment out begin
//      std::cout << "Found SiStripMatchedRecHit in " << detname << " from detid " << detId.rawId()
//                << " subdet = " << subdetId
//                << " layer = " << layerNumber
//                << " Stereo = " << stereo
//                << std::endl;
//      std::cout << "Rechit global x/y/z/r : "
//                << geomDet->surface().toGlobal(iterRecHit->localPosition()).x() << " " 
//                << geomDet->surface().toGlobal(iterRecHit->localPosition()).y() << " " 
//                << geomDet->surface().toGlobal(iterRecHit->localPosition()).z() << " " 
//                << geomDet->surface().toGlobal(iterRecHit->localPosition()).perp() << std::endl;
//comment out end
      unsigned int subid = detId.subdetId();
      fillSRecHit(subid, iterRecHit, geomDet);
      fillEvt(e);
      striptree_->Fill();
      init();
    } // end of rechit loop
  } // end of loop test on rechit size
            
} // end analyze function

void StdHitNtuplizer::fillSRecHit(const int subid, 
                                   SiStripRecHit2DCollection::const_iterator pixeliter,
                                   const GeomDet* theGeom)
{
  LocalPoint lp = pixeliter->localPosition();
  LocalError le = pixeliter->localPositionError();

  striprecHit_.x = lp.x();
  striprecHit_.y = lp.y();
  striprecHit_.xx = le.xx();
  striprecHit_.xy = le.xy();
  striprecHit_.yy = le.yy();
  GlobalPoint GP = theGeom->surface().toGlobal(pixeliter->localPosition());
  striprecHit_.gx = GP.x();
  striprecHit_.gy = GP.y();
  striprecHit_.gz = GP.z();
  striprecHit_.subid = subid;
}
void StdHitNtuplizer::fillSRecHit(const int subid, 
                                   SiStripMatchedRecHit2DCollection::const_iterator pixeliter,
                                   const GeomDet* theGeom)
{
  LocalPoint lp = pixeliter->localPosition();
  LocalError le = pixeliter->localPositionError();

  striprecHit_.x = lp.x();
  striprecHit_.y = lp.y();
  striprecHit_.xx = le.xx();
  striprecHit_.xy = le.xy();
  striprecHit_.yy = le.yy();
  GlobalPoint GP = theGeom->surface().toGlobal(pixeliter->localPosition());
  striprecHit_.gx = GP.x();
  striprecHit_.gy = GP.y();
  striprecHit_.gz = GP.z();
  striprecHit_.subid = subid;
}
void StdHitNtuplizer::fillSRecHit(const int subid, 
                                   SiTrackerGSRecHit2DCollection::const_iterator pixeliter,
                                   const GeomDet* theGeom)
{
  LocalPoint lp = pixeliter->localPosition();
  LocalError le = pixeliter->localPositionError();

  striprecHit_.x = lp.x();
  striprecHit_.y = lp.y();
  striprecHit_.xx = le.xx();
  striprecHit_.xy = le.xy();
  striprecHit_.yy = le.yy();
  //MeasurementPoint mp = topol->measurementPosition(LocalPoint(striprecHit_.x, striprecHit_.y));
  //striprecHit_.row = mp.x();
  //striprecHit_.col = mp.y();
  GlobalPoint GP = theGeom->surface().toGlobal(pixeliter->localPosition());
  striprecHit_.gx = GP.x();
  striprecHit_.gy = GP.y();
  striprecHit_.gz = GP.z();
  striprecHit_.subid = subid;
}
void StdHitNtuplizer::fillPRecHit(const int subid, 
                                  const int layer_num,
                                  SiPixelRecHitCollection::const_iterator pixeliter,
                                  const int num_simhit,
                                  std::vector<PSimHit>::const_iterator closest_simhit,
                                  const GeomDet* PixGeom)
{
  LocalPoint lp = pixeliter->localPosition();
  LocalError le = pixeliter->localPositionError();

  recHit_.x = lp.x();
  recHit_.y = lp.y();
  recHit_.xx = le.xx();
  recHit_.xy = le.xy();
  recHit_.yy = le.yy();
  //MeasurementPoint mp = topol->measurementPosition(LocalPoint(recHit_.x, recHit_.y));
  //recHit_.row = mp.x();
  //recHit_.col = mp.y();
  GlobalPoint GP = PixGeom->surface().toGlobal(pixeliter->localPosition());
  recHit_.gx = GP.x();
  recHit_.gy = GP.y();
  recHit_.gz = GP.z();
  recHit_.subid = subid;
  recHit_.layer = layer_num;
  recHit_.nsimhit = num_simhit;
  //std::cout << "num_simhit = " << num_simhit << std::endl;
  if(num_simhit > 0) {
    float sim_x1 = (*closest_simhit).entryPoint().x();
    float sim_x2 = (*closest_simhit).exitPoint().x();
    recHit_.hx = 0.5*(sim_x1+sim_x2);
    float sim_y1 = (*closest_simhit).entryPoint().y();
    float sim_y2 = (*closest_simhit).exitPoint().y();
    recHit_.hy = 0.5*(sim_y1+sim_y2);
    //std::cout << "num_simhit x, y = " << 0.5*(sim_x1+sim_x2) << " " << 0.5*(sim_y1+sim_y2) << std::endl;
  }
/*
       std::cout << "Found RecHit in " << subid
                 << " global x/y/z : "
                 << PixGeom->surface().toGlobal(pixeliter->localPosition()).x() << " " 
                 << PixGeom->surface().toGlobal(pixeliter->localPosition()).y() << " " 
                 << PixGeom->surface().toGlobal(pixeliter->localPosition()).z() << std::endl;
*/
}
void StdHitNtuplizer::fillPRecHit(const int subid, 
                                  trackingRecHit_iterator ih,
                                  const GeomDet* PixGeom)
{
  TrackingRecHit * pixeliter = (*ih)->clone();
  LocalPoint lp = pixeliter->localPosition();
  LocalError le = pixeliter->localPositionError();

  recHit_.x = lp.x();
  recHit_.y = lp.y();
  recHit_.xx = le.xx();
  recHit_.xy = le.xy();
  recHit_.yy = le.yy();
  GlobalPoint GP = PixGeom->surface().toGlobal(pixeliter->localPosition());
  recHit_.gx = GP.x();
  recHit_.gy = GP.y();
  recHit_.gz = GP.z();
  recHit_.subid = subid;
}

void
StdHitNtuplizer::fillEvt(const edm::Event& E)
{
   evt_.run = E.id().run();
   evt_.evtnum = E.id().event();
}

void StdHitNtuplizer::init()
{
  evt_.init();
  recHit_.init();
  striprecHit_.init();
}

void StdHitNtuplizer::evt::init()
{
  int dummy_int = 9999;
  run = dummy_int;
  evtnum = dummy_int;
}

void StdHitNtuplizer::RecHit::init()
{
  float dummy_float = 9999.0;

  x = dummy_float;
  y = dummy_float;
  xx = dummy_float;
  xy = dummy_float; 
  yy = dummy_float;
  row = dummy_float;
  col = dummy_float;
  gx = dummy_float;
  gy = dummy_float;
  gz = dummy_float;
  layer = 0;
  nsimhit = 0;
  hx = dummy_float;
  hy = dummy_float;
  tx = dummy_float;
  ty = dummy_float;
  theta = dummy_float;
  phi = dummy_float;
}

//define this as a plug-in
DEFINE_ANOTHER_FWK_MODULE(StdHitNtuplizer);

