
// Package:          SiPixelErrorEstimation
// Class:            SiPixelErrorEstimation
// Original Author:  Gavril Giurgiu (JHU)
// Created:          Fri May  4 17:48:24 CDT 2007

#include <memory>
#include <string>
#include <iostream>
#include <TMath.h>
#include "CalibTracker/SiPixelErrorEstimation/interface/SiPixelErrorEstimation.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/DetId/interface/DetId.h" 
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include <TTree.h>
#include <TFile.h>

using namespace std;
using namespace edm;

SiPixelErrorEstimation::SiPixelErrorEstimation(const edm::ParameterSet& ps):tfile_(0), ttree_all_hits_(0), ttree_track_hits_(0) 
{
  //Read config file
  outputFile_ = ps.getUntrackedParameter<string>( "outputFile", "SiPixelErrorEstimation_Ntuple.root" );
  src_ = ps.getUntrackedParameter<std::string>( "src", "ctfWithMaterialTracks" );
  checkType_ = ps.getParameter<bool>( "checkType" );
  genType_ = ps.getParameter<int>( "genType" );
  include_trk_hits_ = ps.getParameter<bool>( "include_trk_hits" );
}

SiPixelErrorEstimation::~SiPixelErrorEstimation()
{}

void SiPixelErrorEstimation::beginJob(const edm::EventSetup& es)
{
  int bufsize = 64000;

  if ( include_trk_hits_ )
    {
      //tfile_ = new TFile ("SiPixelErrorEstimation_Ntuple.root" , "RECREATE");
      //const char* tmp_name = outputFile_.c_str();
      tfile_ = new TFile ( outputFile_.c_str() , "RECREATE");
      
      ttree_track_hits_ = new TTree("TrackHitNtuple", "TrackHitNtuple");
      
      ttree_track_hits_->Branch("evt", &evt, "evt/I", bufsize);
      ttree_track_hits_->Branch("run", &run, "run/I", bufsize);
      
      ttree_track_hits_->Branch("subdetId", &subdetId, "subdetId/I", bufsize);
      
      ttree_track_hits_->Branch("layer" , &layer , "layer/I" , bufsize);
      ttree_track_hits_->Branch("ladder", &ladder, "ladder/I", bufsize);
      ttree_track_hits_->Branch("mod"   , &mod   , "mod/I"   , bufsize);
      
      ttree_track_hits_->Branch("side"  , &side  , "side/I"  , bufsize);
      ttree_track_hits_->Branch("disk"  , &disk  , "disk/I"  , bufsize);
      ttree_track_hits_->Branch("blade" , &blade , "blade/I" , bufsize);
      ttree_track_hits_->Branch("panel" , &panel , "panel/I" , bufsize);
      ttree_track_hits_->Branch("plaq"  , &plaq  , "plaq/I"  , bufsize);
      
      ttree_track_hits_->Branch("half"   , &half   , "half/I"   , bufsize);
      ttree_track_hits_->Branch("flipped", &flipped, "flipped/I", bufsize);
      
      ttree_track_hits_->Branch("rechitx", &rechitx, "rechitx/F"    , bufsize);
      ttree_track_hits_->Branch("rechity", &rechity, "rechity/F"    , bufsize);
      ttree_track_hits_->Branch("rechitz", &rechitz, "rechitz/F"    , bufsize);
      
      ttree_track_hits_->Branch("rechiterrx", &rechiterrx, "rechiterrx/F" , bufsize);
      ttree_track_hits_->Branch("rechiterry", &rechiterry, "rechiterry/F" , bufsize);
      
      ttree_track_hits_->Branch("rechitresx", &rechitresx, "rechitresx/F" , bufsize);
      ttree_track_hits_->Branch("rechitresy", &rechitresy, "rechitresy/F" , bufsize);
      
      ttree_track_hits_->Branch("rechitpullx", &rechitpullx, "rechitpullx/F", bufsize);
      ttree_track_hits_->Branch("rechitpully", &rechitpully, "rechitpully/F", bufsize);
      
      ttree_track_hits_->Branch("npix"  , &npix  , "npix/I"  , bufsize);
      ttree_track_hits_->Branch("nxpix" , &nxpix , "nxpix/I" , bufsize);
      ttree_track_hits_->Branch("nypix" , &nypix , "nypix/I" , bufsize);
      ttree_track_hits_->Branch("charge", &charge, "charge/F", bufsize);
      
      ttree_track_hits_->Branch("edgex", &edgex, "edgex/I", bufsize);
      ttree_track_hits_->Branch("edgey", &edgey, "edgey/I", bufsize);
      
      ttree_track_hits_->Branch("bigx", &bigx, "bigx/I", bufsize);
      ttree_track_hits_->Branch("bigy", &bigy, "bigy/I", bufsize);
      
      ttree_track_hits_->Branch("alpha", &alpha, "alpha/F", bufsize);
      ttree_track_hits_->Branch("beta" , &beta , "beta/F" , bufsize);
      
      ttree_track_hits_->Branch("phi", &phi, "phi/F", bufsize);
      ttree_track_hits_->Branch("eta", &eta, "eta/F", bufsize);
      
      ttree_track_hits_->Branch("simhitx", &simhitx, "simhitx/F", bufsize);
      ttree_track_hits_->Branch("simhity", &simhity, "simhity/F", bufsize);
      
      ttree_track_hits_->Branch("nsimhit", &nsimhit, "nsimhit/I", bufsize);
      ttree_track_hits_->Branch("pidhit" , &pidhit , "pidhit/I" , bufsize);
      ttree_track_hits_->Branch("simproc", &simproc, "simproc/I", bufsize);
    } // if ( include_trk_hits_ )

  // ----------------------------------------------------------------------
  
  ttree_all_hits_ = new TTree("AllHitNtuple", "AllHitNtuple");

  ttree_all_hits_->Branch("evt", &evt, "evt/I", bufsize);
  ttree_all_hits_->Branch("run", &run, "run/I", bufsize);

  ttree_all_hits_->Branch("subdetid", &all_subdetid, "subdetid/I", bufsize);
  
  ttree_all_hits_->Branch("layer" , &all_layer , "layer/I" , bufsize);
  ttree_all_hits_->Branch("ladder", &all_ladder, "ladder/I", bufsize);
  ttree_all_hits_->Branch("mod"   , &all_mod   , "mod/I"   , bufsize);
  
  ttree_all_hits_->Branch("side"  , &all_side  , "side/I"  , bufsize);
  ttree_all_hits_->Branch("disk"  , &all_disk  , "disk/I"  , bufsize);
  ttree_all_hits_->Branch("blade" , &all_blade , "blade/I" , bufsize);
  ttree_all_hits_->Branch("panel" , &all_panel , "panel/I" , bufsize);
  ttree_all_hits_->Branch("plaq"  , &all_plaq  , "plaq/I"  , bufsize);

  ttree_all_hits_->Branch("half"   , &all_half   , "half/I"   , bufsize);
  ttree_all_hits_->Branch("flipped", &all_flipped, "flipped/I", bufsize);

  ttree_all_hits_->Branch("cols", &all_cols, "cols/I", bufsize);
  ttree_all_hits_->Branch("rows", &all_rows, "rows/I", bufsize);

  ttree_all_hits_->Branch("rechitx"    , &all_rechitx    , "rechitx/F"    , bufsize);
  ttree_all_hits_->Branch("rechity"    , &all_rechity    , "rechity/F"    , bufsize);
  ttree_all_hits_->Branch("rechitz"    , &all_rechitz    , "rechitz/F"    , bufsize);
 
  ttree_all_hits_->Branch("rechiterrx" , &all_rechiterrx , "rechiterrx/F" , bufsize);
  ttree_all_hits_->Branch("rechiterry" , &all_rechiterry , "rechiterry/F" , bufsize);
  
  ttree_all_hits_->Branch("rechitresx" , &all_rechitresx , "rechitresx/F" , bufsize);
  ttree_all_hits_->Branch("rechitresy" , &all_rechitresy , "rechitresy/F" , bufsize);
  
  ttree_all_hits_->Branch("rechitpullx", &all_rechitpullx, "rechitpullx/F", bufsize);
  ttree_all_hits_->Branch("rechitpully", &all_rechitpully, "rechitpully/F", bufsize);

  ttree_all_hits_->Branch("npix"  , &all_npix  , "npix/I"  , bufsize);
  ttree_all_hits_->Branch("nxpix" , &all_nxpix , "nxpix/I" , bufsize);
  ttree_all_hits_->Branch("nypix" , &all_nypix , "nypix/I" , bufsize);

  ttree_all_hits_->Branch("edgex", &all_edgex, "edgex/I", bufsize);
  ttree_all_hits_->Branch("edgey", &all_edgey, "edgey/I", bufsize);
 
  ttree_all_hits_->Branch("bigx", &all_bigx, "bigx/I", bufsize);
  ttree_all_hits_->Branch("bigy", &all_bigy, "bigy/I", bufsize);
  
  ttree_all_hits_->Branch("alpha", &all_alpha, "alpha/F", bufsize);
  ttree_all_hits_->Branch("beta" , &all_beta , "beta/F" , bufsize);

  ttree_all_hits_->Branch("simhitx", &all_simhitx, "simhitx/F", bufsize);
  ttree_all_hits_->Branch("simhity", &all_simhity, "simhity/F", bufsize);

  ttree_all_hits_->Branch("nsimhit", &all_nsimhit, "nsimhit/I", bufsize);
  ttree_all_hits_->Branch("pidhit" , &all_pidhit , "pidhit/I" , bufsize);
  ttree_all_hits_->Branch("simproc", &all_simproc, "simproc/I", bufsize);

  ttree_all_hits_->Branch("vtxr", &all_vtxr, "vtxr/F", bufsize);
  ttree_all_hits_->Branch("vtxz", &all_vtxz, "vtxz/F", bufsize);

  ttree_all_hits_->Branch("simpx", &all_simpx, "simpx/F", bufsize);
  ttree_all_hits_->Branch("simpy", &all_simpy, "simpy/F", bufsize);
  ttree_all_hits_->Branch("simpz", &all_simpz, "simpz/F", bufsize);

  ttree_all_hits_->Branch("eloss", &all_eloss, "eloss/F", bufsize);
  
  ttree_all_hits_->Branch("simphi", &all_simphi, "simphi/F", bufsize);
  ttree_all_hits_->Branch("simtheta", &all_simtheta, "simtheta/F", bufsize);
  
  ttree_all_hits_->Branch("trkid", &all_trkid, "trkid/I", bufsize);
  
  ttree_all_hits_->Branch("x1", &all_x1, "x1/F", bufsize);
  ttree_all_hits_->Branch("x2", &all_x2, "x2/F", bufsize);
  ttree_all_hits_->Branch("y1", &all_x1, "y1/F", bufsize);
  ttree_all_hits_->Branch("y2", &all_x2, "y2/F", bufsize);
  ttree_all_hits_->Branch("z1", &all_x1, "z1/F", bufsize);
  ttree_all_hits_->Branch("z2", &all_x2, "z2/F", bufsize);

  ttree_all_hits_->Branch("row1", &all_row1, "row1/F", bufsize);
  ttree_all_hits_->Branch("row2", &all_row2, "row2/F", bufsize);
  ttree_all_hits_->Branch("col1", &all_col1, "col1/F", bufsize);
  ttree_all_hits_->Branch("col2", &all_col2, "col2/F", bufsize);

  ttree_all_hits_->Branch("gx1", &all_gx1, "gx1/F", bufsize);
  ttree_all_hits_->Branch("gx2", &all_gx2, "gx2/F", bufsize);
  ttree_all_hits_->Branch("gy1", &all_gx1, "gy1/F", bufsize);
  ttree_all_hits_->Branch("gy2", &all_gx2, "gy2/F", bufsize);
  ttree_all_hits_->Branch("gz1", &all_gx1, "gz1/F", bufsize);
  ttree_all_hits_->Branch("gz2", &all_gx2, "gz2/F", bufsize);
  
  ttree_all_hits_->Branch("simtrketa", &all_simtrketa, "simtrketa/F", bufsize);
  ttree_all_hits_->Branch("simtrkphi", &all_simtrkphi, "simtrkphi/F", bufsize);

  ttree_all_hits_->Branch("clust_row", &all_clust_row, "clust_row/F", bufsize);
  ttree_all_hits_->Branch("clust_col", &all_clust_col, "clust_col/F", bufsize);
  
  ttree_all_hits_->Branch("clust_x", &all_clust_x, "clust_x/F", bufsize);
  ttree_all_hits_->Branch("clust_y", &all_clust_y, "clust_y/F", bufsize);

  ttree_all_hits_->Branch("clust_q", &all_clust_q, "clust_q/F", bufsize);

  ttree_all_hits_->Branch("clust_maxpixcol", &all_clust_maxpixcol, "clust_maxpixcol/I", bufsize);
  ttree_all_hits_->Branch("clust_maxpixrow", &all_clust_maxpixrow, "clust_maxpixrow/I", bufsize);
  ttree_all_hits_->Branch("clust_minpixcol", &all_clust_minpixcol, "clust_minpixcol/I", bufsize);
  ttree_all_hits_->Branch("clust_minpixrow", &all_clust_minpixrow, "clust_minpixrow/I", bufsize);
  
  ttree_all_hits_->Branch("clust_geoid", &all_clust_geoid, "clust_geoid/I", bufsize);
  
  ttree_all_hits_->Branch("clust_alpha", &all_clust_alpha, "clust_alpha/F", bufsize);
  ttree_all_hits_->Branch("clust_beta" , &all_clust_beta , "clust_beta/F" , bufsize);

  ttree_all_hits_->Branch("rowpix", all_pixrow, "row[npix]/F", bufsize);
  ttree_all_hits_->Branch("colpix", all_pixcol, "col[npix]/F", bufsize);
  ttree_all_hits_->Branch("adc", all_pixadc, "adc[npix]/F", bufsize);
  ttree_all_hits_->Branch("xpix", all_pixx, "x[npix]/F", bufsize);
  ttree_all_hits_->Branch("ypix", all_pixy, "y[npix]/F", bufsize);
  ttree_all_hits_->Branch("gxpix", all_pixgx, "gx[npix]/F", bufsize);
  ttree_all_hits_->Branch("gypix", all_pixgy, "gy[npix]/F", bufsize);
  ttree_all_hits_->Branch("gzpix", all_pixgz, "gz[npix]/F", bufsize);
  
  ttree_all_hits_->Branch("hit_probx", &all_hit_probx, "hit_probx/F" , bufsize);
  ttree_all_hits_->Branch("hit_proby", &all_hit_proby, "hit_proby/F" , bufsize);

}

void SiPixelErrorEstimation::endJob() 
{
  tfile_->Write();
  tfile_->Close();
}

void
SiPixelErrorEstimation::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  using namespace edm;
  
  run = e.id().run();
  evt = e.id().event();
  
  if ( evt%1000 == 0 ) 
    cout << "evt = " << evt << endl;

  float math_pi = 3.14159265;
  float radtodeg = 180.0 / math_pi;
    
  DetId detId;

  LocalPoint position;
  LocalError error;
  float mindist = 999999.9;

  std::vector<PSimHit> matched;
  TrackerHitAssociator associate(e);

  edm::ESHandle<TrackerGeometry> pDD;
  es.get<TrackerDigiGeometryRecord> ().get (pDD);
  const TrackerGeometry* tracker = &(* pDD);
        
  // --------------------------------------- all hits -----------------------------------------------------------
  edm::Handle<SiPixelRecHitCollection> recHitColl;
  e.getByLabel( "siPixelRecHits", recHitColl);

  Handle<edm::SimTrackContainer> simtracks;
  e.getByLabel("g4SimHits", simtracks);

  //-----Iterate over detunits
  for (TrackerGeometry::DetContainer::const_iterator it = pDD->dets().begin(); it != pDD->dets().end(); it++) 
    {
      DetId detId = ((*it)->geographicalId());
      
      SiPixelRecHitCollection::range pixelrechitRange = (recHitColl.product())->get(detId);
      SiPixelRecHitCollection::const_iterator pixelrechitRangeIteratorBegin = pixelrechitRange.first;
      SiPixelRecHitCollection::const_iterator pixelrechitRangeIteratorEnd = pixelrechitRange.second;
      SiPixelRecHitCollection::const_iterator pixeliter = pixelrechitRangeIteratorBegin;
      std::vector<PSimHit> matched;
      
      //----Loop over rechits for this detId
      for ( ; pixeliter != pixelrechitRangeIteratorEnd; ++pixeliter) 
	{
	  matched.clear();
	  matched = associate.associateHit(*pixeliter);
	  // only consider rechits that have associated simhit
	  // otherwise cannot determine residiual
	  if ( matched.empty() )
	    {
	      cout << "SiPixelErrorEstimation::analyze: rechits without associated simhit !!!!!!!" 
		   << endl;
	      continue;
	    }
		
	  all_subdetid = -9999;
	  
	  all_layer = -9999;
	  all_ladder = -9999;
	  all_mod = -9999;
	  
	  all_side = -9999;
	  all_disk = -9999;
	  all_blade = -9999;
	  all_panel = -9999;
	  all_plaq = -9999;
	  
	  all_half = -9999;
	  all_flipped = -9999;
	  
	  all_cols = -9999;
	  all_rows = -9999;
	  
	  all_rechitx = -9999;
	  all_rechity = -9999;
	  all_rechitz = -9999;
	  
	  all_simhitx = -9999;
	  all_simhity = -9999;

	  all_rechiterrx = -9999;
	  all_rechiterry = -9999;
	  
	  all_rechitresx = -9999;
	  all_rechitresy = -9999;
	  
	  all_rechitpullx = -9999;
	  all_rechitpully = -9999;
	  
	  all_npix = -9999;
	  all_nxpix = -9999;
	  all_nypix = -9999;
	 	  
	  all_edgex = -9999;
	  all_edgey = -9999;
	  
	  all_bigx = -9999;
	  all_bigy = -9999;
	  
	  all_alpha = -9999;
	  all_beta = -9999;
	  
	  all_simphi = -9999;
	  all_simtheta = -9999;
	  
	  all_simhitx = -9999;
	  all_simhity = -9999;
	  
	  all_nsimhit = -9999;
	  all_pidhit = -9999;
	  all_simproc = -9999;
	  
	  all_vtxr = -9999;
	  all_vtxz = -9999;
	  
	  all_simpx = -9999;
	  all_simpy = -9999;
	  all_simpz = -9999;
	  
	  all_eloss = -9999;
	  	  
	  all_trkid = -9999;
	  
	  all_x1 = -9999;
	  all_x2 = -9999;
	  all_y1 = -9999;
	  all_y2 = -9999;
	  all_z1 = -9999;
	  all_z2 = -9999;
	  
	  all_row1 = -9999;
	  all_row2 = -9999;
	  all_col1 = -9999;
	  all_col2 = -9999;
	  
	  all_gx1 = -9999;
	  all_gx2 = -9999;
	  all_gy1 = -9999;
	  all_gy2 = -9999;
	  all_gz1 = -9999;
	  all_gz2 = -9999;
	  
	  all_simtrketa = -9999;
	  all_simtrkphi = -9999;
	  
	  all_clust_row = -9999;
	  all_clust_col = -9999;
	  
	  all_clust_x = -9999;
	  all_clust_y = -9999;
	  
	  all_clust_q = -9999;
	  
	  all_clust_maxpixcol = -9999;
	  all_clust_maxpixrow = -9999;
	  all_clust_minpixcol = -9999;
	  all_clust_minpixrow = -9999;
	  
	  all_clust_geoid = -9999;
	  
	  all_clust_alpha = -9999;
	  all_clust_beta = -9999;
	  
	  /*
	    for (int i=0; i<all_npix; ++i)
	    {
	    all_pixrow[i] = -9999;
	    all_pixcol[i] = -9999;
	    all_pixadc[i] = -9999;
	    all_pixx[i] = -9999;
	    all_pixy[i] = -9999;
	    all_pixgx[i] = -9999;
	    all_pixgy[i] = -9999;
	    all_pixgz[i] = -9999;
	    }
	  */

	  all_hit_probx = -9999;
	  all_hit_proby = -9999;
	  
	  all_nsimhit = (int)matched.size();
	  
	  all_subdetid = (int)detId.subdetId();
	  // only consider rechits in pixel barrel and pixel forward 
	  if ( !(all_subdetid==1 || all_subdetid==2) ) 
	    {
	      cout << "SiPixelErrorEstimation::analyze: Not in a pixel detector !!!!!" << endl; 
	      continue;
	    }

	  const PixelGeomDetUnit* theGeomDet 
	    = dynamic_cast<const PixelGeomDetUnit*> ( tracker->idToDet(detId) );
	  
	  const RectangularPixelTopology* topol = 
	    dynamic_cast<const RectangularPixelTopology*>(&(theGeomDet->specificTopology()));

	  const int maxPixelCol = pixeliter->cluster()->maxPixelCol();
	  const int maxPixelRow = pixeliter->cluster()->maxPixelRow();
	  const int minPixelCol = pixeliter->cluster()->minPixelCol();
	  const int minPixelRow = pixeliter->cluster()->minPixelRow();
	  
	  // check whether the cluster is at the module edge 
	  if ( topol->isItEdgePixelInX( minPixelRow ) || 
	       topol->isItEdgePixelInX( maxPixelRow ) )
	    all_edgex = 1;
	  else 
	    all_edgex = 0;
	  
	  if ( topol->isItEdgePixelInY( minPixelCol ) || 
	       topol->isItEdgePixelInY( maxPixelCol ) )
	    all_edgey = 1;
	  else 
	    all_edgey = 0;
	  
	  // check whether this rechit contains big pixels
	  if ( topol->containsBigPixelInX(minPixelRow, maxPixelRow) )
	    all_bigx = 1;
	  else 
	    all_bigx = 0;
	  
	  if ( topol->containsBigPixelInY(minPixelCol, maxPixelCol) )
	    all_bigy = 1;
	  else 
	    all_bigy = 0;
	  
	  if ( (int)detId.subdetId() == (int)PixelSubdetector::PixelBarrel ) 
	    {
	      PXBDetId bdetid(detId);
	      all_layer = bdetid.layer();
	      all_ladder = bdetid.ladder();
	      all_mod = bdetid.module();
	      
	      int tmp_nrows = theGeomDet->specificTopology().nrows();
	      if ( tmp_nrows == 80 ) 
		all_half = 1;
	      else if ( tmp_nrows == 160 ) 
		all_half = 0;
	      else 
		cout << "-------------------------------------------------- Wrong module size !!!" << endl;
	      
	      float tmp1 = theGeomDet->surface().toGlobal(Local3DPoint(0.,0.,0.)).perp();
	      float tmp2 = theGeomDet->surface().toGlobal(Local3DPoint(0.,0.,1.)).perp();
	      
	      if ( tmp2<tmp1 ) 
		all_flipped = 1;
	      else 
		all_flipped = 0;
	    }
	  else if ( (int)detId.subdetId() == (int)PixelSubdetector::PixelEndcap )
	    {
	      PXFDetId fdetid(detId);
	      all_side  = fdetid.side();
	      all_disk  = fdetid.disk();
	      all_blade = fdetid.blade();
	      all_panel = fdetid.panel();
	      all_plaq  = fdetid.module(); // also known as plaquette
	      
	    } // else if ( detId.subdetId()==PixelSubdetector::PixelEndcap )
	  else std::cout << "We are not in the pixel detector" << (int)detId.subdetId() << endl;
	  
	  all_cols = theGeomDet->specificTopology().ncolumns();
	  all_rows = theGeomDet->specificTopology().nrows();
	  	  
	  LocalPoint lp = pixeliter->localPosition();
	  // gavril: change this name
	  all_rechitx = lp.x();
	  all_rechity = lp.y();
	  all_rechitz = lp.z();
	  
	  LocalError le = pixeliter->localPositionError();
	  all_rechiterrx = sqrt( le.xx() );
	  all_rechiterry = sqrt( le.yy() );

	  bool found_hit_from_generated_particle = false;
	  
	  //---Loop over sim hits, fill closest
	  float closest_dist = 99999.9;
	  std::vector<PSimHit>::const_iterator closest_simhit = matched.begin();
	  
	  for (std::vector<PSimHit>::const_iterator m = matched.begin(); m < matched.end(); m++) 
	    {
	      if ( checkType_ )
		{
		  int pid = (*m).particleType();
		  if ( abs(pid) != genType_ )
		    continue;
		} 
	      
	      float simhitx = 0.5 * ( (*m).entryPoint().x() + (*m).exitPoint().x() );
	      float simhity = 0.5 * ( (*m).entryPoint().y() + (*m).exitPoint().y() );
	      
	      float x_res = simhitx - rechitx;
	      float y_res = simhity - rechity;
		  
	      float dist = sqrt( x_res*x_res + y_res*y_res );		  
	      
	      if ( dist < closest_dist ) 
		{
		  closest_dist = dist;
		  closest_simhit = m;
		  found_hit_from_generated_particle = true;
		} 
	    } // end sim hit loop
	  
	  // If this recHit does not have any simHit with the same particleType as the particles generated
	  // ignore it as most probably comes from delta rays.
	  if ( checkType_ && !found_hit_from_generated_particle )
	    continue; 

	  all_x1 = (*closest_simhit).entryPoint().x(); // width (row index, in col direction)
	  all_y1 = (*closest_simhit).entryPoint().y(); // length (col index, in row direction) 
	  all_z1 = (*closest_simhit).entryPoint().z(); 
	  all_x2 = (*closest_simhit).exitPoint().x();
	  all_y2 = (*closest_simhit).exitPoint().y();
	  all_z2 = (*closest_simhit).exitPoint().z();
	  GlobalPoint GP1 = 
	    theGeomDet->surface().toGlobal( Local3DPoint( (*closest_simhit).entryPoint().x(),
							  (*closest_simhit).entryPoint().y(),
							  (*closest_simhit).entryPoint().z() ) );
	  GlobalPoint GP2 = 
	    theGeomDet->surface().toGlobal (Local3DPoint( (*closest_simhit).exitPoint().x(),
							  (*closest_simhit).exitPoint().y(),
							  (*closest_simhit).exitPoint().z() ) );
	  all_gx1 = GP1.x();
	  all_gx2 = GP2.x();
	  all_gy1 = GP1.y();
	  all_gy2 = GP2.y();
	  all_gz1 = GP1.z();
	  all_gz2 = GP2.z();
	  
	  MeasurementPoint mp1 =
	    topol->measurementPosition( LocalPoint( (*closest_simhit).entryPoint().x(),
						    (*closest_simhit).entryPoint().y(),
						    (*closest_simhit).entryPoint().z() ) );
	  MeasurementPoint mp2 =
	    topol->measurementPosition( LocalPoint( (*closest_simhit).exitPoint().x(),
						    (*closest_simhit).exitPoint().y(), 
						    (*closest_simhit).exitPoint().z() ) );
	  all_row1 = mp1.x();
	  all_col1 = mp1.y();
	  all_row2 = mp2.x();
	  all_col2 = mp2.y();
	  
	  all_simhitx = 0.5*(all_x1+all_x2);  
	  all_simhity = 0.5*(all_y1+all_y2);  
	  
	  all_rechitresx = all_rechitx - all_simhitx;
	  all_rechitresy = all_rechity - all_simhity;

	  all_rechitpullx = all_rechitresx / all_rechiterrx;
	  all_rechitpully = all_rechitresy / all_rechiterry;
	  
	  edm::Ref<edm::DetSetVector<SiPixelCluster>, SiPixelCluster> const& clust = pixeliter->cluster();
	  
	  all_npix = clust->size();
	  all_nxpix = clust->sizeX();
	  all_nypix = clust->sizeY();

	  all_clust_row = clust->x();
	  all_clust_col = clust->y();
	  
	  LocalPoint lp2 = topol->localPosition( MeasurementPoint( all_clust_row, all_clust_col ) );
	  all_clust_x = lp2.x();
	  all_clust_y = lp2.y();

	  all_clust_q = clust->charge();

	  all_clust_maxpixcol = clust->maxPixelCol();
	  all_clust_maxpixrow = clust->maxPixelRow();
	  all_clust_minpixcol = clust->minPixelCol();
	  all_clust_minpixrow = clust->minPixelRow();
	  
	  all_clust_geoid = clust->geographicalId();
  
	  all_simpx  = (*closest_simhit).momentumAtEntry().x();
	  all_simpy  = (*closest_simhit).momentumAtEntry().y();
	  all_simpz  = (*closest_simhit).momentumAtEntry().z();
	  all_eloss = (*closest_simhit).energyLoss();
	  all_simphi   = (*closest_simhit).phiAtEntry();
	  all_simtheta = (*closest_simhit).thetaAtEntry();
	  all_pidhit = (*closest_simhit).particleType();
	  all_trkid = (*closest_simhit).trackId();
	  
	  //--- Fill alpha and beta -- more useful for exploring the residuals...
	  all_beta  = atan2(all_simpz, all_simpy);
	  all_alpha = atan2(all_simpz, all_simpx);
	  
	  all_simproc = (int)closest_simhit->processType();
	  
	  const edm::SimTrackContainer& trks = *(simtracks.product());
	  SimTrackContainer::const_iterator trksiter;
	  for (trksiter = trks.begin(); trksiter != trks.end(); trksiter++) 
	    if ( (int)trksiter->trackId() == all_trkid ) 
	      {
		all_simtrketa = trksiter->momentum().eta();
		all_simtrkphi = trksiter->momentum().phi();
	      }
	  
	  all_vtxz = theGeomDet->surface().position().z();
	  all_vtxr = theGeomDet->surface().position().perp();
	  
	  //computeAnglesFromDetPosition(clust, 
	  //		       theGeomDet, 
	  //		       all_clust_alpha, all_clust_beta )

	  const std::vector<SiPixelCluster::Pixel>& pixvector = clust->pixels();
	  for ( int i=0;  i<(int)pixvector.size(); ++i)
	    {
	      SiPixelCluster::Pixel holdpix = pixvector[i];
	      all_pixrow[i] = holdpix.x;
	      all_pixcol[i] = holdpix.y;
	      all_pixadc[i] = holdpix.adc;
	      LocalPoint lp = topol->localPosition(MeasurementPoint(holdpix.x, holdpix.y));
	      all_pixx[i] = lp.x();
	      all_pixy[i]= lp.y();
	      GlobalPoint GP =  theGeomDet->surface().toGlobal(Local3DPoint(lp.x(),lp.y(),lp.z()));
	      all_pixgx[i] = GP.x();	
	      all_pixgy[i]= GP.y();
	      all_pixgz[i]= GP.z();
	    }
	  
	  ttree_all_hits_->Fill();
	  
	} // for ( ; pixeliter != pixelrechitRangeIteratorEnd; ++pixeliter) 
    } // for (TrackerGeometry::DetContainer::const_iterator it = pDD->dets().begin(); it != pDD->dets().end(); it++) 
 
  // ------------------------------------------------ all hits ---------------------------------------------------------------
  

  // ------------------------------------------------ track hits only -------------------------------------------------------- 
  
  if ( include_trk_hits_ )
    {
      // Get tracks
      edm::Handle<reco::TrackCollection> trackCollection;
      e.getByLabel(src_, trackCollection);
      const reco::TrackCollection *tracks = trackCollection.product();
      reco::TrackCollection::const_iterator tciter;
      
      if ( tracks->size() > 0 )
	{
	  // Loop on tracks
	  for ( tciter=tracks->begin(); tciter!=tracks->end(); ++tciter)
	    {
	      // First loop on hits: find matched hits
	      for ( trackingRecHit_iterator it = tciter->recHitsBegin(); it != tciter->recHitsEnd(); ++it) 
		{
		  const TrackingRecHit &thit = **it;
		  // Is it a matched hit?
		  const SiPixelRecHit* matchedhit = dynamic_cast<const SiPixelRecHit*>(&thit);
		  
		  if ( matchedhit ) 
		    {
		      rechitx = -9999.9;
		      rechity = -9999.9;
		      rechitz = -9999.9;
		      rechiterrx = -9999.9;
		      rechiterry = -9999.9;		      
		      rechitresx = -9999.9;
		      rechitresy = -9999.9;
		      rechitpullx = -9999.9;
		      rechitpully = -9999.9;
		      
		      npix = -9999;
		      nxpix = -9999;
		      nypix = -9999;
		      charge = -9999.9;
		      
		      edgex = -9999;
		      edgey = -9999;
		      
		      bigx = -9999;
		      bigy = -9999;
		      
		      alpha = -9999.9;
		      beta  = -9999.9;
		      
		      phi = -9999.9;
		      eta = -9999.9;
		      
		      subdetId = -9999;
		      
		      layer  = -9999; 
		      ladder = -9999; 
		      mod    = -9999; 
		      side   = -9999;  
		      disk   = -9999;  
		      blade  = -9999; 
		      panel  = -9999; 
		      plaq   = -9999; 
		      
		      half = -9999;
		      flipped = -9999;
		      
		      nsimhit = -9999;
		      pidhit  = -9999;
		      simproc = -9999;
		      
		      simhitx = -9999.9;
		      simhity = -9999.9;
		      
		      position = (*it)->localPosition();
		      error = (*it)->localPositionError();
		      
		      rechitx = position.x();
		      rechity = position.y();
		      rechitz = position.z();
		      rechiterrx = sqrt(error.xx());
		      rechiterry = sqrt(error.yy());
		      
		      npix = matchedhit->cluster()->size();
		      nxpix = matchedhit->cluster()->sizeX();
		      nypix = matchedhit->cluster()->sizeY();
		      charge = matchedhit->cluster()->charge();
		      
		      //Association of the rechit to the simhit
		      matched.clear();
		      matched = associate.associateHit(*matchedhit);
		      
		      nsimhit = (int)matched.size();
		      
		      if ( !matched.empty() ) 
			{
			  mindist = 999999.9;
			  float distx, disty, dist;
			  bool found_hit_from_generated_particle = false;
			  
			  int n_assoc_muon = 0;
			  
			  vector<PSimHit>::const_iterator closestit = matched.begin();
			  for (vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); ++m)
			    {
			      if ( checkType_ )
				{ // only consider associated simhits with the generated pid (muons)
				  int pid = (*m).particleType();
				  if ( abs(pid) != genType_ )
				    continue;
				}
			      
			      float simhitx = 0.5 * ( (*m).entryPoint().x() + (*m).exitPoint().x() );
			      float simhity = 0.5 * ( (*m).entryPoint().y() + (*m).exitPoint().y() );
			      
			      distx = fabs(rechitx - simhitx);
			      disty = fabs(rechity - simhity);
			      dist = sqrt( distx*distx + disty*disty );
			      
			      if ( dist < mindist )
				{
				  n_assoc_muon++;
				  
				  mindist = dist;
				  closestit = m;
				  found_hit_from_generated_particle = true;
				}
			    } // for (vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); m++)
			  
			  // This recHit does not have any simHit with the same particleType as the particles generated
			  // Ignore it as most probably come from delta rays.
			  if ( checkType_ && !found_hit_from_generated_particle )
			    continue; 
			  
			  if ( n_assoc_muon > 1 )
			    {
			      cout << " ----- This is not good: n_assoc_muon = " << n_assoc_muon << endl;
			      cout << "evt = " << evt << endl;
			    }
			  
			  pidhit = (*closestit).particleType();
			  simproc = (int)(*closestit).processType();
			  
			  simhitx = 0.5*( (*closestit).entryPoint().x() + (*closestit).exitPoint().x() );
			  simhity = 0.5*( (*closestit).entryPoint().y() + (*closestit).exitPoint().y() );
			  
			  rechitresx = rechitx - simhitx;
			  rechitresy = rechity - simhity;
			  rechitpullx = ( rechitx - simhitx ) / sqrt(error.xx());
			  rechitpully = ( rechity - simhity ) / sqrt(error.yy());
			  
			  float simhitpx = (*closestit).momentumAtEntry().x();
			  float simhitpy = (*closestit).momentumAtEntry().y();
			  float simhitpz = (*closestit).momentumAtEntry().z();
			  
			  beta  = atan2(simhitpz, simhitpy) * radtodeg;
			  alpha = atan2(simhitpz, simhitpx) * radtodeg;
			  
			  //beta  = fabs(atan2(simhitpz, simhitpy)) * radtodeg;
			  //alpha = fabs(atan2(simhitpz, simhitpx)) * radtodeg;
			  
			  phi = tciter->momentum().phi() / math_pi*180.0;
			  eta = tciter->momentum().eta();
			  
			  detId = (*it)->geographicalId();
			  
			  const int maxPixelCol = (*matchedhit).cluster()->maxPixelCol();
			  const int maxPixelRow = (*matchedhit).cluster()->maxPixelRow();
			  const int minPixelCol = (*matchedhit).cluster()->minPixelCol();
			  const int minPixelRow = (*matchedhit).cluster()->minPixelRow();
			  
			  const PixelGeomDetUnit* theGeomDet =
			    dynamic_cast<const PixelGeomDetUnit*> ((*tracker).idToDet(detId) );
			  
			  const RectangularPixelTopology* theTopol = 
			    dynamic_cast<const RectangularPixelTopology*>(&(theGeomDet->specificTopology()));
			  
			  // check whether the cluster is at the module edge 
			  if ( theTopol->isItEdgePixelInX( minPixelRow ) || 
			       theTopol->isItEdgePixelInX( maxPixelRow ) )
			    edgex = 1;
			  else 
			    edgex = 0;
			  
			  if ( theTopol->isItEdgePixelInY( minPixelCol ) || 
			       theTopol->isItEdgePixelInY( maxPixelCol ) )
			    edgey = 1;
			  else 
			    edgey = 0;
			  
			  // check whether this rechit contains big pixels
			  if ( theTopol->containsBigPixelInX(minPixelRow, maxPixelRow) )
			    bigx = 1;
			  else 
			    bigx = 0;
			  
			  if ( theTopol->containsBigPixelInY(minPixelCol, maxPixelCol) )
			    bigy = 1;
			  else 
			    bigy = 0;
			  
			  subdetId = (int)detId.subdetId();
			  
			  if ( (int)detId.subdetId() == (int)PixelSubdetector::PixelBarrel ) 
			    {
			      //const PixelGeomDetUnit* theGeomDet 
			      //= dynamic_cast<const PixelGeomDetUnit*> ( tracker->idToDet(detId) );
			      //const RectangularPixelTopology * topol 
			      //= dynamic_cast<const RectangularPixelTopology*>(&(theGeomDet->specificTopology()));
			      
			      int tmp_nrows = theGeomDet->specificTopology().nrows();
			      if ( tmp_nrows == 80 ) 
				half = 1;
			      else if ( tmp_nrows == 160 ) 
				half = 0;
			      else 
				cout << "-------------------------------------------------- Wrong module size !!!" << endl;
			      
			      float tmp1 = theGeomDet->surface().toGlobal(Local3DPoint(0.,0.,0.)).perp();
			      float tmp2 = theGeomDet->surface().toGlobal(Local3DPoint(0.,0.,1.)).perp();
			      
			      if ( tmp2<tmp1 ) 
				flipped = 1;
			      else 
				flipped = 0;
			      
			      PXBDetId  bdetid(detId);
			      layer  = bdetid.layer();   // Layer: 1,2,3.
			      ladder = bdetid.ladder();  // Ladder: 1-20, 32, 44. 
			      mod   = bdetid.module();  // Mod: 1-8.
			    }			  
			  else if ( (int)detId.subdetId() == (int)PixelSubdetector::PixelEndcap )
			    {
			      PXFDetId fdetid(detId);
			      side  = fdetid.side();
			      disk  = fdetid.disk();
			      blade = fdetid.blade();
			      panel = fdetid.panel();
			      plaq  = fdetid.module(); // also known as plaquette
			      
			    } // else if ( detId.subdetId()==PixelSubdetector::PixelEndcap )
			  //else std::cout << "We are not in the pixel detector. detId.subdetId() = " << (int)detId.subdetId() << endl;
			  
			  ttree_track_hits_->Fill();
			  
			} // if ( !matched.empty() )
		      else
			cout << "---------------- RecHit with no associated SimHit !!! -------------------------- " << endl;
		      
		    } // if ( matchedhit )
		  
		} // end of loop on hits
	      
	    } //end of loop on track 
      
	} // tracks > 0.
     
    } // if ( include_trk_hits_ )

  // ----------------------------------------------- track hits only -----------------------------------------------------------
  
}

void SiPixelErrorEstimation::
computeAnglesFromDetPosition(const SiPixelCluster & cl, 
			     const GeomDetUnit    & det, 
			     float& alpha, float& beta )
{
  //--- This is a new det unit, so cache it
  const PixelGeomDetUnit* theDet = dynamic_cast<const PixelGeomDetUnit*>( &det );
  if (! theDet) 
    {
      cout << "---------------------------------------------- Not a pixel detector !!!!!!!!!!!!!!" << endl;
      assert(0);
    }

  const RectangularPixelTopology* theTopol = 
    dynamic_cast<const RectangularPixelTopology*>(&(theDet->specificTopology()));

  // get cluster center of gravity (of charge)
  float xcenter = cl.x();
  float ycenter = cl.y();
  
  // get the cluster position in local coordinates (cm) 
  LocalPoint lp = theTopol->localPosition( MeasurementPoint(xcenter, ycenter) );
  //float lp_mod = sqrt( lp.x()*lp.x() + lp.y()*lp.y() + lp.z()*lp.z() );

  // get the cluster position in global coordinates (cm)
  GlobalPoint gp = theDet->surface().toGlobal( lp );
  float gp_mod = sqrt( gp.x()*gp.x() + gp.y()*gp.y() + gp.z()*gp.z() );

  // normalize
  float gpx = gp.x()/gp_mod;
  float gpy = gp.y()/gp_mod;
  float gpz = gp.z()/gp_mod;

  // make a global vector out of the global point; this vector will point from the 
  // origin of the detector to the cluster
  GlobalVector gv(gpx, gpy, gpz);

  // make local unit vector along local X axis
  const Local3DVector lvx(1.0, 0.0, 0.0);

  // get the unit X vector in global coordinates/
  GlobalVector gvx = theDet->surface().toGlobal( lvx );

  // make local unit vector along local Y axis
  const Local3DVector lvy(0.0, 1.0, 0.0);

  // get the unit Y vector in global coordinates
  GlobalVector gvy = theDet->surface().toGlobal( lvy );
   
  // make local unit vector along local Z axis
  const Local3DVector lvz(0.0, 0.0, 1.0);

  // get the unit Z vector in global coordinates
  GlobalVector gvz = theDet->surface().toGlobal( lvz );
    
  // calculate the components of gv (the unit vector pointing to the cluster) 
  // in the local coordinate system given by the basis {gvx, gvy, gvz}
  // note that both gv and the basis {gvx, gvy, gvz} are given in global coordinates
  float gv_dot_gvx = gv.x()*gvx.x() + gv.y()*gvx.y() + gv.z()*gvx.z();
  float gv_dot_gvy = gv.x()*gvy.x() + gv.y()*gvy.y() + gv.z()*gvy.z();
  float gv_dot_gvz = gv.x()*gvz.x() + gv.y()*gvz.y() + gv.z()*gvz.z();

  // calculate angles
  alpha = atan2( gv_dot_gvz, gv_dot_gvx );
  beta  = atan2( gv_dot_gvz, gv_dot_gvy );

  // calculate cotalpha and cotbeta
  //   cotalpha_ = 1.0/tan(alpha_);
  //   cotbeta_  = 1.0/tan(beta_ );
  // or like this
  //cotalpha_ = gv_dot_gvx / gv_dot_gvz;
  //cotbeta_  = gv_dot_gvy / gv_dot_gvz;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiPixelErrorEstimation);
