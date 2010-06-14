/**
 *  Class: ChamberSegmentUtility
 *
 *  Description:
 *  utility class for the dynamical truncation algorithm
 *
 *  $Date: 2010/05/10 14:23:50 $
 *  $Revision: 1.1 $
 *
 *  Authors :
 *  D. Pagano & G. Bruno - UCL Louvain
 *
 **/

#include "RecoMuon/GlobalTrackingTools/interface/ChamberSegmentUtility.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include <iostream>
#include <map>
#include <vector>
#include <iostream>


ChamberSegmentUtility::ChamberSegmentUtility(const edm::Event& Event, const edm::EventSetup& Setup)
{
  Setup.get<MuonGeometryRecord>().get(cscGeometry);
  Event.getByLabel("cscSegments", CSCSegments);
  Setup.get<MuonGeometryRecord>().get(dtGeom);
  Event.getByLabel("dt4DSegments", all4DSegments);
}



vector<CSCSegment> ChamberSegmentUtility::getCSCSegmentsInChamber(CSCDetId sel)
{
  unsigned int index = 0;
  for (  CSCSegmentCollection::id_iterator chamberId = CSCSegments->id_begin();
	 chamberId != CSCSegments->id_end(); ++chamberId, ++index ) {
    if ((*chamberId).chamber() != sel.chamber()) continue;
    CSCSegmentCollection::range  range = CSCSegments->get((*chamberId));
    for (CSCSegmentCollection::const_iterator segment = range.first;
	 segment!=range.second; ++segment) {
      cscseg.push_back(*segment);
    }
  }
  return cscseg;
}



vector<DTRecSegment4D> ChamberSegmentUtility::getDTSegmentsInChamber(DTChamberId sel)
{
  DTRecSegment4DCollection::id_iterator chamberIdIt;
  for (chamberIdIt = all4DSegments->id_begin();
       chamberIdIt != all4DSegments->id_end();
       ++chamberIdIt){
    if (*chamberIdIt != sel) continue;
    DTRecSegment4DCollection::range  range = all4DSegments->get((*chamberIdIt));
    for (DTRecSegment4DCollection::const_iterator segment = range.first;
         segment!=range.second; ++segment){
      dtseg.push_back(*segment);
    }
  }
  return dtseg;
}



vector<CSCRecHit2D> ChamberSegmentUtility::getCSCRHmap(CSCSegment selected)
{
  vector<CSCRecHit2D> allchRH;
  unsigned int index = 0;
  for (  CSCSegmentCollection::id_iterator chamberId = CSCSegments->id_begin();
         chamberId != CSCSegments->id_end(); ++chamberId, ++index ) {
    CSCSegmentCollection::range  range = CSCSegments->get((*chamberId));
    for (CSCSegmentCollection::const_iterator segment = range.first;
         segment!=range.second; ++segment) {
      if((*segment).parameters() == selected.parameters()) {
	allchRH = (*segment).specificRecHits();
      }
    }  
  }
  return allchRH;
}

  

vector<DTRecHit1D> ChamberSegmentUtility::getDTRHmap(DTRecSegment4D selected)
{
  vector<DTRecHit1D> allchRH;
  phiSegRH.clear();
  zSegRH.clear();
  DTRecSegment4DCollection::id_iterator chamberIdIt;
  for (chamberIdIt = all4DSegments->id_begin();
       chamberIdIt != all4DSegments->id_end();
       ++chamberIdIt){
    DTRecSegment4DCollection::range  range = all4DSegments->get((*chamberIdIt));
    for (DTRecSegment4DCollection::const_iterator segment = range.first;
         segment!=range.second; ++segment){
      if((*segment).parameters() == selected.parameters()) {
	if((*segment).hasPhi()){
	  const DTChamberRecSegment2D* phiSeg = (*segment).phiSegment();
	  phiSegRH = phiSeg->specificRecHits();
	}
	if((*segment).hasZed()){
	  const DTSLRecSegment2D* zSeg = (*segment).zSegment();
	  zSegRH = zSeg->specificRecHits();
	}
	for (vector<DTRecHit1D>::const_iterator itphi = phiSegRH.begin(); itphi != phiSegRH.end(); itphi++) allchRH.push_back(*itphi);
	for (vector<DTRecHit1D>::iterator itz = zSegRH.begin(); itz < zSegRH.end(); itz++) allchRH.push_back(*itz);
      }
    }
  }
  return allchRH;
}



