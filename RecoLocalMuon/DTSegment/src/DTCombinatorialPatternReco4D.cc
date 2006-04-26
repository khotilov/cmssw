/** \file
 *
 * $Date: 2006/04/26 09:55:08 $
 * $Revision: 1.6 $
 * \author Stefano Lacaprara - INFN Legnaro <stefano.lacaprara@pd.infn.it>
 * \author Riccardo Bellan - INFN TO <riccardo.bellan@cern.ch>
 */

#include "RecoLocalMuon/DTSegment/src/DTRecSegment2DBaseAlgo.h"
#include "RecoLocalMuon/DTSegment/src/DTRecSegment2DAlgoFactory.h"

#include "RecoLocalMuon/DTSegment/src/DTSegmentUpdator.h"

// For the 2D reco I use thei reconstructor!
#include "RecoLocalMuon/DTSegment/src/DTCombinatorialPatternReco.h"

#include "RecoLocalMuon/DTSegment/src/DTCombinatorialPatternReco4D.h"
#include "DataFormats/MuonDetId/interface/DTDetIdAccessor.h"
#include "DataFormats/DTRecHit/interface/DTRecHit1DPair.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"


using namespace std;
using namespace edm;
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment2DPhi.h"
#include "RecoLocalMuon/DTSegment/src/DTSegmentCand.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

DTCombinatorialPatternReco4D::DTCombinatorialPatternReco4D(const ParameterSet& pset):
  DTRecSegment4DBaseAlgo(pset), theAlgoName("DTCombinatorialPatternReco4D"){

  // debug parameter
  debug = pset.getUntrackedParameter<bool>("debug");
  
  // the updator
  theUpdator = new DTSegmentUpdator(pset);

  // the input type. 
  // If true the instructions in setDTRecSegment2DContainer will be schipped and the 
  // theta segment will be recomputed from the 1D rechits
  // If false the theta segment will be taken from the Event. Caveat: in this case the
  // event must contain the 2D segments!
  allDTRecHits = pset.getParameter<bool>("AllDTRecHits");
 
  // Get the concrete 2D-segments reconstruction algo from the factory
  // For the 2D reco I use this reconstructor!
  the2DAlgo = new DTCombinatorialPatternReco(pset.getParameter<ParameterSet>("Reco2DAlgoConfig"));
  
  //   string theReco2DAlgoName = pset.getParameter<string>("Reco2DAlgoName");
  //   cout << "the Reco2D AlgoName is " << theReco2DAlgoName << endl;
  //   the2DAlgo = DTRecSegment2DAlgoFactory::get()->create(theReco2DAlgoName,
  // 						       pset.getParameter<ParameterSet>("Reco2DAlgoConfig"));
}

void DTCombinatorialPatternReco4D::setES(const EventSetup& setup){
  setup.get<MuonGeometryRecord>().get(theDTGeometry);
  the2DAlgo->setES(setup);
  theUpdator->setES(setup);
}

void DTCombinatorialPatternReco4D::setChamber(const DTChamberId &chId){
  // Set the chamber
  theChamber = theDTGeometry->chamber(chId); 
}

void DTCombinatorialPatternReco4D::setDTRecHit1DContainer(Handle<DTRecHitCollection> all1DHits){
  theHitsFromPhi1.clear();
  theHitsFromPhi2.clear();
  theHitsFromTheta.clear();

  // FIXME!!!
  DTRecHitCollection::range rangeHitsFromPhi1 = all1DHits->get(DTLayerId(theChamber->id(),1,1), DTSuperLayerIdComparator());
  DTRecHitCollection::range rangeHitsFromPhi2 = all1DHits->get(DTLayerId(theChamber->id(),3,1), DTSuperLayerIdComparator());
  //

  vector<DTRecHit1DPair> hitsFromPhi1(rangeHitsFromPhi1.first,rangeHitsFromPhi1.second);
  vector<DTRecHit1DPair> hitsFromPhi2(rangeHitsFromPhi2.first,rangeHitsFromPhi2.second);
  if(debug)
    cout<< "Number of DTRecHit1DPair in the SL 1 (Phi 1): " << hitsFromPhi1.size()<<endl
	<< "Number of DTRecHit1DPair in the SL 3 (Phi 2): " << hitsFromPhi2.size()<<endl;
  
  theHitsFromPhi1 = hitsFromPhi1;
  theHitsFromPhi2 = hitsFromPhi2;

  if(allDTRecHits){
    // FIXME!!!
    DTRecHitCollection::range rangeHitsFromTheta = all1DHits->get(DTLayerId(theChamber->id(),2,1), DTSuperLayerIdComparator());
    //
    
    vector<DTRecHit1DPair> hitsFromTheta(rangeHitsFromTheta.first,rangeHitsFromTheta.second);
    if(debug)
      cout<< "Number of DTRecHit1DPair in the SL 2 (Theta): " << hitsFromTheta.size()<<endl;
    theHitsFromTheta = hitsFromTheta;
  }

}

void DTCombinatorialPatternReco4D::setDTRecSegment2DContainer(Handle<DTRecSegment2DCollection> all2DSegments){
  theSegments2DTheta.clear();

  if(!allDTRecHits){

    //Extract the DTRecSegment2DCollection range for the theta SL
    DTRecSegment2DCollection::range rangeTheta = all2DSegments->get(DTDetIdAccessor::bySuperLayer(DTSuperLayerId(theChamber->id(),2)));
    
    // Fill the DTRecSegment2D container for the theta SL
    vector<DTRecSegment2D> segments2DTheta(rangeTheta.first,rangeTheta.second);
    
    if(debug)
      cout << "Number of 2D-segments in the second SL (Theta)" << segments2DTheta.size() << endl;
    theSegments2DTheta = segments2DTheta;
  }

}

  
OwnVector<DTRecSegment4D>
DTCombinatorialPatternReco4D::reconstruct(){

  OwnVector<DTRecSegment4D> result;
  
  if (debug) cout << "Segments in " << theChamber->id() << endl;

  vector<DTSegmentCand*> resultPhi = buildPhiSuperSegmentsCandidates();
  
  if (debug) cout << "There are " << resultPhi.size() << " Phi cand" << endl;
  
  if (allDTRecHits){
    // take the theta SL of this chamber
    const DTSuperLayer* sl = theChamber->superLayer(2);
    // sl points to 0 if the theta SL was not found
    if(sl){
      // reconstruct the theta segments
      OwnVector<DTRecSegment2D> thetaSegs = the2DAlgo->reconstruct(sl, theHitsFromTheta);
      vector<DTRecSegment2D> segments2DTheta(thetaSegs.begin(),thetaSegs.end());
      theSegments2DTheta = segments2DTheta;
    }
  }
  
  bool hasZed=false;
  
  // has this chamber the Z-superlayer?
  if (theSegments2DTheta.size()){
    hasZed = theSegments2DTheta.size()>0;
    if (debug) cout << "There are " << theSegments2DTheta.size() << " Theta cand" << endl;
  } else {
    if (debug) cout << "No Theta SL" << endl;
  }

  // Now I want to build the concrete DTRecSegment4D.
  if (resultPhi.size()) {
    for (vector<DTSegmentCand*>::const_iterator phi=resultPhi.begin();
         phi!=resultPhi.end(); ++phi) {
      
      //FIXME, check the converter and change its name
      DTRecSegment2DPhi* superPhi = (*phi)->convert(theChamber);
      
      theUpdator->update(superPhi);
      
      
      // << start
      if (hasZed) {

	// Create all the 4D-segment combining the Z view with the Phi one
	// loop over the Z segments
	for(vector<DTRecSegment2D>::const_iterator zed = theSegments2DTheta.begin();
	    zed != theSegments2DTheta.end(); ++zed){
	  
	  // Important!!
	  DTSuperLayerId ZedSegSLId(zed->geographicalId().rawId());
	  
	  const LocalPoint posZInCh  = theChamber->toLocal( theChamber->superLayer(ZedSegSLId)->toGlobal(zed->localPosition() )) ;
	  const LocalVector dirZInCh = theChamber->toLocal( theChamber->superLayer(ZedSegSLId)->toGlobal(zed->localDirection() )) ;
	
          DTRecSegment4D* newSeg = new DTRecSegment4D(*superPhi,*zed,posZInCh,dirZInCh);
	  //<<
  
          /// 4d segment: I have the pos along the wire => further update!
          theUpdator->update(newSeg);
          if (debug) cout << "Created a 4D seg " << endl;
	  result.push_back(newSeg);
        }
      } else {
        // Only phi
        DTRecSegment4D* newSeg = new DTRecSegment4D(*superPhi);
	
        if (debug) cout << "Created a 4D segment using only the 2D Phi segment" << endl;
	result.push_back(newSeg);
      }
    }
  } else { 
    // DTRecSegment4D from zed projection only (unlikely, not so useful, but...)
    if (hasZed) {
      for(vector<DTRecSegment2D>::const_iterator zed = theSegments2DTheta.begin();
	  zed != theSegments2DTheta.end(); ++zed){
        
	// Important!!
	DTSuperLayerId ZedSegSLId(zed->geographicalId().rawId());
	  
	const LocalPoint posZInCh  = theChamber->toLocal( theChamber->superLayer(ZedSegSLId)->toGlobal(zed->localPosition() )) ;
	const LocalVector dirZInCh = theChamber->toLocal( theChamber->superLayer(ZedSegSLId)->toGlobal(zed->localDirection() )) ;
	
        DTRecSegment4D* newSeg = new DTRecSegment4D( *zed,posZInCh,dirZInCh);
	// <<
	
        if (debug) cout << "Created a 4D segment using only the 2D Theta segment" << endl;
	result.push_back(newSeg);
      }
    }
  }
  // finally delete the candidates!
  for (vector<DTSegmentCand*>::iterator phi=resultPhi.begin();
       phi!=resultPhi.end(); ++phi) delete *phi;
  
  return result;
}



vector<DTSegmentCand*> DTCombinatorialPatternReco4D::buildPhiSuperSegmentsCandidates(){
  
  DTSuperLayerId slId;

  if(theHitsFromPhi1.size())
    slId = theHitsFromPhi1.front().wireId().superlayerId();
  else
    if(theHitsFromPhi2.size())
      slId = theHitsFromPhi2.front().wireId().superlayerId();
    else{
      if(debug) cout<<"DTCombinatorialPatternReco4D::buildPhiSuperSegmentsCandidates: "
		    <<"No Hits in the two Phi SL"<<endl;
      return vector<DTSegmentCand*>();
    }

  const DTSuperLayer *sl = theDTGeometry->superLayer(slId);
  
  vector<DTHitPairForFit*> pairPhi1 = the2DAlgo->initHits(sl,theHitsFromPhi1);
  // same sl!! Since the fit will be in the sl phi 1!
  vector<DTHitPairForFit*> pairPhi2 = the2DAlgo->initHits(sl,theHitsFromPhi2);
  // copy the pairPhi2 in the pairPhi1 vector 
  copy(pairPhi2.begin(),pairPhi2.end(),back_inserter(pairPhi1));

  // Build the segment candidate
  return the2DAlgo->buildSegments(sl,pairPhi1);
}
