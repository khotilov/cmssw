#include "RecoTracker/Loopers/interface/Methods.h"


#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include <fstream>

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"

#include "SimG4Core/Physics/interface/ProcessTypeEnumerator.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TAxis.h"
#include "TMath.h"
#include "TH2F.h"
#include "TGraph.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/Jet.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"

#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "RecoTracker/TkTrackingRegions/interface/GlobalTrackingRegion.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h"
#include <iostream>

void LooperClusterRemoverMethod::FractionOfTruth::run(edm::Event& iEvent, const edm::EventSetup& iSetup,
						      LooperClusterRemover::products &prod_)
{

  ProcessTypeEnumerator translateG4;
    
  //implement one method
  edm::Handle<SiPixelRecHitCollection> pixelHits;
  iEvent.getByLabel(pixelRecHits_,pixelHits);

  //associator_ = new TrackerHitAssociator(iEvent,associatorConf_);
  TrackerHitAssociator  * associator = new TrackerHitAssociator(iEvent);

  typedef std::map<unsigned int, TrackingParticleCollection::const_iterator> ReverseTPCollection;
  ReverseTPCollection reverseCollection;
  edm::Handle<TrackingParticleCollection> mcTruth;
  iEvent.getByLabel(edm::InputTag("mergedtruth","MergedTrackTruth"),mcTruth);
  for (TrackingParticleCollection::const_iterator tpIt=mcTruth->begin();
       tpIt!=mcTruth->end();++tpIt){
    if (tpIt->pt() > 0.5) continue;
    //    std::cout<<" tp with pt: "<<tpIt->pt();
    for (TrackingParticle::g4t_iterator iG4=tpIt->g4Track_begin();
	 iG4!=tpIt->g4Track_end();++iG4){
      //      std::cout<<" G4: "<<iG4->trackId();
      reverseCollection[iG4->trackId()]=tpIt;
    }
    //    std::cout<<"id :"<<tpIt->pdgId()<<std::endl;
  }


  for ( SiPixelRecHitCollection::const_iterator dPxIt=pixelHits->begin();
        dPxIt!=pixelHits->end();++dPxIt) {
    for (SiPixelRecHitCollection::DetSet::const_iterator pxIt=dPxIt->begin();
         pxIt!=dPxIt->end();++pxIt) {
      //determine whether this a pixel from a true looper
      std::vector<PSimHit> simHits = associator->associateHit(*pxIt);
      //      std::cout<<" for :"<<dPxIt->id()<<" a hit has "<<simHits.size()<<" sim hit associated"<<std::endl;
      bool rejected=false;
      for (uint ips=0;ips!=simHits.size();++ips){
	//	std::cout<<"[P] pType: "<<simHits[ips].particleType()<<", pId: "<<simHits[ips].trackId() <<", process: "<<simHits[ips].processType()<<", event: "<<simHits[ips].eventId().rawId()<<std::endl;
	ReverseTPCollection::iterator tpIt=reverseCollection.find(simHits[ips].trackId());
	if (tpIt!=reverseCollection.end()){
	  //	  std::cout<<"\t belongs to tp with pt: "<<tpIt->second->pt()<<" and pz: "<< tpIt->second->pz()<<std::endl;
	  double Rprod=tpIt->second->parentVertex()->position().Pt();
	  //	  std::cout<<" production vertex: "<<tpIt->second->parentVertex()->position()<<" radius: "<<Rprod<<std::endl;
	  if (tpIt->second->pt() < 0.200 && Rprod<2){
	    rejected=true;
	    edm::LogError("FractionOfTruth")<<"reject [P] pType: "<<simHits[ips].particleType()<<", pId: "<<simHits[ips].trackId() <<", process: "
					    << simHits[ips].processType()<<"::"<<translateG4.processName(simHits[ips].processType())
					    <<", event: "<<simHits[ips].eventId().rawId()
					    <<" production vertex: "<<tpIt->second->parentVertex()->position()<<" radius: "<<Rprod<<std::endl;
	  }
	}
	//	if (simHits[ips].processType() == 2) rejected=true;
	//	if (abs(simHits[ips].particleType()) == 211) rejected=true;
      }
      if (rejected) prod_.collectedPixels[pxIt->cluster().key()]=true;
    }
  }

  edm::Handle<SiStripMatchedRecHit2DCollection> matchedHits;
  iEvent.getByLabel(stripRecHits_,matchedHits);
  
  for (SiStripMatchedRecHit2DCollection::const_iterator dStIt=matchedHits->begin();
       dStIt!=matchedHits->end();++dStIt){
    for (SiStripMatchedRecHit2DCollection::DetSet::const_iterator stIt=dStIt->begin();
	 stIt!=dStIt->end();++stIt){
      std::vector<PSimHit> simHits = associator->associateHit(*stIt);
      bool rejected=false; 
      for (uint ips=0;ips!=simHits.size();++ips){
	//	std::cout<<"[M] pType: "<<simHits[ips].particleType()<<", pId: "<<simHits[ips].trackId() 
	//	<< simHits[ips].processType()<<"::"<<translateG4.processName(simHits[ips].processType())
	//		     <<", event: "<<simHits[ips].eventId().rawId()<<std::endl;	    
	ReverseTPCollection::iterator tpIt=reverseCollection.find(simHits[ips].trackId());
	if (tpIt!=reverseCollection.end()){
	  double Rprod=tpIt->second->parentVertex()->position().Pt();
	  //	  std::cout<<" production vertex: "<<tpIt->second->parentVertex()->position()<<" radius: "<<Rprod<<std::endl;
	  if (tpIt->second->pt() < 0.200 && Rprod<2) {
	    rejected=true;
	    edm::LogError("FractionOfTruth")<<"reject [M] pType: "<<simHits[ips].particleType()<<", pId: "<<simHits[ips].trackId() 
					    << simHits[ips].processType()<<"::"<<translateG4.processName(simHits[ips].processType())
					    <<", event: "<<simHits[ips].eventId().rawId()
					    <<" production vertex: "<<tpIt->second->parentVertex()->position()<<" radius: "<<Rprod<<std::endl;
	  }
	}
      }
      if (rejected) {
	prod_.collectedStrips[stIt->stereoHit()->cluster().key()]=true;
	prod_.collectedStrips[stIt->monoHit()->cluster().key()]=true;
      }
    }
  }  
  
  delete associator;  
}


void LooperClusterRemoverMethod::EveryNMethod::run(edm::Event& iEvent, const edm::EventSetup& iSetup,
						   LooperClusterRemover::products &prod_)
{
  //implement one method
  edm::Handle<SiPixelRecHitCollection> pixelHits;
  edm::Handle<SiStripMatchedRecHit2DCollection> matchedHits;
  //  edm::Handle<SiStripRecHit2DCollection> rphiHits;

  //get the products
  iEvent.getByLabel(pixelRecHits_,pixelHits);
  iEvent.getByLabel(stripRecHits_,matchedHits);

  unsigned int countMe=0;
  for ( SiPixelRecHitCollection::const_iterator dPxIt=pixelHits->begin();
	dPxIt!=pixelHits->end();++dPxIt)  {
    for (SiPixelRecHitCollection::DetSet::const_iterator pxIt=dPxIt->begin();
	 pxIt!=dPxIt->end();++pxIt) {
      if (countMe++==everyNPixel_){
	countMe=0;
	//fill in a ref to mask
	prod_.collectedPixels[pxIt->cluster().key()]=true;
      }
    }
  }

  for (SiStripMatchedRecHit2DCollection::const_iterator dStIt=matchedHits->begin();
       dStIt!=matchedHits->end();++dStIt){
    for (SiStripMatchedRecHit2DCollection::DetSet::const_iterator stIt=dStIt->begin();
	 stIt!=dStIt->end();++stIt){
      if (countMe++==everyNPixel_){
	countMe=0;
	prod_.collectedStrips[stIt->stereoHit()->cluster().key()]=true; 
	prod_.collectedStrips[stIt->monoHit()->cluster().key()]=true;
      }
    }
  }

}

struct classcomp {
  bool operator() ( const GlobalPoint & p1, const GlobalPoint & p2) const{
    if (p1.x() > p2.x()) return true;
    else{
      if (p1.y() > p2.y()) return true;
      else{
	if (p1.z() > p2.z()) return true;
      }
    }
    return false;
  }
};

LooperClusterRemoverMethod::ReadFileIn::ReadFileIn(const edm::ParameterSet& iConfig)
  :  file_(iConfig.getParameter<std::string>("maskFile").c_str(),fstream::in)
{
  pixelRecHits_=iConfig.getParameter<edm::InputTag>("pixelRecHits");
  stripRecHits_=iConfig.getParameter<edm::InputTag>("stripRecHits");
  epsilon_=iConfig.getParameter<double>("epsilon");
  withDetId_=iConfig.getParameter<bool>("withDetId");
}

void LooperClusterRemoverMethod::ReadFileIn::run(edm::Event& iEvent,  const edm::EventSetup& iSetup,
						 LooperClusterRemover::products &prod_)
{
  unsigned int nNext=0;
  file_>>nNext;
  //  std::set<GlobalPoint,classcomp> points;
  uint detid;
  double x,y,z;
  std::cout<<"to read: "<<nNext<<std::endl;
  std::vector<GlobalPoint> points(nNext);
  
  for (unsigned int iRead=0;iRead!=nNext;++iRead){
    if (withDetId_)file_>>detid;
    file_>>x;
    file_>>y;
    file_>>z;
    //    std::cout<<x<<" "<<y<<" "<<z<<std::endl;
    //    points.insert(GlobalPoint(x,y,z));
    points[iRead]=GlobalPoint(x,y,z);
  }
  std::vector<uint> counts(points.size(),0);
  
  TrackerHitAssociator * associator= 0;
  if (!iEvent.isRealData()) associator=new TrackerHitAssociator(iEvent);
  
  //implement one method
  edm::Handle<SiPixelRecHitCollection> pixelHits;
  edm::Handle<SiStripMatchedRecHit2DCollection> matchedHits;
  //  edm::Handle<SiStripRecHit2DCollection> rphiHits;

  //get the products
  iEvent.getByLabel(pixelRecHits_,pixelHits);
  iEvent.getByLabel(stripRecHits_,matchedHits);

  edm::ESHandle<TransientTrackingRecHitBuilder> builder;
  iSetup.get<TransientRecHitRecord>().get("WithTrackAngle",builder);

  for ( SiPixelRecHitCollection::const_iterator dPxIt=pixelHits->begin();
	dPxIt!=pixelHits->end();++dPxIt)  {
    for (SiPixelRecHitCollection::DetSet::const_iterator pxIt=dPxIt->begin();
	 pxIt!=dPxIt->end();++pxIt) {

      //get the transient point for position
      TransientTrackingRecHit::RecHitPointer tHit = builder->build(&*pxIt);
      GlobalPoint tPoint=tHit->globalPosition();
      //      std::cout<<"one pixel at: "<<tPoint<<std::endl;
      //find something in the list to be removed
      float closest=100;
      bool masked=false;
      for (uint iR=0;iR!=points.size();++iR){
	float mag=(points[iR]-tPoint).mag();
	if (mag<closest) closest=mag;
	if (mag < epsilon_) 
	  { 
	    prod_.collectedPixels[pxIt->cluster().key()]=true;
	    if (associator){
	      std::vector<PSimHit> simHits = associator->associateHit(*pxIt);
	      std::stringstream text;
	      text<<" kill a pixel at: "<<dPxIt->id()<<", "<<tPoint<<" using "<<points[iR];
	      for (uint ips=0;ips!=simHits.size();++ips)	      
		text<<"\n\t[P] pType: "<<simHits[ips].particleType()<<", pId: "<<simHits[ips].trackId() <<", process: "<<simHits[ips].processType()<<", event: "<<simHits[ips].eventId().rawId();
	      edm::LogError("ReadFileIn")<<text.str();
	    }
	    counts[iR]++;
	    masked=true;
	    break;
	  }
      }
      if (!masked && closest!=100) edm::LogError("ReadFileIn")<<"closest point to pixel "<<tPoint<<" was at:"<<closest;
    }
  }
  

  for (SiStripMatchedRecHit2DCollection::const_iterator dStIt=matchedHits->begin();
       dStIt!=matchedHits->end();++dStIt){
    for (SiStripMatchedRecHit2DCollection::DetSet::const_iterator stIt=dStIt->begin();
	 stIt!=dStIt->end();++stIt){
      TransientTrackingRecHit::RecHitPointer tHit = builder->build(&*stIt);
      GlobalPoint tPoint=tHit->globalPosition();
      //      std::cout<<"one matched at: "<<tPoint<<std::endl;
      //find something in the list to be removed
      float closest=100;
      bool masked=false;
      for (uint iR=0;iR!=points.size();++iR){     
	float mag=(points[iR]-tPoint).mag();
	//	float mag=(points[iR]-tPoint).perp();
	if (mag<closest) closest=mag;
	if (mag < epsilon_)                     
	  {                  
	    prod_.collectedStrips[stIt->stereoHit()->cluster().key()]=true; 
	    prod_.collectedStrips[stIt->monoHit()->cluster().key()]=true;
	    if (associator){
	      std::stringstream text;
	      text<<" kill a strip at: "<<dStIt->id()<<"="<<stIt->stereoHit()->geographicalId().rawId()<<"+"<<stIt->monoHit()->geographicalId().rawId()<<" "<<tPoint<<" using "<<points[iR];
	      std::vector<PSimHit> simHits = associator->associateHit(*stIt);
	      for (uint ips=0;ips!=simHits.size();++ips)	 
		text<<"\n\t[M] pType: "<<simHits[ips].particleType()<<", pId: "<<simHits[ips].trackId() <<", process: "<<simHits[ips].processType()<<", event: "<<simHits[ips].eventId().rawId();
	      edm::LogError("ReadFileIn")<<text.str();
	    }
	    counts[iR]++;
	    masked=true;
	    break;         
	  }
      }
      if (!masked && closest!=100) edm::LogError("ReadFileIn")<<"closest point to strip "<< tPoint <<" was at:"<<closest;
    }
  }  


  if (associator)    delete associator;
  
  std::stringstream text;
  uint notused=0;
  for (uint iR=0;iR!=points.size();++iR){
    if (counts[iR]==0){
      bool foundAnotherPoint=false;
      for (uint iRo=0;iRo!=points.size();++iRo){
	if (counts[iRo]==0) continue; // only in the considered ones
	if (iRo==iR) continue; // only not this point
	float m =(points[iRo]-points[iR]).mag();
	if (m<epsilon_){
	  foundAnotherPoint=true;
	  break;
	}
      }
      if (!foundAnotherPoint){
	text<<"point not used:"<<points[iR]<<std::endl;
	notused++;
      }
    }
  }
  text<<notused<<" not used points from file";
  edm::LogError("ReadFileIn")<<text.str();
}


#include "DataDumper.h"


void LooperClusterRemoverMethod::LooperMethod::run(edm::Event& iEvent, const edm::EventSetup& iSetup,
						   LooperClusterRemover::products &prod_)
{

  //cheater
  TrackerHitAssociator * associator=0;
  if (!iEvent.isRealData())  associator= new TrackerHitAssociator(iEvent);

  std::vector<fastRecHit> fastHits;

  //implement one method
  edm::Handle<SiPixelRecHitCollection> pixelHits;
  edm::Handle<SiStripMatchedRecHit2DCollection> matchedHits;
  edm::Handle<SiStripRecHit2DCollection> rphiHits;
  edm::Handle<SiStripRecHit2DCollection> stereoHits;

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel("offlineBeamSpot",beamSpot);
  
  //  GlobalPoint bs(0,0,0);
  GlobalPoint bs(beamSpot->position().x(),
		 beamSpot->position().y(),
		 beamSpot->position().z());		 

  
  //get the products
  iEvent.getByLabel(pixelRecHits_,pixelHits);
  iEvent.getByLabel(stripRecHits_,matchedHits);
  if (useUnmatched_){
    iEvent.getByLabel(rphiRecHits_,rphiHits);
    iEvent.getByLabel(stereoRecHits_,stereoHits);
  }

  edm::ESHandle<TransientTrackingRecHitBuilder> builder;
  iSetup.get<TransientRecHitRecord>().get("WithTrackAngle",builder);


  edm::ESHandle<MagneticField> magField;
  iSetup.get<IdealMagneticFieldRecord>().get(magField);

  if (useUnmatched_)
    fastHits.reserve(pixelHits->dataSize()+matchedHits->dataSize()+rphiHits->size()+stereoHits->size());
  else
    fastHits.reserve(pixelHits->dataSize()+matchedHits->dataSize());

  for ( SiPixelRecHitCollection::const_iterator dPxIt=pixelHits->begin();
	dPxIt!=pixelHits->end();++dPxIt)  {
    for (SiPixelRecHitCollection::DetSet::const_iterator pxIt=dPxIt->begin();
	 pxIt!=dPxIt->end();++pxIt) {
      fastHits.push_back(fastRecHit(builder->build(&*pxIt),bs));
    }
  }
  
		   
  for (SiStripMatchedRecHit2DCollection::const_iterator dStIt=matchedHits->begin();
       dStIt!=matchedHits->end();++dStIt){
    for (SiStripMatchedRecHit2DCollection::DetSet::const_iterator stIt=dStIt->begin();
	 stIt!=dStIt->end();++stIt){
      fastHits.push_back(fastRecHit(builder->build(&*stIt),bs));
    }
  }  

  if (useUnmatched_){
  for (SiStripRecHit2DCollection::const_iterator dStIt=rphiHits->begin();
       dStIt!=rphiHits->end();++dStIt){
    for (SiStripRecHit2DCollection::DetSet::const_iterator stIt=dStIt->begin(); 
	 stIt!=dStIt->end();++stIt){
      fastHits.push_back(fastRecHit(builder->build(&*stIt),bs,true));
    }
  }

  for (SiStripRecHit2DCollection::const_iterator dStIt=stereoHits->begin();
       dStIt!=stereoHits->end();++dStIt){
    for (SiStripRecHit2DCollection::DetSet::const_iterator stIt=dStIt->begin(); 
	 stIt!=dStIt->end();++stIt){
      fastHits.push_back(fastRecHit(builder->build(&*stIt),bs,true));
    }
  }
  }

  LogDebug("LooperMethod")<<"collected: "<< fastHits.size() <<" hits"<<std::endl;
  //initialize the aggregator
  DataDumper collector(collectorConf_);
  
  //aggregate hits in the phi, 1/R circle plane
  for (uint iH=0;iH!=fastHits.size();++iH){
    collector.collect(fastHits[iH]);
  }

  LogDebug("LooperMethod")<<"aggregation ended"<<std::endl;


  //make the peaks
  collector.makePeaks();

  TrajectorySeedCollection  seedCollection;
  if (makeTC_){
    prod_.tcOut->reserve(collector.peaks_.size());
    seedCollection.reserve(collector.peaks_.size());
  }

  unsigned int nMasked=0;
  std::ofstream positionsInHelix("positionsInHelix.txt");
    
  for(std::vector<aCell*>::iterator iPeak=collector.peak_begin();
      iPeak!=collector.peak_end();++iPeak){
    aCell * peak=*iPeak;

    LogDebug("LooperMethod")<<" A peak cell has: "<<peak->count()<<" elements"<<std::endl;
    edm::OwnVector<TrackingRecHit> recHits;
    SeedingHitSet forSeedCreator;

    bool goodToMask=true;

    positionsInHelix<<" a looper with: "<<peak->count()<<" hits"<<std::endl;
    for (uint iH=0;iH!=peak->count();++iH){
      positionsInHelix<<peak->elements_[iH]->hit_->globalPosition()<<std::endl;
    }

    if (makeTC_){
      goodToMask=false; //until we made a TC out of it
      recHits.reserve(peak->count());
      for (uint iH=0;iH!=peak->count();++iH){
	recHits.push_back( peak->elements_[iH]->hit_->hit()->clone());
	//put only two hits to create the seed
	if (forSeedCreator.size()<2) forSeedCreator.add(peak->elements_[iH]->hit_);
      }

      //from the first and second hits, get an estimate of the  ~z origine for the seed creator
      GlobalTrackingRegion region (0.050, /* pt min */
				   bs, /*zExtrapolation,*/
				   0.5, /* origine radius */
				   5 /* orignine half length */
				   );
      // make a state from the helix state on the surface of the first hit.
      const TrajectorySeed * newSeed = aCreator->trajectorySeed(seedCollection,
								forSeedCreator,
								region,iSetup);
      if (newSeed){
	const PTrajectoryStateOnDet & state=newSeed->startingState();
	prod_.tcOut->push_back(TrackCandidate(recHits,
					      *newSeed,
					      state));
	goodToMask=true;
      }
      else{
	edm::LogWarning("NoSeed")<<"no seed could be made from the hit set, so no TC added";
      }
    }//insert track candidate
    
    // copy the hits in the given order (already in increasing z |z| lowest first)
    if (goodToMask || maskWithNoTC_){
    std::stringstream text;
    text<<peak->printElements()<<"\n";
    for (uint iH=0;iH!=peak->count();++iH){
      //do the masking
      uint subdetId = DetId(peak->elements_[iH]->id_).subdetId();
      if (subdetId==PixelSubdetector::PixelBarrel || subdetId==PixelSubdetector::PixelEndcap)
	{
	  LogDebug("LooperMethod")<<" in the pixel case"<<std::endl;
	  const SiPixelRecHit * pH=static_cast<const SiPixelRecHit *>(peak->elements_[iH]->hit_->hit());
	  LogTrace("LooperMethod")<<"actively masking:" <<pH<<std::endl;
	  if (!prod_.collectedPixels[pH->cluster().key()])	  nMasked++;
	  prod_.collectedPixels[pH->cluster().key()]=true;
	}//pixel case
      else{
	LogDebug("LooperMethod")<<" in the strip case"<<std::endl;      
	//so far so good
	if (peak->elements_[iH]->single_){
	  const SiStripRecHit2D * rH=static_cast<const SiStripRecHit2D  *>(peak->elements_[iH]->hit_->hit());
	  /*if (!rH){
	    edm::LogError("LooperMethod")<<" not casting back to a 1d rechit. probably projected on the way"<<std::endl;
	    assert(0==1);
	    continue;
	    }*/ //was verified that the static_cast <- dynamic_cast was always functionning
	  LogTrace("LooperMethod")<<"actively masking:" <<rH<<std::endl;
	  if (!prod_.collectedStrips[rH->cluster().key()]) nMasked++;
	  prod_.collectedStrips[rH->cluster().key()]=true;
	}
	else{
	  const SiStripMatchedRecHit2D * mH=dynamic_cast<const SiStripMatchedRecHit2D *>(peak->elements_[iH]->hit_->hit());
	  /*	  if (!mH){
	    edm::LogError("LooperMethod")<<" not casting back to a 2d rechit. probably projected on the way"<<std::endl;
	    assert(0==1);
	    continue;
	    }*/ //was verified that the static_cast <- dynamic_cast was always functionning
	  LogTrace("LooperMethod")<<"actively masking:" <<mH<<std::endl;
	  if (!prod_.collectedStrips[mH->stereoHit()->cluster().key()]) nMasked++;
	  if (!prod_.collectedStrips[mH->monoHit()->cluster().key()]) nMasked++;
	  prod_.collectedStrips[mH->stereoHit()->cluster().key()]=true;
	  prod_.collectedStrips[mH->monoHit()->cluster().key()]=true;
	}
      }//strip case


      //cheat, just for display
      if (associator){
	std::vector<PSimHit> simHits = associator->associateHit(*peak->elements_[iH]->hit_->hit());
	text<<" Kill a hit (by using it in the track candidate) at: "<<peak->elements_[iH]->hit_->geographicalId().rawId();
	for (uint ips=0;ips!=simHits.size();++ips)
	  text<<"\n\t[P] pType: "<<simHits[ips].particleType()
	      <<", pId: "<<simHits[ips].trackId() 
	      <<", process: "<<simHits[ips].processType()
	      <<", event: "<<simHits[ips].eventId().rawId();
	text<<"\n";
      }
      LogDebug("LooperMethod")<<text.str();
    }
    }//good to mask

  }//loop on peaks
  
  if (associator)    delete associator;

  if (makeTC_)
    edm::LogError("LooperMethod")<<collector.nPeaks() <<" loopers identified, leading to "<<prod_.tcOut->size()<<" track candidates, containing " <<nMasked<<" hits to be masked, with duplicates, which does not matter"<<std::endl;
  else
    edm::LogError("LooperMethod")<<collector.nPeaks() <<" loopers identified, leading to "<<nMasked<<" hits to be masked, with duplicates, which does not matter"<<std::endl;

  positionsInHelix.close();
}

bool LooperClusterRemoverMethod::PerJet::inZone(edm::Handle<edm::View<reco::Candidate> > & h_jets,
						TransientTrackingRecHit::RecHitPointer & h){
  bool inIt=false;
  for (unsigned int iJet =0; iJet < h_jets->size(); ++iJet)
    {
      const reco::Candidate & myJet = (*h_jets)[iJet];
      float dR=reco::deltaR(myJet.eta(),myJet.phi(),
			    h->globalPosition().eta(),h->globalPosition().phi());
      if (dR<coneSize_) {
	inIt=true;
	break;
      }
    }
  return inIt;
}
void LooperClusterRemoverMethod::PerJet:: run(edm::Event&iEvent, const edm::EventSetup&iSetup,
					      LooperClusterRemover::products &prod_)
{
  //retrieve previous mask
  
  //retrieve the jets
  edm::Handle<edm::View<reco::Candidate> > h_jets;
  iEvent.getByLabel(jets_, h_jets);

  edm::Handle<SiPixelRecHitCollection> pixelHits;
  edm::Handle<SiStripMatchedRecHit2DCollection> matchedHits;
  //get the products
  bool doPix=iEvent.getByLabel(pixelRecHits_,pixelHits);
  bool doStrip=iEvent.getByLabel(stripRecHits_,matchedHits);


  edm::Handle<LooperClusterRemover::PixelMaskContainer> oldPxlMask;
  edm::Handle<LooperClusterRemover::StripMaskContainer> oldStrMask;

  edm::ESHandle<TransientTrackingRecHitBuilder> builder;
  iSetup.get<TransientRecHitRecord>().get("WithTrackAngle",builder);


  //register what is the previous masking
  iEvent.getByLabel(previousMask_,oldPxlMask);
  oldPxlMask->copyMaskTo(prod_.collectedPixels);
  iEvent.getByLabel(previousMask_,oldStrMask);
  oldStrMask->copyMaskTo(prod_.collectedStrips);
  

  if (doPix){
  //  uint seen=0,reject=0;
  //increment with those not in the jet zone
  for ( SiPixelRecHitCollection::const_iterator dPxIt=pixelHits->begin();
	dPxIt!=pixelHits->end();++dPxIt)  {
    for (SiPixelRecHitCollection::DetSet::const_iterator pxIt=dPxIt->begin();
	 pxIt!=dPxIt->end();++pxIt) {
      TransientTrackingRecHit::RecHitPointer h=builder->build(&*pxIt); 
      //      seen++;
      if (!inZone(h_jets,h)){
	//mask it
	prod_.collectedPixels[pxIt->cluster().key()]=true;
	//	reject++;
      }
    }
  }
  //  std::cout<<" from "<<seen<<" reject "<<reject<<std::endl;
  }
  

  if (doStrip){  
  //  seen=0;reject=0;	
  //increment with those not in the jet zone
  for (SiStripMatchedRecHit2DCollection::const_iterator dStIt=matchedHits->begin();
       dStIt!=matchedHits->end();++dStIt){
    for (SiStripMatchedRecHit2DCollection::DetSet::const_iterator stIt=dStIt->begin();
	 stIt!=dStIt->end();++stIt){
      TransientTrackingRecHit::RecHitPointer h=builder->build(&*stIt);
      //      seen++;
      if (!inZone(h_jets,h)){
	//mask it        
	//	reject++;
	prod_.collectedStrips[stIt->stereoHit()->cluster().key()]=true;
	prod_.collectedStrips[stIt->monoHit()->cluster().key()]=true;
      }
    }
  }

  //  std::cout<<" from "<<seen<<" reject "<<reject<<std::endl;
  }

}
