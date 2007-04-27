#include "SimMuon/MCTruth/interface/MuonTruth.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"


MuonTruth::MuonTruth()
: theSimTrackContainer(0),
  theDigiSimLinks(0),
  theSimHitMap("MuonCSCHits")
{
}

void MuonTruth::eventSetup(const edm::Event & event)
{
  edm::Handle<edm::SimTrackContainer> simTrackCollection;
  event.getByLabel("g4SimHits", simTrackCollection);
  theSimTrackContainer = simTrackCollection.product();

  edm::Handle<DigiSimLinks> digiSimLinks;
  edm::InputTag linksTag("muonCSCDigis" , "MuonCSCStripDigiSimLinks");
  event.getByLabel(linksTag, digiSimLinks);
  theDigiSimLinks = digiSimLinks.product();

  theSimHitMap.fill(event);
}



float MuonTruth::muonFraction()
{
  if(theChargeMap.size() == 0) return 0.;

  float muonCharge = 0.;
  for(std::map<int, float>::const_iterator chargeMapItr = theChargeMap.begin();
      chargeMapItr != theChargeMap.end(); ++chargeMapItr)
  {
    if( abs(particleType(chargeMapItr->first)) == 13)
    {
      muonCharge += chargeMapItr->second;
    }
  }

  return muonCharge / theTotalCharge;
}


std::vector<const PSimHit *> MuonTruth::simHits()
{
  std::vector<const PSimHit *> result;
  for(std::map<int, float>::const_iterator chargeMapItr = theChargeMap.begin();
      chargeMapItr != theChargeMap.end(); ++chargeMapItr)
  {
    std::vector<const PSimHit *> trackHits = hitsFromSimTrack(chargeMapItr->first);
    result.insert(result.end(), trackHits.begin(), trackHits.end());
  }
  return result;
}


std::vector<const PSimHit *> MuonTruth::muonHits()
{
  std::vector<const PSimHit *> result;
  std::vector<const PSimHit *> allHits = simHits();
  std::vector<const PSimHit *>::const_iterator hitItr = allHits.begin(), lastHit = allHits.end();

  for( ; hitItr != lastHit; ++hitItr)
  {
    if(abs((**hitItr).particleType()) == 13)
    {
      result.push_back(*hitItr);
    }
  }
  return result;
}



std::vector<const PSimHit *> MuonTruth::hitsFromSimTrack(int index) const
{
  std::vector<const PSimHit *> result;
  edm::PSimHitContainer hits = theSimHitMap.hits(theDetId);
  edm::PSimHitContainer::const_iterator hitItr = hits.begin(), lastHit = hits.end();

  for( ; hitItr != lastHit; ++hitItr)
  {
    int hitTrack = hitItr->trackId();
    if(hitTrack == index) 
    {
      result.push_back(&*hitItr);
    }
  }
  return result;
}


int MuonTruth::particleType(int simTrack) const
{
  int result = 0;
  std::vector<const PSimHit *> hits = hitsFromSimTrack(simTrack);
  if(!hits.empty())
  {
    result = hits[0]->particleType();
  }
  return result;
}



void MuonTruth::analyze(const CSCRecHit2D & recHit)
{
  theChargeMap.clear();
  theTotalCharge = 0.;
  theDetId = recHit.geographicalId().rawId();

  int nchannels = recHit.channels().size();
  CSCRecHit2D::ADCContainer adcContainer = recHit.adcs();
  for(int idigi = 0; idigi < nchannels; ++idigi)
  {
    int channel = recHit.channels()[idigi];
    float weight = adcContainer[idigi];


    DigiSimLinks::const_iterator layerLinks = theDigiSimLinks->find(theDetId);

    if(layerLinks != theDigiSimLinks->end())
    {
      addChannel(*layerLinks, channel, weight);
    }
  }
}


void MuonTruth::analyze(const CSCStripDigi & stripDigi, int layerId)
{
  theDetId = layerId;
  theChargeMap.clear();
  theTotalCharge = 0.;

  DigiSimLinks::const_iterator layerLinks = theDigiSimLinks->find(theDetId);
  if(layerLinks != theDigiSimLinks->end())
  {
    addChannel(*layerLinks, stripDigi.getStrip(), 1.);
  }
}


void MuonTruth::addChannel(const LayerLinks &layerLinks, int channel, float weight)
{
  LayerLinks::const_iterator linkItr = layerLinks.begin(), 
                             lastLayerLink = layerLinks.end();

  for( ; linkItr != lastLayerLink; ++linkItr)
  {
    int linkChannel = linkItr->channel();
    if(linkChannel == channel)
    {
      float charge = linkItr->fraction() * weight;
      theTotalCharge += charge;
      // see if it's in the map
      int simTrack = linkItr->SimTrackId();
      std::map<int, float>::const_iterator chargeMapItr = theChargeMap.find( simTrack );
      if(chargeMapItr == theChargeMap.end())
      {
        theChargeMap[simTrack] = charge;
      }
      else
      {
        theChargeMap[simTrack] += charge;
      }
    }
  }
}
    

