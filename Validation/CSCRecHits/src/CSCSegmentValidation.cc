#include "Validation/CSCRecHits/src/CSCSegmentValidation.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include <algorithm>



CSCSegmentValidation::CSCSegmentValidation(DaqMonitorBEInterface* dbe, const edm::InputTag & inputTag)
: CSCBaseValidation(dbe, inputTag),
  theLayerHitsPerChamber(),
  theChamberSegmentMap(),
  theShowerThreshold(10),
  theNPerEventPlot( dbe_->book1D("CSCSegmentsPerEvent", "Number of CSC segments per event", 100, 0, 50) ),
  theNRecHitsPlot( dbe_->book1D("CSCRecHitsPerSegment", "Number of CSC rec hits per segment" , 8, 0, 7) ),
  theNPerChamberTypePlot( dbe_->book1D("CSCSegmentsPerChamberType", "Number of CSC segments per chamber type", 11, 0, 10) ),
  theTypePlot4HitsNoShower( dbe_->book1D("CSCSegments4HitsNoShower", "", 100, 0, 10) ),
  theTypePlot4HitsNoShowerSeg( dbe_->book1D("CSCSegments4HitsNoShowerSeg", "", 100, 0, 10) ),
  theTypePlot4HitsShower( dbe_->book1D("CSCSegments4HitsShower", "", 100, 0, 10) ),
  theTypePlot4HitsShowerSeg( dbe_->book1D("CSCSegments4HitsShowerSeg", "", 100, 0, 10) ),
  theTypePlot5HitsNoShower( dbe_->book1D("CSCSegments5HitsNoShower", "", 100, 0, 10) ),
  theTypePlot5HitsNoShowerSeg( dbe_->book1D("CSCSegments5HitsNoShowerSeg", "", 100, 0, 10) ),
  theTypePlot5HitsShower( dbe_->book1D("CSCSegments5HitsShower", "", 100, 0, 10) ),
  theTypePlot5HitsShowerSeg( dbe_->book1D("CSCSegments5HitsShowerSeg", "", 100, 0, 10) ),
  theTypePlot6HitsNoShower( dbe_->book1D("CSCSegments6HitsNoShower", "", 100, 0, 10) ),
  theTypePlot6HitsNoShowerSeg( dbe_->book1D("CSCSegments6HitsNoShowerSeg", "", 100, 0, 10) ),
  theTypePlot6HitsShower( dbe_->book1D("CSCSegments6HitsShower", "", 100, 0, 10) ),
  theTypePlot6HitsShowerSeg( dbe_->book1D("CSCSegments6HitsShowerSeg", "", 100, 0, 10) )
{
   dbe_->setCurrentFolder("CSCRecHitTask");

   for(int i = 0; i < 10; ++i)
  {
    char title1[200], title2[200];
    sprintf(title1, "CSCSegmentResolution%d", i+1);
    sprintf(title2, "CSCSegmentPull%d", i+1);
    theResolutionPlots[i] = dbe_->book1D(title1, title1, 100, -1, 1);
    thePullPlots[i] = dbe_->book1D(title2, title2, 100, -1, 1);
  }
}

void CSCSegmentValidation::analyze(const edm::Event&e, const edm::EventSetup& eventSetup)
{
  // get the collection of CSCRecHsegmentItrD
  edm::Handle<CSCSegmentCollection> hRecHits;
  e.getByLabel(theInputTag, hRecHits);
  const CSCSegmentCollection * cscRecHits = hRecHits.product();

  theChamberSegmentMap.clear();
  unsigned nPerEvent = 0;
  for(CSCSegmentCollection::const_iterator segmentItr = cscRecHits->begin(); 
      segmentItr != cscRecHits->end(); segmentItr++) 
  {
    ++nPerEvent;
    int detId = segmentItr->geographicalId().rawId();
    int chamberType = whatChamberType(detId);

    theNRecHitsPlot->Fill(segmentItr->nRecHits());
    theNPerChamberTypePlot->Fill(chamberType);
    theChamberSegmentMap[detId].push_back(*segmentItr);
  }

  theNPerEventPlot->Fill(nPerEvent);

  fillLayerHitsPerChamber();
  fillEfficiencyPlots();
}


void CSCSegmentValidation::fillEfficiencyPlots()
{
    // now plot efficiency by looping over all chambers with hits
  for(ChamberHitMap::const_iterator mapItr = theLayerHitsPerChamber.begin(),
      mapEnd = theLayerHitsPerChamber.end();
      mapItr != mapEnd;
      ++mapItr)
  {
    int chamberId = mapItr->first;
    int nHitsInChamber = mapItr->second.size();
    bool isShower = (nHitsInChamber > theShowerThreshold);
    bool hasSeg = hasSegment(chamberId);
    int chamberType = whatChamberType(chamberId);
    // find how many layers were hit in this chamber
    std::vector<int> v = mapItr->second;
    std::sort(v.begin(), v.end());
    // maybe can just count
    v.erase(std::unique(v.begin(), v.end()), v.end());
    int nLayersHit = v.size();

std::cout << "SEGMENT VAL NLAYERS HIT " << nLayersHit << " chambertype " << chamberType << " HAS SEG " << hasSeg << " isshower "<< isShower << " nhitsinchamber " << nHitsInChamber << std::endl;
    if(nLayersHit == 4)
    {

      if(isShower) theTypePlot4HitsShower->Fill(chamberType);
      else         theTypePlot4HitsNoShower->Fill(chamberType);

      if(hasSeg) 
      {
        if(isShower) theTypePlot4HitsShowerSeg->Fill(chamberType);
        else         theTypePlot4HitsNoShowerSeg->Fill(chamberType);
      } 
    }
  }
}

bool CSCSegmentValidation::hasSegment(int chamberId) const
{
  return (theChamberSegmentMap.find(chamberId) != theChamberSegmentMap.end());
}


int CSCSegmentValidation::whatChamberType(int detId)
{
  CSCDetId cscDetId(detId);
  return CSCChamberSpecs::whatChamberType(cscDetId.station(), cscDetId.ring());
}


void CSCSegmentValidation::plotResolution(const PSimHit & simHit, const CSCSegment & segment,
                                         const CSCLayer * layer, int chamberType)
{
  GlobalPoint simHitPos = layer->toGlobal(simHit.localPosition());
  GlobalPoint segmentPos = layer->toGlobal(segment.localPosition());

  double dphi = segmentPos.phi() - simHitPos.phi();
  double rdphi = segmentPos.perp() * dphi;
  theResolutionPlots[chamberType-1]->Fill( rdphi );
  thePullPlots[chamberType-1]->Fill( rdphi/segment.localPositionError().xx() );
}


void CSCSegmentValidation::fillLayerHitsPerChamber()
{
  theLayerHitsPerChamber.clear();
  std::vector<int> layersHit = theSimHitMap->detsWithHits();
  for(std::vector<int>::const_iterator layerItr = layersHit.begin(), 
      layersHitEnd = layersHit.end();
      layerItr != layersHitEnd;
      ++layerItr)
  {
    CSCDetId layerId(*layerItr);
    CSCDetId chamberId = layerId.chamberId();
    int nhits = theSimHitMap->hits(*layerItr).size();
    // multiple entries, so we can see showers
    for(int i = 0; i < nhits; ++i) {
      theLayerHitsPerChamber[chamberId.rawId()].push_back(*layerItr);
    }
  }

}



