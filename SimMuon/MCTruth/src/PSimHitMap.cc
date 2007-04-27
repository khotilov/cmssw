#include "SimMuon/MCTruth/interface/PSimHitMap.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"

PSimHitMap::PSimHitMap(const std::string & collectionName)
:  theCollectionName(collectionName),
   theMap(),
   theEmptyContainer()
{
}


void PSimHitMap::fill(const edm::Event & e)
{
  theMap.clear();
  edm::Handle<CrossingFrame> cf;
  e.getByType(cf);
  
  MixCollection<PSimHit> simHits(cf.product(), theCollectionName);

  // arrange the hits by detUnit
  for(MixCollection<PSimHit>::MixItr hitItr = simHits.begin();
      hitItr != simHits.end(); ++hitItr)
  {
    theMap[hitItr->detUnitId()].push_back(*hitItr);
  }


}


const edm::PSimHitContainer & PSimHitMap::hits(int detId) const
{
  std::map<int, edm::PSimHitContainer>::const_iterator mapItr
    = theMap.find(detId);
  if(mapItr != theMap.end())
  {
    return mapItr->second;
  }
  else
  {
    return theEmptyContainer;
  }
}


std::vector<int> PSimHitMap::detsWithHits() const 
{
  std::vector<int> result;
  result.reserve(theMap.size());
  for(std::map<int, edm::PSimHitContainer>::const_iterator mapItr = theMap.begin(),
      mapEnd = theMap.end();
      mapItr != mapEnd;
      ++mapItr)
  {
    result.push_back(mapItr->first);
  }
  return result;
} 
     

