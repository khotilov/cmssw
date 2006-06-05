#ifndef TrackerRecHit2D_CLASSES_H
#define TrackerRecHit2D_CLASSES_H

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DLocalPosCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DMatchedLocalPosCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include <vector>

namespace {
  namespace {    
    SiStripRecHit2DLocalPos a1;
    SiStripRecHit2DMatchedLocalPos a2;
    SiPixelRecHit b1;
    DetId a3;

    edm::ClonePolicy<SiStripRecHit2DLocalPos>  a4;
    edm::ClonePolicy<SiStripRecHit2DMatchedLocalPos > a5;
    edm::ClonePolicy<SiPixelRecHit> b2;

    edm::OwnVector<SiStripRecHit2DLocalPos,
      edm::ClonePolicy<SiStripRecHit2DLocalPos> > a6;
    edm::OwnVector<SiStripRecHit2DMatchedLocalPos,
      edm::ClonePolicy<SiStripRecHit2DMatchedLocalPos> > a7;
    edm::OwnVector<SiPixelRecHit,
      edm::ClonePolicy<SiPixelRecHit> > b3;
    
    edm::Wrapper< edm::RangeMap<DetId,
      edm::OwnVector<SiStripRecHit2DLocalPos,
      edm::ClonePolicy<SiStripRecHit2DLocalPos> >, 
      edm::ClonePolicy<SiStripRecHit2DLocalPos> > >    siStripRecHit2DLocalPosCollectionWrapper;
    
    edm::Wrapper< edm::RangeMap<DetId,
      edm::OwnVector<SiStripRecHit2DMatchedLocalPos,
      edm::ClonePolicy<SiStripRecHit2DMatchedLocalPos> >, 
      edm::ClonePolicy<SiStripRecHit2DMatchedLocalPos> > > siStripRecHit2DMatchedLocalPosCollectionWrapper;

    edm::Wrapper< edm::RangeMap<DetId,
      edm::OwnVector<SiPixelRecHit,
      edm::ClonePolicy<SiPixelRecHit> >, 
      edm::ClonePolicy<SiPixelRecHit> > >  siPixelRecHitCollectionWrapper;

    edm::Ref<edm::DetSetVector<SiPixelCluster>, 
      SiPixelCluster, 
      edm::refhelper::FindForDetSetVector<SiPixelCluster> > siPixelClusterReference;

    edm::Ref<edm::DetSetVector<SiStripCluster>, 
      SiStripCluster, 
      edm::refhelper::FindForDetSetVector<SiStripCluster> > siStripClusterReference;

    edm::RefBase<std::pair<unsigned int,unsigned int> > p1;
    edm::RefItem<std::pair<unsigned int,unsigned int> > p2;

  }
}

#endif // SISTRIPRECHIT_CLASSES_H
