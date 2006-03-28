#ifndef ROADSEARCHCLOUD_CLASSES_H
#define ROADSEARCHCLOUD_CLASSES_H

#include "DataFormats/RoadSearchCloud/interface/RoadSearchCloudCollection.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include <vector>

namespace {
  namespace {
    edm::Wrapper<RoadSearchCloudCollection> roadSearchCloudCollectionWrapper;
    edm::Ref<std::vector<RoadSearchCloud>, RoadSearchCloud> roadSearchCloudRef;
    edm::RefVector<std::vector<RoadSearchCloud>, RoadSearchCloud> roadSearchCloudRefVector;

  }
}

#endif // ROADSEARCHCLOUD_CLASSES_H
