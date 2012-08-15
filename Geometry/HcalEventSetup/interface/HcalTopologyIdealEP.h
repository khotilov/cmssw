#ifndef GEOMETRY_HCALEVENTSETUP_HCALTOPOLOGYIDEALEP_H
#define GEOMETRY_HCALEVENTSETUP_HCALTOPOLOGYIDEALEP_H 1


// system include files
#include <memory>
#include "boost/shared_ptr.hpp"

// user include files
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"

//
// class decleration
//

class HcalTopologyIdealEP : public edm::ESProducer {

public:
  HcalTopologyIdealEP(const edm::ParameterSet&);
  ~HcalTopologyIdealEP();

  typedef boost::shared_ptr<HcalTopology> ReturnType;

  ReturnType produce(const IdealGeometryRecord&);

private:
  // ----------member data ---------------------------
  std::string m_restrictions;
  HcalTopology::Mode  m_mode;
  int m_maxDepthHB;
  int m_maxDepthHE;
    
  // can be specified in the config
  struct RingSegmentation {
    int ring;
    std::vector<int> segmentation;
  };
  std::vector<RingSegmentation> m_segmentation;
};
#endif
