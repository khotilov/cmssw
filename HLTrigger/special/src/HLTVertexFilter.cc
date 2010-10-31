// -*- C++ -*-
//
// Package:    HLTVertexFilter
// Class:      HLTVertexFilter
// 
/**\class HLTVertexFilter HLTVertexFilter.cc

 Description: HLTFilter to accept events with at least a given number of vertices

 Implementation:
     This class implements an HLTFilter to select events with at least
     a certain number of vertices matching some selection criteria.
*/
// Original Author:  Andrea Bocci
//         Created:  Tue Apr 20 12:34:27 CEST 2010


#include <cmath>
#include <math.h>
#include <boost/foreach.hpp>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"

//
// class declaration
//
class HLTVertexFilter : public HLTFilter {
public:
  explicit HLTVertexFilter(const edm::ParameterSet & config);
  ~HLTVertexFilter();
    
private:
  virtual 
  bool filter(edm::Event & event, const edm::EventSetup & setup);

  edm::InputTag m_inputTag;     // input vertex collection
  double        m_minNDoF;      // minimum vertex NDoF
  double        m_maxChi2;      // maximum vertex chi2
  double        m_maxD0;        // maximum transverse distance from the beam
  double        m_maxZ;         // maximum longitudinal distance nominal center of the detector
  unsigned int  m_minVertices;

};

//
// constructors and destructor
//
HLTVertexFilter::HLTVertexFilter(const edm::ParameterSet& config) :
  m_inputTag(config.getParameter<edm::InputTag>("inputTag")),
  m_minNDoF(config.getParameter<double>("minNDoF")),
  m_maxChi2(config.getParameter<double>("maxChi2")),
  m_maxD0(config.getParameter<double>("maxD0")),
  m_maxZ(config.getParameter<double>("maxZ")),
  m_minVertices(config.getParameter<unsigned int>("minVertices"))
{
  // register your products
  produces<trigger::TriggerFilterObjectWithRefs>();
}


HLTVertexFilter::~HLTVertexFilter()
{
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HLTVertexFilter::filter(edm::Event &  event, edm::EventSetup const & setup) {

  // The filter object
  std::auto_ptr<trigger::TriggerFilterObjectWithRefs> filterobject (new trigger::TriggerFilterObjectWithRefs(path(),module()));

  // get hold of collection of objects
  edm::Handle<reco::VertexCollection> vertices;
  event.getByLabel(m_inputTag, vertices);

  unsigned int n = 0;
  if (vertices.isValid()) {
    BOOST_FOREACH(const reco::Vertex & vertex, * vertices) {
      if (vertex.isValid()
          and not vertex.isFake()
          and (vertex.chi2() <= m_maxChi2)
          and (vertex.ndof() >= m_minNDoF)
          and (std::abs(vertex.z()) <= m_maxZ)
          and (::hypot(vertex.x(), vertex.y()) <= m_maxD0)
          )
        ++n;
    }
  }

  // put filter object into the Event
  event.put(filterobject);

  // filter decision
  return (n >= m_minVertices);
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTVertexFilter);
