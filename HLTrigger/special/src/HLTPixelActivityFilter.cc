#include "HLTrigger/HLTcore/interface/HLTFilter.h"

//
// class declaration
//

class HLTPixelActivityFilter : public HLTFilter {
public:
  explicit HLTPixelActivityFilter(const edm::ParameterSet&);
  ~HLTPixelActivityFilter();

private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);

  edm::InputTag inputTag_;          // input tag identifying product containing pixel clusters
  bool          saveTag_;           // whether to save this tag
  unsigned int  min_clusters_;      // minimum number of clusters

};

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

//
// constructors and destructor
//
 
HLTPixelActivityFilter::HLTPixelActivityFilter(const edm::ParameterSet& config) :
  inputTag_     (config.getParameter<edm::InputTag>("inputTag")),
  saveTag_      (config.getUntrackedParameter<bool>("saveTag", false)),
  min_clusters_ (config.getParameter<unsigned int>("minClusters"))
{
  LogDebug("") << "Using the " << inputTag_ << " input collection";
  LogDebug("") << "Requesting " << min_clusters_ << " clusters";

  // register your products
  produces<trigger::TriggerFilterObjectWithRefs>();
}

HLTPixelActivityFilter::~HLTPixelActivityFilter()
{
}

//
// member functions
//

// ------------ method called to produce the data  ------------
bool HLTPixelActivityFilter::filter(edm::Event& event, const edm::EventSetup& iSetup)
{
  // All HLT filters must create and fill an HLT filter object,
  // recording any reconstructed physics objects satisfying (or not)
  // this HLT filter, and place it in the Event.

  // The filter object
  std::auto_ptr<trigger::TriggerFilterObjectWithRefs> filterobject (new trigger::TriggerFilterObjectWithRefs(path(),module()));
  if (saveTag_) filterobject->addCollectionTag(inputTag_);

  // get hold of products from Event
  edm::Handle<edmNew::DetSetVector<SiPixelCluster> > clusterColl;
  event.getByLabel(inputTag_, clusterColl);

  unsigned int clusterSize = clusterColl->size();
  LogDebug("") << "Number of clusters accepted: " << clusterSize;
  bool accept = (clusterSize >= min_clusters_);

  // put filter object into the Event
  event.put(filterobject);

  // return with final filter decision
  return accept;
}

// define as a framework module
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTPixelActivityFilter);
