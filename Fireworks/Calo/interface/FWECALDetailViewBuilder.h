#include "Rtypes.h"
#include <map>
#include <vector>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

namespace fwlite {
   class Event;
}

class DetIdToMatrix;
class TEveCaloDataVec;
class TEveCaloLego;

// Less than operator for sorting clusters according to eta
class superClsterEtaLess : public std::binary_function<const reco::CaloCluster&, const reco::CaloCluster&, bool>
{
public:
   bool operator()(const reco::CaloCluster &lhs, const reco::CaloCluster &rhs)
   {
      return ( lhs.eta() < rhs.eta()) ;
   }
};

// builder class for ecal detail view
class FWECALDetailViewBuilder {

public:

   // construct an ecal detail view builder
   // the arguments are the event, a pointer to geometry object,
   // the eta and phi position to show,
   // the half width of the region (in indices, e.g. iEta) and
   // the default color for the hits.
   FWECALDetailViewBuilder(const fwlite::Event *event, const DetIdToMatrix* geom,
                           float eta, float phi, int size = 50,
                           Color_t defaultColor = kMagenta)
      : m_event(event), m_geom(geom),
        m_eta(eta), m_phi(phi), m_size(size),
        m_defaultColor(defaultColor){
   }

   // draw the ecal information with the preset colors
   // (if any colors have been preset)
   TEveCaloLego* build();

   // set colors of some predefined detids
   void setColor(Color_t color, const std::vector<DetId> &detIds);

   // show superclusters using two alternating colors
   // to make adjacent clusters visible
   void showSuperClusters(Color_t color1, Color_t color2);

   // show a specific supercluster in a specific color
   void showSuperCluster(const reco::SuperCluster &cluster, Color_t color);

private:

   // fill data
   void fillData(const EcalRecHitCollection *hits,
                 TEveCaloDataVec *data);
   const fwlite::Event *m_event;                                // the event
   const DetIdToMatrix* m_geom;                                 // the geometry
   float m_eta;                                                 // eta position view centred on
   float m_phi;                                                 // phi position view centred on
   int m_size;                                                  // view half width in number of crystals
   Color_t m_defaultColor;                                      // default color for crystals

   // for keeping track of what det id goes in what slice
   std::map<DetId, int> m_detIdsToColor;

   // for keeping track of the colors to use for each slice
   std::vector<Color_t> m_colors;

   // sorting function to sort super clusters by eta.
   static bool superClusterEtaLess(const reco::CaloCluster &lhs, const reco::CaloCluster &rhs)
   {
      return ( lhs.eta() < rhs.eta());
   }

};
