#include "Fireworks/Core/interface/FWSimpleProxyBuilderTemplate.h"
#include "Fireworks/Core/interface/FWEventItem.h"
#include "Fireworks/Core/interface/DetIdToMatrix.h"
#include "Fireworks/Calo/interface/CaloUtils.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

class FWCaloClusterProxyBuilder : public FWSimpleProxyBuilderTemplate<reco::CaloCluster>
{
public:
   FWCaloClusterProxyBuilder(void) {}  
   virtual ~FWCaloClusterProxyBuilder(void) {}

   REGISTER_PROXYBUILDER_METHODS();

private:
   FWCaloClusterProxyBuilder(const FWCaloClusterProxyBuilder&); 			// stop default
   const FWCaloClusterProxyBuilder& operator=(const FWCaloClusterProxyBuilder&); 	// stop default

   void build(const reco::CaloCluster& iData, unsigned int iIndex, TEveElement& oItemHolder, const FWViewContext*);
};

void
FWCaloClusterProxyBuilder::build(const reco::CaloCluster& iData, unsigned int iIndex, TEveElement& oItemHolder, const FWViewContext*) 
{
   std::vector<std::pair<DetId, float> > clusterDetIds = iData.hitsAndFractions ();
   Float_t scale = 10.0; 	// FIXME: The scale should be taken form somewhere else
   
   for(std::vector<std::pair<DetId, float> >::iterator id = clusterDetIds.begin (), idend = clusterDetIds.end ();
       id != idend; ++id)
   {
      std::vector<TEveVector> corners = item()->getGeom()->getPoints((*id).first);
      if( corners.empty() ) {
	 continue;
      }

      fireworks::drawEnergyTower3D(corners, (*id).second * scale, &oItemHolder, this);
   }
}

REGISTER_FWPROXYBUILDER(FWCaloClusterProxyBuilder, reco::CaloCluster, "Calo Cluster", FWViewType::kISpyBit );
