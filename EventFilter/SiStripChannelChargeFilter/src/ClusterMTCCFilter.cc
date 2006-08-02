#include "EventFilter/SiStripChannelChargeFilter/interface/ClusterMTCCFilter.h"
#include "FWCore/Framework/interface/Handle.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


namespace cms
{

ClusterMTCCFilter::ClusterMTCCFilter(const edm::ParameterSet& ps){
   ChargeThresholdTIB=ps.getParameter<int>("ChargeThresholdTIB");
   ChargeThresholdTOB=ps.getParameter<int>("ChargeThresholdTOB");
   ChargeThresholdTEC=ps.getParameter<int>("ChargeThresholdTEC");
   MinClustersDiffComponents=ps.getParameter<int>("MinClustersDiffComponents");
   clusterProducer = ps.getParameter<std::string>("ClusterProducer");
   //
   produces <int>();
   produces <unsigned int >();
   produces < std::map<uint,std::pair<SiStripCluster,uint32_t> > >();
}

bool ClusterMTCCFilter::filter(edm::Event & e, edm::EventSetup const& c) {
  //get SiStripCluster
  edm::Handle< edm::DetSetVector<SiStripCluster> > h;
  e.getByLabel(clusterProducer,h);

  //
  unsigned int sum_of_cluster_charges=0;
  clusters_in_subcomponents.clear();
  // first find all clusters that are over the threshold
  // uint generalized_layer:  31 = TIB1, 32 = TIB2, 33 = TIB3, 50 = TOB, 60 = TEC
  for (edm::DetSetVector<SiStripCluster>::const_iterator it=h->begin();it!=h->end();it++) {
    for(std::vector<SiStripCluster>::const_iterator vit=(it->data).begin(); vit!=(it->data).end(); vit++){
      // calculate sum of amplitudes
      unsigned int amplclus=0;
      for(std::vector<short>::const_iterator ia=vit->amplitudes().begin(); ia!=vit->amplitudes().end(); ia++) {
        if ((*ia)>0) amplclus+=(*ia); // why should this be negative?
      }
      sum_of_cluster_charges += amplclus;
      DetId thedetId = DetId(it->detId());
      uint generalized_layer = 0;
      // apply different thresholds for TIB/TOB/TEC
      if ( ( thedetId.subdetId()==StripSubdetector::TIB && amplclus>ChargeThresholdTIB )
        || ( thedetId.subdetId()==StripSubdetector::TOB && amplclus>ChargeThresholdTOB )
        || ( thedetId.subdetId()==StripSubdetector::TEC && amplclus>ChargeThresholdTEC )
        ){
        if(thedetId.subdetId()==StripSubdetector::TIB){
           TIBDetId ptib = TIBDetId(thedetId.rawId());
           generalized_layer = 10*thedetId.subdetId() + ptib.layer() + ptib.stereo();
	   if (ptib.layer()==2){
	     generalized_layer++;
	     if (ptib.glued()) std::cout<<"WRONGGGG"<<std::endl;
	   }

        }else{
          generalized_layer = 10*thedetId.subdetId();
	  if(thedetId.subdetId()==StripSubdetector::TOB){
	    TOBDetId ptob = TOBDetId(thedetId.rawId());
	    generalized_layer += ptob.layer();
	  }
        }
        clusters_in_subcomponents.insert( std::make_pair(generalized_layer,  std::make_pair(*vit,it->detId())) );
      }
    }
  }

  bool decision=false; // default value, only accept if set true in this loop
  uint nr_of_subcomps_with_clusters=0;
  if( clusters_in_subcomponents.count(31)>0 ) nr_of_subcomps_with_clusters++; // TIB1
  if( clusters_in_subcomponents.count(32)>0 ) nr_of_subcomps_with_clusters++; // TIB2
  if( clusters_in_subcomponents.count(33)>0 ) nr_of_subcomps_with_clusters++; // TIB3
  if( clusters_in_subcomponents.count(51)>0 ) nr_of_subcomps_with_clusters++; // TOB1
  if( clusters_in_subcomponents.count(52)>0 ) nr_of_subcomps_with_clusters++; // TOB2
  if( clusters_in_subcomponents.count(60)>0 ) nr_of_subcomps_with_clusters++; // TEC
  if(
     nr_of_subcomps_with_clusters >= MinClustersDiffComponents // more than 'MinClustersDiffComponents' components have at least 1 cluster
     ) {
      decision = true; // accept event
//      std::cout<<"Nr of layers with clusters: "<<nr_of_subcomps_with_clusters<<" "<< clusters_in_subcomponents.count(31)<<" "<< clusters_in_subcomponents.count(32)<<" "<< clusters_in_subcomponents.count(33)<<" "<< clusters_in_subcomponents.count(51)<<" "<< clusters_in_subcomponents.count(51)<<" "<< clusters_in_subcomponents.count(60)<<std::endl;
  }

  std::auto_ptr< int > output_decision( new int(decision) );
  e.put(output_decision);

  std::auto_ptr< unsigned int > output_sumofcharges( new unsigned int(sum_of_cluster_charges) );
  e.put(output_sumofcharges);

  std::auto_ptr< std::map< uint, std::pair<SiStripCluster,uint32_t> >  > output_clusters(new std::map< uint, std::pair<SiStripCluster,uint32_t> > (clusters_in_subcomponents));
  e.put(output_clusters);

  return decision;

}
}

