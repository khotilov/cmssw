// -*- C++ -*-
//
// Package:    SiStripMonitorCluster
// Class:      SiStripMonitorHLT
// 
//class SiStripMonitorHLT SiStripMonitorHLT.cc DQM/SiStripMonitorCluster/src/SiStripMonitorHLT.cc
#include <vector>

#include <numeric>
#include <iostream>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"

#include "DQM/SiStripMonitorCluster/interface/SiStripMonitorHLT.h"
#include "DQMServices/Core/interface/DQMStore.h"


SiStripMonitorHLT::SiStripMonitorHLT(const edm::ParameterSet& iConfig)
{
  HLTDirectory="HLTResults";
  dqmStore_  = edm::Service<DQMStore>().operator->();
  conf_ = iConfig;
}


void SiStripMonitorHLT::beginJob(){

  dqmStore_->setCurrentFolder(HLTDirectory);
  std::string HLTProducer = conf_.getParameter<std::string>("HLTProducer");
  HLTDecision = dqmStore_->book1D(HLTProducer+"_HLTDecision", HLTProducer+"HLTDecision", 2, -0.5, 1.5);
  // all
  SumOfClusterCharges_all = dqmStore_->book1D("SumOfClusterCharges_all", "SumOfClusterCharges_all", 50, 0, 2000);
  ChargeOfEachClusterTIB_all = dqmStore_->book1D("ChargeOfEachClusterTIB_all", "ChargeOfEachClusterTIB_all", 400, -0.5, 400.5);
  ChargeOfEachClusterTOB_all = dqmStore_->book1D("ChargeOfEachClusterTOB_all", "ChargeOfEachClusterTOB_all", 400, -0.5, 400.5);
  ChargeOfEachClusterTEC_all = dqmStore_->book1D("ChargeOfEachClusterTEC_all", "ChargeOfEachClusterTEC_all", 400, -0.5, 400.5);
  NumberOfClustersAboveThreshold_all = dqmStore_->book1D("NumberOfClustersAboveThreshold_all", "NumberOfClustersAboveThreshold_all", 30, 30.5, 60.5);
  // 31 = TIB2, 32 = TIB2, 33 = TIB3, 51 = TOB1, 52=TOB2, 60 = TEC
  // accepted from HLT
  SumOfClusterCharges_hlt = dqmStore_->book1D("SumOfClusterCharges_hlt", "SumOfClusterCharges_hlt", 50, 0, 2000);
  ChargeOfEachClusterTIB_hlt = dqmStore_->book1D("ChargeOfEachClusterTIB_hlt", "ChargeOfEachClusterTIB_hlt", 400, -0.5, 400.5);
  ChargeOfEachClusterTOB_hlt = dqmStore_->book1D("ChargeOfEachClusterTOB_hlt", "ChargeOfEachClusterTOB_hlt", 400, -0.5, 400.5);
  ChargeOfEachClusterTEC_hlt = dqmStore_->book1D("ChargeOfEachClusterTEC_hlt", "ChargeOfEachClusterTEC_hlt", 400, -0.5, 400.5);
  NumberOfClustersAboveThreshold_hlt = dqmStore_->book1D("NumberOfClustersAboveThreshold_hlt", "NumberOfClustersAboveThreshold_hlt", 30, 30.5, 60.5);
}

void SiStripMonitorHLT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // get from event
  std::string HLTProducer = conf_.getParameter<std::string>("HLTProducer");
  edm::Handle<int> filter_decision; iEvent.getByLabel(HLTProducer, "", filter_decision); // filter decision
  edm::Handle<uint> sum_of_clustch; iEvent.getByLabel(HLTProducer, "", sum_of_clustch); // sum of cluster charges
  // first element of pair: layer: TIB1, ...., TEC; second element: nr of clusters above threshold
  edm::Handle<std::map<uint,std::vector<SiStripCluster> > > clusters_in_subcomponents;
  if(HLTProducer=="ClusterMTCCFilter") iEvent.getByLabel(HLTProducer, "", clusters_in_subcomponents);

  // trigger decision
  HLTDecision->Fill(*filter_decision);

  // sum of charges of clusters
  SumOfClusterCharges_all->Fill(*sum_of_clustch);
  if(*filter_decision) SumOfClusterCharges_hlt->Fill(*sum_of_clustch);

  //clusters in different layers
  if(HLTProducer=="ClusterMTCCFilter"){
    // loop over layers ("subcomponents")
    for(std::map<uint,std::vector<SiStripCluster> >::const_iterator it = clusters_in_subcomponents->begin(); it != clusters_in_subcomponents->end(); it++){
      int generalized_layer = it->first;
      std::vector<SiStripCluster> theclusters = it->second;
      NumberOfClustersAboveThreshold_all->Fill( generalized_layer, theclusters.size() ); // number of clusters in this generalized layer
      if(*filter_decision) NumberOfClustersAboveThreshold_hlt->Fill( generalized_layer, theclusters.size() );
      //loop over clusters (and detids)
      for(std::vector<SiStripCluster>::const_iterator icluster = theclusters.begin(); icluster != theclusters.end(); icluster++){
        // calculate sum of amplitudes
        unsigned int amplclus=0;
        for(std::vector<uint8_t>::const_iterator ia=icluster->amplitudes().begin(); ia!=icluster->amplitudes().end(); ia++) {
          if ((*ia)>0) amplclus+=(*ia); // why should this be negative?
        }
        if(generalized_layer==31 || generalized_layer==32 || generalized_layer==33){ // you can also ask the detid here whether is TIB
          ChargeOfEachClusterTIB_all->Fill(amplclus,1.);
          if(*filter_decision) ChargeOfEachClusterTIB_hlt->Fill(amplclus,1.);
        }
        if(generalized_layer==51 || generalized_layer==52){
          ChargeOfEachClusterTOB_all->Fill(amplclus,1.);
          if(*filter_decision) ChargeOfEachClusterTOB_hlt->Fill(amplclus,1.);
        }
        if(generalized_layer==60 ){
          ChargeOfEachClusterTEC_all->Fill(amplclus,1.);
          if(*filter_decision) ChargeOfEachClusterTEC_hlt->Fill(amplclus,1.);
        }
      }
    }
  }
}

void SiStripMonitorHLT::endJob(void){
  edm::LogInfo("DQM|SiStripMonitorHLT")<<"Events rejected/accepted "<<HLTDecision->getBinContent(1)<<"/"<<HLTDecision->getBinContent(2);
  bool outputMEsInRootFile = conf_.getParameter<bool>("OutputMEsInRootFile");
  std::string outputFileName = conf_.getParameter<std::string>("OutputFileName");
  if(outputMEsInRootFile){
    dqmStore_->save(outputFileName);
  }

  // delete MEs
//  LogInfo("SiStripTkDQM|SiStripMonitorHLT")<<"pwd="<<dqmStore_->pwd();
////  std::string folder_to_delete = dqmStore_->monitorDirName + "/" + HLTDirectory;
//  dqmStore_->cd();
//  std::string folder_to_delete = HLTDirectory;
//  LogInfo("SiStripTkDQM|SiStripMonitorHLT")<<" Removing whole directory "<<folder_to_delete;
//  dqmStore_->rmdir(folder_to_delete);

}

