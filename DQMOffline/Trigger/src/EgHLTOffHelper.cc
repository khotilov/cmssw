#include "DQMOffline/Trigger/interface/EgHLTOffHelper.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

void EgHLTOffHelper::setup(const edm::ParameterSet& conf)
{
  barrelShapeAssocProd_ = conf.getParameter<edm::InputTag>("barrelClusterShapeAssociation");
  endcapShapeAssocProd_ = conf.getParameter<edm::InputTag>("endcapClusterShapeAssociation");

  cuts_.setHighNrgy();
  tagCuts_.setHighNrgy();
  probeCuts_.setPreSel();
}


//this function coverts GsfElectrons to a format which is actually useful to me
void EgHLTOffHelper::fillEgHLTOffEleVec(edm::Handle<reco::GsfElectronCollection> gsfElectrons,std::vector<EgHLTOffEle>& egHLTOffEles)
{
  egHLTOffEles.clear();
  egHLTOffEles.reserve(gsfElectrons->size());
  for(reco::PixelMatchGsfElectronCollection::const_iterator gsfIter=gsfElectrons->begin(); gsfIter!=gsfElectrons->end();++gsfIter){
    //for now use dummy isolation data
    EgHLTOffEle::IsolData isolData;
    isolData.nrTrks=1;
    isolData.ptTrks=1.6;
    isolData.em= 0.42;
    isolData.had=0.42;
    
    //get cluster shape and we're done construction
    const reco::ClusterShape* clusShape = getClusterShape(&*gsfIter);
    egHLTOffEles.push_back(EgHLTOffEle(*gsfIter,clusShape,isolData));
    
    //now we would like to set the cut results
    EgHLTOffEle& ele =  egHLTOffEles.back();
    ele.setTagCutCode(tagCuts_.getCutCode(ele));
    ele.setProbeCutCode(probeCuts_.getCutCode(ele));
    ele.setCutCode(cuts_.getCutCode(ele));
      

  }//end loop over gsf electron collection

}

//so I do this rather needlessly expensive operation once an event instead of once an electron
void EgHLTOffHelper::getHandles(const edm::Event& event)
{
  try{
    event.getByLabel(barrelShapeAssocProd_, clusterShapeHandleBarrel_) ;
  }catch(...){} //the worlds most pointless try, catch pair, damn you CMSSW framework,  DAMN you
  
  try{
    event.getByLabel(endcapShapeAssocProd_, clusterShapeHandleEndcap_) ;
  }catch(...){}

  if (!clusterShapeHandleBarrel_.isValid()) {
    edm::LogError ("EgHLTOffHelper") << "Can't get ECAL barrel Cluster Shape Collection" ; 
  }
  if (!clusterShapeHandleEndcap_.isValid()) {
    edm::LogError ("EgHLTOffHelper") << "Can't get ECAL endcap Cluster Shape Collection" ; 
  }


}



//ripped of from the electronIDAlgo (there must be a better way, I *cannot* believe that there isnt a better way)
//I've made some minor mods for speed and robustness (it could still be faster though)
//I'm sorry for the pain you are about to go through
//in summary it determines where the electron is barrel or endcap and if the clusterShape association map handle is valid for it
//it then looks in the map for the electrons seed cluster and if found, returns a pointer to the shape
//a complication arrises as electrons which are endcap may be classified as in the barrel-endcap gap and therefore have classification 40
//and therefore be labeled barrel someplaces (like here) and endcap others
const reco::ClusterShape* EgHLTOffHelper::getClusterShape(const reco::GsfElectron* electron)
{
  // Find entry in map corresponding to seed BasicCluster of SuperCluster
  reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;
 
  if ( electron->classification() < 100 && clusterShapeHandleBarrel_.isValid() ) {
    const reco::BasicClusterShapeAssociationCollection& barrelClShpMap = *clusterShapeHandleBarrel_;
    reco::SuperClusterRef sclusRef = electron->get<reco::SuperClusterRef> () ;
    seedShpItr = barrelClShpMap.find ( sclusRef->seed () ) ;
    if (seedShpItr!=barrelClShpMap.end()) return &*(seedShpItr->val);
    else if (clusterShapeHandleEndcap_.isValid()){
      const reco::BasicClusterShapeAssociationCollection& endcapClShpMap = *clusterShapeHandleEndcap_;
      seedShpItr = endcapClShpMap.find ( sclusRef->seed ()) ;
      if(seedShpItr!=endcapClShpMap.end()) return &*(seedShpItr->val);
    }//end check of valid endcap cluster shape in barrel section
  } else if(electron->classification()>=100 && clusterShapeHandleEndcap_.isValid()) {
    const reco::BasicClusterShapeAssociationCollection& endcapClShpMap = *clusterShapeHandleEndcap_;
    reco::SuperClusterRef sclusRef = electron->get<reco::SuperClusterRef> () ;
    seedShpItr = endcapClShpMap.find ( sclusRef->seed () ) ;
    if(seedShpItr!=endcapClShpMap.end()) return &*(seedShpItr->val); 
  }//end check of endcap electron with valid shape map
  
  return NULL;
}

