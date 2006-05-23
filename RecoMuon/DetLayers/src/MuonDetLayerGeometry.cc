/** \file
 *
 *  $Date: 2006/05/03 15:20:09 $
 *  $Revision: 1.4 $
 *  \author N. Amapane - CERN
 */

#include <RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h>

#include <FWCore/Utilities/interface/Exception.h>
#include <TrackingTools/DetLayers/interface/DetLayer.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/MuonDetId/interface/DTChamberId.h>

using namespace std;

MuonDetLayerGeometry::MuonDetLayerGeometry() {}

MuonDetLayerGeometry::~MuonDetLayerGeometry(){}

void MuonDetLayerGeometry::addCSCLayers(pair<vector<DetLayer*>, vector<DetLayer*> > csclayers) {
    
    vector<DetLayer*>::const_iterator it;
    for(it=csclayers.first.begin(); it!=csclayers.first.end(); it++) {
        cscLayers_fw.push_back(*it);
        cscLayers_all.push_back(*it);
        allForward.push_back(*it);
        allEndcap.push_back(*it);
        allDetLayers.push_back(*it);

	detLayersMap[ makeDetLayerId(*it) ] = *it;
    }
    
    for(it=csclayers.second.begin(); it!=csclayers.second.end(); it++) {
        cscLayers_bk.push_back(*it);
        cscLayers_all.push_back(*it);
        allBackward.push_back(*it);
        allEndcap.push_back(*it);
        allDetLayers.push_back(*it);

	detLayersMap[ makeDetLayerId(*it) ] = *it;
    }    
}    

/*
void MuonDetLayerGeometry::addRPCLayers(pair<vector<DetLayer*>, vector<DetLayer*> > csclayers) {
    
    cscLayers_fw = csclayers.first;
    cscLayers_bg = csclayers.second;
}    
*/
void MuonDetLayerGeometry::addDTLayers(vector<DetLayer*> dtlayers) {

    vector<DetLayer*>::const_iterator it;
    for(it=dtlayers.begin(); it!=dtlayers.end(); it++) {
        dtLayers.push_back(*it);
        allBarrel.push_back(*it);
        allDetLayers.push_back(*it);

	detLayersMap[ makeDetLayerId(*it) ] = *it;
    }
}    

DetId MuonDetLayerGeometry::makeDetLayerId(DetLayer* detLayer){

  if(detLayer->module() ==  csc){
    CSCDetId id( detLayer->basicComponents().front()->geographicalId().rawId() ) ;
    return CSCDetId(id.endcap(),id.station(),0,0,0);
  }
  else if(detLayer->module() == dt){
    DTChamberId id( detLayer->basicComponents().front()->geographicalId().rawId() ) ;
    return  DTChamberId(0,id.station(),0);
  }
  else throw cms::Exception("InvalidModuleIdentification"); // << detLayer->module();
}


const vector<DetLayer*>& 
MuonDetLayerGeometry::allDTLayers() const {    
    return dtLayers; 
}

const vector<DetLayer*>&
MuonDetLayerGeometry::allCSCLayers() const {
    return cscLayers_all;
}


const vector<DetLayer*>&
MuonDetLayerGeometry::forwardCSCLayers() const {
    return cscLayers_fw;
}


const vector<DetLayer*>& 
MuonDetLayerGeometry::backwardCSCLayers() const {
    return cscLayers_bk;
}


const vector<DetLayer*>& 
MuonDetLayerGeometry::allRPCLayers() const {
    return rpcLayers_all;    
}


const vector<DetLayer*>& 
MuonDetLayerGeometry::barrelRPCLayers() const {
    return rpcLayers_barrel; 
}


const vector<DetLayer*>& 
MuonDetLayerGeometry::endcapRPCLayers() const {
    return rpcLayers_endcap;    
}


const vector<DetLayer*>& 
MuonDetLayerGeometry::forwardRPCLayers() const {
     return rpcLayers_fw; 
}


const vector<DetLayer*>& 
MuonDetLayerGeometry::backwardRPCLayers() const {
    return rpcLayers_bk; 
}


const vector<DetLayer*> 
MuonDetLayerGeometry::allLayers() const {
    return allDetLayers;    
}    


const vector<DetLayer*> 
MuonDetLayerGeometry::allBarrelLayers() const {
    return allBarrel;    
}    

const vector<DetLayer*> 
MuonDetLayerGeometry::allEndcapLayers() const {
    return allEndcap;    
}    


const vector<DetLayer*> 
MuonDetLayerGeometry::allForwardLayers() const {
    return allForward;    
}    


const vector<DetLayer*> 
MuonDetLayerGeometry::allBackwardLayers() const {
    return allBackward;    
}    

DetLayer* MuonDetLayerGeometry::idToLayer(DetId detId){
  
  if(detId.subdetId() == MuonSubdetId::CSC){
    CSCDetId cscId( detId.rawId() );
    CSCDetId id(cscId.endcap(),cscId.station(),0,0,0);
    return detLayersMap[ DetId(id.rawId()) ]; 
  }
  else if (detId.subdetId() == MuonSubdetId::DT){
    DTChamberId dtId( detId.rawId() );
    DTChamberId id(0,dtId.station(),0);
    return detLayersMap[ DetId(id.rawId()) ]; 
  }
  else throw cms::Exception("InvalidSubdetId")<< detId.subdetId();
}


