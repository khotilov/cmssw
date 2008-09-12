#include <FWCore/MessageLogger/interface/MessageLogger.h>

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"

#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "DQM/SiStripMonitorTrack/interface/SiStripMonitorTrack.h"

#include "DQM/SiStripCommon/interface/SiStripHistoId.h"
#include "TMath.h"

static const uint16_t _NUM_SISTRIP_SUBDET_ = 4;
static TString SubDet[_NUM_SISTRIP_SUBDET_]={"TIB","TID","TOB","TEC"};
static std::string flags[2] = {"OnTrack","OffTrack"};

SiStripMonitorTrack::SiStripMonitorTrack(const edm::ParameterSet& conf):
  dbe(edm::Service<DQMStore>().operator->()),
  conf_(conf),
  Cluster_src_( conf.getParameter<edm::InputTag>( "Cluster_src" ) ),
  Mod_On_(conf.getParameter<bool>("Mod_On")),
  Trend_On_(conf.getParameter<bool>("Trend_On")),
  OffHisto_On_(conf.getParameter<bool>("OffHisto_On")),
  RawDigis_On_(conf.getParameter<bool>("RawDigis_On")),
  CCAnalysis_On_(conf.getParameter<bool>("CCAnalysis_On")),
  folder_organizer(), tracksCollection_in_EventTree(true),
  firstEvent(-1)
{
  for(int i=0;i<4;++i) for(int j=0;j<2;++j) NClus[i][j]=0;
  if(OffHisto_On_){
    off_Flag = 2;
  }else{
    off_Flag = 1;
  }
}

//------------------------------------------------------------------------
SiStripMonitorTrack::~SiStripMonitorTrack() { }

//------------------------------------------------------------------------
//void SiStripMonitorTrack::beginJob(const edm::EventSetup& es)

void SiStripMonitorTrack::beginRun(const edm::Run& run,const edm::EventSetup& es)
{
  //get geom 
  es.get<TrackerDigiGeometryRecord>().get( tkgeom );
  edm::LogInfo("SiStripMonitorTrack") << "[SiStripMonitorTrack::beginRun] There are "<<tkgeom->detUnits().size() <<" detectors instantiated in the geometry" << std::endl;  
  es.get<SiStripDetCablingRcd>().get( SiStripDetCabling_ );

  book();
}

//------------------------------------------------------------------------
void SiStripMonitorTrack::endJob(void)
{
  if(conf_.getParameter<bool>("OutputMEsInRootFile")){
    dbe->showDirStructure();
    dbe->save(conf_.getParameter<std::string>("OutputFileName"));
  }
}

// ------------ method called to produce the data  ------------
void SiStripMonitorTrack::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  
  tracksCollection_in_EventTree=true;
  trackAssociatorCollection_in_EventTree=true;
  
  //initialization of global quantities
  edm::LogInfo("SiStripMonitorTrack") << "[SiStripMonitorTrack::analyse]  " << "Run " << e.id().run() << " Event " << e.id().event() << std::endl;
  runNb   = e.id().run();
  eventNb = e.id().event();
  vPSiStripCluster.clear();
  countOn=0;
  countOff=0;
  
  e.getByLabel( Cluster_src_, dsv_SiStripCluster); 
  
  // track input  
  std::string TrackProducer = conf_.getParameter<std::string>("TrackProducer");
  std::string TrackLabel = conf_.getParameter<std::string>("TrackLabel");
  
  e.getByLabel(TrackProducer, TrackLabel, trackCollection);//takes the track collection
 
  if (trackCollection.isValid()){
  }else{
    edm::LogError("SiStripMonitorTrack")<<" Track Collection is not valid !! " << TrackLabel<<std::endl;
    tracksCollection_in_EventTree=false;
  }
  
  // get RawDigi collection
  if(RawDigis_On_) {
    std::string digiProducer = conf_.getParameter<std::string>("RawDigiProducer");
    std::string    digiLabel = conf_.getParameter<std::string>("RawDigiLabel");
    e.getByLabel(digiProducer,digiLabel,dsv_SiStripRawDigi);
    if( ! dsv_SiStripRawDigi.failedToGet() ){
      // digi check
      edm::LogInfo("SiStripMonitorTrack") << "[SiStripMonitorTrack::analyse]  " << "Digi " << digiLabel
					  << " found " << dsv_SiStripRawDigi->size() << "\n";
      /* digi check
	 for( edm::DetSetVector<SiStripRawDigi>::const_iterator iDigi = dsv_SiStripRawDigi->begin();
	 iDigi != dsv_SiStripRawDigi->end(); ++iDigi)
	 edm::LogInfo("SiStripMonitorTrack") << "\n\t Digi in detid " << iDigi->detId() << ": " << iDigi->size();
      */
    }else{
      LogDebug("SiStripMonitorTrack") << "RawDigis not found " << std::endl;
    }
  } else {
    LogDebug("SiStripMonitorTrack") << "RawDigis Off " << std::endl;
  }
  
  // trajectory input
  e.getByLabel(TrackProducer, TrackLabel, TrajectoryCollection);
  e.getByLabel(TrackProducer, TrackLabel, TItkAssociatorCollection);
  if( TItkAssociatorCollection.isValid()){
  }else{
    edm::LogError("SiStripMonitorTrack")<<"Association not found "<<std::endl;
    trackAssociatorCollection_in_EventTree=false;
  }
  
  //Perform track study
  if (tracksCollection_in_EventTree && trackAssociatorCollection_in_EventTree) trackStudy(es);
  
  //Perform Cluster Study (irrespectively to tracks)

  if(OffHisto_On_){
    if (dsv_SiStripCluster.isValid()){
      AllClusters(es);//analyzes the off Track Clusters
    }else{
      edm::LogError("SiStripMonitorTrack")<< "ClusterCollection is not valid!!" << std::endl;
    }
  }

  //Summary Counts of clusters
  std::map<TString, MonitorElement*>::iterator iME;
  std::map<TString, ModMEs>::iterator          iModME ;
  std::map<TString, LayerMEs>::iterator        iLayerME;
  
  for (int j=0;j<off_Flag;++j){ // loop over ontrack, offtrack
    int nTot=0;
    for (int i=0;i<4;++i){ // loop over TIB, TID, TOB, TEC
      name=flags[j]+"_in_"+SubDet[i];      
      iLayerME = LayerMEsMap.find(name);
      if(iLayerME!=LayerMEsMap.end()) {
	if(flags[j]=="OnTrack" && NClus[i][j]){
	  fillME(   iLayerME->second.nClusters,      NClus[i][j]);
	}else if(flags[j]=="OffTrack"){
	  fillME(   iLayerME->second.nClusters,      NClus[i][j]);
	}
	if(Trend_On_)
	  fillTrend(iLayerME->second.nClustersTrend, NClus[i][j]);
      }
      nTot+=NClus[i][j];
      NClus[i][j]=0;
    } // loop over TIB, TID, TOB, TEC
    
    name=flags[j]+"_TotalNumberOfClusters";
    iME = MEMap.find(name);
    if(iME!=MEMap.end() && nTot) iME->second->Fill(nTot);
    if(Trend_On_){
      iME = MEMap.find(name+"Trend");
      if(iME!=MEMap.end()) fillTrend(iME->second,nTot);
    }
  } // loop over ontrack, offtrack
  
}

//------------------------------------------------------------------------  
void SiStripMonitorTrack::book() 
{
  
  std::vector<uint32_t> vdetId_;
  SiStripDetCabling_->addActiveDetectorsRawIds(vdetId_);
  //Histos for each detector, layer and module
  for (std::vector<uint32_t>::const_iterator detid_iter=vdetId_.begin();detid_iter!=vdetId_.end();detid_iter++){  //loop on all the active detid
    uint32_t detid = *detid_iter;
    
    if (detid < 1){
      edm::LogError("SiStripMonitorTrack")<< "[" <<__PRETTY_FUNCTION__ << "] invalid detid " << detid<< std::endl;
      continue;
    }

    // set the DQM directory
    std::string MEFolderName = conf_.getParameter<std::string>("FolderName");    
    dbe->setCurrentFolder(MEFolderName);
    
    for (int j=0;j<off_Flag;j++) { // Loop on onTrack, offTrack
      name=flags[j]+"_TotalNumberOfClusters";
      if(MEMap.find(name)==MEMap.end()) {
	if(flags[j] == "OnTrack"){
	  MEMap[name]=bookME1D("TH1nClustersOn", name.Data());
	}else{
	  MEMap[name]=bookME1D("TH1nClustersOff", name.Data());
	}
	if(Trend_On_){
	  name+="Trend";
	  if(flags[j] == "OnTrack"){
	    MEMap[name]=bookMETrend("TH1nClustersOn", name.Data());
	  }else{
	    MEMap[name]=bookMETrend("TH1nClustersOff", name.Data());
	  }
	}      
      }
    }// End Loop on onTrack, offTrack
    
    // 	LogTrace("SiStripMonitorTrack") << " Detid " << *detid 
    //					<< " SubDet " << GetSubDetAndLayer(*detid).first 
    //					<< " Layer "  << GetSubDetAndLayer(*detid).second;

    // book Layer plots      
    if (DetectedLayers.find(folder_organizer.GetSubDetAndLayer(detid)) == DetectedLayers.end()){
      
      DetectedLayers[folder_organizer.GetSubDetAndLayer(detid)]=true;
      for (int j=0;j<off_Flag;j++){     
	folder_organizer.setLayerFolder(*detid_iter,folder_organizer.GetSubDetAndLayer(*detid_iter).second); 
	bookTrendMEs("layer",folder_organizer.GetSubDetAndLayer(*detid_iter).second,*detid_iter,flags[j]);
      }
    }
    
    if(Mod_On_){
      //    book module plots
      folder_organizer.setDetectorFolder(*detid_iter);
      bookModMEs("det",*detid_iter);
    }
    DetectedLayers[folder_organizer.GetSubDetAndLayer(*detid_iter)] |= (DetectedLayers.find(folder_organizer.GetSubDetAndLayer(*detid_iter)) == DetectedLayers.end());
    //      }
  }//end loop on detectors detid
  
  //  book SubDet plots
  for (std::map<std::pair<std::string,int32_t>,bool>::const_iterator iter=DetectedLayers.begin(); iter!=DetectedLayers.end();iter++){
    for (int j=0;j<off_Flag;j++){ // Loop on onTrack, offTrack
      folder_organizer.setDetectorFolder(0);
      dbe->cd("SiStrip/MechanicalView/"+iter->first.first);
      name=flags[j]+"_in_"+iter->first.first; 
      bookSubDetMEs(name,flags[j]); //for subdets
    }//end loop on onTrack,offTrack
  }  

}
  
//--------------------------------------------------------------------------------
void SiStripMonitorTrack::bookModMEs(TString name, uint32_t id)//Histograms at MODULE level
{
  SiStripHistoId hidmanager;
  std::string hid = hidmanager.createHistoId("",name.Data(),id);
  std::map<TString, ModMEs>::iterator iModME  = ModMEsMap.find(TString(hid));
  if(iModME==ModMEsMap.end()){
    ModMEs theModMEs; 
    // Cluster Width
    theModMEs.ClusterWidth=bookME1D("TH1ClusterWidth", hidmanager.createHistoId("ClusterWidth_OnTrack",name.Data(),id).c_str()); 
    dbe->tag(theModMEs.ClusterWidth,id); 
    // Cluster Charge
    theModMEs.ClusterCharge=bookME1D("TH1ClusterCharge", hidmanager.createHistoId("ClusterCharge_OnTrack",name.Data(),id).c_str());
    dbe->tag(theModMEs.ClusterCharge,id); 
    // Cluster StoN
    theModMEs.ClusterStoN=bookME1D("TH1ClusterStoN", hidmanager.createHistoId("ClusterStoN_OnTrack",name.Data(),id).c_str());
    dbe->tag(theModMEs.ClusterStoN,id); 
    // Cluster Charge Corrected
    theModMEs.ClusterChargeCorr=bookME1D("TH1ClusterChargeCorr", hidmanager.createHistoId("ClusterChargeCorr_OnTrack",name.Data(),id).c_str());
    dbe->tag(theModMEs.ClusterChargeCorr,id); 
    // Cluster StoN Corrected
    theModMEs.ClusterStoNCorr=bookME1D("TH1ClusterStoNCorr", hidmanager.createHistoId("ClusterStoNCorr_OnTrack",name.Data(),id).c_str());
    dbe->tag(theModMEs.ClusterStoNCorr,id); 
    // Cluster Position
    short total_nr_strips = SiStripDetCabling_->nApvPairs(id) * 2 * 128;
    theModMEs.ClusterPos=dbe->book1D(hidmanager.createHistoId("ClusterPosition_OnTrack",name.Data(),id).c_str(),hidmanager.createHistoId("ClusterPosition_OnTrack",name.Data(),id).c_str(),total_nr_strips,0.5,total_nr_strips+0.5);
    dbe->tag(theModMEs.ClusterPos,id); 
    // Cluster PGV
    theModMEs.ClusterPGV=bookMEProfile("TProfileClusterPGV", hidmanager.createHistoId("PGV_OnTrack",name.Data(),id).c_str()); 
    dbe->tag(theModMEs.ClusterPGV,id); 
    // Symmetric Eta function Eta=(L+R)/(2C) (Capacitive Coupling)
    //    theModMEs.ClusterSymmEtaCC=bookME1D("TH1ClusterSymmEtaCC", hidmanager.createHistoId("ClusterSymmEtaCC",name.Data(),id).c_str()); 
    //    dbe->tag(theModMEs.ClusterSymmEtaCC,id); 
    //bookeeping
    ModMEsMap[hid]=theModMEs;
  }
}

void SiStripMonitorTrack::bookTrendMEs(TString name,int32_t layer,uint32_t id,std::string flag)//Histograms and Trends at LAYER LEVEL
{
  char rest[1024];
  int subdetid = ((id>>25)&0x7);
  if(       subdetid==3 ){
  // ---------------------------  TIB  --------------------------- //
    TIBDetId tib1 = TIBDetId(id);
    sprintf(rest,"TIB__layer__%d",tib1.layer());
  }else if( subdetid==4){
  // ---------------------------  TID  --------------------------- //
    TIDDetId tid1 = TIDDetId(id);
    sprintf(rest,"TID__side__%d__wheel__%d",tid1.side(),tid1.wheel());
  }else if( subdetid==5){
  // ---------------------------  TOB  --------------------------- //
    TOBDetId tob1 = TOBDetId(id);
    sprintf(rest,"TOB__layer__%d",tob1.layer());
  }else if( subdetid==6){
  // ---------------------------  TEC  --------------------------- //
    TECDetId tec1 = TECDetId(id);
    sprintf(rest,"TEC__side__%d__wheel__%d",tec1.side(),tec1.wheel());
  }else{
  // ---------------------------  ???  --------------------------- //
    edm::LogError("SiStripTkDQM|WrongInput")<<"no such subdetector type :"<<subdetid<<" no folder set!"<<std::endl;
    return;
  }

  SiStripHistoId hidmanager;
  std::string hid = hidmanager.createHistoLayer("",name.Data(),rest,flag);
  std::map<TString, LayerMEs>::iterator iLayerME  = LayerMEsMap.find(TString(hid));
  if(iLayerME==LayerMEsMap.end()){
    LayerMEs theLayerMEs; 
 
    LogDebug("SiStripMonitorTrack") << "Booking " << rest << flag << std::endl;   
    // Cluster Width
    theLayerMEs.ClusterWidth=bookME1D("TH1ClusterWidth", hidmanager.createHistoLayer("Summary_ClusterWidth",name.Data(),rest,flag).c_str()); 
    dbe->tag(theLayerMEs.ClusterWidth,layer); 
    
    // Cluster Noise
    theLayerMEs.ClusterNoise=bookME1D("TH1ClusterNoise", hidmanager.createHistoLayer("Summary_ClusterNoise",name.Data(),rest,flag).c_str()); 
    dbe->tag(theLayerMEs.ClusterNoise,layer); 
    
    // Cluster Charge
    theLayerMEs.ClusterCharge=bookME1D("TH1ClusterCharge", hidmanager.createHistoLayer("Summary_ClusterCharge",name.Data(),rest,flag).c_str());
    dbe->tag(theLayerMEs.ClusterCharge,layer);
    
    // Cluster StoN
    theLayerMEs.ClusterStoN=bookME1D("TH1ClusterStoN", hidmanager.createHistoLayer("Summary_ClusterStoN",name.Data(),rest,flag).c_str());
    dbe->tag(theLayerMEs.ClusterStoN,layer); 

    // Trends
    if(Trend_On_){
      // Cluster Width
      theLayerMEs.ClusterWidthTrend=bookMETrend("TH1ClusterWidth", hidmanager.createHistoLayer("Trend_ClusterWidth",name.Data(),rest,flag).c_str()); 
      dbe->tag(theLayerMEs.ClusterWidthTrend,layer); 
      // Cluster Noise
      theLayerMEs.ClusterNoiseTrend=bookMETrend("TH1ClusterNoise", hidmanager.createHistoLayer("Trend_ClusterNoise",name.Data(),rest,flag).c_str()); 
      dbe->tag(theLayerMEs.ClusterNoiseTrend,layer); 
      // Cluster Charge
      theLayerMEs.ClusterChargeTrend=bookMETrend("TH1ClusterCharge", hidmanager.createHistoLayer("Trend_ClusterCharge",name.Data(),rest,flag).c_str());
      dbe->tag(theLayerMEs.ClusterChargeTrend,layer); 
      // Cluster StoN
      theLayerMEs.ClusterStoNTrend=bookMETrend("TH1ClusterStoN", hidmanager.createHistoLayer("Trend_ClusterStoN",name.Data(),rest,flag).c_str());
      dbe->tag(theLayerMEs.ClusterStoNTrend,layer); 
    }
    
    if(flag=="OnTrack"){
      // Cluster Charge Corrected
      theLayerMEs.ClusterChargeCorr=bookME1D("TH1ClusterChargeCorr", hidmanager.createHistoLayer("Summary_ClusterChargeCorr",name.Data(),rest,flag).c_str());
      dbe->tag(theLayerMEs.ClusterChargeCorr,layer); 
      // Cluster StoN Corrected
      theLayerMEs.ClusterStoNCorr=bookME1D("TH1ClusterStoNCorr", hidmanager.createHistoLayer("Summary_ClusterStoNCorr",name.Data(),rest,flag).c_str());
      dbe->tag(theLayerMEs.ClusterStoNCorr,layer); 
      // Symmetric Eta function Eta=(L+R)/(2C) (Capacitive Coupling)
      theLayerMEs.ClusterSymmEtaCC=bookME1D("TH1ClusterSymmEtaCC", hidmanager.createHistoLayer("Summary_ClusterSymmEtaCC",name.Data(),rest,flag).c_str()); 
      dbe->tag(theLayerMEs.ClusterSymmEtaCC,layer); 
      
      // Histograms booked and filled only if Charge Coupling Analysis is selected
      if(CCAnalysis_On_) {
	// Cluster Width (perpendicular tracks, Capacitive Coupling analysis)
	theLayerMEs.ClusterWidthCC=bookME1D("TH1ClusterWidthCC", hidmanager.createHistoLayer("Summary_ClusterWidthCC",name.Data(),rest,flag).c_str()); 
	dbe->tag(theLayerMEs.ClusterWidthCC,layer); 
	// Charge Coupling Estimator x=Eta/(1+2xEta) (Capacitive Coupling analysis)
	theLayerMEs.ClusterEstimatorCC=bookME1D("TH1ClusterEstimatorCC", hidmanager.createHistoLayer("Summary_ClusterEstimatorCC",name.Data(),rest,flag).c_str()); 
	dbe->tag(theLayerMEs.ClusterEstimatorCC,layer); 
      }
      
      if(Trend_On_){
	// Cluster Charge Corrected
	theLayerMEs.ClusterChargeCorrTrend=bookMETrend("TH1ClusterChargeCorr", hidmanager.createHistoLayer("Trend_ClusterChargeCorr",name.Data(),rest,flag).c_str());
	dbe->tag(theLayerMEs.ClusterChargeCorrTrend,layer); 
	
	// Cluster StoN Corrected
	theLayerMEs.ClusterStoNCorrTrend=bookMETrend("TH1ClusterStoNCorr", hidmanager.createHistoLayer("Trend_ClusterStoNCorr",name.Data(),rest,flag).c_str());
	dbe->tag(theLayerMEs.ClusterStoNCorrTrend,layer); 
	
	// Symmetric Eta function Eta=(L+R)/(2C) (Capacitive Coupling)
	theLayerMEs.ClusterSymmEtaCCTrend=bookMETrend("TH1ClusterSymmEtaCC", hidmanager.createHistoLayer("Trend_ClusterSymmEtaCC",name.Data(),rest,flag).c_str()); 
	dbe->tag(theLayerMEs.ClusterSymmEtaCCTrend,layer); 
      }
      
    }
    
    //Cluster Position
    short total_nr_strips = SiStripDetCabling_->nApvPairs(id) * 2 * 128; 
    theLayerMEs.ClusterPos= dbe->book1D(hidmanager.createHistoLayer("Summary_ClusterPosition",name.Data(),rest,flag).c_str(),hidmanager.createHistoLayer("Summary_ClusterPosition",name.Data(),rest,flag).c_str(),total_nr_strips, 0.5,total_nr_strips+0.5);
    dbe->tag(theLayerMEs.ClusterPos,layer); 
    
    //bookeeping
    LayerMEsMap[hid]=theLayerMEs;
  }
  
}

void SiStripMonitorTrack::bookSubDetMEs(TString name,TString flag)//Histograms at SubDet level
{
  std::map<TString, LayerMEs>::iterator iLayerME  = LayerMEsMap.find(name);
  char completeName[1024];
  if(iLayerME==LayerMEsMap.end()){
    LayerMEs theLayerMEs;
    
    // TotalNumber of Cluster 

    if (flag=="OnTrack"){
      sprintf(completeName,"Summary_TotalNumberOfClusters_%s",name.Data());
      theLayerMEs.nClusters=bookME1D("TH1nClustersOn", completeName);
      theLayerMEs.nClusters->getTH1()->StatOverflows(kTRUE);
    }else{
      sprintf(completeName,"Summary_TotalNumberOfClusters_%s",name.Data());
      theLayerMEs.nClusters=bookME1D("TH1nClustersOff", completeName);
      theLayerMEs.nClusters->getTH1()->StatOverflows(kTRUE);
    }
    
    // Cluster Width
    sprintf(completeName,"Summary_ClusterWidth_%s",name.Data());
    theLayerMEs.ClusterWidth=bookME1D("TH1ClusterWidth", completeName);
    
    // Cluster Noise
    sprintf(completeName,"Summary_ClusterNoise_%s",name.Data());
    theLayerMEs.ClusterNoise=bookME1D("TH1ClusterNoise", completeName);
    
    // Cluster Charge
    sprintf(completeName,"Summary_ClusterCharge_%s",name.Data());
    theLayerMEs.ClusterCharge=bookME1D("TH1ClusterCharge", completeName);
    
    // Cluster StoN
    sprintf(completeName,"Summary_ClusterStoN_%s",name.Data());
    theLayerMEs.ClusterStoN=bookME1D("TH1ClusterStoN", completeName);


    if(Trend_On_){
      if (flag=="OnTrack"){
	// TotalNumber of Cluster 
	sprintf(completeName,"Trend_TotalNumberOfClusters_%s",name.Data());
	theLayerMEs.nClustersTrend=bookMETrend("TH1nClustersOn", completeName);
      }else{
	sprintf(completeName,"Trend_TotalNumberOfClusters_%s",name.Data());
	theLayerMEs.nClustersTrend=bookMETrend("TH1nClustersOff", completeName);
      }
      // Cluster Width
      sprintf(completeName,"Trend_ClusterWidth_%s",name.Data());
      theLayerMEs.ClusterWidthTrend=bookMETrend("TH1ClusterWidth", completeName);
      // Cluster Noise
      sprintf(completeName,"Trend_ClusterNoise_%s",name.Data());
      theLayerMEs.ClusterNoiseTrend=bookMETrend("TH1ClusterNoise", completeName);
      // Cluster Charge
      sprintf(completeName,"Trend_ClusterCharge_%s",name.Data());
      theLayerMEs.ClusterChargeTrend=bookMETrend("TH1ClusterCharge", completeName);
      // Cluster StoN
      sprintf(completeName,"Trend_ClusterStoN_%s",name.Data());
      theLayerMEs.ClusterStoNTrend=bookMETrend("TH1ClusterStoN", completeName); 
    }

    if (flag=="OnTrack"){
      //Cluster StoNCorr
      sprintf(completeName,"Summary_ClusterStoNCorr_%s",name.Data());
      theLayerMEs.ClusterStoNCorr=bookME1D("TH1ClusterStoNCorr", completeName);
      
      // Cluster ChargeCorr
      sprintf(completeName,"Summary_ClusterChargeCorr_%s",name.Data());
      theLayerMEs.ClusterChargeCorr=bookME1D("TH1ClusterChargeCorr", completeName);

      // Symmetric Eta function Eta=(L+R)/(2C) (Capacitive Coupling)
      sprintf(completeName,"Summary_ClusterSymmEtaCC_%s",name.Data());
      theLayerMEs.ClusterSymmEtaCC=bookME1D("TH1ClusterSymmEtaCC", completeName);
         
      // Histograms booked and filled only if Charge Coupling Analysis is selected
      if(CCAnalysis_On_) {
	// Cluster Width (perpendicular tracks, Capacitive Coupling analysis)
	sprintf(completeName,"Summary_ClusterWidthCC_%s",name.Data());
	theLayerMEs.ClusterWidthCC=bookME1D("TH1ClusterWidthCC", completeName);
	// Charge Coupling Estimator x=Eta/(1+2xEta) (Capacitive Coupling analysis)
	sprintf(completeName,"Summary_ClusterEstimatorCC_%s",name.Data());
	theLayerMEs.ClusterEstimatorCC=bookME1D("TH1ClusterEstimatorCC", completeName);
      }
      
      if(Trend_On_){ 
	// Cluster StoNCorr
	sprintf(completeName,"Trend_ClusterStoNCorr_%s",name.Data());
	theLayerMEs.ClusterStoNCorrTrend=bookMETrend("TH1ClusterStoNCorr", completeName);     
	// Cluster ChargeCorr
	sprintf(completeName,"Trend_ClusterChargeCorr_%s",name.Data());
	theLayerMEs.ClusterChargeCorrTrend=bookMETrend("TH1ClusterChargeCorr", completeName);
	// Cluster Eta function
	sprintf(completeName,"Trend_ClusterSymmEtaCC_%s",name.Data());
	theLayerMEs.ClusterSymmEtaCCTrend=bookMETrend("TH1ClusterSymmEtaCC", completeName);


      }
    }
    
    //bookeeping
    LayerMEsMap[name]=theLayerMEs;
  }
}
//--------------------------------------------------------------------------------

MonitorElement* SiStripMonitorTrack::bookME1D(const char* ParameterSetLabel, const char* HistoName)
{
  Parameters =  conf_.getParameter<edm::ParameterSet>(ParameterSetLabel);
  return dbe->book1D(HistoName,HistoName,
		         Parameters.getParameter<int32_t>("Nbinx"),
		         Parameters.getParameter<double>("xmin"),
		         Parameters.getParameter<double>("xmax")
		    );
}

//--------------------------------------------------------------------------------
MonitorElement* SiStripMonitorTrack::bookME2D(const char* ParameterSetLabel, const char* HistoName)
{
  Parameters =  conf_.getParameter<edm::ParameterSet>(ParameterSetLabel);
  return dbe->book2D(HistoName,HistoName,
		     Parameters.getParameter<int32_t>("Nbinx"),
		     Parameters.getParameter<double>("xmin"),
		     Parameters.getParameter<double>("xmax"),
		     Parameters.getParameter<int32_t>("Nbiny"),
		     Parameters.getParameter<double>("ymin"),
		     Parameters.getParameter<double>("ymax")
		     );
}

//--------------------------------------------------------------------------------
MonitorElement* SiStripMonitorTrack::bookME3D(const char* ParameterSetLabel, const char* HistoName)
{
  Parameters =  conf_.getParameter<edm::ParameterSet>(ParameterSetLabel);
  return dbe->book3D(HistoName,HistoName,
		     Parameters.getParameter<int32_t>("Nbinx"),
		     Parameters.getParameter<double>("xmin"),
		     Parameters.getParameter<double>("xmax"),
		     Parameters.getParameter<int32_t>("Nbiny"),
		     Parameters.getParameter<double>("ymin"),
		     Parameters.getParameter<double>("ymax"),
		     Parameters.getParameter<int32_t>("Nbinz"),
		     Parameters.getParameter<double>("zmin"),
		     Parameters.getParameter<double>("zmax")
		     );
}

//--------------------------------------------------------------------------------
MonitorElement* SiStripMonitorTrack::bookMEProfile(const char* ParameterSetLabel, const char* HistoName)
{
    Parameters =  conf_.getParameter<edm::ParameterSet>(ParameterSetLabel);
    return dbe->bookProfile(HistoName,HistoName,
                            Parameters.getParameter<int32_t>("Nbinx"),
                            Parameters.getParameter<double>("xmin"),
                            Parameters.getParameter<double>("xmax"),
                            Parameters.getParameter<int32_t>("Nbiny"),
                            Parameters.getParameter<double>("ymin"),
                            Parameters.getParameter<double>("ymax"),
                            "" );
}

//--------------------------------------------------------------------------------
MonitorElement* SiStripMonitorTrack::bookMETrend(const char* ParameterSetLabel, const char* HistoName)
{
  Parameters =  conf_.getParameter<edm::ParameterSet>(ParameterSetLabel);
  edm::ParameterSet ParametersTrend =  conf_.getParameter<edm::ParameterSet>("Trending");
  MonitorElement* me = dbe->bookProfile(HistoName,HistoName,
					ParametersTrend.getParameter<int32_t>("Nbins"),
					0,
					ParametersTrend.getParameter<int32_t>("Nbins"),
					100, //that parameter should not be there !?
					Parameters.getParameter<double>("xmin"),
					Parameters.getParameter<double>("xmax"),
					"" );
  if(!me) return me;
  char buffer[256];
  sprintf(buffer,"EventId/%d",ParametersTrend.getParameter<int32_t>("Steps"));
  me->setAxisTitle(std::string(buffer),1);
  return me;
}

//------------------------------------------------------------------------------------------
void SiStripMonitorTrack::trackStudy(const edm::EventSetup& es)
{

  const reco::TrackCollection tC = *(trackCollection.product());
  int i=0;
  std::vector<TrajectoryMeasurement> measurements;
  for(TrajTrackAssociationCollection::const_iterator it =  TItkAssociatorCollection->begin();it !=  TItkAssociatorCollection->end(); ++it){
    const edm::Ref<std::vector<Trajectory> > traj_iterator = it->key;  
    // Trajectory Map, extract Trajectory for this track
    reco::TrackRef trackref = it->val;
    LogTrace("SiStripMonitorTrack")
      << "Track number "<< i+1 
      << "\n\tmomentum: " << trackref->momentum()
      << "\n\tPT: " << trackref->pt()
      << "\n\tvertex: " << trackref->vertex()
      << "\n\timpact parameter: " << trackref->d0()
      << "\n\tcharge: " << trackref->charge()
      << "\n\tnormalizedChi2: " << trackref->normalizedChi2() 
      <<"\n\tFrom EXTRA : "
      <<"\n\t\touter PT "<< trackref->outerPt()<<std::endl;
    i++;

    measurements =traj_iterator->measurements();
    std::vector<TrajectoryMeasurement>::iterator traj_mes_iterator;
    int nhit=0;
    for(traj_mes_iterator=measurements.begin();traj_mes_iterator!=measurements.end();traj_mes_iterator++){//loop on measurements
      //trajectory local direction and position on detector
      LocalPoint  stateposition;
      LocalVector statedirection;
      
      TrajectoryStateOnSurface  updatedtsos=traj_mes_iterator->updatedState();
      ConstRecHitPointer ttrh=traj_mes_iterator->recHit();
      if (!ttrh->isValid()) {continue;}
      
      std::stringstream ss;
      
      nhit++;
      
      const ProjectedSiStripRecHit2D* phit=dynamic_cast<const ProjectedSiStripRecHit2D*>( ttrh->hit() );
      const SiStripMatchedRecHit2D* matchedhit=dynamic_cast<const SiStripMatchedRecHit2D*>( ttrh->hit() );
      const SiStripRecHit2D* hit=dynamic_cast<const SiStripRecHit2D*>( ttrh->hit() );	
      
      RecHitType type=Single;

      if(matchedhit){
	LogTrace("SiStripMonitorTrack")<<"\nMatched recHit found"<< std::endl;
	type=Matched;
	
	GluedGeomDet * gdet=(GluedGeomDet *)tkgeom->idToDet(matchedhit->geographicalId());
	GlobalVector gtrkdirup=gdet->toGlobal(updatedtsos.localMomentum());	    
	//mono side
	const GeomDetUnit * monodet=gdet->monoDet();
	statedirection=monodet->toLocal(gtrkdirup);
	if(statedirection.mag() != 0)	  RecHitInfo(matchedhit->monoHit(),statedirection,trackref,es);
	//stereo side
	const GeomDetUnit * stereodet=gdet->stereoDet();
	statedirection=stereodet->toLocal(gtrkdirup);
	if(statedirection.mag() != 0)	  RecHitInfo(matchedhit->stereoHit(),statedirection,trackref,es);
	ss<<"\nLocalMomentum (stereo): " <<  statedirection;
      }
      else if(phit){
	LogTrace("SiStripMonitorTrack")<<"\nProjected recHit found"<< std::endl;
	type=Projected;
	GluedGeomDet * gdet=(GluedGeomDet *)tkgeom->idToDet(phit->geographicalId());
	
	GlobalVector gtrkdirup=gdet->toGlobal(updatedtsos.localMomentum());
	const SiStripRecHit2D&  originalhit=phit->originalHit();
	const GeomDetUnit * det;
	if(!StripSubdetector(originalhit.geographicalId().rawId()).stereo()){
	  //mono side
	  LogTrace("SiStripMonitorTrack")<<"\nProjected recHit found  MONO"<< std::endl;
	  det=gdet->monoDet();
	  statedirection=det->toLocal(gtrkdirup);
	  if(statedirection.mag() != 0) RecHitInfo(&(phit->originalHit()),statedirection,trackref,es);
	}
	else{
	  LogTrace("SiStripMonitorTrack")<<"\nProjected recHit found STEREO"<< std::endl;
	  //stereo side
	  det=gdet->stereoDet();
	  statedirection=det->toLocal(gtrkdirup);
	  if(statedirection.mag() != 0) RecHitInfo(&(phit->originalHit()),statedirection,trackref,es);
	}
      }else {
	if(hit!=0){
	  ss<<"\nSingle recHit found"<< std::endl;	  
	  statedirection=updatedtsos.localMomentum();
	  if(statedirection.mag() != 0) RecHitInfo(hit,statedirection,trackref,es);
	}
      }
      ss <<"LocalMomentum: "<<statedirection
	 << "\nLocal x-z plane angle: "<<atan2(statedirection.x(),statedirection.z());	      
      LogTrace("SiStripMonitorTrack") <<ss.str() << std::endl;
    }
    
  }
}

  void SiStripMonitorTrack::RecHitInfo(const SiStripRecHit2D* tkrecHit, LocalVector LV,reco::TrackRef track_ref, const edm::EventSetup& es){
    
    if(!tkrecHit->isValid()){
      LogTrace("SiStripMonitorTrack") <<"\t\t Invalid Hit " << std::endl;
      return;  
    }
    
    const uint32_t& detid = tkrecHit->geographicalId().rawId();
    if (find(ModulesToBeExcluded_.begin(),ModulesToBeExcluded_.end(),detid)!=ModulesToBeExcluded_.end()){
      LogTrace("SiStripMonitorTrack") << "Modules Excluded" << std::endl;
      return;
    }
    
    LogTrace("SiStripMonitorTrack")
      <<"\n\t\tRecHit on det "<<tkrecHit->geographicalId().rawId()
      <<"\n\t\tRecHit in LP "<<tkrecHit->localPosition()
      <<"\n\t\tRecHit in GP "<<tkgeom->idToDet(tkrecHit->geographicalId())->surface().toGlobal(tkrecHit->localPosition()) 
      <<"\n\t\tRecHit trackLocal vector "<<LV.x() << " " << LV.y() << " " << LV.z() <<std::endl; 
    
    //Get SiStripCluster from SiStripRecHit
    if ( tkrecHit != NULL ){
      LogTrace("SiStripMonitorTrack") << "GOOD hit" << std::endl;
      const SiStripCluster* SiStripCluster_ = &*(tkrecHit->cluster());
      SiStripClusterInfo* SiStripClusterInfo_ = new SiStripClusterInfo(detid,*SiStripCluster_,es);
            
      if ( clusterInfos(SiStripClusterInfo_,detid,"OnTrack", LV ) ) {
	vPSiStripCluster.push_back(SiStripCluster_);
	countOn++;
      }
      delete SiStripClusterInfo_; 
      //}
    }else{
     edm::LogError("SiStripMonitorTrack") << "NULL hit" << std::endl;
    }	  
  }

//------------------------------------------------------------------------

void SiStripMonitorTrack::AllClusters( const edm::EventSetup& es)
{

  //Loop on Dets
  for ( edmNew::DetSetVector<SiStripCluster>::const_iterator DSViter=dsv_SiStripCluster->begin(); DSViter!=dsv_SiStripCluster->end();DSViter++){
    uint32_t detid=DSViter->id();
    if (find(ModulesToBeExcluded_.begin(),ModulesToBeExcluded_.end(),detid)!=ModulesToBeExcluded_.end()) continue;
    //Loop on Clusters
    edm::LogInfo("SiStripMonitorTrack") << "on detid "<< detid << " N Cluster= " << DSViter->size();
    edmNew::DetSet<SiStripCluster>::const_iterator ClusIter = DSViter->begin();
    for(; ClusIter!=DSViter->end(); ClusIter++) {
      SiStripClusterInfo* SiStripClusterInfo_= new SiStripClusterInfo(detid,*ClusIter,es);
	LogDebug("SiStripMonitorTrack") << "ClusIter " << &*ClusIter << "\t " 
	                                << std::find(vPSiStripCluster.begin(),vPSiStripCluster.end(),&*ClusIter)-vPSiStripCluster.begin();
	if (std::find(vPSiStripCluster.begin(),vPSiStripCluster.end(),&*ClusIter) == vPSiStripCluster.end()){
	  if ( clusterInfos(SiStripClusterInfo_,detid,"OffTrack",LV) ) {
	    countOff++;
	  }
	}
	delete SiStripClusterInfo_; 
    }
  }
}

//------------------------------------------------------------------------
bool SiStripMonitorTrack::clusterInfos(SiStripClusterInfo* cluster, const uint32_t& detid,std::string flag, const LocalVector LV)
{
  LogTrace("SiStripMonitorTrack") << "\n["<<__PRETTY_FUNCTION__<<"]" << std::endl;
  //folder_organizer.setDetectorFolder(0);
  if (cluster==0) return false;
  // if one imposes a cut on the clusters, apply it
  const  edm::ParameterSet ps = conf_.getParameter<edm::ParameterSet>("ClusterConditions");
  if( ps.getParameter<bool>("On") &&
      (cluster->getSignalOverNoise() < ps.getParameter<double>("minStoN") ||
       cluster->getSignalOverNoise() > ps.getParameter<double>("maxStoN") ||
       cluster->getWidth() < ps.getParameter<double>("minWidth") ||
       cluster->getWidth() > ps.getParameter<double>("maxWidth")                    )) return false;
  // start of the analysis
  
  int SubDet_enum = StripSubdetector(detid).subdetId()-3;
  int iflag =0;
  if      (flag=="OnTrack")  iflag=0;
  else if (flag=="OffTrack") iflag=1;
  NClus[SubDet_enum][iflag]++;
  std::stringstream ss;
  //  const_cast<SiStripClusterInfo*>(cluster)->print(ss);
  LogTrace("SiStripMonitorTrack") << "\n["<<__PRETTY_FUNCTION__<<"]\n" << ss.str() << std::endl;
  
  float cosRZ = -2;
  LogTrace("SiStripMonitorTrack")<< "\n\tLV " << LV.x() << " " << LV.y() << " " << LV.z() << " " << LV.mag() << std::endl;
  if (LV.mag()!=0){
    cosRZ= fabs(LV.z())/LV.mag();
    LogTrace("SiStripMonitorTrack")<< "\n\t cosRZ " << cosRZ << std::endl;
  }
  std::string name;
  
  // Filling SubDet Plots (on Track + off Track)
  std::pair<std::string,int32_t> SubDetAndLayer = folder_organizer.GetSubDetAndLayer(detid);
  name=flag+"_in_"+SubDetAndLayer.first;
  fillTrendMEs(cluster,name,cosRZ,flag);
  fillCapacitiveCouplingMEs(cluster,name,cosRZ,flag);
  
  char rest[1024];
  int subdetid = ((detid>>25)&0x7);
  if(       subdetid==3 ){
    // ---------------------------  TIB  --------------------------- //
    TIBDetId tib1 = TIBDetId(detid);
     sprintf(rest,"TIB__layer__%d",tib1.layer());
   }else if( subdetid==4){
     // ---------------------------  TID  --------------------------- //
     TIDDetId tid1 = TIDDetId(detid);
     sprintf(rest,"TID__side__%d__wheel__%d",tid1.side(),tid1.wheel());
   }else if( subdetid==5){
     // ---------------------------  TOB  --------------------------- //
     TOBDetId tob1 = TOBDetId(detid);
     sprintf(rest,"TOB__layer__%d",tob1.layer());
   }else if( subdetid==6){
     // ---------------------------  TEC  --------------------------- //
     TECDetId tec1 = TECDetId(detid);
     sprintf(rest,"TEC__side__%d__wheel__%d",tec1.side(),tec1.wheel());
   }else{
     // ---------------------------  ???  --------------------------- //
     edm::LogError("SiStripTkDQM|WrongInput")<<"no such subdetector type :"<<subdetid<<" no folder set!"<<std::endl;
     return 0;
   }
   
   SiStripHistoId hidmanager1;
   
   // Filling Layer Plots
   name= hidmanager1.createHistoLayer("","layer",rest,flag);
   fillTrendMEs(cluster,name,cosRZ,flag);
   fillCapacitiveCouplingMEs(cluster,name,cosRZ,flag);
   
   // Module plots filled only for onTrack Clusters
   if(Mod_On_){
     if(flag=="OnTrack"){
       SiStripHistoId hidmanager2;
       name =hidmanager2.createHistoId("","det",detid);
       fillModMEs(cluster,name,cosRZ); 
     }
   }
     return true;
   }

//--------------------------------------------------------------------------------
void SiStripMonitorTrack::fillTrend(MonitorElement* me ,float value)
{
  if(!me) return;
  //check the origin and check options
  int option = conf_.getParameter<edm::ParameterSet>("Trending").getParameter<int32_t>("UpdateMode");
  if(firstEvent==-1) firstEvent = eventNb;
  int CurrentStep = atoi(me->getAxisTitle(1).c_str()+8);
  int firstEventUsed = firstEvent;
  int presentOverflow = (int)me->getBinEntries(me->getNbinsX()+1);
  if(option==2) firstEventUsed += CurrentStep * int(me->getBinEntries(me->getNbinsX()+1));
  else if(option==3) firstEventUsed += CurrentStep * int(me->getBinEntries(me->getNbinsX()+1)) * me->getNbinsX();
  //fill
  me->Fill((eventNb-firstEventUsed)/CurrentStep,value);
  if(eventNb-firstEvent<1) return;
  // check if we reached the end
  if(presentOverflow == me->getBinEntries(me->getNbinsX()+1)) return;
  switch(option) {
  case 1:
    {
      // mode 1: rebin and change X scale
      int NbinsX = me->getNbinsX();
      float entries = 0.;
      float content = 0.;
      float error = 0.;
      int bin = 1;
      int totEntries = int(me->getEntries());
      for(;bin<=NbinsX/2;++bin) {
	content = (me->getBinContent(2*bin-1) + me->getBinContent(2*bin))/2.; 
	error   = pow((me->getBinError(2*bin-1)*me->getBinError(2*bin-1)) + (me->getBinError(2*bin)*me->getBinError(2*bin)),0.5)/2.; 
	entries = me->getBinEntries(2*bin-1) + me->getBinEntries(2*bin);
	me->setBinContent(bin,content*entries);
	me->setBinError(bin,error);
	me->setBinEntries(bin,entries);
      }
      for(;bin<=NbinsX+1;++bin) {
	me->setBinContent(bin,0);
	me->setBinError(bin,0);
	me->setBinEntries(bin,0); 
      }
      me->setEntries(totEntries);
      char buffer[256];
      sprintf(buffer,"EventId/%d",CurrentStep*2);
      me->setAxisTitle(std::string(buffer),1);
      break;
    }
  case 2:
    {
      // mode 2: slide
      int bin=1;
      int NbinsX = me->getNbinsX();
      for(;bin<=NbinsX;++bin) {
	me->setBinContent(bin,me->getBinContent(bin+1)*me->getBinEntries(bin+1));
	me->setBinError(bin,me->getBinError(bin+1));
	me->setBinEntries(bin,me->getBinEntries(bin+1));
      }
      break;
    }
  case 3:
    {
      // mode 3: reset
      int NbinsX = me->getNbinsX();
      for(int bin=0;bin<=NbinsX;++bin) {
	me->setBinContent(bin,0);
	me->setBinError(bin,0);
	me->setBinEntries(bin,0); 
      }
      break;
    }
  }
}

//--------------------------------------------------------------------------------
void SiStripMonitorTrack::fillModMEs(SiStripClusterInfo* cluster,TString name,float cos)
{
  std::map<TString, ModMEs>::iterator iModME  = ModMEsMap.find(name);
  if(iModME!=ModMEsMap.end()){
    fillME(iModME->second.ClusterStoN  ,cluster->getSignalOverNoise());
    fillME(iModME->second.ClusterStoNCorr ,cluster->getSignalOverNoise()*cos);
    fillME(iModME->second.ClusterCharge,cluster->getCharge());
    fillME(iModME->second.ClusterChargeCorr,cluster->getCharge()*cos);
    fillME(iModME->second.ClusterWidth ,cluster->getWidth());
    fillME(iModME->second.ClusterPos   ,cluster->getPosition());
    
    //fill the PGV histo
    float PGVmax = cluster->getMaxCharge();
    int PGVposCounter = cluster->getFirstStrip() - cluster->getMaxPosition();
    for (int i= int(conf_.getParameter<edm::ParameterSet>("TProfileClusterPGV").getParameter<double>("xmin"));i<PGVposCounter;++i)
      fillME(iModME->second.ClusterPGV, i,0.);
    for (std::vector<uint8_t>::const_iterator it=cluster->getStripAmplitudes().begin();it<cluster->getStripAmplitudes().end();++it) {
      fillME(iModME->second.ClusterPGV, PGVposCounter++,(*it)/PGVmax);
    }
    for (int i= PGVposCounter;i<int(conf_.getParameter<edm::ParameterSet>("TProfileClusterPGV").getParameter<double>("xmax"));++i)
      fillME(iModME->second.ClusterPGV, i,0.);
    //end fill the PGV histo
  }
}

//------------------------------------------------------------------------
void SiStripMonitorTrack::fillTrendMEs(SiStripClusterInfo* cluster,std::string name,float cos, std::string flag)
{ 
  std::map<TString, LayerMEs>::iterator iLayerME  = LayerMEsMap.find(name);
  if(iLayerME!=LayerMEsMap.end()){
    if(flag=="OnTrack"){
      fillME(iLayerME->second.ClusterStoNCorr,(cluster->getSignalOverNoise())*cos);
      fillTrend(iLayerME->second.ClusterStoNCorrTrend,(cluster->getSignalOverNoise())*cos);
      fillME(iLayerME->second.ClusterChargeCorr,cluster->getCharge()*cos);
      fillTrend(iLayerME->second.ClusterChargeCorrTrend,cluster->getCharge()*cos);
    }
    fillME(iLayerME->second.ClusterStoN  ,cluster->getSignalOverNoise());
    fillTrend(iLayerME->second.ClusterStoNTrend,cluster->getSignalOverNoise());
    fillME(iLayerME->second.ClusterCharge,cluster->getCharge());
    fillTrend(iLayerME->second.ClusterChargeTrend,cluster->getCharge());
    fillME(iLayerME->second.ClusterNoise ,cluster->getNoise());
    fillTrend(iLayerME->second.ClusterNoiseTrend,cluster->getNoise());
    fillME(iLayerME->second.ClusterWidth ,cluster->getWidth());
    fillTrend(iLayerME->second.ClusterWidthTrend,cluster->getWidth());
    fillME(iLayerME->second.ClusterPos   ,cluster->getPosition());
  }
}

//------------------------------------------------------------------------
void SiStripMonitorTrack::fillCapacitiveCouplingMEs(SiStripClusterInfo* cluster,std::string name, float cos, std::string flag) {
  std::map<TString, LayerMEs>::iterator iLayerME  = LayerMEsMap.find(name);
  if(iLayerME!=LayerMEsMap.end()){
    // Capacitive Coupling analysis
    if( cos > 0.9 ) { // perpendicular track
      LogTrace("SiStripMonitorTrack") << "\t\t Perpendicular Track, cluster center " << cluster->getMaxPosition() << std::endl;
      
      // calculate symmetric eta function with only the Central + First Left / First Right strips
      float chargeCentral = 0.;
      std::pair< float,float > chargeLeftRight;
      //
      if(RawDigis_On_) {
	// with RawDigi
	if(!dsv_SiStripRawDigi.failedToGet()) {
	  if((*(dsv_SiStripRawDigi.product())).find(cluster->getDetId()) != (*(dsv_SiStripRawDigi.product())).end() ) { // ...if they exist
	    const edm::DetSet<SiStripRawDigi>    ds_SiStripRawDigi = (*(dsv_SiStripRawDigi.product()))[cluster->getDetId()];
	    edmNew::DetSet<SiStripCluster> ds_SiStripCluster = (*(dsv_SiStripCluster.product()))[cluster->getDetId()];
	    LogDebug("SiStripMonitorTrack") << "RawDigis found " << ds_SiStripRawDigi.size() << std::endl;
	    std::vector<float> chargesCLR = cluster->getRawChargeCLR(ds_SiStripRawDigi,ds_SiStripCluster,std::string("VirginRaw"));
	    chargeCentral = chargesCLR[0];
	    chargeLeftRight.first = chargesCLR[1];
	    chargeLeftRight.second = chargesCLR[2];
	  } // no SiStripRawDigi at all
	}
      } else {
	LogDebug("SiStripMonitorTrack") << "RawDigis " << RawDigis_On_ << " take info from cluster strips" << std::endl;
	// with cluster strips only
	chargeCentral = cluster->getMaxCharge();
	chargeLeftRight = cluster->getChargeLRFirstNeighbour();
      }
      float symmetricEta = SymEta(chargeCentral,chargeLeftRight.first,chargeLeftRight.second);
      float ccEstimator  = symmetricEta / ( 1 + 2 * symmetricEta );
      
      // Summary
      LogTrace("SiStripMonitorTrack")    
	<<"\n\t\t Cluster with multiplicity " << cluster->getWidth()
	<< " in det " << cluster->getDetId()
	<< ": Eta Function from Charge"
	<<"\n\t\t\t Signal = " << cluster->getCharge()
	<<"\n\t\t\t Noise  = " << cluster->getNoise()
	<<"\n\t\t\t S/N    = " << cluster->getSignalOverNoise()
	<<"\n\t\t\t  Cluster Central C = " << cluster->getMaxCharge()
	<<"\n\t\t\t  Cluster    Left L = " << cluster->getChargeLRFirstNeighbour().first
	<<"\n\t\t\t  Cluster   Right R = " << cluster->getChargeLRFirstNeighbour().second
	<<"\n\t\t\t  RawDigi Central C = " << chargeCentral
	<<"\n\t\t\t  RawDigi    Left L = " << chargeLeftRight.first
	<<"\n\t\t\t  RawDigi   Right R = " << chargeLeftRight.second
	<<"\n\t\t\t     Max = " << cluster->getMaxCharge()
	<<"\n\t\t\t     Pos = " << cluster->getMaxPosition()
	<<"\n\t\t\t     1st = " << cluster->getFirstStrip()
	<<"\n\t\t\t Symmetric Eta = " << symmetricEta
	<<"\n\t\t\t CC Estimator  = " << ccEstimator
	<< std::endl;
      for (std::vector<uint8_t>::const_iterator it=cluster->getStripAmplitudes().begin();it<cluster->getStripAmplitudes().end();++it) {
	LogTrace("SiStripMonitorTrack")    
	  <<"\t\t\t       " << (short)(*it) << std::endl;
      }
      //
      
      // fill monitor elements
      fillME(iLayerME->second.ClusterSymmEtaCC , symmetricEta);
      fillTrend(iLayerME->second.ClusterSymmEtaCCTrend, symmetricEta);
      // Histograms booked and filled only if Charge Coupling Analysis is selected
      if(CCAnalysis_On_) {
	fillME(iLayerME->second.ClusterWidthCC, cluster->getWidth() );
	fillME(iLayerME->second.ClusterEstimatorCC , ccEstimator );
      }
      //
      
    } // perpendicular track
  }
}

//------------------------------------------------------------------------
// Symmetric Eta function
float SiStripMonitorTrack::SymEta( float clusterCentralCharge, float clusterLeftCharge, float clusterRightCharge){
  float symeta = ( clusterLeftCharge + clusterRightCharge )/( 2 * clusterCentralCharge );
  return symeta;
}
