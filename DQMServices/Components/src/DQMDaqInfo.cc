#include "DQMServices/Components/src/DQMDaqInfo.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"

DQMDaqInfo::DQMDaqInfo(const edm::ParameterSet& iConfig)  
{
  
   
}


DQMDaqInfo::~DQMDaqInfo()
{  
}


void DQMDaqInfo::beginLuminosityBlock(const edm::LuminosityBlock& lumiBlock, const  edm::EventSetup& iSetup){
  
  
    
  edm::eventsetup::EventSetupRecordKey recordKey(edm::eventsetup::EventSetupRecordKey::TypeTag::findType("RunInfoRcd"));
  
  if( 0 != iSetup.find( recordKey ) ) {
    
    edm::ESHandle<RunInfo> sumFED;
    iSetup.get<RunInfoRcd>().get(sumFED);    
  
  
    //const RunInfo* summaryFED=sumFED.product();
  
    std::vector<int> FedsInIds= sumFED->m_fed_in;   
    
      
    std::pair<int,int> PixelRange   = FEDNumbering::getSiPixelFEDIds();;
    std::pair<int,int> TrackerRange = FEDNumbering::getSiStripFEDIds();
    std::pair<int,int> CSCRange     = FEDNumbering::getCSCFEDIds();
    std::pair<int,int> RPCRange     = FEDNumbering::getRPCFEDIds();
    std::pair<int,int> DTRange      = FEDNumbering::getDTFEDIds();
    std::pair<int,int> HcalRange    = FEDNumbering::getHcalFEDIds();

    std::pair<int,int> ECALBarrRange(610,645);
    std::pair<int,int> ECALEndcapRangeLow(601,609);
    std::pair<int,int> ECALEndcapRangeHigh(646,654);
  
    std::pair<int,int> L1TRange    = FEDNumbering::getTriggerGTPFEDIds();
  
    float  FedCount[9]={0., 0., 0., 0., 0., 0., 0., 0., 0.};
    float  NumberOfFeds[9];
    NumberOfFeds[Pixel]   = PixelRange.second-PixelRange.first +1;
    NumberOfFeds[SiStrip] = TrackerRange.second-TrackerRange.first +1;
    NumberOfFeds[CSC]     = CSCRange.second-CSCRange.first  +1;
    NumberOfFeds[RPC]     = RPCRange.second-RPCRange.first  +1;
    NumberOfFeds[DT]      = DTRange.second-DTRange.first +1;
    NumberOfFeds[Hcal]    = HcalRange.second-HcalRange.first +1;

    NumberOfFeds[EcalBarrel]    = ECALBarrRange.second-ECALBarrRange.first ;
    NumberOfFeds[EcalEndcap]    = (ECALEndcapRangeLow.second-ECALEndcapRangeLow.first)+(ECALEndcapRangeHigh.second-ECALEndcapRangeHigh.first) ;

    NumberOfFeds[L1T]    = L1TRange.second-L1TRange.first +1;
    
    for(unsigned int fedItr=0;fedItr<FedsInIds.size(); ++fedItr) {
      int fedID=FedsInIds[fedItr];
      

      if(fedID>=PixelRange.first   &  fedID<=PixelRange.second)    ++FedCount[Pixel]  ;
      if(fedID>=TrackerRange.first &  fedID<=TrackerRange.second)  ++FedCount[SiStrip];
      if(fedID>=CSCRange.first     &  fedID<=CSCRange.second)      ++FedCount[CSC]    ;
      if(fedID>=RPCRange.first     &  fedID<=RPCRange.second)      ++FedCount[RPC]    ;
      if(fedID>=DTRange.first      &  fedID<=DTRange.second)       ++FedCount[DT]     ;
      if(fedID>=HcalRange.first    &  fedID<=HcalRange.second)     ++FedCount[Hcal]	;       
      if(fedID>=ECALBarrRange.first    &  fedID<=ECALBarrRange.second)     ++FedCount[EcalBarrel]   ;      
      if((fedID>=ECALEndcapRangeLow.first & fedID<=ECALEndcapRangeLow.second)
	 ||(fedID>=ECALEndcapRangeHigh.first & fedID<=ECALEndcapRangeHigh.second)) ++FedCount[EcalEndcap]   ;
      if(fedID>=L1TRange.first    &  fedID<=L1TRange.second)     ++FedCount[L1T]   ;
    
    }   
    
    for(int detIndex=0; detIndex<9; ++detIndex) { 
      DaqFraction[detIndex]->Fill( FedCount[detIndex]/NumberOfFeds[detIndex]);
    }
    


  }else{    
  
    for(int detIndex=0; detIndex<9; ++detIndex)  DaqFraction[detIndex]->Fill(-1);               
    return; 
  }
  
 
}


void DQMDaqInfo::endLuminosityBlock(const edm::LuminosityBlock&  lumiBlock, const  edm::EventSetup& iSetup){
}


void 
DQMDaqInfo::beginJob(const edm::EventSetup& iSetup)
{
  dbe_ = 0;
  dbe_ = edm::Service<DQMStore>().operator->();
  
  std::string commonFolder = "/EventInfo/DAQContents";  
  std::string subsystFolder;
  std::string curentFolder;
  
  subsystFolder="Pixel";  
  curentFolder= subsystFolder+commonFolder;
  dbe_->setCurrentFolder(curentFolder.c_str());
  DaqFraction[Pixel]   = dbe_->bookFloat("PixelDaqFraction");
  

  subsystFolder="SiStrip";  
  curentFolder=subsystFolder+commonFolder;
  dbe_->setCurrentFolder(curentFolder.c_str());
  DaqFraction[SiStrip]    = dbe_->bookFloat("SiStripDaqFraction");
  
  subsystFolder="RPC";  
  curentFolder=subsystFolder+commonFolder;
  dbe_->setCurrentFolder(curentFolder.c_str());
  DaqFraction[RPC]        = dbe_->bookFloat("RPCDaqFraction");
  
  subsystFolder="CSC";  
  curentFolder=subsystFolder+commonFolder;
  dbe_->setCurrentFolder(curentFolder.c_str());
  DaqFraction[CSC]       = dbe_->bookFloat("CSCDaqFraction");

  subsystFolder="DT";  
  curentFolder=subsystFolder+commonFolder;
  dbe_->setCurrentFolder(curentFolder.c_str());
  DaqFraction[DT]         = dbe_->bookFloat("DTDaqFraction");

  subsystFolder="Hcal";  
  curentFolder=subsystFolder+commonFolder;
  dbe_->setCurrentFolder(curentFolder.c_str());
  DaqFraction[Hcal]       = dbe_->bookFloat("HcalDaqFraction");

  subsystFolder="EcalBarrel";  
  curentFolder=subsystFolder+commonFolder;
  dbe_->setCurrentFolder(curentFolder.c_str());
  DaqFraction[EcalBarrel]       = dbe_->bookFloat("EcalBarrDaqFraction");

  subsystFolder="EcalEndcap";  
  curentFolder=subsystFolder+commonFolder;
  dbe_->setCurrentFolder(curentFolder.c_str());
  DaqFraction[EcalEndcap]       = dbe_->bookFloat("EcalEndDaqFraction");

  subsystFolder="L1T";  
  curentFolder=subsystFolder+commonFolder;
  dbe_->setCurrentFolder(curentFolder.c_str());
  DaqFraction[L1T]       = dbe_->bookFloat("L1TDaqFraction");


}


void 
DQMDaqInfo::endJob() {
}



void
DQMDaqInfo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
 

}


//define this as a plug-in
//DEFINE_FWK_MODULE(DQMDaqInfo);
