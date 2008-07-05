#include "DQM/HcalMonitorTasks/interface/HcalTrigPrimMonitor.h"
#include "DQMServices/Core/interface/DQMStore.h"

HcalTrigPrimMonitor::HcalTrigPrimMonitor() {
  ievt_=0;
  occThresh_=0;
}

HcalTrigPrimMonitor::~HcalTrigPrimMonitor() {
}

void HcalTrigPrimMonitor::reset(){}

void HcalTrigPrimMonitor::clearME(){

  if(m_dbe){
    m_dbe->setCurrentFolder(baseFolder_);
    m_dbe->removeContents();
    meEVT_= 0;
  }

}

void HcalTrigPrimMonitor::setup(const edm::ParameterSet& ps, DQMStore* dbe){
  HcalBaseMonitor::setup(ps,dbe);
  baseFolder_ = rootFolder_+"TrigPrimMonitor";

  etaMax_ = ps.getUntrackedParameter<double>("MaxEta", 41.5);
  etaMin_ = ps.getUntrackedParameter<double>("MinEta", -41.5);
  etaBins_ = (int)(etaMax_ - etaMin_);
  cout << "TrigPrim eta min/max set to " << etaMin_ << "/" << etaMax_ << endl;
  
  phiMax_ = ps.getUntrackedParameter<double>("MaxPhi", 73);
  phiMin_ = ps.getUntrackedParameter<double>("MinPhi", 0);
  phiBins_ = (int)(phiMax_ - phiMin_);
  cout << "TrigPrim phi min/max set to " << phiMin_ << "/" << phiMax_ << endl;
  
  occThresh_ = ps.getUntrackedParameter<double>("TPOccThresh", 1.0);
  cout << "TrigPrim occupancy threshold set to " << occThresh_ << endl;

  TPThresh_ = ps.getUntrackedParameter<double>("TPThreshold", 1.0);
  cout << "TrigPrim threshold set to " << TPThresh_ << endl;

  TPdigi_ = ps.getUntrackedParameter<int>("TPdigiTS", 1);
  cout << "TP digi set to " << TPdigi_ << endl;

  ADCdigi_ = ps.getUntrackedParameter<int>("ADCdigiTS", 3);
  cout << "ADC dgisi set to " << ADCdigi_ << endl;


  ievt_=0;
  
  if ( m_dbe !=NULL ) {    

    char* type;
//    char name[128];
    m_dbe->setCurrentFolder(baseFolder_);

//ZZ Expert Plots
    m_dbe->setCurrentFolder(baseFolder_ + "/ZZ Expert Plots/ZZ DQM Expert Plots");
    type = "TrigPrim Event Number";
    meEVT_ = m_dbe->bookInt(type);

// 00 TP Occupancy
    m_dbe->setCurrentFolder(baseFolder_);
    type = "00 TP Occupancy";
    TPOcc_          = m_dbe->book2D(type,type,etaBins_,etaMin_,etaMax_,phiBins_,phiMin_,phiMax_);

//Timing Plots
    m_dbe->setCurrentFolder(baseFolder_+"/Timing Plots");
    type = "TP Size";
    tpSize_ = m_dbe->book1D(type,type,20,-0.5,19.5);
    type = "TP Timing";
    TPTiming_       = m_dbe->book1D(type,type,10,0,10);
    type = "TP Timing (Top wedges)";
    TPTimingTop_    = m_dbe->book1D(type,type,10,0,10);
    type = "TP Timing (Bottom wedges)";
    TPTimingBot_    = m_dbe->book1D(type,type,10,0,10);
    type = "TS with max ADC";
    TS_MAX_         = m_dbe->book1D(type,type,10,0,10);
//Energy Plots
    m_dbe->setCurrentFolder(baseFolder_+"/Energy Plots");
    type = "# TP Digis";
    tpCount_ = m_dbe->book1D(type,type,500,-0.5,4999.5);
    type = "# TP Digis over Threshold";
    tpCountThr_ = m_dbe->book1D(type,type,100,-0.5,999.5);
    type = "ADC spectrum positive TP";
    TP_ADC_         = m_dbe->book1D(type,type,200,-0.5,199.5);
    type = "Max ADC in TP";
    MAX_ADC_        = m_dbe->book1D(type,type,20,-0.5,19.5);
    type = "Full TP Spectrum";
    tpSpectrumAll_ = m_dbe->book1D(type,type,200,-0.5,199.5);
    type = "TP ET Sum";
    tpETSumAll_ = m_dbe->book1D(type,type,200,-0.5,199.5);
    type = "TP SOI ET";
    tpSOI_ET_ = m_dbe->book1D(type,type,100,-0.5,99.5);
    
    m_dbe->setCurrentFolder(baseFolder_+"/Energy Plots/TP Spectra by TS");
    for (int i=0; i<10; i++) {
      type = "TP Spectrum sample ";
      std::stringstream samp;
      std::string teststr;
      samp << i;
      samp >> teststr;
      teststr = type + teststr;
      tpSpectrum_[i]= m_dbe->book1D(teststr,teststr,100,-0.5,99.5);      
    }

//Electronics Plots
    m_dbe->setCurrentFolder(baseFolder_+"/Electronics Plots");
    type = "TP vs Digi";
    TPvsDigi_       = m_dbe->book2D(type,type,128,0,128,200,0,200);
    TPvsDigi_->setAxisTitle("lin ADC digi",1);
    TPvsDigi_->setAxisTitle("TP digi",2);
    type = "TrigPrim VME Occupancy Map";
    OCC_ELEC_VME = m_dbe->book2D(type,type,40,-0.25,19.75,18,-0.5,17.5);
    type = "TrigPrim Spigot Occupancy Map";
    OCC_ELEC_DCC = m_dbe->book2D(type,type,HcalDCCHeader::SPIGOT_COUNT,-0.5,HcalDCCHeader::SPIGOT_COUNT-0.5,36,-0.5,35.5);
    type = "TrigPrim VME Energy Map";
    EN_ELEC_VME = m_dbe->book2D(type,type,40,-0.25,19.75,18,-0.5,17.5);
    type = "TrigPrim Spigot Energy Map";
    EN_ELEC_DCC = m_dbe->book2D(type,type,HcalDCCHeader::SPIGOT_COUNT,-0.5,HcalDCCHeader::SPIGOT_COUNT-0.5,36,-0.5,35.5);

//Geometry Plots
    m_dbe->setCurrentFolder(baseFolder_+"/Geometry Plots");
    type = "TrigPrim Eta Occupancy Map";
    OCC_ETA = m_dbe->book1D(type,type,etaBins_,etaMin_,etaMax_);
    type = "TrigPrim Phi Occupancy Map";
    OCC_PHI = m_dbe->book1D(type,type,phiBins_,phiMin_,phiMax_);
    type = "TrigPrim Geo Occupancy Map";
    OCC_MAP_GEO = m_dbe->book2D(type,type,etaBins_,etaMin_,etaMax_,phiBins_,phiMin_,phiMax_);
    type = "TrigPrim Geo Threshold Map";
    OCC_MAP_THR = m_dbe->book2D(type,type,etaBins_,etaMin_,etaMax_,phiBins_,phiMin_,phiMax_);
    type = "TrigPrim Eta Energy Map";
    EN_ETA = m_dbe->book1D(type,type,etaBins_,etaMin_,etaMax_);
    type = "TrigPrim Phi Energy Map";
    EN_PHI = m_dbe->book1D(type,type,phiBins_,phiMin_,phiMax_);
    type = "TrigPrim Geo Energy Map";
    EN_MAP_GEO = m_dbe->book2D(type,type,etaBins_,etaMin_,etaMax_,phiBins_,phiMin_,phiMax_);
  
    meEVT_->Fill(ievt_);

  }

  return;
}

void HcalTrigPrimMonitor::processEvent(const HBHERecHitCollection& hbHits, 
				       const HORecHitCollection& hoHits, 
				       const HFRecHitCollection& hfHits,
				       const HBHEDigiCollection& hbhedigi,
				       const HODigiCollection& hodigi,
				       const HFDigiCollection& hfdigi,
                                       const HcalTrigPrimDigiCollection& tpDigis,
				       const HcalElectronicsMap& emap
				       ){
  

  if(!m_dbe) { 
    printf("HcalTrigPrimMonitor::processEvent   DQMStore not instantiated!!!\n");  
    return; 
  }

  ievt_++;
  meEVT_->Fill(ievt_);

  tpCount_->Fill(tpDigis.size()*1.0);  // number of TPGs collected per event
  
 float data[10];
  ClearEvent();
  try{
         edm::Handle<HcalTrigPrimDigiCollection> Tr_hbhe;
     

    int TPGsOverThreshold = 0;
    // for(HBHEDigiCollection::const_iterator j=hbhedigi.begin(); j!=hbhedigi.end(); j++){
    // HBHEDataFrame digi = (const HBHEDataFrame)(*j);

       for (HcalTrigPrimDigiCollection::const_iterator j=tpDigis.begin(); j!=tpDigis.end(); j++){
    const HcalTriggerPrimitiveDigi digi = (const HcalTriggerPrimitiveDigi)(*j);

// for(HcalTrigPrimDigiCollection::const_iterator j=Tr_hbhe->begin();j!=Tr_hbhe->end();j++){
  //        HcalTriggerPrimitiveDigi digi = (const HcalTriggerPrimitiveDigi)(*j);

      // find corresponding rechit and digis
      HcalTrigTowerDetId tpid=digi.id();
      HcalElectronicsId eid = emap.lookupTrigger(tpid);

      tpSOI_ET_->Fill(digi.SOI_compressedEt());      
      if(digi.SOI_compressedEt()>0 || true){	
	tpSize_->Fill(digi.size());	
	OCC_ETA->Fill(tpid.ieta());
	OCC_PHI->Fill(tpid.iphi());
	OCC_MAP_GEO->Fill(tpid.ieta(), tpid.iphi());

	EN_ETA->Fill(tpid.ieta(),digi.SOI_compressedEt());
	EN_PHI->Fill(tpid.iphi(),digi.SOI_compressedEt());
	EN_MAP_GEO->Fill(tpid.ieta(), tpid.iphi(),digi.SOI_compressedEt());
	
	float slotnum = eid.htrSlot() + 0.5*eid.htrTopBottom();	
	OCC_ELEC_VME->Fill(slotnum,eid.readoutVMECrateId());
	OCC_ELEC_DCC->Fill(eid.spigot(),eid.dccid());
	EN_ELEC_VME->Fill(slotnum,eid.readoutVMECrateId(),digi.SOI_compressedEt());
	EN_ELEC_DCC->Fill(eid.spigot(),eid.dccid(),digi.SOI_compressedEt());
	double etSum = 0;
	bool threshCond = false;
	//	printf("\nSampling\n");
	for (int j=0; j<digi.size(); j++) {
	  //	  printf("Sample %d\n",j);
	  float compressedEt = digi.sample(j).compressedEt();
	  //	  float compressedEt =1;
	  tpSpectrum_[j]->Fill(compressedEt);
	  tpSpectrumAll_->Fill(compressedEt);
	  etSum += compressedEt;
	  if (compressedEt>occThresh_) threshCond = true;
	}
	tpETSumAll_->Fill(etSum);
	
	if (threshCond){
	  OCC_MAP_THR->Fill(tpid.ieta(),tpid.iphi());  // which ieta and iphi positions the TPGs have for overThreshold cut
	  
	  TPGsOverThreshold++;
	}
      }
  ///////
   for (int i=0; i<digi.size(); i++) {
        data[i]=digi.sample(i).compressedEt();
	if(digi.sample(i).compressedEt()>TPThresh_){
	    tpSpectrum_[i]->Fill(digi.sample(i).compressedEt());
	    tpSpectrumAll_->Fill(digi.sample(i).compressedEt());
	    TPTiming_->Fill(i);
	    if(digi.id().iphi()>1  && digi.id().iphi()<36)  TPTimingTop_->Fill(i);
	    if(digi.id().iphi()>37 && digi.id().iphi()<72)  TPTimingBot_->Fill(i);
	    TPOcc_->Fill(digi.id().ieta(),digi.id().iphi());
	}
      }
    set_tp(digi.id().ieta(),digi.id().iphi(),1,data);
    /////////
    }
    tpCountThr_->Fill(TPGsOverThreshold*1.0);  // number of TPGs collected per event
  } catch (...) {    
    printf("HcalTrigPrimMonitor:  no tp digis\n");
  }
 try{
    for(HBHEDigiCollection::const_iterator j=hbhedigi.begin(); j!=hbhedigi.end(); j++){
       HBHEDataFrame digi = (const HBHEDataFrame)(*j);
       for(int i=0; i<digi.size(); i++) data[i]=digi.sample(i).adc();
       set_adc(digi.id().ieta(),digi.id().iphi(),digi.id().depth(),data);
    }
  } catch (...) {    
    printf("HcalTrigPrimMonitor:  no hcal digis\n");
  }
    // Correlation plots...
  int eta,phi;
  for(eta=-16;eta<=16;eta++) for(phi=1;phi<=72;phi++){
    for(int i=1;i<10;i++){
      int j1=(int)get_adc(eta,phi,1)[ADCdigi_];
      float tmp11 = (TrigMonAdc2fc[j1]+0.5);
      int j2=(int)get_adc(eta,phi,1)[ADCdigi_+1];
      float tmp21 = (TrigMonAdc2fc[j2]+0.5);
      int j3=(int)get_adc(eta,phi,2)[ADCdigi_];
      float tmp12 = (TrigMonAdc2fc[j3]+0.5);
      int j4=(int)get_adc(eta,phi,2)[ADCdigi_+1];
      float tmp22 = (TrigMonAdc2fc[j4]+0.5);
      if(IsSet_adc(eta,phi,1) && IsSet_tp(eta,phi,1)){
        if(get_tp(eta,phi)[TPdigi_]>TPThresh_){ 
          TPvsDigi_->Fill(tmp11+tmp21+tmp12+tmp22,get_tp(eta,phi)[TPdigi_]);

          float Energy=0;
          int TS = 0;
          for(int j=0;j<10;j++){
            if (get_adc(eta,phi,1)[j]>Energy){
              Energy=get_adc(eta,phi,1)[j];
              TS = j;
            }
          }
          MAX_ADC_->Fill(Energy);
          TS_MAX_->Fill(TS);
          //This may need to continue?
          Energy=0; for(int j=0;j<10;j++) Energy+=get_adc(eta,phi,1)[j]; 
	  TP_ADC_->Fill(Energy);	       	       
	}
      }
    }	   
  }

  return;
}
