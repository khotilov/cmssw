#include "DQM/HcalMonitorTasks/interface/HcalDigiMonitor.h"

HcalDigiMonitor::HcalDigiMonitor() {
  occThresh_ = 1;
  ievt_=0;
}

HcalDigiMonitor::~HcalDigiMonitor() {}

static bool bitUpset(int last, int now){
  if(last ==-1) return false;
  int v = last+1; if(v==4) v=0;
  if(v==now) return false;
  return true;
}

namespace HcalDigiMap{
  template<class Digi>
  inline void fillErrors(const Digi& digi, MonitorElement* mapG, MonitorElement* mapE){
    if(digiErr(digi)){
      mapG->Fill(digi.id().ieta(),digi.id().iphi());
      mapE->Fill(digi.elecId().readoutVMECrateId(),digi.elecId().htrSlot());
    }
    return;
  }

  template<class Digi>
  inline void fillOccupancy(const Digi& digi, MonitorElement* mapG, MonitorElement* mapE, float thr){
    if(digiOccupied(digi,thr)){
      mapG->Fill(digi.id().ieta(),digi.id().iphi());
      mapE->Fill(digi.elecId().readoutVMECrateId(),digi.elecId().htrSlot());
    }
    return;
  }

  template<class Digi>
  static bool digiErr(const Digi& digi){
    int last = -1;
    for (int i=0; i<digi.size(); i++) { 
      if(bitUpset(last,digi.sample(i).capid())) return true;
      if(digi.sample(i).er()) return true;
    }
    return false;
  }

  template<class Digi>
  static bool digiOccupied(const Digi& digi, float thr){
    for (int i=0; i<digi.size(); i++) { 
      if(digi.sample(i).adc()>thr) return true;
    }
    return false;
  }

}

void HcalDigiMonitor::setup(const edm::ParameterSet& ps, DaqMonitorBEInterface* dbe){
  HcalBaseMonitor::setup(ps,dbe);
  
  occThresh_ = ps.getUntrackedParameter<int>("DigiOccThresh", 10);
  cout << "Digi occupancy threshold set to " << occThresh_ << endl;
  ievt_=0;

  if ( m_dbe ) {

    m_dbe->setCurrentFolder("HcalMonitor/DigiMonitor");        
    meEVT_ = m_dbe->bookInt("Digi Task Event Number");    
    meEVT_->Fill(ievt_);

    m_dbe->setCurrentFolder("HcalMonitor/DigiMonitor/HBHE");
    hbHists.DIGI_NUM =  m_dbe->book1D("HBHE # of Digis","HBHE # of Digis",100,0,200);
    hbHists.DIGI_SIZE =  m_dbe->book1D("HBHE Digi Size","HBHE Digi Size",50,0,50);
    hbHists.DIGI_PRESAMPLE =  m_dbe->book1D("HBHE Digi Presamples","HBHE Digi Presamples",50,0,50);
    hbHists.QIE_CAPID =  m_dbe->book1D("HBHE QIE Cap-ID","HBHE QIE Cap-ID",6,-0.5,5.5);
    hbHists.QIE_ADC = m_dbe->book1D("HBHE QIE ADC Value","HBHE QIE ADC Value",100,0,200);
    hbHists.QIE_DV = m_dbe->book1D("HBHE QIE Data Value","HBHE QIE Data Value",2,-0.5,1.5);
    hbHists.ERR_MAP_GEO = m_dbe->book2D("HBHE Digi Geo Error Map","HBHE Digi Geo Error Map",59,-29.5,29.5,74,0,73);
    hbHists.ERR_MAP_ELEC = m_dbe->book2D("HBHE Digi Elec Error Map","HBHE Digi Elec Error Map",20,0,20,20,0,20);
    hbHists.OCC_MAP_GEO = m_dbe->book2D("HBHE Digi Geo Occupancy Map","HBHE Digi Geo Occupancy Map",59,-29.5,29.5,74,0,73);
    hbHists.OCC_MAP_ELEC = m_dbe->book2D("HBHE Digi Elec Occupancy Map","HBHE Digi Elec Occupancy Map",20,0,20,20,0,20);

    m_dbe->setCurrentFolder("HcalMonitor/DigiMonitor/HF");
    hfHists.DIGI_NUM =  m_dbe->book1D("HF # of Digis","HF # of Digis",100,0,200);
    hfHists.DIGI_SIZE =  m_dbe->book1D("HF Digi Size","HF Digi Size",50,0,50);
    hfHists.DIGI_PRESAMPLE =  m_dbe->book1D("HF Digi Presamples","HF Digi Presamples",50,0,50);
    hfHists.QIE_CAPID =  m_dbe->book1D("HF QIE Cap-ID","HF QIE Cap-ID",6,-0.5,5.5);
    hfHists.QIE_ADC = m_dbe->book1D("HF QIE ADC Value","HF QIE ADC Value",100,0,200);
    hfHists.QIE_DV = m_dbe->book1D("HF QIE Data Value","HF QIE Data Value",2,-0.5,1.5);
    hfHists.ERR_MAP_GEO = m_dbe->book2D("HF Digi Geo Error Map","HF Digi Geo Error Map",59,-29.5,29.5,74,0,73);
    hfHists.ERR_MAP_ELEC = m_dbe->book2D("HF Digi Elec Error Map","HF Digi Elec Error Map",20,0,20,20,0,20);
    hfHists.OCC_MAP_GEO = m_dbe->book2D("HF Digi Geo Occupancy Map","HF Digi Geo Occupancy Map",59,-29.5,29.5,74,0,73);
    hfHists.OCC_MAP_ELEC = m_dbe->book2D("HF Digi Elec Occupancy Map","HF Digi Elec Occupancy Map",20,0,20,20,0,20);

    m_dbe->setCurrentFolder("HcalMonitor/DigiMonitor/HO");
    hoHists.DIGI_NUM =  m_dbe->book1D("HO # of Digis","HO # of Digis",100,0,200);
    hoHists.DIGI_SIZE =  m_dbe->book1D("HO Digi Size","HO Digi Size",50,0,50);
    hoHists.DIGI_PRESAMPLE =  m_dbe->book1D("HO Digi Presamples","HO Digi Presamples",50,0,50);
    hoHists.QIE_CAPID =  m_dbe->book1D("HO QIE Cap-ID","HO QIE Cap-ID",6,-0.5,5.5);
    hoHists.QIE_ADC = m_dbe->book1D("HO QIE ADC Value","HO QIE ADC Value",100,0,200);
    hoHists.QIE_DV = m_dbe->book1D("HO QIE Data Value","HO QIE Data Value",2,-0.5,1.5);
    hoHists.ERR_MAP_GEO = m_dbe->book2D("HO Digi Geo Error Map","HO Digi Geo Error Map",59,-29.5,29.5,74,0,73);
    hoHists.ERR_MAP_ELEC = m_dbe->book2D("HO Digi Elec Error Map","HO Digi Elec Error Map",20,0,20,20,0,20);
    hoHists.OCC_MAP_GEO = m_dbe->book2D("HO Digi Geo Occupancy Map","HO Digi Geo Occupancy Map",59,-29.5,29.5,74,0,73);
    hoHists.OCC_MAP_ELEC = m_dbe->book2D("HO Digi Elec Occupancy Map","HO Digi Elec Occupancy Map",20,0,20,20,0,20);
    
}

  return;
}

void HcalDigiMonitor::processEvent(const HBHEDigiCollection& hbhe,
				   const HODigiCollection& ho,
				   const HFDigiCollection& hf)
{

  if(!m_dbe) { printf("HcalDigiMonitor::processEvent   DaqMonitorBEInterface not instantiated!!!\n");  return; }

  ievt_++;
  meEVT_->Fill(ievt_);

  try{
    hbHists.DIGI_NUM->Fill(hbhe.size());
    for (HBHEDigiCollection::const_iterator j=hbhe.begin(); j!=hbhe.end(); j++){
      const HBHEDataFrame digi = (const HBHEDataFrame)(*j);	
      HcalDigiMap::fillErrors<HBHEDataFrame>(digi,hbHists.ERR_MAP_GEO,hbHists.ERR_MAP_ELEC);	  
      HcalDigiMap::fillOccupancy<HBHEDataFrame>(digi,hbHists.OCC_MAP_GEO,hbHists.OCC_MAP_ELEC,occThresh_);	  
      hbHists.DIGI_SIZE->Fill(digi.size());
      hbHists.DIGI_PRESAMPLE->Fill(digi.presamples());
      int last = -1;
      //    printf("hb/he digi crate: %d, %d\n",digi.elecId().readoutVMECrateId(),digi.elecId().htrSlot());
      for (int i=0; i<digi.size(); i++) {	    
	hbHists.QIE_CAPID->Fill(digi.sample(i).capid());
	hbHists.QIE_ADC->Fill(digi.sample(i).adc());
	hbHists.QIE_CAPID->Fill(5,bitUpset(last,digi.sample(i).capid()));
	last = digi.sample(i).capid();
	hbHists.QIE_DV->Fill(0,digi.sample(i).dv());
	hbHists.QIE_DV->Fill(1,digi.sample(i).er());
      }    
    }
  } catch (...) {

    printf("HcalDigiMonitor::processEvent  No HBHE Digis.\n");
  }
  
  try{
     hoHists.DIGI_NUM->Fill(ho.size());
    for (HODigiCollection::const_iterator j=ho.begin(); j!=ho.end(); j++){
      const HODataFrame digi = (const HODataFrame)(*j);	
      HcalDigiMap::fillErrors<HODataFrame>(digi,hoHists.ERR_MAP_GEO,hoHists.ERR_MAP_ELEC);  
      HcalDigiMap::fillOccupancy<HODataFrame>(digi,hoHists.OCC_MAP_GEO,hoHists.OCC_MAP_ELEC,occThresh_);	  
      hoHists.DIGI_SIZE->Fill(digi.size());
      hoHists.DIGI_PRESAMPLE->Fill(digi.presamples());
      //     printf("ho digi crate: %d, %d\n",digi.elecId().readoutVMECrateId(),digi.elecId().htrSlot());
     int last = -1;
      for (int i=0; i<digi.size(); i++) {	    
	hoHists.QIE_CAPID->Fill(digi.sample(i).capid());
	hoHists.QIE_ADC->Fill(digi.sample(i).adc());
	hoHists.QIE_CAPID->Fill(5,bitUpset(last,digi.sample(i).capid()));
	last = digi.sample(i).capid();
	hoHists.QIE_DV->Fill(0,digi.sample(i).dv());
	hoHists.QIE_DV->Fill(1,digi.sample(i).er());
      }    
    }    
  } catch (...) {
    cout << "HcalDigiMonitor::processEvent  No HO Digis." << endl;
  }
  
  try{
    hfHists.DIGI_NUM->Fill(hf.size());
    for (HFDigiCollection::const_iterator j=hf.begin(); j!=hf.end(); j++){
      const HFDataFrame digi = (const HFDataFrame)(*j);	
      HcalDigiMap::fillErrors<HFDataFrame>(digi,hfHists.ERR_MAP_GEO,hfHists.ERR_MAP_ELEC); 
      HcalDigiMap::fillOccupancy<HFDataFrame>(digi,hfHists.OCC_MAP_GEO,hfHists.OCC_MAP_ELEC,occThresh_);	  
      hfHists.DIGI_SIZE->Fill(digi.size());
      hfHists.DIGI_PRESAMPLE->Fill(digi.presamples());
      //    printf("hf digi crate: %d, %d\n",digi.elecId().readoutVMECrateId(),digi.elecId().htrSlot());
      int last = -1;
      for (int i=0; i<digi.size(); i++) {	    
	hfHists.QIE_CAPID->Fill(digi.sample(i).capid());
	hfHists.QIE_ADC->Fill(digi.sample(i).adc());
	hfHists.QIE_CAPID->Fill(5,bitUpset(last,digi.sample(i).capid()));
	last = digi.sample(i).capid();
	hfHists.QIE_DV->Fill(0,digi.sample(i).dv());
	hfHists.QIE_DV->Fill(1,digi.sample(i).er());
      }
    }
  } catch (...) {
    cout << "HcalDigiMonitor::processEvent  No HF Digis." << endl;
  }

}
