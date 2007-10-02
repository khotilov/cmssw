#include <DQM/HcalMonitorModule/src/HcalMonitorModule.h>

/*
 * \file HcalMonitorModule.cc
 * 
 * $Date: 2007/08/09 00:33:06 $
 * $Revision: 1.36 $
 * \author W Fisher
 *
*/

HcalMonitorModule::HcalMonitorModule(const edm::ParameterSet& ps){

  printf("\n\nHcal Monitor Module\n\n");

  // verbosity switch
  m_verbose = ps.getUntrackedParameter<bool>("verbose", false);

  if(m_verbose) cout << "HcalMonitorModule: constructor...." << endl;

  m_logFile.open("HcalMonitorModule.log");

  m_dbe = NULL;
  if ( ps.getUntrackedParameter<bool>("DaqMonitorBEInterface", false)){
    m_dbe = edm::Service<DaqMonitorBEInterface>().operator->();
    m_dbe->setVerbose(0);
  }
  
  m_monitorDaemon = false;
  if ( ps.getUntrackedParameter<bool>("MonitorDaemon", false)){
    edm::Service<MonitorDaemon> daemon;
    daemon.operator->();
    m_monitorDaemon = true;
  }

  m_outputFile = ps.getUntrackedParameter<string>("outputFile", "");
  if ( m_outputFile.size() != 0 ) 
    cout << "Hcal Monitoring histograms will be saved to " << m_outputFile.c_str() << endl;    
  else m_outputFile = "DQM_Hcal_";


  bool disable = ps.getUntrackedParameter<bool>("disableROOToutput", false);
  if(disable) m_outputFile="";

  inputLabelDigi_        = ps.getParameter<edm::InputTag>("digiLabel");
  inputLabelRecHitHBHE_  = ps.getParameter<edm::InputTag>("hbheRecHitLabel");
  inputLabelRecHitHF_    = ps.getParameter<edm::InputTag>("hfRecHitLabel");
  inputLabelRecHitHO_    = ps.getParameter<edm::InputTag>("hoRecHitLabel");

  m_runNum = 0; m_meStatus=0;
  m_meRunNum=0; m_meRunType=0;
  m_meEvtNum=0; m_meEvtMask=0;
  
  if ( m_dbe != NULL ){
    m_dbe->setCurrentFolder("HcalMonitor");
    m_meStatus  = m_dbe->bookInt("STATUS");
    m_meRunNum  = m_dbe->bookInt("RUN NUMBER");
    m_meRunType = m_dbe->bookInt("RUN TYPE");
    m_meEvtNum  = m_dbe->bookInt("EVT NUMBER");
    m_meEvtMask = m_dbe->bookInt("EVT MASK");
    m_meBeamE   = m_dbe->bookInt("BEAM ENERGY");
    
    m_meTrigger = m_dbe->book1D("TB Trigger Type","TB Trigger Type",6,0,5);

    m_meStatus->Fill(-1);
    m_meRunNum->Fill(0);
    m_meRunType->Fill(-1);
    m_meEvtNum->Fill(-1);
    m_meEvtMask->Fill(-1);
    m_meBeamE->Fill(-1);
  }

  m_evtSel = new HcalMonitorSelector(ps);
  m_digiMon = NULL; m_dfMon = NULL; 
  m_rhMon = NULL;   m_pedMon = NULL; 
  m_ledMon = NULL;  m_mtccMon = NULL;
  m_hotMon = NULL;  m_tempAnalysis = NULL;
  m_commisMon = NULL; m_tpMon = NULL;

  if ( ps.getUntrackedParameter<bool>("DataFormatMonitor", false) ) {
    if(m_verbose) cout << "HcalMonitorModule: DataFormat monitor flag is on...." << endl;
    m_dfMon = new HcalDataFormatMonitor();
    m_dfMon->setup(ps, m_dbe);
  }

  if ( ps.getUntrackedParameter<bool>("DigiMonitor", false) ) {
     if(m_verbose) cout << "HcalMonitorModule: Digi monitor flag is on...." << endl;
   m_digiMon = new HcalDigiMonitor();
    m_digiMon->setup(ps, m_dbe);
  }

  if ( ps.getUntrackedParameter<bool>("RecHitMonitor", false) ) {
    if(m_verbose) cout << "HcalMonitorModule: RecHit monitor flag is on...." << endl;
    m_rhMon = new HcalRecHitMonitor();
    m_rhMon->setup(ps, m_dbe);
  }

  if ( ps.getUntrackedParameter<bool>("TrigPrimMonitor", false) ) {
    if(m_verbose) cout << "HcalMonitorModule: TrigPrim monitor flag is on...." << endl;
    m_tpMon = new HcalTrigPrimMonitor();
    m_tpMon->setup(ps, m_dbe);
  }

  if ( ps.getUntrackedParameter<bool>("PedestalMonitor", false) ) {
    if(m_verbose) cout << "HcalMonitorModule: Pedestal monitor flag is on...." << endl;
    m_pedMon = new HcalPedestalMonitor();
    m_pedMon->setup(ps, m_dbe);
  }

  if ( ps.getUntrackedParameter<bool>("LEDMonitor", false) ) {
    if(m_verbose) cout << "HcalMonitorModule: LED monitor flag is on...." << endl;
    m_ledMon = new HcalLEDMonitor();
    m_ledMon->setup(ps, m_dbe);
  }

  if ( ps.getUntrackedParameter<bool>("MTCCMonitor", false) ) {
    if(m_verbose) cout << "HcalMonitorModule: MTCC monitor flag is on...." << endl;
    m_mtccMon = new HcalMTCCMonitor();
    m_mtccMon->setup(ps, m_dbe);
  }

  if ( ps.getUntrackedParameter<bool>("HotCellMonitor", false) ) {
    if(m_verbose) cout << "HcalMonitorModule: Hot Cell monitor flag is on...." << endl;
    m_hotMon = new HcalHotCellMonitor();
    m_hotMon->setup(ps, m_dbe);
  }

  if ( ps.getUntrackedParameter<bool>("CommissioningMonitor", false) ) {
    if(m_verbose) cout << "HcalMonitorModule: Commissioning monitor flag is on...." << endl;

    m_commisMon = new HcalCommissioningMonitor();
    m_commisMon->setup(ps, m_dbe);
  }

  if ( ps.getUntrackedParameter<bool>("HcalAnalysis", false) ) {
    if(m_verbose) cout << "HcalMonitorModule: Hcal Analysis flag is on...." << endl;

    m_tempAnalysis = new HcalTemplateAnalysis();
    m_tempAnalysis->setup(ps);
  }

  offline_ = ps.getUntrackedParameter<bool>("OffLine", false);

}

HcalMonitorModule::~HcalMonitorModule(){
  
  if(m_verbose) printf("HcalMonitorModule: Destructor.....");

  if ( offline_ ) sleep(15); 

  if (m_dbe && !offline_){    
    if(m_digiMon!=NULL) {  m_digiMon->clearME();}
    if(m_dfMon!=NULL)   {  m_dfMon->clearME();}
    if(m_pedMon!=NULL)  {  m_pedMon->clearME();}
    if(m_ledMon!=NULL)  {  m_ledMon->clearME();}
    if(m_hotMon!=NULL)  {  m_hotMon->clearME();}
    if(m_commisMon!=NULL) {  m_commisMon->clearME();}
    if(m_mtccMon!=NULL) {  m_mtccMon->clearME();}
    if(m_rhMon!=NULL)   {  m_rhMon->clearME();}
    if(m_tpMon!=NULL)   {  m_tpMon->clearME();}
    
    m_dbe->setCurrentFolder("HcalMonitor");
    m_dbe->removeContents();
  }
  
  if(m_digiMon!=NULL) { delete m_digiMon; m_digiMon=NULL; }
  if(m_dfMon!=NULL) { delete m_dfMon; m_dfMon=NULL; }
  if(m_pedMon!=NULL) { delete m_pedMon; m_pedMon=NULL; }
  if(m_ledMon!=NULL) { delete m_ledMon; m_ledMon=NULL; }
  if(m_hotMon!=NULL) { delete m_hotMon; m_hotMon=NULL; }
  if(m_commisMon!=NULL) { delete m_commisMon; m_commisMon=NULL; }
  if(m_mtccMon!=NULL) { delete m_mtccMon; m_mtccMon=NULL; }
  if(m_rhMon!=NULL) { delete m_rhMon; m_rhMon=NULL; }
  if(m_tpMon!=NULL) { delete m_tpMon; m_tpMon=NULL; }
  if(m_tempAnalysis!=NULL) { delete m_tempAnalysis; m_tempAnalysis=NULL; }
  delete m_evtSel;

  m_logFile.close();
}

void HcalMonitorModule::beginJob(const edm::EventSetup& c){
  m_ievt = 0;
  
  if(m_verbose) cout << "HcalMonitorModule: begin job...." << endl;
  
  if ( m_meStatus ) m_meStatus->Fill(0);
  if ( m_meEvtNum ) m_meEvtNum->Fill(m_ievt);
  
  // get the hcal mapping
  edm::ESHandle<HcalDbService> pSetup;
  c.get<HcalDbRecord>().get( pSetup );
  m_readoutMap=pSetup->getHcalMapping();
  
  // get conditions
  c.get<HcalDbRecord>().get(m_conditions);

  return;
}

void HcalMonitorModule::endJob(void) {

  if(m_verbose) cout << "HcalMonitorModule: end job...." << endl;  
  cout << "HcalMonitorModule::endJob, analyzed " << m_ievt << " events" << endl;
  
  if ( m_meStatus ) m_meStatus->Fill(2);
  if ( m_meRunNum ) m_meRunNum->Fill(m_runNum);
  if ( m_meEvtNum ) m_meEvtNum->Fill(m_ievt);

  if(m_rhMon!=NULL) m_rhMon->done();
  if(m_tpMon!=NULL) m_tpMon->done();
  if(m_digiMon!=NULL) m_digiMon->done();
  if(m_dfMon!=NULL) m_dfMon->done();
  if(m_pedMon!=NULL) m_pedMon->done();
  if(m_ledMon!=NULL) m_ledMon->done();
  if(m_hotMon!=NULL) m_hotMon->done();
  if(m_commisMon!=NULL) m_commisMon->done();
  if(m_mtccMon!=NULL) m_mtccMon->done();
  if(m_tempAnalysis!=NULL) m_tempAnalysis->done();

  char tmp[150]; bool update = true;
  for ( unsigned int i = 0; i < m_outputFile.size(); i++ ) {
    if ( m_outputFile.substr(i, 5) == ".root" )  {
      update = false;
    }
  }
  string saver = m_outputFile;
  if(update){
    sprintf(tmp,"%09d.root", m_runNum);
    saver = m_outputFile+tmp;
  }
  if ( m_outputFile.size() != 0  && m_dbe ) m_dbe->save(saver);

  return;
}

void HcalMonitorModule::analyze(const edm::Event& e, const edm::EventSetup& eventSetup){

  if(m_verbose) cout << "HcalMonitorModule: analyze...." << endl;

  // Do default setup...
  m_ievt++;
  int evtMask=DO_HCAL_DIGIMON|DO_HCAL_DFMON|DO_HCAL_RECHITMON|DO_HCAL_PED_CALIBMON;

  int trigMask=0;
  if(m_mtccMon==NULL){
    m_evtSel->processEvent(e);
    evtMask = m_evtSel->getEventMask();
    trigMask =  m_evtSel->getTriggerMask();

    if(trigMask&0x01) m_meTrigger->Fill(1);
    if(trigMask&0x02) m_meTrigger->Fill(2);
    if(trigMask&0x04) m_meTrigger->Fill(3);
    if(trigMask&0x08) m_meTrigger->Fill(4);
    if(trigMask&0x10) m_meTrigger->Fill(5);
  }

  edm::EventID id_ = e.id();
  m_runNum = (int)(id_.run());

  if ( m_dbe ){ 
    m_meStatus->Fill(1);
    m_meRunNum->Fill(m_runNum);
    m_meEvtNum->Fill(m_ievt);
    m_meEvtMask->Fill(evtMask);
  }
  
  bool rawOK_    = true;
  bool digiOK_   = true;
  bool rechitOK_ = true;
  bool ltcOK_    = true;
  bool tpdOK_    = true;

  // try to get raw data and unpacker report
  edm::Handle<FEDRawDataCollection> rawraw;  
  try{e.getByType(rawraw);} catch(exception& ex){ rawOK_=false;};           
  edm::Handle<HcalUnpackerReport> report;  
  try{e.getByType(report);} catch(exception& ex){ rawOK_=false;};

  // try to get digis
  edm::Handle<HBHEDigiCollection> hbhe_digi;
  edm::Handle<HODigiCollection> ho_digi;
  edm::Handle<HFDigiCollection> hf_digi;
  edm::Handle<HcalTrigPrimDigiCollection> tp_digi;
  try{e.getByLabel(inputLabelDigi_,hbhe_digi);} catch(exception& ex){digiOK_=false;};
  try{e.getByLabel(inputLabelDigi_,hf_digi);} catch(exception& ex){digiOK_=false;};
  try{e.getByLabel(inputLabelDigi_,ho_digi);} catch(exception& ex){digiOK_=false;};
  try{e.getByLabel(inputLabelDigi_,tp_digi);} catch(exception& ex){tpdOK_=false;};

  // try to get rechits
  edm::Handle<HBHERecHitCollection> hb_hits;
  edm::Handle<HORecHitCollection> ho_hits;
  edm::Handle<HFRecHitCollection> hf_hits;
  try{e.getByLabel(inputLabelRecHitHBHE_,hb_hits);} catch(exception& ex){rechitOK_ = false;}; 
  try{e.getByLabel(inputLabelRecHitHO_,ho_hits);} catch(exception& ex){rechitOK_ = false;}; 
  try{e.getByLabel(inputLabelRecHitHF_,hf_hits);} catch(exception& ex){rechitOK_ = false;}; 

  // try to get trigger info
  edm::Handle<LTCDigiCollection> ltc;
  try{ e.getByType(ltc); } catch(exception& ex){ltcOK_=false;};         


  /// Run the configured tasks

  // Data Format monitor task
  if((m_dfMon != NULL) && (evtMask&DO_HCAL_DFMON) && rawOK_) 
    m_dfMon->processEvent(*rawraw,*report,*m_readoutMap);

  // Digi monitor task
  if((m_digiMon!=NULL) && (evtMask&DO_HCAL_DIGIMON) && digiOK_) 
    m_digiMon->processEvent(*hbhe_digi,*ho_digi,*hf_digi,*m_conditions);

  // Pedestal monitor task
  if((m_pedMon!=NULL) && (evtMask&DO_HCAL_PED_CALIBMON) && digiOK_) 
    m_pedMon->processEvent(*hbhe_digi,*ho_digi,*hf_digi,*m_conditions);

  // LED monitor task
  //  if((m_ledMon!=NULL) && (evtMask&DO_HCAL_LED_CALIBMON)) m_ledMon->processEvent(*hbhe_digi,*ho_digi,*hf_digi);
  if(m_ledMon!=NULL && digiOK_) 
    m_ledMon->processEvent(*hbhe_digi,*ho_digi,*hf_digi,*m_conditions, *report);
  
  // Rec Hit monitor task
  if((m_rhMon != NULL) && (evtMask&DO_HCAL_RECHITMON) && rechitOK_) 
    m_rhMon->processEvent(*hb_hits,*ho_hits,*hf_hits);

  // Rec Hit monitor task
  if((m_tpMon != NULL) && (evtMask&DO_HCAL_RECHITMON) && rechitOK_ && digiOK_) 
    m_tpMon->processEvent(*hb_hits,*ho_hits,*hf_hits,*tp_digi);

  // Hot Cell monitor task
  if((m_hotMon != NULL) && (evtMask&DO_HCAL_RECHITMON) && rechitOK_) 
    m_hotMon->processEvent(*hb_hits,*ho_hits,*hf_hits);

  // Old MTCC monitor task
  if(m_mtccMon != NULL && digiOK_ && ltcOK_) m_mtccMon->processEvent(*hbhe_digi,*ho_digi, *ltc,*m_conditions);


  if(m_commisMon != NULL && digiOK_ && ltcOK_ && rechitOK_) m_commisMon->processEvent(*hbhe_digi,*ho_digi, *hf_digi,
						    *hb_hits,*ho_hits,*hf_hits,
						    *ltc,*m_conditions);

  if(m_tempAnalysis != NULL && digiOK_ && ltcOK_ && rechitOK_) 
    m_tempAnalysis->processEvent(*hbhe_digi,*ho_digi, *hf_digi,
				 *hb_hits,*ho_hits,*hf_hits,
				 *ltc,*m_conditions);



  if(m_ievt%1000 == 0)
    cout << "HcalMonitorModule: analyzed " << m_ievt << " events" << endl;

  if(m_verbose){
    cout << "HcalMonitorModule: analyzed " << m_ievt << " events" << endl;
    cout << "    RAW Data==> " << rawOK_<< endl;
    cout << "    Digis   ==> " << digiOK_<< endl;
    cout << "    RecHits ==> " << rechitOK_<< endl;
    cout << "    LTCdigi ==> " << ltcOK_<< endl;
    cout << "    TPdigis ==> " << tpdOK_<< endl;    
  }


  if(offline_) usleep(30);

  return;
}

#include "FWCore/Framework/interface/MakerMacros.h"
#include <DQM/HcalMonitorModule/src/HcalMonitorModule.h>

DEFINE_FWK_MODULE(HcalMonitorModule);
