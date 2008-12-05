#include <DQM/HcalMonitorModule/src/HcalMonitorModule.h>
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"

/*
 * \file HcalMonitorModule.cc
 * 
 * $Date: 2008/12/05 13:09:04 $
 * $Revision: 1.97 $
 * \author W Fisher
 *
*/

//--------------------------------------------------------
HcalMonitorModule::HcalMonitorModule(const edm::ParameterSet& ps){

  irun_=0; ilumisec_=0; ievent_=0; itime_=0;
  actonLS_=false;
  meStatus_=0;  meRunType_=0;
  meEvtMask_=0; meFEDS_=0;
  meLatency_=0; meQuality_=0;
  fedsListed_ = false;
  digiMon_ = NULL;   dfMon_ = NULL;
  diTask_ = NULL;
  rhMon_ = NULL;     pedMon_ = NULL; 
  ledMon_ = NULL;    mtccMon_ = NULL;
  hotMon_ = NULL;    tempAnalysis_ = NULL;
  deadMon_ = NULL;   tpMon_ = NULL;
  ctMon_ = NULL;     beamMon_ = NULL;
  laserMon_ = NULL;
  expertMon_ = NULL;

  // initialize hcal quality object
  

  // All subdetectors assumed out of the run by default
  HBpresent_=0;
  HEpresent_=0;
  HOpresent_=0;
  HFpresent_=0;

  inputLabelDigi_        = ps.getParameter<edm::InputTag>("digiLabel");
  inputLabelRecHitHBHE_  = ps.getParameter<edm::InputTag>("hbheRecHitLabel");
  inputLabelRecHitHF_    = ps.getParameter<edm::InputTag>("hfRecHitLabel");
  inputLabelRecHitHO_    = ps.getParameter<edm::InputTag>("hoRecHitLabel");
  inputLabelRecHitZDC_   = ps.getParameter<edm::InputTag>("zdcRecHitLabel");
  inputLabelCaloTower_   = ps.getParameter<edm::InputTag>("caloTowerLabel");
  inputLabelLaser_       = ps.getParameter<edm::InputTag>("hcalLaserLabel");

  checkHB_=ps.getUntrackedParameter<bool>("checkHB", 1); 
  checkHE_=ps.getUntrackedParameter<bool>("checkHE", 1);  
  checkHO_=ps.getUntrackedParameter<bool>("checkHO", 1);  
  checkHF_=ps.getUntrackedParameter<bool>("checkHF", 1);   

  evtSel_ = new HcalMonitorSelector(ps);
  
  dbe_ = Service<DQMStore>().operator->();
  
  debug_ = ps.getUntrackedParameter<int>("debug", 0);
  
  showTiming_ = ps.getUntrackedParameter<bool>("showTiming", false);
  dump2database_   = ps.getUntrackedParameter<bool>("dump2database",false); // dumps output to database file

  // Valgrind complained when the test was simply:  if ( ps.getUntrackedParameter<bool>("DataFormatMonitor", false))
  // try assigning value to bool first?
  bool taskOn = ps.getUntrackedParameter<bool>("DataFormatMonitor", false);
  if (taskOn) {
    if(debug_>0) cout << "HcalMonitorModule: DataFormat monitor flag is on...." << endl;
    dfMon_ = new HcalDataFormatMonitor();
    dfMon_->setup(ps, dbe_);
  }

  taskOn = ps.getUntrackedParameter<bool>("DataIntegrityTask", false); 
  if (taskOn ) 
    {
      if (debug_>0) cout <<"HcalMonitorModule: DataIntegrity monitor flag is on...."<<endl;
      diTask_ = new HcalDataIntegrityTask();
      diTask_->setup(ps, dbe_);
    }

  if ( ps.getUntrackedParameter<bool>("DigiMonitor", false) ) {
    if(debug_>0) cout << "HcalMonitorModule: Digi monitor flag is on...." << endl;
    digiMon_ = new HcalDigiMonitor();
    digiMon_->setup(ps, dbe_);
  }
  
  if ( ps.getUntrackedParameter<bool>("RecHitMonitor", false) ) {
    if(debug_>0) cout << "HcalMonitorModule: RecHit monitor flag is on...." << endl;
    rhMon_ = new HcalRecHitMonitor();
    rhMon_->setup(ps, dbe_);
  }
  
  if ( ps.getUntrackedParameter<bool>("PedestalMonitor", false) ) {
    if(debug_>0) cout << "HcalMonitorModule: Pedestal monitor flag is on...." << endl;
    pedMon_ = new HcalPedestalMonitor();
    pedMon_->setup(ps, dbe_);
  }
  
  if ( ps.getUntrackedParameter<bool>("LEDMonitor", false) ) {
    if(debug_>0) cout << "HcalMonitorModule: LED monitor flag is on...." << endl;
    ledMon_ = new HcalLEDMonitor();
    ledMon_->setup(ps, dbe_);
  }
  
  if ( ps.getUntrackedParameter<bool>("LaserMonitor", false) ) {
    if(debug_>0) cout << "HcalMonitorModule: Laser monitor flag is on...." << endl;
    laserMon_ = new HcalLaserMonitor();
    laserMon_->setup(ps, dbe_);
  }

  if ( ps.getUntrackedParameter<bool>("MTCCMonitor", false) ) {
    if(debug_>0) cout << "HcalMonitorModule: MTCC monitor flag is on...." << endl;
    mtccMon_ = new HcalMTCCMonitor();
    mtccMon_->setup(ps, dbe_);
  }
  
  if ( ps.getUntrackedParameter<bool>("HotCellMonitor", false) ) {
    if(debug_>0) cout << "HcalMonitorModule: Hot Cell monitor flag is on...." << endl;
    hotMon_ = new HcalHotCellMonitor();
    hotMon_->setup(ps, dbe_);
  }
  
  if ( ps.getUntrackedParameter<bool>("DeadCellMonitor", false) ) {
    if(debug_>0) cout << "HcalMonitorModule: Dead Cell monitor flag is on...." << endl;
    deadMon_ = new HcalDeadCellMonitor();
    deadMon_->setup(ps, dbe_);
  }

  if ( ps.getUntrackedParameter<bool>("TrigPrimMonitor", false) ) { 	 
    if(debug_>0) cout << "HcalMonitorModule: TrigPrim monitor flag is on...." << endl; 	 
    tpMon_ = new HcalTrigPrimMonitor(); 	 
    tpMon_->setup(ps, dbe_); 	 
  }  

  if (ps.getUntrackedParameter<bool>("CaloTowerMonitor",false)){
    if(debug_>0) cout << "HcalMonitorModule: CaloTower monitor flag is on...." << endl; 	 
    ctMon_ = new HcalCaloTowerMonitor(); 	 
    ctMon_->setup(ps, dbe_); 	 
  }  

  if (ps.getUntrackedParameter<bool>("BeamMonitor",false)){
    if(debug_>0) cout << "HcalMonitorModule: Beam monitor flag is on...."<<endl;
    beamMon_ = new HcalBeamMonitor();
    beamMon_->setup(ps, dbe_);
  }

  if (ps.getUntrackedParameter<bool>("ExpertMonitor",false)){
    if(debug_>0) cout << "HcalMonitorModule: Expert monitor flag is on...."<<endl;
    expertMon_ = new HcalExpertMonitor();
    expertMon_->setup(ps, dbe_);
  }

  

  if ( ps.getUntrackedParameter<bool>("HcalAnalysis", false) ) {
    if(debug_>0) cout << "HcalMonitorModule: Hcal Analysis flag is on...." << endl;
    tempAnalysis_ = new HcalTemplateAnalysis();
    tempAnalysis_->setup(ps);
  }
  

  // set parameters   
  prescaleEvt_ = ps.getUntrackedParameter<int>("diagnosticPrescaleEvt", -1);
  if(debug_>1) cout << "===>HcalMonitor event prescale = " << prescaleEvt_ << " event(s)"<< endl;

  prescaleLS_ = ps.getUntrackedParameter<int>("diagnosticPrescaleLS", -1);
  if(debug_>1) cout << "===>HcalMonitor lumi section prescale = " << prescaleLS_ << " lumi section(s)"<< endl;
  if (prescaleLS_>0) actonLS_=true;

  prescaleUpdate_ = ps.getUntrackedParameter<int>("diagnosticPrescaleUpdate", -1);
  if(debug_>1) cout << "===>HcalMonitor update prescale = " << prescaleUpdate_ << " update(s)"<< endl;

  prescaleTime_ = ps.getUntrackedParameter<int>("diagnosticPrescaleTime", -1);
  if(debug_>1) cout << "===>HcalMonitor time prescale = " << prescaleTime_ << " minute(s)"<< endl;
  
  // Base folder for the contents of this job
  string subsystemname = ps.getUntrackedParameter<string>("subSystemFolder", "Hcal") ;
  if(debug_>0) cout << "===>HcalMonitor name = " << subsystemname << endl;
  rootFolder_ = subsystemname + "/";
  
  gettimeofday(&psTime_.updateTV,NULL);
  /// get time in milliseconds, convert to minutes
  psTime_.updateTime = (psTime_.updateTV.tv_sec*1000.0+psTime_.updateTV.tv_usec/1000.0);
  psTime_.updateTime /= 1000.0;
  psTime_.elapsedTime=0;
  psTime_.vetoTime=psTime_.updateTime;
}

//--------------------------------------------------------
HcalMonitorModule::~HcalMonitorModule(){
  
// if (dbe_){    
//   if(digiMon_!=NULL)   {  digiMon_->clearME();}
//   if(dfMon_!=NULL)     {  dfMon_->clearME();}
//   if(diTask_!=NULL)    {  diTask_->clearME();}
//   if(pedMon_!=NULL)    {  pedMon_->clearME();}
//   if(ledMon_!=NULL)    {  ledMon_->clearME();}
//   if(laserMon_!=NULL)  {  laserMon_->clearME();}
//   if(hotMon_!=NULL)    {  hotMon_->clearME();}
//   if(deadMon_!=NULL)   {  deadMon_->clearME();}
//   if(mtccMon_!=NULL)   {  mtccMon_->clearME();}
//   if(rhMon_!=NULL)     {  rhMon_->clearME();}
//   
//   dbe_->setCurrentFolder(rootFolder_);
//   dbe_->removeContents();
// }
//
//  if(digiMon_!=NULL) { delete digiMon_;  digiMon_=NULL; }
//  if(dfMon_!=NULL) { delete dfMon_;     dfMon_=NULL; }
//  if(diTask_!=NULL) { delete diTask_;   diTask_=NULL; }
//  if(pedMon_!=NULL) { delete pedMon_;   pedMon_=NULL; }
//  if(ledMon_!=NULL) { delete ledMon_;   ledMon_=NULL; }
//  if(laserMon_!=NULL) { delete laserMon_;   laserMon_=NULL; }
//  if(hotMon_!=NULL) { delete hotMon_;   hotMon_=NULL; }
//  if(deadMon_!=NULL) { delete deadMon_; deadMon_=NULL; }
//  if(mtccMon_!=NULL) { delete mtccMon_; mtccMon_=NULL; }
//  if(rhMon_!=NULL) { delete rhMon_;     rhMon_=NULL; }
//  if(tempAnalysis_!=NULL) { delete tempAnalysis_; tempAnalysis_=NULL; }
//  delete evtSel_; evtSel_ = NULL;
//
} //void HcalMonitorModule::~HcalMonitorModule()

//--------------------------------------------------------
void HcalMonitorModule::beginJob(const edm::EventSetup& c){
  ievt_ = 0;
  
  ievt_pre_=0;

  if ( dbe_ != NULL ){
    dbe_->setCurrentFolder(rootFolder_+"DQM Job Status" );
    meStatus_  = dbe_->bookInt("STATUS");
    meRunType_ = dbe_->bookInt("RUN TYPE");
    meEvtMask_ = dbe_->bookInt("EVT MASK");
    meFEDS_    = dbe_->book1D("FEDs Unpacked","FEDs Unpacked",100,700,799);
    // process latency was (200,0,1), but that gave overflows
    meLatency_ = dbe_->book1D("Process Latency","Process Latency",2000,0,10);
    meQuality_ = dbe_->book1D("Quality Status","Quality Status",100,0,1);
    // Store whether or not subdetectors are present
    meHB_ = dbe_->bookInt("HBpresent");
    meHE_ = dbe_->bookInt("HEpresent");
    meHO_ = dbe_->bookInt("HOpresent");
    meHF_ = dbe_->bookInt("HFpresent");
    
    meStatus_->Fill(0);
    meRunType_->Fill(-1);
    meEvtMask_->Fill(-1);

    // Should fill with 0 to start
    meHB_->Fill(HBpresent_);
    meHE_->Fill(HEpresent_);
    meHO_->Fill(HOpresent_);
    meHF_->Fill(HFpresent_);
  }

  edm::ESHandle<HcalDbService> pSetup;
  c.get<HcalDbRecord>().get( pSetup );

  readoutMap_=pSetup->getHcalMapping();
  DetId detid_;
  HcalDetId hcaldetid_; 

  // Build a map of readout hardware unit to calorimeter channel
  std::vector <HcalElectronicsId> AllElIds = readoutMap_->allElectronicsIdPrecision();
  int dccid;
  pair <int,int> dcc_spgt;
  // by looping over all precision (non-trigger) items.
  for (std::vector <HcalElectronicsId>::iterator eid = AllElIds.begin();
       eid != AllElIds.end();
       eid++) {

    //Get the HcalDetId from the HcalElectronicsId
    detid_ = readoutMap_->lookup(*eid);
    

    // NULL if illegal; ignore
    if (!detid_.null()) {
      try {
	hcaldetid_ = HcalDetId(detid_);

	dccid = eid->dccid();
	dcc_spgt = pair <int,int> (dccid, eid->spigot());
      
	thisDCC = DCCtoCell.find(dccid);
	thisHTR = HTRtoCell.find(dcc_spgt);
      
	// If this DCC has no entries, make this its first one.
	if (thisDCC == DCCtoCell.end()) {
	  std::vector <HcalDetId> tempv;
	  tempv.push_back(hcaldetid_);
	  pair <int, std::vector<HcalDetId> > thispair;
	  thispair = pair <int, std::vector<HcalDetId> > (dccid,tempv);
	  DCCtoCell.insert(thispair); 
	}
	else {
	  thisDCC->second.push_back(hcaldetid_);
	}
      
	// If this HTR has no entries, make this its first one.
	if (thisHTR == HTRtoCell.end()) {
	  std::vector <HcalDetId> tempv;
	  tempv.push_back(hcaldetid_);
	  pair < pair <int,int>, std::vector<HcalDetId> > thispair;
	  thispair = pair <pair <int,int>, std::vector<HcalDetId> > (dcc_spgt,tempv);
	  HTRtoCell.insert(thispair); 
	}
	else {
	  thisHTR->second.push_back(hcaldetid_);	
	}

      } catch (...) {
      }
    } // fi (!detid_.null()) 
  } 
  if (dfMon_) {
    dfMon_->smuggleMaps(DCCtoCell, HTRtoCell);
  }

  //get conditions
  c.get<HcalDbRecord>().get(conditions_);

  // fill reference pedestals with database values
  // Need to repeat this so many times?  Just do it once? And then we can be smarter about the whole fC/ADC thing?
  if (pedMon_!=NULL)
    pedMon_->fillDBValues(*conditions_);
  if (deadMon_!=NULL)
    deadMon_->createMaps(*conditions_);
  if (hotMon_!=NULL)
    hotMon_->createMaps(*conditions_);


  edm::ESHandle<HcalChannelQuality> p;
  c.get<HcalChannelQualityRcd>().get(p);
  chanquality_= new HcalChannelQuality(*p.product());
  return;
} // HcalMonitorModule::beginJob(...)

//--------------------------------------------------------
void HcalMonitorModule::beginRun(const edm::Run& run, const edm::EventSetup& c) {
  fedsListed_ = false;

  // I think we want to reset these at 0 at the start of each run
  HBpresent_ = 0;
  HEpresent_ = 0;
  HOpresent_ = 0;
  HFpresent_ = 0;

  // Should fill with 0 to start
  meHB_->Fill(HBpresent_);
  meHE_->Fill(HEpresent_);
  meHO_->Fill(HOpresent_);
  meHF_->Fill(HFpresent_);
  reset();
}

//--------------------------------------------------------
void HcalMonitorModule::beginLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
     const edm::EventSetup& context) {
  
  if(actonLS_ && !prescale()){
    // do scheduled tasks...
  }
}


//--------------------------------------------------------
void HcalMonitorModule::endLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
					   const edm::EventSetup& context) {
  if(actonLS_ && !prescale()){
    // do scheduled tasks...
  }
}

//--------------------------------------------------------
void HcalMonitorModule::endRun(const edm::Run& r, const edm::EventSetup& context)
{
  if (debug_>1)  
    cout <<"HcalMonitorModule::endRun(...) "<<endl;
  // Do final pedestal histogram filling
  if (pedMon_!=NULL)
    pedMon_->fillPedestalHistos();
  if (deadMon_!=NULL)
    deadMon_->fillDeadHistosAtEndRun();

  return;
    }


//--------------------------------------------------------
void HcalMonitorModule::endJob(void) {
  
  if ( meStatus_ ) meStatus_->Fill(2);

  if(rhMon_!=NULL) rhMon_->done();
  if(digiMon_!=NULL) digiMon_->done();
  if(dfMon_!=NULL) dfMon_->done();
  if(diTask_!=NULL) diTask_->done();
  if(pedMon_!=NULL) pedMon_->done();
  if(ledMon_!=NULL) ledMon_->done();
  if(laserMon_!=NULL) laserMon_->done();
  if(hotMon_!=NULL) hotMon_->done(myquality_);
  if(deadMon_!=NULL) deadMon_->done(myquality_);
  if(mtccMon_!=NULL) mtccMon_->done();
  if (tpMon_!=NULL) tpMon_->done();
  if (ctMon_!=NULL) ctMon_->done();
  if (beamMon_!=NULL) beamMon_->done();
  if (expertMon_!=NULL) expertMon_->done();
  if(tempAnalysis_!=NULL) tempAnalysis_->done();

  if (dump2database_)
    {
      std::vector<DetId> mydetids = chanquality_->getAllChannels();
      HcalChannelQuality* newChanQual = new HcalChannelQuality();
      for (unsigned int i=0;i<mydetids.size();++i)
	{
	  if (mydetids[i].det()!=4) continue; // not hcal
	  //HcalDetId id(mydetids[i]);
	  HcalDetId id=mydetids[i];
	  // get original channel status item
	  const HcalChannelStatus* origstatus=chanquality_->getValues(mydetids[i]);
	  // make copy of status
	  HcalChannelStatus* mystatus=new HcalChannelStatus(origstatus->rawId(),origstatus->getValue());
	  if (myquality_.find(id)!=myquality_.end())
	    {
	      // Set bit 1 for cells which aren't present 	 
	      if ((id.subdet()==HcalBarrel &&!HBpresent_) || 	 
		  (id.subdet()==HcalEndcap &&!HEpresent_) || 	 
		  (id.subdet()==HcalOuter  &&!HOpresent_) || 	 
		  (id.subdet()==HcalForward&&!HFpresent_)) 	 
		{ 	 
		  mystatus->setBit(1); 	 
		} 	 
	      // Only perform these checks if bit 0 not set?
	      // check dead cells
	      if ((myquality_[id]>>5)&0x1)
		  mystatus->setBit(5);
	      else
		mystatus->unsetBit(5);
	      // check hot cells
	      if ((myquality_[id]>>6)&0x1)
		mystatus->setBit(6);
	      else
		mystatus->unsetBit(6);
	    } // if (myquality.find_...)
	  newChanQual->addValues(*mystatus);
	} // for (unsigned int i=0;...)
      // Now dump out to text file
      std::ostringstream file;
      file <<"HcalDQMstatus_"<<irun_<<".txt";
      std::ofstream outStream(file.str().c_str());
      HcalDbASCIIIO::dumpObject (outStream, (*newChanQual));
      /*
      std::ofstream dumb("orig.txt");
      HcalDbASCIIIO::dumpObject (dumb,(*chanquality_));
      */
    } // if (dump2databse_)
  return;
}

//--------------------------------------------------------
void HcalMonitorModule::reset(){

  if(rhMon_!=NULL)   rhMon_->reset();
  if(digiMon_!=NULL) digiMon_->reset();
  if(dfMon_!=NULL)   dfMon_->reset();
  if(diTask_!=NULL)  diTask_->reset();
  if(pedMon_!=NULL)  pedMon_->reset();
  if(ledMon_!=NULL)  ledMon_->reset();
  if(laserMon_!=NULL)  laserMon_->reset();
  if(hotMon_!=NULL)  hotMon_->reset();
  if(deadMon_!=NULL)  deadMon_->reset();
  if(mtccMon_!=NULL)   mtccMon_->reset();
  if(tempAnalysis_!=NULL) tempAnalysis_->reset();
  if(tpMon_!=NULL) tpMon_->reset();
  if(ctMon_!=NULL) ctMon_->reset();
  if(beamMon_!=NULL) beamMon_->reset();
  if(expertMon_!=NULL) expertMon_->reset();
}

//--------------------------------------------------------
void HcalMonitorModule::analyze(const edm::Event& e, const edm::EventSetup& eventSetup){

  // environment datamembers
  irun_     = e.id().run();
  ilumisec_ = e.luminosityBlock();
  ievent_   = e.id().event();
  itime_    = e.time().value();

  if (debug_>1) cout << "HcalMonitorModule: evts: "<< nevt_ << ", run: " << irun_ << ", LS: " << ilumisec_ << ", evt: " << ievent_ << ", time: " << itime_ << endl <<"\t counter = "<<ievt_pre_<<"\t total count = "<<ievt_<<endl; 

  // skip this event if we're prescaling...
  ievt_pre_++; // need to increment counter before calling prescale
  if(prescale()) return;

  meLatency_->Fill(psTime_.elapsedTime);

  // Do default setup...
  ievt_++;

  int evtMask=DO_HCAL_DIGIMON|DO_HCAL_DFMON|DO_HCAL_RECHITMON|DO_HCAL_PED_CALIBMON|DO_HCAL_LED_CALIBMON|DO_HCAL_LASER_CALIBMON; // add in DO_HCAL_TPMON, DO_HCAL_CTMON?  (in HcalMonitorSelector.h)

  //  int trigMask=0;
  if(mtccMon_==NULL){
    evtSel_->processEvent(e);
    evtMask = evtSel_->getEventMask();
    //    trigMask =  evtSel_->getTriggerMask();
  }
  if ( dbe_ ){ 
    meStatus_->Fill(1);
    meEvtMask_->Fill(evtMask);
  }
  
  ///See if our products are in the event...
  bool rawOK_    = true;
  bool digiOK_   = true;
  bool rechitOK_ = true;
  bool zdchitOK_ = true;
  bool trigOK_   = false;
  bool tpdOK_    = true;
  bool calotowerOK_ = true;
  bool laserOK_  = true;

  // try to get raw data and unpacker report
  edm::Handle<FEDRawDataCollection> rawraw;  

  try{
    e.getByType(rawraw);
  }
  catch(...)
    {
      rawOK_=false;
    }
  if (rawOK_&&!rawraw.isValid()) {
    rawOK_=false;
  }

  edm::Handle<HcalUnpackerReport> report;  
  try{
    e.getByType(report);
  }
  catch(...)
    {
      rawOK_=false;
    }
  if (rawOK_&&!report.isValid()) {
    rawOK_=false;
  }
  else 
    {
      if(!fedsListed_){
	const std::vector<int> feds =  (*report).getFedsUnpacked();    
	for(unsigned int f=0; f<feds.size(); f++){
	  meFEDS_->Fill(feds[f]);    
	}
	fedsListed_ = true;
      }
    }

  // try to get digis
  edm::Handle<HBHEDigiCollection> hbhe_digi;
  edm::Handle<HODigiCollection> ho_digi;
  edm::Handle<HFDigiCollection> hf_digi;
  edm::Handle<HcalTrigPrimDigiCollection> tp_digi;
  edm::Handle<HcalLaserDigi> laser_digi;

  try 
    {
      e.getByLabel(inputLabelDigi_,hbhe_digi);
    }
  catch(...)
    {
      digiOK_=false;
    }
  if (digiOK_&&!hbhe_digi.isValid()) {

    digiOK_=false;
  }

  try{
  e.getByLabel(inputLabelDigi_,hf_digi);
  }
  catch(...)
    {digiOK_=false;}
  if (digiOK_&&!hf_digi.isValid()) {
    digiOK_=false;
  }

  try
    {e.getByLabel(inputLabelDigi_,ho_digi);}
  catch(...)
    {digiOK_=false;}
  if (digiOK_&&!ho_digi.isValid()) {
    digiOK_=false;
  }

  // check which Subdetectors are on by seeing which are reading out FED data
  // Assume subdetectors aren't present, unless we explicitly find otherwise
  if ((checkHB_ && HBpresent_==0) ||
      (checkHE_ && HEpresent_==0) ||
      (checkHO_ && HOpresent_==0) ||
      (checkHF_ && HFpresent_==0))
  
    CheckSubdetectorStatus(*rawraw,*report,*readoutMap_,*hbhe_digi, *ho_digi, *hf_digi);
    
  // Case where all subdetectors have no raw data -- skip event
  if ((checkHB_ && HBpresent_==0) &&
      (checkHE_ && HEpresent_==0) &&
      (checkHO_ && HOpresent_==0) &&
      (checkHF_ && HFpresent_==0))
    {
      if (debug_>1) cout <<"<HcalMonitorModule::analyze>  No HCAL raw data found for event "<<ievt_<<endl;
      return;
    }

  try{
    e.getByLabel(inputLabelDigi_,tp_digi);
  }
  catch(...)
    {tpdOK_=false;}

  if (tpdOK_ && !tp_digi.isValid()) {
    tpdOK_=false;
  }
  try{
  e.getByLabel(inputLabelLaser_,laser_digi);
  }
  catch(...)
    {laserOK_=false;}
  if (laserOK_&&!laser_digi.isValid()) {
    laserOK_=false;
  }

  // try to get rechits
  edm::Handle<HBHERecHitCollection> hb_hits;
  edm::Handle<HORecHitCollection> ho_hits;
  edm::Handle<HFRecHitCollection> hf_hits;
  edm::Handle<ZDCRecHitCollection> zdc_hits;
  edm::Handle<CaloTowerCollection> calotowers;

  try{
  e.getByLabel(inputLabelRecHitHBHE_,hb_hits);
  }
  catch(...)
    {rechitOK_=false;}
  
  if (rechitOK_&&!hb_hits.isValid()) {
    rechitOK_ = false;
  }
  try{
  e.getByLabel(inputLabelRecHitHO_,ho_hits);
  }
  catch(...)
    {rechitOK_=false;}
  if (rechitOK_&&!ho_hits.isValid()) {
    rechitOK_ = false;
  }
  try{
    e.getByLabel(inputLabelRecHitHF_,hf_hits);
  }
  catch(...)
    {rechitOK_=false;}
  if (rechitOK_&&!hf_hits.isValid()) {
    rechitOK_ = false;
  }
  
  try{
    e.getByLabel(inputLabelRecHitZDC_,zdc_hits);
  }
  catch(...)
    {zdchitOK_=false;}
  if (zdchitOK_&&!zdc_hits.isValid()) {
    zdchitOK_ = false;
    //cout <<"CANNOT GET ZDC HITS!!!!"<<endl;
    //cout <<"input label = "<<inputLabelRecHitZDC_<<endl;
  }

  // try to get calotowers 
  if (ctMon_!=NULL)
    {
      try{
      e.getByLabel(inputLabelCaloTower_,calotowers);
      }
      catch(...)
	{calotowerOK_=false;}
      if(calotowerOK_&&!calotowers.isValid()){
	calotowerOK_=false;
      }
    }
  else
    calotowerOK_=false;

  // Run the configured tasks, protect against missing products

  // Data Format monitor task
  
  if (showTiming_)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

  if((dfMon_ != NULL) && (evtMask&DO_HCAL_DFMON) && rawOK_) 
    dfMon_->processEvent(*rawraw,*report,*readoutMap_);
  if (showTiming_)
    {
      cpu_timer.stop();
      if (dfMon_ !=NULL) cout <<"TIMER:: DATAFORMAT MONITOR ->"<<cpu_timer.cpuTime()<<endl;
      cpu_timer.reset(); cpu_timer.start();
    }

  if ((diTask_ != NULL) && (evtMask&DO_HCAL_DFMON) && rawOK_)
    diTask_->processEvent(*rawraw,*report,*readoutMap_);
  if (showTiming_)
    {
      cpu_timer.stop();
      if (diTask_ !=NULL) cout <<"TIMER:: DATA INTEGRITY TASK ->"<<cpu_timer.cpuTime()<<endl;
      cpu_timer.reset(); cpu_timer.start();
    }

  // Digi monitor task
  if((digiMon_!=NULL) && (evtMask&DO_HCAL_DIGIMON) && digiOK_) 
    {
      digiMon_->setSubDetectors(HBpresent_,HEpresent_, HOpresent_, HFpresent_);
      digiMon_->processEvent(*hbhe_digi,*ho_digi,*hf_digi,
			     *conditions_,*report);
    }
  if (showTiming_)
    {
      cpu_timer.stop();
      if (digiMon_ != NULL) cout <<"TIMER:: DIGI MONITOR ->"<<cpu_timer.cpuTime()<<endl;
      cpu_timer.reset(); cpu_timer.start();
    }
  // Pedestal monitor task
  if((pedMon_!=NULL) && (evtMask&DO_HCAL_PED_CALIBMON) && digiOK_) 
    pedMon_->processEvent(*hbhe_digi,*ho_digi,*hf_digi,*conditions_);
  if (showTiming_)
    {
      cpu_timer.stop();
      if (pedMon_!=NULL) cout <<"TIMER:: PEDESTAL MONITOR ->"<<cpu_timer.cpuTime()<<endl;
      cpu_timer.reset(); cpu_timer.start();
    }

  // LED monitor task
  if((ledMon_!=NULL) && (evtMask&DO_HCAL_LED_CALIBMON) && digiOK_)
    ledMon_->processEvent(*hbhe_digi,*ho_digi,*hf_digi,*conditions_);
  if (showTiming_)
    {
      cpu_timer.stop();
      if (ledMon_!=NULL) cout <<"TIMER:: LED MONITOR ->"<<cpu_timer.cpuTime()<<endl;
      cpu_timer.reset(); cpu_timer.start();
    }

  // Laser monitor task
  if((laserMon_!=NULL) && (evtMask&DO_HCAL_LASER_CALIBMON) && digiOK_ && laserOK_)
    laserMon_->processEvent(*hbhe_digi,*ho_digi,*hf_digi,*laser_digi,*conditions_);
  if (showTiming_)
    {
      cpu_timer.stop();
      if (laserMon_!=NULL) cout <<"TIMER:: LASER MONITOR ->"<<cpu_timer.cpuTime()<<endl;
      cpu_timer.reset(); cpu_timer.start();
    }

  // Rec Hit monitor task
  if((rhMon_ != NULL) && (evtMask&DO_HCAL_RECHITMON) && rechitOK_) 
    {
    rhMon_->processEvent(*hb_hits,*ho_hits,*hf_hits);
    // This lets us process rec hits regardless of ZDC status.
    // But is ZDC is okay, we'll make rec hit plots for that as well.
    if (zdchitOK_)
      {
	if (debug_>1) cout <<"PROCESSING ZDC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	//rhMon_->processZDC(*zdc_hits);
      }
    }
  if (showTiming_)
    {
      cpu_timer.stop();
      if (rhMon_!=NULL) cout <<"TIMER:: RECHIT MONITOR ->"<<cpu_timer.cpuTime()<<endl;
      cpu_timer.reset(); cpu_timer.start();
    }
  
  // Beam Monitor task
  if ((beamMon_ != NULL) && (evtMask&DO_HCAL_RECHITMON) && rechitOK_)
    {
      beamMon_->processEvent(*hb_hits,*ho_hits,*hf_hits,*hf_digi);
    }
  if (showTiming_)
    {
      cpu_timer.stop();
      if (beamMon_!=NULL) cout <<"TIMER:: BEAM MONITOR ->"<<cpu_timer.cpuTime( \
)<<endl;
      cpu_timer.reset(); cpu_timer.start();
    }

  // Hot Cell monitor task
  if((hotMon_ != NULL) && (evtMask&DO_HCAL_RECHITMON) && rechitOK_) 
    {
      hotMon_->processEvent(*hb_hits,*ho_hits,*hf_hits, 
			    *hbhe_digi,*ho_digi,*hf_digi,*conditions_);
      //hotMon_->setSubDetectors(HBpresent_,HEpresent_, HOpresent_, HFpresent_);
    }
  if (showTiming_)
    {
      cpu_timer.stop();
      if (hotMon_!=NULL) cout <<"TIMER:: HOTCELL MONITOR ->"<<cpu_timer.cpuTime()<<endl;
      cpu_timer.reset(); cpu_timer.start();
    }
  // Dead Cell monitor task -- may end up using both rec hits and digis?
  if((deadMon_ != NULL) && (evtMask&DO_HCAL_RECHITMON) && rechitOK_ && digiOK_) 
    {
      //deadMon_->setSubDetectors(HBpresent_,HEpresent_, HOpresent_, HFpresent_);
      deadMon_->processEvent(*hb_hits,*ho_hits,*hf_hits,
			     *hbhe_digi,*ho_digi,*hf_digi,*conditions_); 
    }
  if (showTiming_)
    {
      cpu_timer.stop();
      if (deadMon_!=NULL) cout <<"TIMER:: DEADCELL MONITOR ->"<<cpu_timer.cpuTime()<<endl;
      cpu_timer.reset(); cpu_timer.start();
    }

  // CalotowerMonitor
  if ((ctMon_ !=NULL) )
    ctMon_->processEvent(*calotowers);

  if (showTiming_)
    {
      cpu_timer.stop();
      if (ctMon_ !=NULL) cout <<"TIMER:: CALOTOWER MONITOR ->"<<cpu_timer.cpuTime()<<endl;
      cpu_timer.reset(); cpu_timer.start();
    }


  // Trigger Primitive monitor task -- may end up using both rec hits and digis?
  if((tpMon_ != NULL) && rechitOK_ && digiOK_ && tpdOK_) 
    tpMon_->processEvent(*hb_hits,*ho_hits,*hf_hits,
			 *hbhe_digi,*ho_digi,*hf_digi,*tp_digi, *readoutMap_);			     
  

  if (showTiming_)
    {
      cpu_timer.stop();
      if (tpMon_!=NULL) cout <<"TIMER:: TRIGGERPRIMITIVE MONITOR ->"<<cpu_timer.cpuTime()<<endl;
    }

  // Expert monitor plots
  if (expertMon_ != NULL) 
    {
      expertMon_->processEvent(*hb_hits,*ho_hits,*hf_hits,
			       *hbhe_digi,*ho_digi,*hf_digi,
			       *tp_digi,
			       *rawraw,*report,*readoutMap_);
    }
  if (showTiming_)
    {
      cpu_timer.stop();
      if (expertMon_!=NULL) cout <<"TIMER:: EXPERT MONITOR ->"<<cpu_timer.cpuTime()<<endl;
      cpu_timer.reset(); cpu_timer.start();
    }


  if(debug_>0 && ievt_%1000 == 0)
    cout << "HcalMonitorModule: processed " << ievt_ << " events" << endl;

  if(debug_>1)
    {
      cout << "HcalMonitorModule: processed " << ievt_ << " events" << endl;
      cout << "    RAW Data   ==> " << rawOK_<< endl;
      cout << "    Digis      ==> " << digiOK_<< endl;
      cout << "    RecHits    ==> " << rechitOK_<< endl;
      cout << "    TrigRec    ==> " << trigOK_<< endl;
      cout << "    TPdigis    ==> " << tpdOK_<< endl;    
      cout << "    CaloTower  ==> " << calotowerOK_ <<endl;
      cout << "    LaserDigis ==> " << laserOK_ << endl;
    }
  
  return;
}

//--------------------------------------------------------
bool HcalMonitorModule::prescale()
{
  ///Return true if this event should be skipped according to the prescale condition...
  ///    Accommodate a logical "OR" of the possible tests
  if (debug_>0) cout <<"HcalMonitorModule::prescale"<<endl;
  
  gettimeofday(&psTime_.updateTV,NULL);
  double time = (psTime_.updateTV.tv_sec*1000.0+psTime_.updateTV.tv_usec/1000.0);
  time/= (1000.0); ///in seconds
  psTime_.elapsedTime = time - psTime_.updateTime;
  psTime_.updateTime = time;
  //First determine if we care...
  bool evtPS =    prescaleEvt_>0;
  bool lsPS =     prescaleLS_>0;
  bool timePS =   prescaleTime_>0;
  bool updatePS = prescaleUpdate_>0;

  // If no prescales are set, keep the event
  if(!evtPS && !lsPS && !timePS && !updatePS)
    {
      return false;
    }
  //check each instance
  if(lsPS && (ilumisec_%prescaleLS_)!=0) lsPS = false; //LS veto
  //if(evtPS && (ievent_%prescaleEvt_)!=0) evtPS = false; //evt # veto
  // we can't just call (ievent_%prescaleEvt_) because ievent values not consecutive
  if (evtPS && (ievt_pre_%prescaleEvt_)!=0) evtPS = false;
  if(timePS)
    {
      double elapsed = (psTime_.updateTime - psTime_.vetoTime)/60.0;
      if(elapsed<prescaleTime_){
	timePS = false;  //timestamp veto
	psTime_.vetoTime = psTime_.updateTime;
      }
    } //if (timePS)

  //  if(prescaleUpdate_>0 && (nupdates_%prescaleUpdate_)==0) updatePS=false; ///need to define what "updates" means
  
  if (debug_>1) 
    {
      cout<<"HcalMonitorModule::prescale  evt: "<<ievent_<<"/"<<evtPS<<", ";
      cout <<"ls: "<<ilumisec_<<"/"<<lsPS<<",";
      cout <<"time: "<<psTime_.updateTime - psTime_.vetoTime<<"/"<<timePS<<endl;
    }  
  // if any criteria wants to keep the event, do so
  if(evtPS || lsPS || timePS) return false; //FIXME updatePS left out for now
  return true;
} // HcalMonitorModule::prescale(...)


void HcalMonitorModule::CheckSubdetectorStatus(const FEDRawDataCollection& rawraw, 
					       const HcalUnpackerReport& report, 
					       const HcalElectronicsMap& emap,
					       const HBHEDigiCollection& hbhedigi,
					       const HODigiCollection& hodigi,
					       const HFDigiCollection& hfdigi
					       //const ZDCDigiCollection& zdcdigi,

					       )
{
  vector<int> fedUnpackList;
  for (int i=FEDNumbering::getHcalFEDIds().first; 
       i<=FEDNumbering::getHcalFEDIds().second; 
       i++) 
    {
      fedUnpackList.push_back(i);
    }
  
  for (vector<int>::const_iterator i=fedUnpackList.begin();
       i!=fedUnpackList.end(); 
       ++i) 
    {
      const FEDRawData& fed = rawraw.FEDData(*i);
      if (fed.size()<12) continue; // Was 16. How do such tiny events even get here?
      
      // get the DCC header
      const HcalDCCHeader* dccHeader=(const HcalDCCHeader*)(fed.data());
      if (!dccHeader) return;
      int dccid=dccHeader->getSourceId();
      // check for HF
      if (dccid>717 && dccid<724)
	{
	  if (HFpresent_==0 && hfdigi.size()>0)
	    {
	      HFpresent_ = 1;
	      meHF_->Fill(HFpresent_);
	    }
	  continue;
	}

      // check for HO
      if (dccid>723)
	{
	  if (HOpresent_==0 && hodigi.size()>0)
	    {
	      HOpresent_ = 1;
	      meHO_->Fill(HOpresent_);
	    }
	  continue;
	}
      
      // Looking at HB and HE is more complicated, since they're combined into HBHE
      // walk through the HTR data...
      HcalHTRData htr;  
      for (int spigot=0; spigot<HcalDCCHeader::SPIGOT_COUNT; spigot++) {    
	if (!dccHeader->getSpigotPresent(spigot)) continue;
	
	// Load the given decoder with the pointer and length from this spigot.
	dccHeader->getSpigotData(spigot,htr, fed.size()); 
	
	// check min length, correct wordcount, empty event, or total length if histo event.
	if (!htr.check()) continue;
	if (htr.isHistogramEvent()) continue;
	
	int firstFED =  FEDNumbering::getHcalFEDIds().first; 
	
	// Tease out HB and HE, which share HTRs in HBHE
	for(int fchan=0; fchan<3; ++fchan) //0,1,2 are valid
	  {
	    for(int fib=1; fib<9; ++fib) //1...8 are valid
	      {
		HcalElectronicsId eid(fchan,fib,spigot,dccid-firstFED);
		eid.setHTR(htr.readoutVMECrateId(),
			   htr.htrSlot(),htr.htrTopBottom());
		DetId did=emap.lookup(eid);
		if (!did.null()) 
		  {
		    
		    switch (((HcalSubdetector)did.subdetId()))
		      {
		      case (HcalBarrel): 
			{
			  if (HBpresent_==0)
			    {
			      HBpresent_ = 1;
			      meHB_->Fill(HBpresent_);
			    }
			} break; // case (HcalBarrel)
		      case (HcalEndcap): 
			{
			  if (HEpresent_==0)
			    {
			      HEpresent_ = 1;
			      meHE_->Fill(HEpresent_);
			    }
			} break; // case (HcalEndcap)
		      case (HcalOuter): 
			{ // shouldn't reach these last two cases
			  if (HOpresent_==0)
			    {
			      {
				HOpresent_ = 1;
				meHO_->Fill(HOpresent_);
				return;
			      }
			    } 
			} break; // case (HcalOuter)
		      case (HcalForward): 
			{
			  if (HFpresent_==0)
			    {
			      meHF_->Fill(HFpresent_);
			      HFpresent_ = 1;
			    }
			} break; //case (HcalForward)
		      default: break;
		      } // switch ((HcalSubdetector...)
		  } // if (!did.null())
	      } // for (int fib=0;...)
	  } // for (int fchan = 0;...)
	
      } // for (int spigot=0;...)
    } //  for (vector<int>::const_iterator i=fedUnpackList.begin();
  return;
} // void HcalMonitorModule::CheckSubdetectorStatus(...)

#include "FWCore/Framework/interface/MakerMacros.h"
#include <DQM/HcalMonitorModule/src/HcalMonitorModule.h>
#include "DQMServices/Core/interface/DQMStore.h"

DEFINE_FWK_MODULE(HcalMonitorModule);
