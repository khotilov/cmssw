#include "EcalTPGParamBuilder.h"

#include "CalibCalorimetry/EcalTPGTools/plugins/EcalTPGDBApp.h"


#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

#include "SimCalorimetry/EcalSimAlgos/interface/EcalSimParameterMap.h"

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <time.h>

#include <TF1.h>
#include <iomanip>
#include <fstream>


double oneOverEtResolEt(double *x, double *par) { 
  double Et = x[0] ;
  if (Et<1e-6) return 1./par[1] ; // to avoid division by 0.
  double resolEt_overEt = sqrt( (par[0]/sqrt(Et))*(par[0]/sqrt(Et)) + (par[1]/Et)*(par[1]/Et) + par[2]*par[2] ) ;
  return 1./(Et*resolEt_overEt) ;
}

EcalTPGParamBuilder::EcalTPGParamBuilder(edm::ParameterSet const& pSet)
  : xtal_LSB_EB_(0), xtal_LSB_EE_(0), nSample_(5), complement2_(7)
{
  std::cout<<"here we are in EcalTPGParamBuilder::EcalTPGParamBuilder"<<endl;

  ped_conf_id_=0;
  lin_conf_id_=0;
  lut_conf_id_=0;
  wei_conf_id_=0;
  fgr_conf_id_=0;
  sli_conf_id_=0;
  bxt_conf_id_=0;
  btt_conf_id_=0;
  tag_="";
  version_=0;

  readFromDB_ = pSet.getParameter<bool>("readFromDB") ;
  writeToDB_  = pSet.getParameter<bool>("writeToDB") ;
  DBEE_ = pSet.getParameter<bool>("allowDBEE") ;
  string DBsid    = pSet.getParameter<std::string>("DBsid") ;
  string DBuser   = pSet.getParameter<std::string>("DBuser") ;
  string DBpass   = pSet.getParameter<std::string>("DBpass") ;
  uint32_t DBport = pSet.getParameter<unsigned int>("DBport") ;
  DBrunNb_        = pSet.getParameter<unsigned int>("DBrunNb") ;

  tag_   = pSet.getParameter<std::string>("TPGtag") ;
  version_ = pSet.getParameter<unsigned int>("TPGversion") ;

  std::cout << "data will be saved with tag and version="<< tag_<< ".version"<<version_<< endl;

  std::cout << "DB RUN NB="<< DBrunNb_<< endl;
 
  if (readFromDB_ || writeToDB_) {
    try {
      cout << "Warning: using the DB is not yet implemented " <<endl ;
      db_ = new EcalTPGDBApp(DBsid, DBuser, DBpass) ;
    } catch (exception &e) {
      cout << "ERROR:  " << e.what() << endl;
    } catch (...) {
      cout << "Unknown error caught" << endl;
    }
  }

  writeToFiles_ =  pSet.getParameter<bool>("writeToFiles") ;
  if (writeToFiles_) {
    std::string outFile = pSet.getParameter<std::string>("outFile") ;
    out_file_ = new std::ofstream(outFile.c_str(), std::ios::out) ;  
    geomFile_   = new std::ofstream("geomFile.txt", std::ios::out) ;  
  }



  useTransverseEnergy_ = pSet.getParameter<bool>("useTransverseEnergy") ;
  
  Et_sat_EB_ = pSet.getParameter<double>("Et_sat_EB") ;
  Et_sat_EE_ = pSet.getParameter<double>("Et_sat_EE") ;
  sliding_ = pSet.getParameter<unsigned int>("sliding") ;
  sampleMax_ = pSet.getParameter<unsigned int>("weight_sampleMax") ;

  forcedPedestalValue_ = pSet.getParameter<int>("forcedPedestalValue") ;
  forceEtaSlice_ = pSet.getParameter<bool>("forceEtaSlice") ;
    
  LUT_option_ = pSet.getParameter<std::string>("LUT_option") ;
  LUT_threshold_EB_ = pSet.getParameter<double>("LUT_threshold_EB") ;
  LUT_threshold_EE_ = pSet.getParameter<double>("LUT_threshold_EE") ;
  LUT_stochastic_EB_ = pSet.getParameter<double>("LUT_stochastic_EB") ;
  LUT_noise_EB_ =pSet.getParameter<double>("LUT_noise_EB") ;
  LUT_constant_EB_ =pSet.getParameter<double>("LUT_constant_EB") ;
  LUT_stochastic_EE_ = pSet.getParameter<double>("LUT_stochastic_EE") ;
  LUT_noise_EE_ =pSet.getParameter<double>("LUT_noise_EE") ;
  LUT_constant_EE_ =pSet.getParameter<double>("LUT_constant_EE") ;

  TTF_lowThreshold_EB_ = pSet.getParameter<double>("TTF_lowThreshold_EB") ;
  TTF_highThreshold_EB_ = pSet.getParameter<double>("TTF_highThreshold_EB") ;
  TTF_lowThreshold_EE_ = pSet.getParameter<double>("TTF_lowThreshold_EE") ;
  TTF_highThreshold_EE_ = pSet.getParameter<double>("TTF_highThreshold_EE") ;

  FG_lowThreshold_EB_ = pSet.getParameter<double>("FG_lowThreshold_EB") ;
  FG_highThreshold_EB_ = pSet.getParameter<double>("FG_highThreshold_EB") ;
  FG_lowRatio_EB_ = pSet.getParameter<double>("FG_lowRatio_EB") ;
  FG_highRatio_EB_ = pSet.getParameter<double>("FG_highRatio_EB") ;
  FG_lut_EB_ = pSet.getParameter<unsigned int>("FG_lut_EB") ;
  FG_Threshold_EE_ = pSet.getParameter<double>("FG_Threshold_EE") ;
  FG_lut_strip_EE_ = pSet.getParameter<unsigned int>("FG_lut_strip_EE") ;
  FG_lut_tower_EE_ = pSet.getParameter<unsigned int>("FG_lut_tower_EE") ;

  std::cout<<"here we are in EcalTPGParamBuilder::EcalTPGParamBuilder done"<<endl;


}

EcalTPGParamBuilder::~EcalTPGParamBuilder()
{ 
  if (writeToFiles_) {
    (*out_file_ )<<"EOF"<<std::endl ;
    out_file_->close() ;
    delete out_file_ ;
  }
}


bool EcalTPGParamBuilder::checkIfOK(     EcalPedestals::Item item) {

  bool result=true;
  if( item.mean_x1 <150. || item.mean_x1 >250) result=false;
  if( item.mean_x6 <150. || item.mean_x6 >250) result=false;
  if( item.mean_x12 <150. || item.mean_x12 >250) result=false;
  if( item.rms_x1 <0 || item.rms_x1 > 2) result=false;
  if( item.rms_x6 <0 || item.rms_x1 > 3) result=false;
  if( item.rms_x12 <0 || item.rms_x1 > 5) result=false;
  return result; 

}

void EcalTPGParamBuilder::analyze(const edm::Event& evt, const edm::EventSetup& evtSetup) 
{
  using namespace edm;
  using namespace std;

  std::cout<< "we are in analyze"<< endl; 

  // geometry
  ESHandle<CaloGeometry> theGeometry;
  ESHandle<CaloSubdetectorGeometry> theEndcapGeometry_handle, theBarrelGeometry_handle;
  evtSetup.get<CaloGeometryRecord>().get( theGeometry );
  evtSetup.get<EcalEndcapGeometryRecord>().get("EcalEndcap",theEndcapGeometry_handle);
  evtSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel",theBarrelGeometry_handle);
  evtSetup.get<IdealGeometryRecord>().get(eTTmap_);
  theEndcapGeometry_ = &(*theEndcapGeometry_handle);
  theBarrelGeometry_ = &(*theBarrelGeometry_handle);

  // electronics mapping
  ESHandle< EcalElectronicsMapping > ecalmapping;
  evtSetup.get< EcalMappingRcd >().get(ecalmapping);
  theMapping_ = ecalmapping.product();


  ////////////////////////////
  // Initialization section //
  ////////////////////////////
  list<uint32_t> towerListEB ;
  list<uint32_t> stripListEB ;
  list<uint32_t> towerListEE ;
  list<uint32_t> stripListEE ;
  list<uint32_t>::iterator itList ;


  std::cout <<"we get the pedestals from offline DB"<<endl;

  // Pedestals
  ESHandle<EcalPedestals> pedHandle;
  evtSetup.get<EcalPedestalsRcd>().get( pedHandle );
  const EcalPedestalsMap & pedMap = pedHandle.product()->getMap() ;

   

  // we copy the last valid record to a temporary object peds
  EcalPedestals* peds = new EcalPedestals();
  for(int iEta=-EBDetId::MAX_IETA; iEta<=EBDetId::MAX_IETA ;++iEta) {
    if(iEta==0) continue;
    for(int iPhi=EBDetId::MIN_IPHI; iPhi<=EBDetId::MAX_IPHI; ++iPhi) {
      // make an EBDetId since we need EBDetId::rawId() to be used as the key for the pedestals
      if (EBDetId::validDetId(iEta,iPhi))
	{
	  EBDetId ebdetid(iEta,iPhi);

	  EcalPedestals::Item aped= *(pedMap.find(ebdetid));

	  // here I copy the last valid value in the peds object
	  EcalPedestals::Item item;
	  item.mean_x1  = aped.mean_x1;
	  item.rms_x1   = aped.rms_x1;
	  item.mean_x6  = aped.mean_x6;
	  item.rms_x6   = aped.rms_x6;
	  item.mean_x12 = aped.mean_x12;
	  item.rms_x12  = aped.rms_x12;

	  peds->insert(std::make_pair(ebdetid.rawId(),item));
	}
    }
  }

// I comment this one 
// because the best pedestals are already in the offline DB
 
/*  
std::cout <<"we get the pedestals from online DB"<<endl;
  map<EcalLogicID, MonPedestalsDat> pedMapDB ;
  int iovId = 0 ;
  std::cout << "DB RUN NB="<< DBrunNb_<< endl;

  if (readFromDB_) {
    iovId = db_->readFromCondDB_Pedestals(pedMapDB, DBrunNb_) ;

    typedef map<EcalLogicID, MonPedestalsDat>::const_iterator CImon;
    EcalLogicID ecid_xt;
    MonPedestalsDat  rd_ped;
          
    for (CImon p = pedMapDB.begin(); p != pedMapDB.end(); p++) {
      ecid_xt = p->first;
      rd_ped  = p->second;
      int sm_num=ecid_xt.getID1();
      int xt_num=ecid_xt.getID2(); 
      
      EcalPedestals::Item item;
      item.mean_x1  =rd_ped.getPedMeanG1() ;
      item.rms_x1   =rd_ped.getPedRMSG1();
      item.mean_x6  =rd_ped.getPedMeanG6();
      item.rms_x6   =rd_ped.getPedRMSG6() ;
      item.mean_x12 =rd_ped.getPedMeanG12();
      item.rms_x12  =rd_ped.getPedRMSG12();
      
      EBDetId ebdetid(sm_num,xt_num,EBDetId::SMCRYSTALMODE);
      // here we change in the peds object only the values that are available in the online DB 
      // otherwise we keep the old value
      if(checkIfOK(item)) peds->insert(std::make_pair(ebdetid.rawId(),item));
    }
  }
  // now peds is complete 

  */


  const EcalPedestalsMap & pedMapNew = peds->getMap() ;


  std::cout <<"we get the intercalib from offline DB"<<endl;
  // Intercalib constants
  ESHandle<EcalIntercalibConstants> pIntercalib ;
  evtSetup.get<EcalIntercalibConstantsRcd>().get(pIntercalib) ;
  const EcalIntercalibConstants * intercalib = pIntercalib.product() ;
  const EcalIntercalibConstantMap & calibMap = intercalib->getMap() ;


  std::cout <<"we get the gain ratios from offline DB"<<endl;
  // Gain Ratios
  ESHandle<EcalGainRatios> pRatio;
  evtSetup.get<EcalGainRatiosRcd>().get(pRatio);
  const EcalGainRatioMap & gainMap = pRatio.product()->getMap();
  

  std::cout <<"we get the ADC to GEV from offline DB"<<endl;
  // ADCtoGeV
  ESHandle<EcalADCToGeVConstant> pADCToGeV ;
  evtSetup.get<EcalADCToGeVConstantRcd>().get(pADCToGeV) ;
  const EcalADCToGeVConstant * ADCToGeV = pADCToGeV.product() ;
  xtal_LSB_EB_ = ADCToGeV->getEBValue() ;
  xtal_LSB_EE_ = ADCToGeV->getEEValue() ;
  std::cout<<"xtal_LSB_EB_ = "<<xtal_LSB_EB_<<std::endl ;
  std::cout<<"xtal_LSB_EE_ = "<<xtal_LSB_EE_<<std::endl ;

  
  vector<EcalLogicID> my_EcalLogicId;
  vector<EcalLogicID> my_TTEcalLogicId;
  vector<EcalLogicID> my_StripEcalLogicId;
  EcalLogicID my_EcalLogicId_EB;
  EcalLogicID my_EcalLogicId_EE;
  if (writeToDB_ || readFromDB_){
    std::cout<<"going to get the ecal logic id set"<< endl;

    my_EcalLogicId_EB = db_->getEcalLogicID( "EB",EcalLogicID::NULLID,EcalLogicID::NULLID,EcalLogicID::NULLID,"EB");
    my_EcalLogicId_EE = db_->getEcalLogicID( "EE",EcalLogicID::NULLID,EcalLogicID::NULLID,EcalLogicID::NULLID,"EE");

    my_EcalLogicId = db_->getEcalLogicIDSetOrdered( "EB_crystal_number",
						    1, 36,
						    1, 1700,
						    EcalLogicID::NULLID,EcalLogicID::NULLID,
						    "EB_crystal_number",12 );
    my_TTEcalLogicId = db_->getEcalLogicIDSetOrdered( "EB_trigger_tower",
						    1, 36,
						    1, 68,
						    EcalLogicID::NULLID,EcalLogicID::NULLID,
						    "EB_trigger_tower",12 );
    my_StripEcalLogicId = db_->getEcalLogicIDSetOrdered( "EB_VFE",   1, 36,   1, 68,   1,5 ,  "EB_VFE",12 );
    std::cout<<"got the 3 ecal barrel logic id set"<< endl;

  }

  /////////////////////////////////////////
  // Compute linearization coeff section //
  /////////////////////////////////////////

  map<EcalLogicID, FEConfigPedDat> pedset ;
  map<EcalLogicID, FEConfigLinDat> linset ;
  map<EcalLogicID, FEConfigParamDat> linparamset ;

  map<int, linStruc> linEtaSlice ;

  // loop on EB xtals
  if (writeToFiles_) (*out_file_)<<"COMMENT ====== barrel crystals ====== "<<std::endl ;
  const std::vector<DetId>& ebCells = theBarrelGeometry_->getValidDetIds(DetId::Ecal, EcalBarrel);
  std::cout <<" number of valid ebcells "<<ebCells.size()<<std::endl;

  // special case of eta slices
  for (vector<DetId>::const_iterator it = ebCells.begin(); it != ebCells.end(); ++it) {
    EBDetId id(*it) ;
    double theta = theBarrelGeometry_->getGeometry(id)->getPosition().theta() ;
    if (!useTransverseEnergy_) theta = acos(0.) ;
    const EcalTrigTowerDetId towid= id.tower();
    towerListEB.push_back(towid.rawId()) ;
    const EcalTriggerElectronicsId elId = theMapping_->getTriggerElectronicsId(id) ;
    int dccNb = theMapping_->DCCid(towid) ;
    int tccNb = theMapping_->TCCid(towid) ;
    int towerInTCC = theMapping_->iTT(towid) ; // from 1 to 68 (EB)
    int stripInTower = elId.pseudoStripId() ;  // from 1 to 5
    int xtalInStrip = elId.channelId() ;       // from 1 to 5

    if (tccNb == 37 && stripInTower == 3 && xtalInStrip == 3 && (towerInTCC-1)%4==0) {
      int etaSlice = (towerInTCC-1)/4+1 ;
      coeffStruc coeff ;
      getCoeff(coeff, calibMap, id.rawId()) ;
      getCoeff(coeff, gainMap, id.rawId()) ;
      getCoeff(coeff, pedMapNew, id.rawId()) ;
      linStruc lin ;
      for (int i=0 ; i<3 ; i++) {
	int mult, shift ;
	bool ok = computeLinearizerParam(theta, coeff.gainRatio_[i], coeff.calibCoeff_, "EB", mult , shift) ;
	if (!ok) std::cout << "unable to compute the parameters for SM="<< id.ism()<<" xt="<< id.ic()<<" " <<dec<<id.rawId()<<std::endl ;  
	else {
	  lin.pedestal_[i] = coeff.pedestals_[i] ;
	  lin.mult_[i] = mult ;
	  lin.shift_[i] = shift ;
	}
      }
      linEtaSlice[etaSlice] = lin ;
    }
  }

  // general case
  for (vector<DetId>::const_iterator it = ebCells.begin(); it != ebCells.end(); ++it) {
    EBDetId id(*it) ;
    double theta = theBarrelGeometry_->getGeometry(id)->getPosition().theta() ;
    if (!useTransverseEnergy_) theta = acos(0.) ;
    const EcalTrigTowerDetId towid= id.tower();
    towerListEB.push_back(towid.rawId()) ;
    const EcalTriggerElectronicsId elId = theMapping_->getTriggerElectronicsId(id) ;
    stripListEB.push_back(elId.rawId() & 0xfffffff8) ;
    int dccNb = theMapping_->DCCid(towid) ;
    int tccNb = theMapping_->TCCid(towid) ;
    int towerInTCC = theMapping_->iTT(towid) ; // from 1 to 68 (EB)
    int stripInTower = elId.pseudoStripId() ;  // from 1 to 5
    int xtalInStrip = elId.channelId() ;       // from 1 to 5
    int etaSlice = (towerInTCC-1)/4+1 ;

    FEConfigPedDat ped ;
    FEConfigLinDat lin ;
    if (writeToFiles_) (*out_file_)<<"CRYSTAL "<<dec<<id.rawId()<<std::endl ;
    //  if (writeToDB_ || readFromDB_) logicId = db_->getEcalLogicID ("EB_crystal_number", id.ism(), id.ic()) ;

    coeffStruc coeff ;
    getCoeff(coeff, calibMap, id.rawId()) ;
    getCoeff(coeff, gainMap, id.rawId()) ;
    getCoeff(coeff, pedMapNew, id.rawId()) ;
    

    // compute and fill linearization parameters

    // case of eta slice
    if (forceEtaSlice_) {
      map<int, linStruc>::const_iterator itLin = linEtaSlice.find(etaSlice);
      if (itLin != linEtaSlice.end()) {
	if (writeToFiles_) {
	  for (int i=0 ; i<3 ; i++) 
	    (*out_file_) << hex <<" 0x"<<itLin->second.pedestal_[i]<<" 0x"<<itLin->second.mult_[i]<<" 0x"<<itLin->second.shift_[i]<<std::endl;
	}
	if (writeToDB_) {
	  for (int i=0 ; i<3 ; i++) {
	    if (i==0)  {ped.setPedMeanG12(itLin->second.pedestal_[i]) ; lin.setMultX12(itLin->second.mult_[i]) ; lin.setShift12(itLin->second.shift_[i]) ; } 
	    if (i==1)  {ped.setPedMeanG6(itLin->second.pedestal_[i]) ; lin.setMultX6(itLin->second.mult_[i]) ; lin.setShift6(itLin->second.shift_[i]) ; } 
	    if (i==2)  {ped.setPedMeanG1(itLin->second.pedestal_[i]) ; lin.setMultX1(itLin->second.mult_[i]) ; lin.setShift1(itLin->second.shift_[i]) ; } 
	  }
	}
      }
    }
    else {
      // general case
      for (int i=0 ; i<3 ; i++) {
	int mult, shift ;
	bool ok = computeLinearizerParam(theta, coeff.gainRatio_[i], coeff.calibCoeff_, "EB", mult , shift) ;
	if (!ok) std::cout << "unable to compute the parameters for SM="<< id.ism()<<" xt="<< id.ic()<<" " <<dec<<id.rawId()<<std::endl ;  
	else {
	  if (writeToFiles_) (*out_file_) << hex <<" 0x"<<coeff.pedestals_[i]<<" 0x"<<mult<<" 0x"<<shift<<std::endl; 
	  if (writeToDB_) {
	    if (i==0)  {ped.setPedMeanG12(coeff.pedestals_[i]) ; lin.setMultX12(mult) ; lin.setShift12(shift) ; } 
	    if (i==1)  {ped.setPedMeanG6(coeff.pedestals_[i]) ; lin.setMultX6(mult) ; lin.setShift6(shift) ; } 
	    if (i==2)  {ped.setPedMeanG1(coeff.pedestals_[i]) ; lin.setMultX1(mult) ; lin.setShift1(shift) ; }
	  }
	}
      }
    }
    if (writeToDB_) {
      int ixtal=(id.ism()-1)*1700+(id.ic()-1);
      EcalLogicID logicId =my_EcalLogicId[ixtal];
      pedset[logicId] = ped ;
      linset[logicId] = lin ;	
    }
  } //ebCells

  if (writeToDB_) {
    // EcalLogicID  of the whole barrel is: my_EcalLogicId_EB
    FEConfigParamDat linparam ;
    linparam.setETSat(Et_sat_EB_); 
    linparam.setTTThreshlow(TTF_lowThreshold_EB_); 
    linparam.setTTThreshhigh(TTF_highThreshold_EB_); 
    linparam.setFGlowthresh(FG_lowThreshold_EB_); 
    linparam.setFGhighthresh(FG_highThreshold_EB_); 
    linparam.setFGlowratio(FG_lowRatio_EB_); 
    linparam.setFGhighratio(FG_highRatio_EB_);
    linparamset[my_EcalLogicId_EB] = linparam ;
  }


  // loop on EE xtals
  if (writeToFiles_) (*out_file_)<<"COMMENT ====== endcap crystals ====== "<<std::endl ;

  
  const std::vector<DetId> & eeCells = theEndcapGeometry_->getValidDetIds(DetId::Ecal, EcalEndcap);
  for (vector<DetId>::const_iterator it = eeCells.begin(); it != eeCells.end(); ++it) {
    EEDetId id(*it);
    double theta = theEndcapGeometry_->getGeometry(id)->getPosition().theta() ;
    if (!useTransverseEnergy_) theta = acos(0.) ;
    const EcalTrigTowerDetId towid= (*eTTmap_).towerOf(id) ;
    towerListEE.push_back(towid.rawId()) ;
    // special case of towers in inner rings of EE
    if (towid.ietaAbs() == 27 || towid.ietaAbs() == 28) {
      EcalTrigTowerDetId additionalTower(towid.zside(), towid.subDet(), towid.ietaAbs(), towid.iphi()+1) ;
      towerListEE.push_back(additionalTower.rawId()) ;
    }
    const EcalTriggerElectronicsId elId = theMapping_->getTriggerElectronicsId(id) ;
    stripListEE.push_back(elId.rawId() & 0xfffffff8) ;
    int dccNb = theMapping_->DCCid(towid) ;
    int tccNb = theMapping_->TCCid(towid) ;
    int towerInTCC = theMapping_->iTT(towid) ;
    int stripInTower = elId.pseudoStripId() ;
    int xtalInStrip = elId.channelId() ;

    EcalLogicID logicId ;
    FEConfigPedDat ped ;
    FEConfigLinDat lin ;
    if (writeToFiles_) (*out_file_)<<"CRYSTAL "<<dec<<id.rawId()<<std::endl ;
    if ((writeToDB_ || readFromDB_) && DBEE_) {
      int iz = id.positiveZ() ;
      if (iz ==0) iz = -1 ;
      logicId = db_->getEcalLogicID ("EE_crystal_number", iz, id.ix(), id.iy()) ;
    }

// comment to can write to the OnlineDB at P5    
    coeffStruc coeff ;
    if (readFromDB_ && DBEE_) {
      getCoeff(coeff, calibMap, id.rawId()) ;
      getCoeff(coeff, gainMap, id.rawId()) ;
      getCoeff(coeff, pedMap, id.rawId()) ;
    }
  
    // compute and fill linearization parameters
    for (int i=0 ; i<3 ; i++) {
      int mult, shift ;
      bool ok = computeLinearizerParam(theta, coeff.gainRatio_[i], coeff.calibCoeff_, "EE", mult , shift) ;
      if (!ok) std::cout << "unable to compute the parameters for "<<dec<<id.rawId()<<std::endl ;  
      else {
	if (writeToFiles_) (*out_file_) << hex <<" 0x"<<coeff.pedestals_[i]<<" 0x"<<mult<<" 0x"<<shift<<std::endl; 
	if (writeToDB_ && DBEE_) {
	  if (i==0)  {ped.setPedMeanG12(coeff.pedestals_[i]) ; lin.setMultX12(mult) ; lin.setShift12(shift) ; } 
	  if (i==1)  {ped.setPedMeanG6(coeff.pedestals_[i]) ; lin.setMultX6(mult) ; lin.setShift6(shift) ; } 
	  if (i==2)  {ped.setPedMeanG1(coeff.pedestals_[i]) ; lin.setMultX1(mult) ; lin.setShift1(shift) ; } 
	}	
      }
    }
    if (writeToDB_ && DBEE_) {
      pedset[logicId] = ped ;
      linset[logicId] = lin ;
    }
  } //eeCells

  if (writeToDB_ ) {
    // EcalLogicID  of the whole barrel is: my_EcalLogicId_EB
    FEConfigParamDat linparam ;
    linparam.setETSat(Et_sat_EE_); 
    linparam.setTTThreshlow(TTF_lowThreshold_EE_); 
    linparam.setTTThreshhigh(TTF_highThreshold_EE_); 
    linparam.setFGlowthresh(FG_Threshold_EE_); 
    linparam.setFGhighthresh(FG_Threshold_EE_); 
    linparamset[my_EcalLogicId_EE] = linparam ;
  }

  std::cout<< "we are in analyze 2"<< endl; 


  if (writeToDB_) {
    ped_conf_id_=db_->writeToConfDB_TPGPedestals(pedset, 1, "from_OfflineDB") ;
    lin_conf_id_=db_->writeToConfDB_TPGLinearCoef(linset,linparamset, 1, "from_CondDB") ;
  }

  /////////////////////////////
  // Compute weights section //
  /////////////////////////////

  // loading reference signal representation
  EcalSimParameterMap parameterMap;  
  EBDetId   barrel(1,1);
  double    phase = parameterMap.simParameters(barrel).timePhase();
  EcalShape shape(phase); 
  std::vector<unsigned int> weights = computeWeights(shape) ;

  if (weights.size() == 5) {
    if (writeToFiles_) {
      (*out_file_) <<std::endl ;
      (*out_file_) <<"WEIGHT 0"<<endl ;
      for (uint sample=0 ; sample<5 ; sample++) (*out_file_) << "0x" <<hex<<weights[sample]<<" " ;
      (*out_file_)<<std::endl ; 
    }
    if (writeToDB_) {
      std::cout<<"going to write the weights "<< endl;
      map<EcalLogicID, FEConfigWeightGroupDat> dataset;
      // we create 1 group
      int NWEIGROUPS =1; 
      for (int ich=0; ich<NWEIGROUPS; ich++){
	FEConfigWeightGroupDat gut;
	gut.setWeightGroupId(ich);
	gut.setWeight0(weights[0] );
	gut.setWeight1(weights[1] );
	gut.setWeight2(weights[2] );
	gut.setWeight3(weights[3] );
	gut.setWeight4(weights[4] );
	EcalLogicID ecid = EcalLogicID( "DUMMY", ich,ich);
	// Fill the dataset
	dataset[ecid] = gut; // we use any logic id but different, because it is in any case ignored... 
      }
      
      // now we store in the DB the correspondence btw channels and groups 
      map<EcalLogicID, FEConfigWeightDat> dataset2;
      // in this case I decide in a stupid way which channel belongs to which group 
      for (int ich=0; ich<my_StripEcalLogicId.size() ; ich++){
	FEConfigWeightDat wut;
	int igroup=0;
	wut.setWeightGroupId(igroup);
	// in case of real groups one has to look for the right logic id 
	// the logic ids are ordered in the vector by SM (EB+ 1,..,18, EB-, 19,..36), TT (1,...68), strip (1,...5)
	// Fill the dataset
	dataset2[my_StripEcalLogicId[ich]] = wut;
      }

      // endcap loop missing ... FIXME 
      //
      //
      //

      // Insert the dataset
      ostringstream wtag;
      wtag.str(""); wtag<<"SimShape_Phase"<<phase<<"_NGroups_"<<NWEIGROUPS;
      std::string weight_tag=wtag.str();
      std::cout<< " weight tag "<<weight_tag<<endl; 
      wei_conf_id_=db_->writeToConfDB_TPGWeight(dataset, dataset2, NWEIGROUPS, weight_tag) ;
      
    }
  }

  /////////////////////////
  // Compute FG section //
  /////////////////////////

  // barrel
  uint lowRatio, highRatio, lowThreshold, highThreshold, lutFG ;
  computeFineGrainEBParameters(lowRatio, highRatio, lowThreshold, highThreshold, lutFG) ;
  if (writeToFiles_) {
    (*out_file_) <<std::endl ;
    (*out_file_) <<"FG 0"<<std::endl ;
    (*out_file_)<<hex<<"0x"<<lowThreshold<<" 0x"<<highThreshold
		  <<" 0x"<<lowRatio<<" 0x"<<highRatio<<" 0x"<<lutFG
		  <<std::endl ;
  }

  // endcap
  uint threshold, lut_strip, lut_tower ;
  computeFineGrainEEParameters(threshold, lut_strip, lut_tower) ; 

  // and here we store the fgr part

  
  if (writeToDB_) {
    std::cout<<"going to write the fgr "<< endl;
      map<EcalLogicID, FEConfigFgrGroupDat> dataset;
      // we create 1 group
      int NFGRGROUPS =1; 
      for (int ich=0; ich<NFGRGROUPS; ich++){
	FEConfigFgrGroupDat gut;
	gut.setFgrGroupId(ich);
	gut.setThreshLow(lowRatio );
	gut.setThreshHigh(highRatio);
	gut.setRatioLow(lowThreshold);
	gut.setRatioHigh(highThreshold);
	gut.setLUTConfId(lutFG);
	EcalLogicID ecid = EcalLogicID( "DUMMY", ich,ich);
	// Fill the dataset
	dataset[ecid] = gut; // we use any logic id but different, because it is in any case ignored... 
      }
      
      // now we store in the DB the correspondence btw channels and groups 
      map<EcalLogicID, FEConfigFgrDat> dataset2;
      // in this case I decide in a stupid way which channel belongs to which group 
      for (int ich=0; ich<my_TTEcalLogicId.size() ; ich++){
	FEConfigFgrDat wut;
	int igroup=0;
	wut.setFgrGroupId(igroup);
	// Fill the dataset
	// the logic ids are ordered by SM (1,...36) and TT (1,...68)  
	// you have to calculate the right index here 
	dataset2[my_TTEcalLogicId[ich]] = wut;
      }

      // endcap loop missing ... FIXME 
      //
      //
      //

      // Insert the dataset
      ostringstream wtag;
      wtag.str(""); wtag<<"FGR_"<<lutFG<<"_NGroups_"<<NFGRGROUPS;
      std::string weight_tag=wtag.str();
      std::cout<< " weight tag "<<weight_tag<<endl; 
      fgr_conf_id_=db_->writeToConfDB_TPGFgr(dataset, dataset2, NFGRGROUPS, weight_tag) ;
  }

  if (writeToDB_) {
    std::cout<<"going to write the sliding "<< endl;
      map<EcalLogicID, FEConfigSlidingDat> dataset;
      // in this case I decide in a stupid way which channel belongs to which group 
      for (int ich=0; ich<my_StripEcalLogicId.size() ; ich++){
	FEConfigSlidingDat wut;
	wut.setSliding(sliding_);
	// Fill the dataset
	// the logic ids are ordered by SM (1,...36) , TT (1,...68) and strip (1..5) 
	// you have to calculate the right index here 
	dataset[my_StripEcalLogicId[ich]] = wut;
      }

      // endcap loop missing ... FIXME 
      //
      //
      //

      // Insert the dataset
      ostringstream wtag;
      wtag.str(""); wtag<<"Sliding_"<<sliding_;
      std::string justatag=wtag.str();
      std::cout<< " sliding tag "<<justatag<<endl;
      int iov_id=0; // just a parameter ... 
      sli_conf_id_=db_->writeToConfDB_TPGSliding(dataset,iov_id, justatag) ;
  }

  



  /////////////////////////
  // Compute LUT section //
  /////////////////////////

  int lut_EB[1024], lut_EE[1024] ;

  // barrel
  computeLUT(lut_EB, "EB") ; 
  if (writeToFiles_) {
    (*out_file_) <<std::endl ;
    (*out_file_) <<"LUT 0"<<std::endl ;
    for (int i=0 ; i<1024 ; i++) (*out_file_)<<"0x"<<hex<<lut_EB[i]<<endl ;
    (*out_file_)<<endl ;
  }
  
  // endcap
  computeLUT(lut_EE, "EE") ;
  // check first if lut_EB and lut_EE are the same
  bool newLUT(false) ;
  for (int i=0 ; i<1024 ; i++) if (lut_EE[i] != lut_EB[i]) newLUT = true ;
  if (newLUT && writeToFiles_) { 
    (*out_file_) <<std::endl ;
    (*out_file_) <<"LUT 1"<<std::endl ;
    for (int i=0 ; i<1024 ; i++) (*out_file_)<<"0x"<<hex<<lut_EE[i]<<endl ;
    (*out_file_)<<endl ;
  }



  if (writeToDB_) {
    map<EcalLogicID, FEConfigLUTGroupDat> dataset;
    // we create 1 LUT group
    int NLUTGROUPS =1; 
    for (int ich=0; ich<NLUTGROUPS; ich++){
      FEConfigLUTGroupDat lut;
      lut.setLUTGroupId(ich);
      for (int i=0; i<1024; i++){
	lut.setLUTValue(i, lut_EB[i] );
      }
      EcalLogicID ecid = EcalLogicID( "DUMMY", ich,ich);
      // Fill the dataset
      dataset[ecid] = lut; // we use any logic id but different, because it is in any case ignored... 
    }

    // now we store in the DB the correspondence btw channels and LUT groups 
    map<EcalLogicID, FEConfigLUTDat> dataset2;
    // in this case I decide in a stupid way which channel belongs to which group 
    for (int ich=0; ich<my_TTEcalLogicId.size() ; ich++){
      FEConfigLUTDat lut;
      int igroup=0;
      lut.setLUTGroupId(igroup);
      // calculate the right TT - in the vector they are ordere by SM and by TT 
      // Fill the dataset
      dataset2[my_TTEcalLogicId[ich]] = lut;
    }

    // endcap loop missing ... FIXME 
    //
    //
    //

    // Insert the dataset
    ostringstream ltag;
    ltag.str(""); ltag<<LUT_option_<<"_NGroups_"<<NLUTGROUPS;
    std::string lut_tag=ltag.str();
    std::cout<< " LUT tag "<<lut_tag<<endl; 
    lut_conf_id_=db_->writeToConfDB_TPGLUT(dataset, dataset2, NLUTGROUPS, lut_tag) ;

  }

  // last we insert the FE_CONFIG_MAIN table 
 if (writeToDB_) {
   
   int conf_id_=db_->writeToConfDB_TPGMain(ped_conf_id_,lin_conf_id_, lut_conf_id_, fgr_conf_id_, 
					sli_conf_id_, wei_conf_id_, bxt_conf_id_, btt_conf_id_, tag_, version_) ;

 }


  ///////////////////////////////////////////////////////////
  // loop on strips and associate them with default values //
  ///////////////////////////////////////////////////////////

  // Barrel
  stripListEB.sort() ;
  stripListEB.unique() ;
  cout<<"Number of EB strips="<<stripListEB.size()<<endl ;
  if (writeToFiles_) {
    (*out_file_) <<std::endl ;
    for (itList = stripListEB.begin(); itList != stripListEB.end(); itList++ ) {
      (*out_file_) <<"STRIP_EB "<<dec<<(*itList)<<endl ;
      (*out_file_) << hex << "0x" <<sliding_<<std::endl ;
      (*out_file_) <<" 0" << std::endl ;
    }
  }

  // Endcap
  stripListEE.sort() ;
  stripListEE.unique() ;
  cout<<"Number of EE strips="<<stripListEE.size()<<endl ;
  if (writeToFiles_) {
    (*out_file_) <<std::endl ;
    for (itList = stripListEE.begin(); itList != stripListEE.end(); itList++ ) {
      (*out_file_) <<"STRIP_EE "<<dec<<(*itList)<<endl ;
      (*out_file_) << hex << "0x" <<sliding_<<std::endl ;
      (*out_file_) <<" 0" << std::endl ;
      (*out_file_)<<hex<<"0x"<<threshold<<" 0x"<<lut_strip<<std::endl ;  
    }
  }


  ///////////////////////////////////////////////////////////
  // loop on towers and associate them with default values //
  ///////////////////////////////////////////////////////////

  // Barrel
  towerListEB.sort() ;
  towerListEB.unique() ;
  cout<<"Number of EB towers="<<towerListEB.size()<<endl ;
  if (writeToFiles_) {
    (*out_file_) <<std::endl ;
    (*geomFile_)<<"BARREL"<<endl ;
    for (itList = towerListEB.begin(); itList != towerListEB.end(); itList++ ) {
      (*out_file_) <<"TOWER_EB "<<dec<<(*itList)<<endl ;
      (*out_file_) <<" 0\n 0\n" ;
      EcalTrigTowerDetId towerId((*itList)) ;
      int dccNb = theMapping_->DCCid(towerId) ;
      int tccNb = theMapping_->TCCid(towerId) ;
      int towerInTCC = theMapping_->iTT(towerId) ;
      (*geomFile_)<<"towerId="<<(*itList)<<" ieta="<<towerId.ietaAbs()<<" iphi="<<towerId.iphi()
		  <<" dccNb="<<dccNb<<" tccNb="<<tccNb<<" towerInTCC="<<towerInTCC<<endl ;
    }
  }

  // Endcap
  towerListEE.sort() ;
  towerListEE.unique() ;
  cout<<"Number of EE towers="<<towerListEE.size()<<endl ;
  if (writeToFiles_) {
    (*out_file_) <<std::endl ;
    (*geomFile_)<<"ENDCAP"<<endl ;
    for (itList = towerListEE.begin(); itList != towerListEE.end(); itList++ ) {
      (*out_file_) <<"TOWER_EE "<<dec<<(*itList)<<endl ;
      if (newLUT) (*out_file_) <<" 1\n" ;
      else  (*out_file_) <<" 0\n" ;
      (*out_file_)<<hex<<"0x"<<lut_tower<<std::endl ;
      EcalTrigTowerDetId towerId((*itList)) ;
      int dccNb = theMapping_->DCCid(towerId) ;
      int tccNb = theMapping_->TCCid(towerId) ;
      int towerInTCC = theMapping_->iTT(towerId) ;
      (*geomFile_)<<"towerId="<<(*itList)<<" ieta="<<towerId.ietaAbs()<<" iphi="<<towerId.iphi()
		  <<" dccNb="<<dccNb<<" tccNb="<<tccNb<<" towerInTCC="<<towerInTCC<<endl ;      
    }
  }

}

void EcalTPGParamBuilder::beginJob(const edm::EventSetup& evtSetup)
{
  using namespace edm;
  using namespace std;

  std::cout<<"we are in beginJob"<<endl;

  create_header() ; 
  std::cout<<"we are in beginJob after create header"<<endl;

  DetId eb(DetId::Ecal,EcalBarrel) ;
  DetId ee(DetId::Ecal,EcalEndcap) ;

  std::cout<<"we are in beginJob after detid"<<endl;

  if (writeToFiles_) {
    (*out_file_)<<"PHYSICS_EB "<<dec<<eb.rawId()<<std::endl ;
    (*out_file_)<<Et_sat_EB_<<" "<<TTF_lowThreshold_EB_<<" "<<TTF_highThreshold_EB_<<std::endl ;
    (*out_file_)<<FG_lowThreshold_EB_<<" "<<FG_highThreshold_EB_<<" "
		  <<FG_lowRatio_EB_<<" "<<FG_highRatio_EB_<<std::endl ;
    (*out_file_) <<std::endl ;

    (*out_file_)<<"PHYSICS_EE "<<dec<<ee.rawId()<<std::endl ;
    (*out_file_)<<Et_sat_EE_<<" "<<TTF_lowThreshold_EE_<<" "<<TTF_highThreshold_EE_<<std::endl ;
    (*out_file_)<<FG_Threshold_EE_<<" "<<-1<<" "
		  <<-1<<" "<<-1<<std::endl ;
    (*out_file_) <<std::endl ;
  }
  std::cout<<"we are in beginJob ending..."<<endl;


}



bool EcalTPGParamBuilder::computeLinearizerParam(double theta, double gainRatio, double calibCoeff, std::string subdet, int & mult , int & shift) 
{
  /*
    Linearization coefficient are determined in order to satisfy:
    tpg(ADC_sat) = 1024
    where: 
    tpg() is a model of the linearized tpg response on 10b 
    ADC_sat is the number of ADC count corresponding the Et_sat, the maximum scale of the transverse energy
    
    Since we have:
    Et_sat = xtal_LSB * ADC_sat * gainRatio * calibCoeff * sin(theta)
    and a simple model of tpg() being given by:
    tpg(X) = [ (X*mult) >> (shift+2) ] >> (sliding+shiftDet) 
    we must satisfy:
    [ (Et_sat/(xtal_LSB * gainRatio * calibCoeff * sin(theta)) * mult) >> (shift+2) ] >> (sliding+shiftDet) = 1024 
    that is:
    mult = 1024/Et_sat * xtal_LSB * gainRatio * calibCoeff * sin(theta) * 2^-(sliding+shiftDet+2) * 2^-shift
    mult = factor * 2^-shift
  */

  // case barrel:
  int shiftDet = 2 ;
  double ratio = xtal_LSB_EB_/Et_sat_EB_ ;
  // case endcap:
  if (subdet=="EE") {
    shiftDet = 0 ;
    ratio = xtal_LSB_EE_/Et_sat_EE_ ;
  }



  double factor = 1024 * ratio * gainRatio * calibCoeff * sin(theta) * (1 << (sliding_ + shiftDet + 2)) ;
  // Let's try first with shift = 0 (trivial solution)
  mult = (int)(factor+0.5) ; 
  for (shift = 0 ; shift<15 ; shift++) {
    if (mult>=128  && mult<256) return true ;
    factor *= 2 ; 
    mult = (int)(factor+0.5) ;
  }
  std::cout << "too bad we did not manage to calculate the factor for calib="<<calibCoeff<<endl;
  return false ;
}

void EcalTPGParamBuilder::create_header() 
{
  if (!writeToFiles_) return ;
  (*out_file_) <<"COMMENT put your comments here"<<std::endl ; 

  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT           physics EB structure"<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;
  (*out_file_) <<"COMMENT  EtSaturation (GeV), ttf_threshold_Low (GeV), ttf_threshold_High (GeV)"<<std::endl ;
  (*out_file_) <<"COMMENT  FG_lowThreshold (GeV), FG_highThreshold (GeV), FG_lowRatio, FG_highRatio"<<std::endl ;
  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;

  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT           physics EE structure"<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;
  (*out_file_) <<"COMMENT  EtSaturation (GeV), ttf_threshold_Low (GeV), ttf_threshold_High (GeV)"<<std::endl ;
  (*out_file_) <<"COMMENT  FG_Threshold (GeV), dummy, dummy, dummy"<<std::endl ;
  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;

  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT           crystal structure (same for EB and EE)"<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;
  (*out_file_) <<"COMMENT  ped, mult, shift [gain12]"<<std::endl ;
  (*out_file_) <<"COMMENT  ped, mult, shift [gain6]"<<std::endl ;
  (*out_file_) <<"COMMENT  ped, mult, shift [gain1]"<<std::endl ;
  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;

  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT           strip EB structure"<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;
  (*out_file_) <<"COMMENT  sliding_window"<<std::endl ;
  (*out_file_) <<"COMMENT  weightGroupId"<<std::endl ;
  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;

  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT           strip EE structure"<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;
  (*out_file_) <<"COMMENT  sliding_window"<<std::endl ;
  (*out_file_) <<"COMMENT  weightGroupId"<<std::endl ;
  (*out_file_) <<"COMMENT  threshold_fg strip_lut_fg"<<std::endl ;
  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;

  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT           tower EB structure"<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;
  (*out_file_) <<"COMMENT  LUTGroupId"<<std::endl ;
  (*out_file_) <<"COMMENT  FgGroupId"<<std::endl ;
  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;

  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT           tower EE structure"<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;
  (*out_file_) <<"COMMENT  LUTGroupId"<<std::endl ;
  (*out_file_) <<"COMMENT  tower_lut_fg"<<std::endl ;
  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;
  
  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT           Weight structure"<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;
  (*out_file_) <<"COMMENT  weightGroupId"<<std::endl ;
  (*out_file_) <<"COMMENT  w0, w1, w2, w3, w4"<<std::endl ;
  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;

  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT           lut structure"<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;
  (*out_file_) <<"COMMENT  LUTGroupId"<<std::endl ;
  (*out_file_) <<"COMMENT  LUT[1-1024]"<<std::endl ;
  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;

  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT           fg EB structure"<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;
  (*out_file_) <<"COMMENT  FgGroupId"<<std::endl ;
  (*out_file_) <<"COMMENT  el, eh, tl, th, lut_fg"<<std::endl ;
  (*out_file_) <<"COMMENT ================================="<<std::endl ;
  (*out_file_) <<"COMMENT"<<std::endl ;

  (*out_file_) <<std::endl ;
}


int EcalTPGParamBuilder::uncodeWeight(double weight, int complement2)
{
  int iweight ;
  uint max = (uint)(pow(2.,complement2)-1) ;
  if (weight>0) iweight=int((1<<6)*weight+0.5) ; // +0.5 for rounding pb
  else iweight= max - int(-weight*(1<<6)+0.5) +1 ;
  iweight = iweight & max ;
  return iweight ;
}

double EcalTPGParamBuilder::uncodeWeight(int iweight, int complement2)
{
  double weight = double(iweight)/pow(2., 6.) ;
  // test if negative weight:
  if ( (iweight & (1<<(complement2-1))) != 0) weight = (double(iweight)-pow(2., complement2))/pow(2., 6.) ;
  return weight ;
}

std::vector<unsigned int> EcalTPGParamBuilder::computeWeights(EcalShape & shape)
{
  double timeMax = shape.computeTimeOfMaximum() - shape.computeT0() ; // timeMax w.r.t begining of pulse
  double max = shape(timeMax) ;

  double sumf = 0. ;
  double sumf2 = 0. ;
  for (uint sample = 0 ; sample<nSample_ ; sample++) {
    double time = timeMax - ((double)sampleMax_-(double)sample)*25. ;
    sumf += shape(time)/max ;
    sumf2 += shape(time)/max * shape(time)/max ;
  }
  double lambda = 1./(sumf2-sumf*sumf/nSample_) ;
  double gamma = -lambda*sumf/nSample_ ;
  double * weight = new double[nSample_] ;
  for (uint sample = 0 ; sample<nSample_ ; sample++) {
    double time = timeMax - ((double)sampleMax_-(double)sample)*25. ;
    weight[sample] = lambda*shape(time)/max + gamma ;
  }

//   double ampl = 0. ;
//   for (uint sample = 0 ; sample<nSample_ ; sample++) {
//     double time = timeMax - ((double)sampleMax_-(double)sample)*25. ;
//     ampl += weight[sample]*shape(time) ;
//   }
//   std::cout<<max<<" "<<ampl<<std::endl ;


  int * iweight = new int[nSample_] ;
  for (uint sample = 0 ; sample<nSample_ ; sample++)   iweight[sample] = uncodeWeight(weight[sample], complement2_) ;

  // Let's check:  
  int isumw  = 0 ;  
  for (uint sample = 0 ; sample<nSample_ ; sample++) isumw  += iweight[sample] ;
  uint imax = (uint)(pow(2.,int(complement2_))-1) ;
  isumw = (isumw & imax ) ;

  // Let's correct for bias if any
  if (isumw != 0) {
    double min = 99. ;
    uint index = 0 ;
    if ( (isumw & (1<<(complement2_-1))) != 0) {
      // add 1:
      for (uint sample = 0 ; sample<nSample_ ; sample++) {
	int new_iweight = iweight[sample]+1 ; 
	double new_weight = uncodeWeight(new_iweight, complement2_) ;
	if (fabs(new_weight-weight[sample])<min) {
	  min = fabs(new_weight-weight[sample]) ;
	  index = sample ;
	}
      }
      iweight[index] ++ ; 
    } else {
      // Sub 1:
      for (uint sample = 0 ; sample<nSample_ ; sample++) {
        int new_iweight = iweight[sample]-1 ;    
	double new_weight = uncodeWeight(new_iweight, complement2_) ;
        if (fabs(new_weight-weight[sample])<min) {
          min = fabs(new_weight-weight[sample]) ;
          index = sample ;
        }
      }
      iweight[index] -- ; 
    } 
  }

  std::vector<unsigned int> theWeights ;
  for (uint sample = 0 ; sample<nSample_ ; sample++) theWeights.push_back(iweight[sample]) ;

  delete weight ;
  delete iweight ;
  return theWeights ;
}

void EcalTPGParamBuilder::computeLUT(int * lut, std::string det) 
{
  double Et_sat = Et_sat_EB_ ;
  double LUT_threshold = LUT_threshold_EB_ ;
  double LUT_stochastic = LUT_stochastic_EB_ ;
  double LUT_noise = LUT_noise_EB_ ;
  double LUT_constant = LUT_constant_EB_ ;
  double TTF_lowThreshold = TTF_lowThreshold_EB_ ;
  double TTF_highThreshold = TTF_highThreshold_EB_ ;
  if (det == "EE") {
    Et_sat = Et_sat_EE_ ;
    LUT_threshold = LUT_threshold_EE_ ;
    LUT_stochastic = LUT_stochastic_EE_ ;
    LUT_noise = LUT_noise_EE_ ;
    LUT_constant = LUT_constant_EE_ ;
    TTF_lowThreshold = TTF_lowThreshold_EE_ ;
    TTF_highThreshold = TTF_highThreshold_EE_ ;
  }

  // initialisation with identity
  for (int i=0 ; i<1024 ; i++) {
    lut[i] = i ;
    if (lut[i]>0xff) lut[i] = 0xff ;
  }

  // case linear LUT
  if (LUT_option_ == "Linear") {
    int mylut = 0 ;
    for (int i=0 ; i<1024 ; i++) {
      lut[i] = mylut ;
      if ((i+1)%4 == 0 ) mylut++ ;
    }
  }

  // case LUT following Ecal resolution
  if (LUT_option_ == "EcalResolution") {
    TF1 * func = new TF1("func",oneOverEtResolEt, 0., Et_sat,3) ;
    func->SetParameters(LUT_stochastic, LUT_noise, LUT_constant) ;
    double norm = func->Integral(0., Et_sat) ;
    for (int i=0 ; i<1024 ; i++) {   
      double Et = i*Et_sat/1024. ;
      lut[i] =  int(0xff*func->Integral(0., Et)/norm + 0.5) ;
    }
  }

  // Now, add TTF thresholds to LUT and apply LUT threshold if needed
  for (int j=0 ; j<1024 ; j++) {
    double Et_GeV = Et_sat/1024*(j+0.5) ;
    if (Et_GeV <= LUT_threshold) lut[j] = 0 ; // LUT threshold
    int ttf = 0x0 ;    
    if (Et_GeV >= TTF_highThreshold) ttf = 3 ;
    if (Et_GeV >= TTF_lowThreshold && Et_GeV < TTF_highThreshold) ttf = 1 ;
    ttf = ttf << 8 ;
    lut[j] += ttf ;
  }

}

void EcalTPGParamBuilder::getCoeff(coeffStruc & coeff, const EcalIntercalibConstantMap & calibMap, uint rawId)
{
  // get current intercalibration coeff
  coeff.calibCoeff_ = 1. ;
  EcalIntercalibConstantMap::const_iterator icalit = calibMap.find(rawId);
  if( icalit != calibMap.end() ) coeff.calibCoeff_ = (*icalit) ;
  else std::cout<<"getCoeff: "<<rawId<<" not found in EcalIntercalibConstantMap"<<std::endl ;
}

void EcalTPGParamBuilder::getCoeff(coeffStruc & coeff, const EcalGainRatioMap & gainMap, uint rawId)
{
  // get current gain ratio
  coeff.gainRatio_[0]  = 1. ;
  coeff.gainRatio_[1]  = 2. ;
  coeff.gainRatio_[2]  = 12. ;
  EcalGainRatioMap::const_iterator gainIter = gainMap.find(rawId);
  if (gainIter != gainMap.end()) {
    const EcalMGPAGainRatio & aGain = (*gainIter) ;
    coeff.gainRatio_[1] = aGain.gain12Over6() ;
    coeff.gainRatio_[2] = aGain.gain6Over1() * aGain.gain12Over6() ;
  }
  else std::cout<<"getCoeff: "<<rawId<<" not found in EcalGainRatioMap"<<std::endl ;
}

void EcalTPGParamBuilder::getCoeff(coeffStruc & coeff, const EcalPedestalsMap & pedMap, uint rawId)
{
  coeff.pedestals_[0] = 0 ;
  coeff.pedestals_[1] = 0 ;
  coeff.pedestals_[2] = 0 ;

  if (forcedPedestalValue_ >= 0) {
    coeff.pedestals_[0] = forcedPedestalValue_ ;
    coeff.pedestals_[1] = forcedPedestalValue_ ;
    coeff.pedestals_[2] = forcedPedestalValue_ ;  
    return ;
  }

  // get current pedestal
  EcalPedestalsMapIterator pedIter = pedMap.find(rawId);
  if (pedIter != pedMap.end()) {
    EcalPedestals::Item aped = (*pedIter);
    coeff.pedestals_[0] = int(aped.mean_x12 + 0.5) ; 
    coeff.pedestals_[1] = int(aped.mean_x6 + 0.5) ;
    coeff.pedestals_[2] = int(aped.mean_x1 + 0.5) ;
  }
  else std::cout<<"getCoeff: "<<rawId<<" not found in EcalPedestalsMap"<<std::endl ;
}

void EcalTPGParamBuilder::getCoeff(coeffStruc & coeff, const map<EcalLogicID, MonPedestalsDat> & pedMap, const EcalLogicID & logicId)
{
  // get current pedestal
  coeff.pedestals_[0] = 0 ;
  coeff.pedestals_[1] = 0 ;
  coeff.pedestals_[2] = 0 ;

  map<EcalLogicID, MonPedestalsDat>::const_iterator it =  pedMap.find(logicId);
  if (it != pedMap.end()) {
    MonPedestalsDat ped = it->second ;
    coeff.pedestals_[0] = int(ped.getPedMeanG12() + 0.5) ; 
    coeff.pedestals_[1] = int(ped.getPedMeanG6() + 0.5) ; 
    coeff.pedestals_[2] = int(ped.getPedMeanG1() + 0.5) ; 
  } 
  else std::cout<<"getCoeff: "<<logicId.getID1()<<", "<<logicId.getID2()<<", "<<logicId.getID3()
		<<" not found in map<EcalLogicID, MonPedestalsDat>"<<std::endl ;
}

void EcalTPGParamBuilder::computeFineGrainEBParameters(uint & lowRatio, uint & highRatio,
						       uint & lowThreshold, uint & highThreshold, uint & lut)
{
  lowRatio = int(0x80*FG_lowRatio_EB_ + 0.5) ;
  if (lowRatio>0x7f) lowRatio = 0x7f ;
  highRatio = int(0x80*FG_highRatio_EB_ + 0.5) ;
  if (highRatio>0x7f) highRatio = 0x7f ;
  
  // lsb at the stage of the FG calculation is:
  double lsb_FG = Et_sat_EB_/1024./4 ;
  lowThreshold = int(FG_lowThreshold_EB_/lsb_FG+0.5) ;
  if (lowThreshold>0xff) lowThreshold = 0xff ;
  highThreshold = int(FG_highThreshold_EB_/lsb_FG+0.5) ;
  if (highThreshold>0xff) highThreshold = 0xff ;

  // FG lut: FGVB response is LUT(adress) where adress is: 
  // bit3: maxof2/ET >= lowRatio, bit2: maxof2/ET >= highRatio, bit1: ET >= lowThreshold, bit0: ET >= highThreshold
  // FGVB =1 if jet-like (veto active), =0 if E.M.-like
  // the condition for jet-like is: ET>Threshold and  maxof2/ET < Ratio (only TT with enough energy are vetoed)

  // With the following lut, what matters is only max(TLow, Thigh) and max(Elow, Ehigh)
  // So, jet-like if maxof2/ettot<max(TLow, Thigh) && ettot >= max(Elow, Ehigh)
  if (FG_lut_EB_ == 0) lut = 0x0888 ; 
  else lut = FG_lut_EB_ ; // let's use the users value (hope he/she knows what he/she does!)
}

void EcalTPGParamBuilder::computeFineGrainEEParameters(uint & threshold, uint & lut_strip, uint & lut_tower) 
{
  // lsb for EE:
  double lsb_FG = Et_sat_EE_/1024. ; // FIXME is it true????
  threshold = int(FG_Threshold_EE_/lsb_FG+0.5) ;
  lut_strip = FG_lut_strip_EE_  ;
  lut_tower = FG_lut_tower_EE_  ;
}

