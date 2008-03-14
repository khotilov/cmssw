#include <stdexcept>
#include <string>
#include "OnlineDB/Oracle/interface/Oracle.h"

#include "OnlineDB/EcalCondDB/interface/ODEcalCycle.h"
#include "OnlineDB/EcalCondDB/interface/ODCCSCycle.h"
#include "OnlineDB/EcalCondDB/interface/ODDCCCycle.h"
#include "OnlineDB/EcalCondDB/interface/ODJBH4Cycle.h"
#include "OnlineDB/EcalCondDB/interface/ODLaserCycle.h"
#include "OnlineDB/EcalCondDB/interface/ODLTCCycle.h"
#include "OnlineDB/EcalCondDB/interface/ODLTSCycle.h"
#include "OnlineDB/EcalCondDB/interface/ODMataqCycle.h"
#include "OnlineDB/EcalCondDB/interface/ODRunConfigCycleInfo.h"
#include "OnlineDB/EcalCondDB/interface/ODScanCycle.h"
#include "OnlineDB/EcalCondDB/interface/ODTCCCycle.h"
#include "OnlineDB/EcalCondDB/interface/ODTTCciCycle.h"

using namespace std;
using namespace oracle::occi;

ODEcalCycle::ODEcalCycle()
{
  m_env = NULL;
  m_conn = NULL;
  m_writeStmt = NULL;
  m_readStmt = NULL;

   m_ID=0;
   clear();
}



ODEcalCycle::~ODEcalCycle()
{
}



void ODEcalCycle::prepareWrite()
  throw(runtime_error)
{

  std::cout<< "ODEcalCycle::prepareWrite(): this is a view writing specific tables  "<< endl;

}



void ODEcalCycle::writeDB()
  throw(runtime_error)
{

  ODRunConfigCycleInfo cyc;
  cyc.setTag(this->getCycleTag());
  cyc.setDescription(this->getCycleDescription());
  cyc.setCycleNumber(this->getCycleNum());
  cyc.setSequenceID(this->getSequenceId());
  // here we insert the cycle and get back the id
  cyc.setConnection(m_env, m_conn);
  cyc.insertConfig();
  int cyc_id=cyc.getId();
  cout << "Cycle inserted with ID "<< cyc_id<< endl;
  m_ID=cyc_id;

  if(getCCSId()!=0){ 
    ODCCSCycle ccs_cycle;
    ccs_cycle.setId(cyc_id);
    ccs_cycle.setCCSConfigurationID(getCCSId());
    ccs_cycle.setConnection(m_env, m_conn);
    ccs_cycle.insertConfig();
    cout << "Inserting CCS-cycle in DB..." << flush;
  } 
  if(getDCCId()!=0){
    ODDCCCycle dcc_cycle;
    dcc_cycle.setId(cyc_id);
    dcc_cycle.setDCCConfigurationID(getDCCId());
    dcc_cycle.setConnection(m_env, m_conn);
    dcc_cycle.insertConfig();
    cout << "Inserting DCC-cycle in DB..." << flush;
  }
  if(getLaserId()!=0){
    ODLaserCycle laser_cycle;
    laser_cycle.setId(cyc_id);
    laser_cycle.setLaserConfigurationID(getLaserId());
    laser_cycle.setConnection(m_env, m_conn);
    laser_cycle.insertConfig();
    cout << "Inserting LASER-cycle in DB..." << flush;
  }
  if(getLTCId()!=0){
    ODLTCCycle ltc_cycle;
    ltc_cycle.setId(cyc_id);
    ltc_cycle.setLTCConfigurationID(getLTCId());
    ltc_cycle.setConnection(m_env, m_conn);
    ltc_cycle.insertConfig();
    cout << "Inserting LTC-cycle in DB..." << flush;
  }
  if(getLTSId()!=0){
    ODLTSCycle lts_cycle;
    lts_cycle.setId(cyc_id);
    lts_cycle.setLTSConfigurationID(getLTSId());
    lts_cycle.setConnection(m_env, m_conn);
    lts_cycle.insertConfig();
    cout << "Inserting LTS-cycle in DB..." << flush;
  }
  if(getTCCId()!=0){
    ODTCCCycle tcc_cycle;
    tcc_cycle.setId(cyc_id);
    tcc_cycle.setTCCConfigurationID(getTCCId());
    tcc_cycle.setConnection(m_env, m_conn);
    tcc_cycle.insertConfig();
    cout << "Inserting TCC-cycle in DB..." << flush;
  }
  if(getTTCCIId()!=0){
    ODTTCciCycle ttcci_cycle;
    ttcci_cycle.setId(cyc_id);
    ttcci_cycle.setTTCciConfigurationID(getTTCCIId());
    ttcci_cycle.setConnection(m_env, m_conn);
    ttcci_cycle.insertConfig();
    cout << "Inserting TTCCI-cycle in DB..." << flush;
  }
  if(getMataqId()!=0){
    ODMataqCycle mataq_cycle;
    mataq_cycle.setId(cyc_id);
    mataq_cycle.setMataqConfigurationID(getMataqId());
    mataq_cycle.setConnection(m_env, m_conn);
    mataq_cycle.insertConfig();
    cout << "Inserting MATAQ-cycle in DB..." << flush;
  }
  if(getJBH4Id()!=0){
    ODJBH4Cycle jbh4_cycle;
    jbh4_cycle.setId(cyc_id);
    jbh4_cycle.setJBH4ConfigurationID(getJBH4Id());
    jbh4_cycle.setConnection(m_env, m_conn);
    jbh4_cycle.insertConfig();
    cout << "Inserting JBH4-cycle in DB..." << flush;
  }
  if(getScanId()!=0){
    ODScanCycle scan_cycle;
    scan_cycle.setId(cyc_id);
    scan_cycle.setScanConfigurationID(getScanId());
    scan_cycle.setConnection(m_env, m_conn);
    scan_cycle.insertConfig();
    cout << "Inserting SCAN-cycle in DB..." << flush;
  }



}


void ODEcalCycle::clear(){
   m_tag="";
   m_version=0;
   m_seq_num=0;
   m_cycle_num=0;
   m_cycle_tag="";
   m_cycle_description="";
   m_ccs=0;
   m_dcc=0;
   m_laser=0;
   m_ltc=0;
   m_lts=0;
   m_tcc=0;
   m_ttcci=0;
   m_mataq=0;
   m_jbh4=0;
   m_scan=0;
   m_seq_id=0;

}

void ODEcalCycle::printout(){

  std::cout << "**** ODEcalCycle ****"<< std::endl;
  std::cout << "**** id: "<< m_ID << std::endl;
  std::cout << "**** tag: "<< m_tag << std::endl;
  std::cout << "**** version: "<< m_version << std::endl;
  std::cout << "**** seq_num: "<< m_seq_num << std::endl;
  std::cout << "**** cycle num: "<< m_cycle_num << std::endl;
  std::cout << "**** cycle tag: "<< m_cycle_tag << std::endl;
  std::cout << "**** cycle description: "<< m_cycle_description<< std::endl;
  std::cout << "**** ccs_id: "<< m_ccs <<std::endl;
  std::cout << "**** dcc_id: "<< m_dcc <<std::endl;
  std::cout << "**** laser: "<< m_laser << std::endl;
  std::cout << "**** ODEcalCycle ****"<< std::endl;
 
  /*   m_ltc=0;
   m_lts=0;
   m_tcc=0;
   m_ttcci=0;
   m_mataq=0;
   m_jbh4=0;
   m_scan=0;
  */
}


void ODEcalCycle::fetchData(ODEcalCycle * result)
  throw(runtime_error)
{
  this->checkConnection();
  result->clear();
  if(result->getId()==0){
    throw(runtime_error("ODEcalCycle::fetchData(): no Id defined for this ODEcalCycle "));
  }

  try {

    m_readStmt->setSQL("SELECT d.tag, d.version, d.sequence_num, d.cycle_num, d.cycle_tag, d.description, d.ccs_configuration_id,    "
		       " d.dcc_configuration_id, d.laser_configuration_id, d.ltc_configuration_id, d.lts_configuration_id, d.tcc_configuration_id,"
		       " d.\"TTCci_CONFIGURATION_ID\" , d.matacq_configuration_id, d.jbh4_configuration_id, d.scan_id " 
		       "FROM ECAL_CYCLE d "
		       " where d.cycle_id = :1 " );
    m_readStmt->setInt(1, result->getId());
    ResultSet* rset = m_readStmt->executeQuery();

    rset->next();

    result->setTag(          rset->getString(1) );
    result->setVersion(      rset->getInt(2) );
    result->setSeqNum(       rset->getInt(3) );
    result->setCycleNum(     rset->getInt(4) );
    result->setCycleTag(     rset->getString(5) );
    result->setCycleDescription( rset->getString(6) );
    result->setCCSId(        rset->getInt(7) );
    result->setDCCId(        rset->getInt(8) );
    result->setLaserId(      rset->getInt(9) );
    result->setLTCId(        rset->getInt(10) );
    result->setLTSId(        rset->getInt(11) );
    result->setTCCId(        rset->getInt(12) );
    result->setTTCCIId(        rset->getInt(13) );
    result->setMataqId(        rset->getInt(14) );
    result->setJBH4Id(        rset->getInt(15) );
    result->setScanId(        rset->getInt(16) );


  } catch (SQLException &e) {
    throw(runtime_error("ODEcalCycle::fetchData():  "+e.getMessage()));
  }
}

