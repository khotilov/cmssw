#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondCore/DBCommon/interface/DbScopedTransaction.h"
#include "CondCore/DBCommon/interface/DbTransaction.h"
#include "CondCore/DBCommon/interface/TagInfo.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondCore/IOVService/interface/IOVService.h"
#include "CondCore/IOVService/interface/IOVEditor.h"
#include "CondCore/IOVService/interface/IOVIterator.h"
//#include "CondCore/IOVService/interface/IOVNames.h"
#include "CondCore/IOVService/interface/IOVSchemaUtility.h"
#include "CondCore/DBCommon/interface/Exception.h"
#include "CondCore/DBOutputService/interface/Exception.h"
//#include "FWCore/Framework/interface/IOVSyncValue.h"

//POOL include
//#include "FileCatalog/IFileCatalog.h"
#include "serviceCallbackToken.h"
#include "CondCore/DBCommon/interface/UserLogInfo.h"
#include "CondCore/DBCommon/interface/IOVInfo.h"

//#include <iostream>
#include <vector>
#include<memory>

namespace {
  std::string dsw("DataWrapper");
}

unsigned int cond::service::GetToken::sizeDSW() {
  return dsw.size();
}

void 
cond::service::PoolDBOutputService::fillRecord( edm::ParameterSet & pset) {
  cond::service::serviceCallbackRecord thisrecord;

  thisrecord.m_idName = pset.getParameter<std::string>("record");
  thisrecord.m_tag = pset.getParameter<std::string>("tag");
  
  thisrecord.m_withWrapper =  
    pset.getUntrackedParameter<bool>("withWrapper", m_withWrapper);
 
  thisrecord.m_freeInsert = 
    pset.getUntrackedParameter<bool>("outOfOrder",m_freeInsert);
  
  thisrecord.m_timetype=cond::findSpecs(pset.getUntrackedParameter< std::string >("timetype",m_timetypestr)).type;

  m_callbacks.insert(std::make_pair(cond::service::serviceCallbackToken::build(thisrecord.m_idName),thisrecord));
 
  if(m_logdbOn){
      cond::UserLogInfo userloginfo;
      m_logheaders.insert(std::make_pair(cond::service::serviceCallbackToken::build(thisrecord.m_idName),userloginfo));
  }
 

}


cond::service::PoolDBOutputService::PoolDBOutputService(const edm::ParameterSet & iConfig,edm::ActivityRegistry & iAR ): 
  m_currentTime( 0 ),
  m_connection(),
  m_session(),
  m_logSession(),
  m_dbstarted( false ),
  m_logdb( 0 ),
  m_logdbOn( false ),
  m_freeInsert(false),
  m_withWrapper(false)
{

  if( iConfig.exists("withWrapper") ){
     m_withWrapper=iConfig.getUntrackedParameter<bool>("withWrapper");
  }  

  if( iConfig.exists("outOfOrder") ){
     m_freeInsert=iConfig.getUntrackedParameter<bool>("outOfOrder");
  }  

  m_timetypestr=iConfig.getUntrackedParameter< std::string >("timetype","runnumber");
  m_timetype=cond::findSpecs( m_timetypestr).type;

  std::string connect=iConfig.getParameter<std::string>("connect");
  std::string logconnect("");
  if( iConfig.exists("logconnect") ){
    logconnect=iConfig.getUntrackedParameter<std::string>("logconnect");
  }  

  edm::ParameterSet connectionPset = iConfig.getParameter<edm::ParameterSet>("DBParameters");
  m_connection.configuration().setParameters( connectionPset );
  m_connection.configure();
  
  m_session = m_connection.createSession();
  std::string blobstreamerName("");
  if( iConfig.exists("BlobStreamerName") ){
    blobstreamerName=iConfig.getUntrackedParameter<std::string>("BlobStreamerName");
    blobstreamerName.insert(0,"COND/Services/");
    m_session.setBlobStreamingService(blobstreamerName);
  }
  m_session.open( connect );
  
  if( !logconnect.empty() ){
    m_logdbOn=true;
    m_logSession = m_connection.createSession();
    m_logSession.open( logconnect );
  }

  typedef std::vector< edm::ParameterSet > Parameters;
  Parameters toPut=iConfig.getParameter<Parameters>("toPut");
  for(Parameters::iterator itToPut = toPut.begin(); itToPut != toPut.end(); ++itToPut)
    fillRecord( *itToPut);


  iAR.watchPreProcessEvent(this,&cond::service::PoolDBOutputService::preEventProcessing);
  iAR.watchPostEndJob(this,&cond::service::PoolDBOutputService::postEndJob);
  iAR.watchPreModule(this,&cond::service::PoolDBOutputService::preModule);
  iAR.watchPostModule(this,&cond::service::PoolDBOutputService::postModule);
  iAR.watchPreBeginLumi(this,&cond::service::PoolDBOutputService::preBeginLumi);
}

cond::DbSession
cond::service::PoolDBOutputService::session() const{
  return m_session;
}

std::string 
cond::service::PoolDBOutputService::tag( const std::string& EventSetupRecordName ){
  return this->lookUpRecord(EventSetupRecordName).m_tag;
}

bool 
cond::service::PoolDBOutputService::isNewTagRequest( const std::string& EventSetupRecordName ){
  cond::service::serviceCallbackRecord& myrecord=this->lookUpRecord(EventSetupRecordName);
  if(!m_dbstarted) this->initDB();
  return myrecord.m_isNewTag;
}


void 
cond::service::PoolDBOutputService::initDB()
{
  if(m_dbstarted) return;
  try{
    cond::DbScopedTransaction transaction( m_session );
    transaction.start(false);
    IOVSchemaUtility ut(m_session);
    ut.create();
    cond::MetaData metadata(m_session);
    for(std::map<size_t,cond::service::serviceCallbackRecord>::iterator it=m_callbacks.begin(); it!=m_callbacks.end(); ++it){
      //std::string iovtoken;
      if( !metadata.hasTag(it->second.m_tag) ){
        it->second.m_iovtoken="";
        it->second.m_isNewTag=true;
      }else{
        it->second.m_iovtoken=metadata.getToken(it->second.m_tag);
        it->second.m_isNewTag=false;
      }
    }
    transaction.commit();    
    //init logdb if required
    if(m_logdbOn){
      m_logdb=new cond::Logger(m_logSession);
      //m_logdb->getWriteLock();
      m_logdb->createLogDBIfNonExist();
      //m_logdb->releaseWriteLock();
    }
  }catch( const std::exception& er ){
    throw cond::Exception( "PoolDBOutputService::initDB "+std::string(er.what()) );
  }
  m_dbstarted=true;
}

void 
cond::service::PoolDBOutputService::postEndJob()
{
  if(m_logdb){
    delete m_logdb;
  }
}

void 
cond::service::PoolDBOutputService::preEventProcessing(const edm::EventID& iEvtid, const edm::Timestamp& iTime)
{
  if( m_timetype == cond::runnumber ){//runnumber
    m_currentTime=iEvtid.run();
  }else if( m_timetype == cond::timestamp ){ //timestamp
    m_currentTime=iTime.value();
  }
}

void
cond::service::PoolDBOutputService::preModule(const edm::ModuleDescription& desc){
}

void 
cond::service::PoolDBOutputService::preBeginLumi(const edm::LuminosityBlockID& iLumiid,  const edm::Timestamp& iTime ){
  if( m_timetype == cond::lumiid ){
    m_currentTime=iLumiid.value();
  }
}

void
cond::service::PoolDBOutputService::postModule(const edm::ModuleDescription& desc){
}

cond::service::PoolDBOutputService::~PoolDBOutputService(){
}

size_t 
cond::service::PoolDBOutputService::callbackToken(const std::string& EventSetupRecordName ) const {
  return cond::service::serviceCallbackToken::build(EventSetupRecordName);
}

cond::Time_t 
cond::service::PoolDBOutputService::endOfTime() const{
  return timeTypeSpecs[m_timetype].endValue;
}

cond::Time_t 
cond::service::PoolDBOutputService::beginOfTime() const{
  return timeTypeSpecs[m_timetype].beginValue;
}

cond::Time_t 
cond::service::PoolDBOutputService::currentTime() const{
  return m_currentTime;
}

void 
cond::service::PoolDBOutputService::createNewIOV( GetToken const & payloadToken, cond::Time_t firstSinceTime, cond::Time_t firstTillTime,const std::string& EventSetupRecordName, bool withlogging){
  cond::service::serviceCallbackRecord& myrecord=this->lookUpRecord(EventSetupRecordName);
  if (!m_dbstarted) this->initDB();
  if(!myrecord.m_isNewTag) throw cond::Exception("PoolDBOutputService::createNewIOV not a new tag");
  std::string iovToken;
  if(withlogging){
    if(!m_logdb) throw cond::Exception("cannot log to non-existing log db");
    m_logdb->getWriteLock();
  }
 
  std::string objToken;
  unsigned int payloadIdx=0;
  try{
    cond::DbScopedTransaction transaction(m_session);
    transaction.start(false);
    
    cond::IOVService iovmanager(m_session);
    std::auto_ptr<cond::IOVEditor> editor(iovmanager.newIOVEditor(""));
    editor->create(myrecord.m_timetype, firstTillTime);
    objToken = payloadToken(m_session,myrecord.m_withWrapper);
    unsigned int payloadIdx=editor->append(firstSinceTime, objToken);
    iovToken=editor->token();
    editor->stamp(cond::userInfo(),false);
    
    cond::MetaData metadata(m_session);

    /*
    MetaDataEntry imetadata;
    imetadata.tagname=myrecord.m_tag;
    imetadata.iovtoken=iovToken;
    imetadata.timetype=m_timetype;
    imetadata.firstsince=firstSinceTime;
    metadata.addMapping(imetadata);
   */
    metadata.addMapping(myrecord.m_tag,iovToken,myrecord.m_timetype);
    transaction.commit();

    m_newtags.push_back( std::make_pair<std::string,std::string>(myrecord.m_tag,iovToken) );
    myrecord.m_iovtoken=iovToken;
    myrecord.m_isNewTag=false;
    if(withlogging){
      std::string destconnect=m_session.connectionString();
      cond::UserLogInfo a=this->lookUpUserLogInfo(EventSetupRecordName);
      m_logdb->logOperationNow(a,destconnect,objToken,myrecord.m_tag,myrecord.timetypestr(),payloadIdx);
    }
  }catch(const std::exception& er){ 
    if(withlogging){
      std::string destconnect=m_session.connectionString();
      cond::UserLogInfo a=this->lookUpUserLogInfo(EventSetupRecordName);
      m_logdb->logFailedOperationNow(a,destconnect,objToken,myrecord.m_tag,myrecord.timetypestr(),payloadIdx,std::string(er.what()));
      m_logdb->releaseWriteLock();
    }
    throw cond::Exception("PoolDBOutputService::createNewIOV "+std::string(er.what()));
  }
  if(withlogging){
    m_logdb->releaseWriteLock();
  }
}


void 
cond::service::PoolDBOutputService::add( bool sinceNotTill, 
					 GetToken const & payloadToken,  
					 cond::Time_t time,
					 const std::string& EventSetupRecordName,
					 bool withlogging) {
  cond::service::serviceCallbackRecord& myrecord=this->lookUpRecord(EventSetupRecordName);
  if (!m_dbstarted) this->initDB();
  if(withlogging){
    if(!m_logdb) throw cond::Exception("cannot log to non-existing log db");
    m_logdb->getWriteLock();
  }

  std::string objToken;
  unsigned int payloadIdx=0;

  try{
    cond::DbScopedTransaction transaction(m_session);
    transaction.start(false);
    objToken = payloadToken(m_session,myrecord.m_withWrapper);
    payloadIdx= sinceNotTill ?
      this->appendIOV(m_session,myrecord,objToken,time) :
      this->insertIOV(m_session,myrecord,objToken,time);
    transaction.commit();
    if(withlogging){
      std::string destconnect=m_session.connectionString();
      cond::UserLogInfo a=this->lookUpUserLogInfo(EventSetupRecordName);
      m_logdb->logOperationNow(a,destconnect,objToken,myrecord.m_tag,myrecord.timetypestr(),payloadIdx);
    }
  }catch(const std::exception& er){
    if(withlogging){
      std::string destconnect=m_session.connectionString();
      cond::UserLogInfo a=this->lookUpUserLogInfo(EventSetupRecordName);
      m_logdb->logFailedOperationNow(a,destconnect,objToken,myrecord.m_tag,myrecord.timetypestr(),payloadIdx,std::string(er.what()));
      m_logdb->releaseWriteLock();
    }
    throw cond::Exception("PoolDBOutputService::add "+std::string(er.what()));
  }
  if(withlogging){
    m_logdb->releaseWriteLock();
  }
}

cond::service::serviceCallbackRecord& 
cond::service::PoolDBOutputService::lookUpRecord(const std::string& EventSetupRecordName){
  size_t callbackToken=this->callbackToken( EventSetupRecordName );
  std::map<size_t,cond::service::serviceCallbackRecord>::iterator it=m_callbacks.find(callbackToken);
  if(it==m_callbacks.end()) throw cond::UnregisteredRecordException(EventSetupRecordName);
  return it->second;
}

cond::UserLogInfo& 
cond::service::PoolDBOutputService::lookUpUserLogInfo(const std::string& EventSetupRecordName){
  size_t callbackToken=this->callbackToken( EventSetupRecordName );
  std::map<size_t,cond::UserLogInfo>::iterator it=m_logheaders.find(callbackToken);
  if(it==m_logheaders.end()) throw cond::UnregisteredRecordException(EventSetupRecordName);
  return it->second;
}


unsigned int 
cond::service::PoolDBOutputService::appendIOV(cond::DbSession& pooldb,
						   cond::service::serviceCallbackRecord& record, 
						   const std::string& payloadToken, 
						   cond::Time_t sinceTime){
  if( record.m_isNewTag ) {
    throw cond::Exception(std::string("PoolDBOutputService::appendIOV: cannot append to non-existing tag ")+record.m_tag );  
  }

  cond::IOVService iovmanager(pooldb);
  std::auto_ptr<cond::IOVEditor> editor(iovmanager.newIOVEditor(record.m_iovtoken));
 
  unsigned int payloadIdx =  record.m_freeInsert ? 
    editor->freeInsert(sinceTime,payloadToken) :
    editor->append(sinceTime,payloadToken);
  editor->stamp(cond::userInfo(),false);

  return payloadIdx;
}

unsigned int
cond::service::PoolDBOutputService::insertIOV( cond::DbSession& pooldb,
					       cond::service::serviceCallbackRecord& record, 
					       const std::string& payloadToken,
					       cond::Time_t tillTime){
  
  if( record.m_isNewTag ) {
    throw cond::Exception(std::string("PoolDBOutputService::insertIOV: cannot append to non-existing tag ")+record.m_tag );  
  }
  
  cond::IOVService iovmanager(pooldb);
  std::auto_ptr<cond::IOVEditor> editor(iovmanager.newIOVEditor(record.m_iovtoken));
  unsigned int payloadIdx=editor->insert(tillTime,payloadToken);
  editor->stamp(cond::userInfo(),false);

  delete editor;    
  return payloadIdx;
}



void
cond::service::PoolDBOutputService::setLogHeaderForRecord(const std::string& EventSetupRecordName,const std::string& dataprovenance,const std::string& usertext)
{
  cond::UserLogInfo& myloginfo=this->lookUpUserLogInfo(EventSetupRecordName);
  myloginfo.provenance=dataprovenance;
  myloginfo.usertext=usertext;
}
const cond::Logger& 
cond::service::PoolDBOutputService::queryLog()const{
  if(!m_logdb) throw cond::Exception("PoolDBOutputService::queryLog ERROR: logging is off");
  return *m_logdb;
}


void 
cond::service::PoolDBOutputService::tagInfo(const std::string& EventSetupRecordName,cond::TagInfo& result ){
  cond::service::serviceCallbackRecord& record=this->lookUpRecord(EventSetupRecordName);
  result.name=record.m_tag;
  result.token=record.m_iovtoken;
  //use ioviterator to find out.
  cond::DbScopedTransaction transaction(m_session);
  transaction.start(true);
  cond::IOVService iovmanager(m_session);
  std::auto_ptr<cond::IOVIterator> iit(iovmanager.newIOVIterator(result.token,cond::IOVService::backwardIter));
  iit->next(); // just to initialize
  result.lastInterval=iit->validity();
  result.lastPayloadToken=iit->payloadToken();
  result.size=iit->size();
 }
