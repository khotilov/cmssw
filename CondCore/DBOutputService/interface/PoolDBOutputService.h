#ifndef CondCore_PoolDBOutputService_h
#define CondCore_PoolDBOutputService_h
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "CondCore/DBCommon/interface/PoolTransaction.h"
#include "CondCore/DBCommon/interface/CoralTransaction.h"
#include "CondCore/DBCommon/interface/TypedRef.h"
#include "CondCore/DBCommon/interface/Time.h"
#include "CondCore/MetaDataService/interface/MetaData.h"
#include "serviceCallbackRecord.h"
#include "Logger.h"
#include "UserLogInfo.h"
#include "TagInfo.h"
#include <string>
#include <map>
//#include <iostream>
//
// Package:     DBOutputService
// Class  :     PoolDBOutputService
// 
/**\class PoolDBOutputService PoolDBOutputService.h CondCore/DBOutputService/interface/PoolDBOutputService.h
   Description: edm service for writing conditions object to DB.  
*/
//
// Author:      Zhen Xie
//
namespace edm{
  class Event;
  class EventSetup;
  class ParameterSet;
}
namespace cond{
  namespace service {
    class serviceCallbackToken;
    /** transaction and data consistency
	for create new tag, 
	start write metadata transaction only if the first pool commit 
	successful;
	for append,start readonly metadata transaction. start pool transaction only if metadata transaction successful.
	
    */
    class PoolDBOutputService{
    public:
      PoolDBOutputService( const edm::ParameterSet & iConfig, 
			   edm::ActivityRegistry & iAR );
      //use these to control connections
      //void  postBeginJob();
      void  postEndJob();
      //
      //use these to control transaction interval
      //
      void preEventProcessing( const edm::EventID & evtID, 
      			       const edm::Timestamp & iTime );
      void preModule(const edm::ModuleDescription& desc);
      void postModule(const edm::ModuleDescription& desc);
      //
      // return the database session in use
      //
      cond::DBSession& session() const;
      //
      // return the database connection handle in use
      //
      cond::Connection& connection() const;
      // 
      std::string tag( const std::string& EventSetupRecordName );
      bool isNewTagRequest( const std::string& EventSetupRecordName );
      const cond::Logger& queryLog() const;
      //
      // insert the payload and its valid till time into the database
      // Note: user looses the ownership of the pointer to the payloadObj
      // The payload object will be stored as well
      // 
      template<typename T>
	void createNewIOV( T* firstPayloadObj, 
			   cond::Time_t firstSinceTime,
			   cond::Time_t firstTillTime,
			   const std::string& EventSetupRecordName,
			   bool withlogging=false){
	cond::service::serviceCallbackRecord& myrecord=this->lookUpRecord(EventSetupRecordName);
	if ( !m_dbstarted ) {
	  this->initDB();
	}
	if(!myrecord.m_isNewTag) throw cond::Exception("PoolDBOutputService::createNewIOV not a new tag");
	//aquire writer lock
	cond::PoolTransaction& pooldb=m_connection->poolTransaction();
	std::string iovToken;
	std::string payloadToken("");
	if(withlogging){
	  m_logdb->getWriteLock();
	}
	std::pair<unsigned int,std::string> result;
	try{
	  pooldb.start(false);
	  cond::TypedRef<T> myPayload(pooldb,firstPayloadObj);
	  myPayload.markWrite(EventSetupRecordName);
	  payloadToken=myPayload.token();
	  result=this->insertIOV(pooldb, myrecord,payloadToken,firstSinceTime,firstTillTime);
	  iovToken=result.second;
	  pooldb.commit();
	  cond::CoralTransaction& coraldb=m_connection->coralTransaction();
	  cond::MetaData metadata(coraldb);
	  coraldb.start(false);
	  MetaDataEntry imetadata;
	  imetadata.tagname=myrecord.m_tag;
	  imetadata.iovtoken=iovToken;
	  imetadata.timetype=m_timetype;
	  imetadata.firstsince=firstSinceTime;
	  //metadata.addMapping(myrecord.m_tag,iovToken,m_timetype);
	  metadata.addMapping(imetadata);
	  coraldb.commit();
	  myrecord.m_isNewTag=false;
	  myrecord.m_iovtoken=iovToken;
	  m_newtags.push_back( std::make_pair<std::string,std::string>(myrecord.m_tag,iovToken) );	  
	  if(withlogging){
	    if(!m_logdb)throw cond::Exception("cannot log to non-existing log db");
	    std::string destconnect=m_connection->connectStr();
	    cond::service::UserLogInfo a=this->lookUpUserLogInfo(EventSetupRecordName);
	    m_logdb->logOperationNow(a,destconnect,payloadToken,myrecord.m_tag,m_timetypestr,result.first);
	  }
	}catch(const std::exception& er){
	  if(withlogging){
	    std::string destconnect=m_connection->connectStr();
	    cond::service::UserLogInfo a=this->lookUpUserLogInfo(EventSetupRecordName);
	    m_logdb->logFailedOperationNow(a,destconnect,payloadToken,myrecord.m_tag,m_timetypestr,result.first,std::string(er.what()));
	    m_logdb->releaseWriteLock();
	  }
	  throw cond::Exception("PoolDBOutputService::createNewIOV "+std::string(er.what()));
	}
	if(withlogging){
	  m_logdb->releaseWriteLock();
	}
      }
      void createNewIOV( const std::string& firstPayloadToken, 
			 cond::Time_t firstSinceTime, 
			 cond::Time_t firstTillTime,
			 const std::string& EventSetupRecordName,
			 bool withlogging=false);
      /*template<typename T>
	void appendTillTime( T* payloadObj, 
			     cond::Time_t tillTime,
			     const std::string& EventSetupRecordName,
			     bool withlogging=false){
	cond::service::serviceCallbackRecord& myrecord=this->lookUpRecord(EventSetupRecordName);
	if (!m_dbstarted) this->initDB();
	cond::PoolTransaction& pooldb=m_connection->poolTransaction();
	if(withlogging){
	  m_logdb->getWriteLock();
	}
	std::string payloadToken("");
	std::pair<unsigned int, std::string> result;
	try{
	  pooldb.start(false);
	  cond::TypedRef<T> myPayload(pooldb,payloadObj);
	  myPayload.markWrite(EventSetupRecordName);
	  payloadToken=myPayload.token();
	  result=this->insertIOV(pooldb,myrecord,payloadToken,tillTime);
	  std::string iovToken=result.second;

	  pooldb.commit();
	  if(withlogging){
	    if(!m_logdb)throw cond::Exception("cannot log to non-existing log db");
	    std::string destconnect=m_connection->connectStr();
	    cond::service::UserLogInfo a=this->lookUpUserLogInfo(EventSetupRecordName);
	    m_logdb->logOperationNow(a,destconnect,payloadToken,myrecord.m_tag,m_timetypestr,result.first);
	  }
	}catch(const std::exception& er){
	  if(withlogging){
	    std::string destconnect=m_connection->connectStr();
	    cond::service::UserLogInfo a=this->lookUpUserLogInfo(EventSetupRecordName);
	    m_logdb->logFailedOperationNow(a,destconnect,payloadToken,myrecord.m_tag,m_timetypestr,result.first,std::string(er.what()));
	    m_logdb->releaseWriteLock();
	  }
	  throw cond::Exception("PoolDBOutputService::appendTillTime "+std::string(er.what()));
	}
	if(withlogging){
	  m_logdb->releaseWriteLock();
	}
      }
      void appendTillTime( const std::string& payloadToken, 
			   cond::Time_t tillTime,
			   const std::string& EventSetupRecordName,
			   bool withlogging=false
			   );
      */
      template<typename T>
	void appendSinceTime( T* payloadObj, 
			      cond::Time_t sinceTime,
			      const std::string& EventSetupRecordName,
			      bool withlogging=false){
	cond::service::serviceCallbackRecord& myrecord=this->lookUpRecord(EventSetupRecordName);
	if (!m_dbstarted) {
	  this->initDB();
	}
	cond::PoolTransaction& pooldb=m_connection->poolTransaction();
	if(withlogging){
	  m_logdb->getWriteLock();
	}
	std::string payloadToken("");
	unsigned int payloadIdx=0;
	try{	  
	  pooldb.start(false);
	  cond::TypedRef<T> myPayload(pooldb,payloadObj);
	  myPayload.markWrite(EventSetupRecordName);
	  payloadToken=myPayload.token();
	  payloadIdx=this->appendIOV(pooldb,myrecord,payloadToken,sinceTime);
	  pooldb.commit();
	  if(withlogging){
	    if(!m_logdb)throw cond::Exception("cannot log to non-existing log db");
	    std::string destconnect=m_connection->connectStr();
	    cond::service::UserLogInfo a=this->lookUpUserLogInfo(EventSetupRecordName);
	    m_logdb->logOperationNow(a,destconnect,payloadToken,myrecord.m_tag,m_timetypestr,payloadIdx);
	  }
	}catch(const std::exception& er){
	  if(withlogging){
	    std::string destconnect=m_connection->connectStr();
	    cond::service::UserLogInfo a=this->lookUpUserLogInfo(EventSetupRecordName);
	    m_logdb->logFailedOperationNow(a,destconnect,payloadToken,myrecord.m_tag,m_timetypestr,payloadIdx,std::string(er.what()));
	    m_logdb->releaseWriteLock();
	  }
	  throw cond::Exception("PoolDBOutputService::appendSinceTime "+std::string(er.what()));
	}
	if(withlogging){
	  m_logdb->releaseWriteLock();
	}
      }
      //
      // Append the payload and its valid sinceTime into the database
      // Note: user looses the ownership of the pointer to the payloadObj
      // Note: the iov index appended to MUST pre-existing and the existing 
      // conditions data are retrieved from EventSetup 
      // 
      void appendSinceTime( const std::string& payloadToken, 
			   cond::Time_t sinceTime,
			   const std::string& EventSetupRecordName,
			   bool withlogging=false);

      //
      // Service time utility callback method 
      // return the infinity value according to the given timetype
      //
      cond::Time_t endOfTime() const;
      //
      // Service time utility callback method 
      // return beginning of time value according to the given timetype
      //
      cond::Time_t beginOfTime() const;
      //
      // Service time utility callback method 
      // return the current conditions time value according to the 
      // given timetype
      //
      cond::Time_t currentTime() const;
      // optional. User can inject additional information into the log associated with a given record
      void setLogHeaderForRecord(const std::string& EventSetupRecordName,
			   const std::string& provenance,
			   const std::string& usertext);
      // 
      // Retrieve tag information of the data
      // 
      void tagInfo(const std::string& EventSetupRecordName,
		   cond::TagInfo& result );
      virtual ~PoolDBOutputService();  
    private:
      void connect();    
      void disconnect();
      void initDB();
      size_t callbackToken(const std::string& EventSetupRecordName ) const ;
      unsigned int appendIOV(cond::PoolTransaction&,
			     cond::service::serviceCallbackRecord& record,
			     const std::string& payloadToken, 
			     cond::Time_t sinceTime);
      /// Returns payload location index and iov token 
      std::pair<unsigned int,std::string> 
	insertIOV(cond::PoolTransaction& pooldb,
		  cond::service::serviceCallbackRecord& record,
		  const std::string& payloadToken, 
		  cond::Time_t firstSinceTime,				    
		  cond::Time_t tillTime);
      //			    const std::string& EventSetupRecordName);
      serviceCallbackRecord& lookUpRecord(const std::string& EventSetupRecordName);
      UserLogInfo& lookUpUserLogInfo(const std::string& EventSetupRecordName);
      
    private:
      cond::TimeType m_timetype; 
      std::string m_timetypestr;
      cond::Time_t m_currentTime;
      cond::DBSession* m_session;
      cond::Connection* m_connection;
      std::map<size_t, cond::service::serviceCallbackRecord> m_callbacks;
      std::vector< std::pair<std::string,std::string> > m_newtags;
      bool m_dbstarted;
      cond::Logger* m_logdb;
      bool m_logdbOn;
      std::map<size_t, cond::service::UserLogInfo> m_logheaders;
      //cond::IOVService* m_iovservice;
      //edm::ParameterSet m_connectionPset;
    };//PoolDBOutputService
  }//ns service
}//ns cond
#endif
