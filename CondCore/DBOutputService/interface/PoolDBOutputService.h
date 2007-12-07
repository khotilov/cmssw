#ifndef CondCore_PoolDBOutputService_h
#define CondCore_PoolDBOutputService_h
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
//#include "CondCore/DBCommon/interface/DBSession.h"
//#include "CondCore/IOVService/interface/IOVService.h"
#include "CondCore/DBCommon/interface/PoolTransaction.h"
#include "CondCore/DBCommon/interface/CoralTransaction.h"
#include "CondCore/DBCommon/interface/TypedRef.h"
#include "CondCore/DBCommon/interface/Time.h"
#include "CondCore/MetaDataService/interface/MetaData.h"
#include "serviceCallbackRecord.h"
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
      std::string tag( const std::string& EventSetupRecordName );
      bool isNewTagRequest( const std::string& EventSetupRecordName );
      //
      // insert the payload and its valid till time into the database
      // Note: user looses the ownership of the pointer to the payloadObj
      // The payload object will be stored as well
      // 
      template<typename T>
	void createNewIOV( T* firstPayloadObj, 
			   cond::Time_t firstTillTime,
			   const std::string& EventSetupRecordName
			   ){
	cond::service::serviceCallbackRecord& myrecord=this->lookUpRecord(EventSetupRecordName);
	if ( !m_dbstarted ) {
	  this->initDB();
	}
	if(!myrecord.m_isNewTag) throw cond::Exception("PoolDBOutputService::createNewIOV not a new tag");
	cond::PoolTransaction& pooldb=m_connection->poolTransaction();
	std::string iovToken;
	try{
	  pooldb.start(false);
	  cond::TypedRef<T> myPayload(pooldb,firstPayloadObj);
	  myPayload.markWrite(EventSetupRecordName);
	  std::string payloadToken=myPayload.token();
	  iovToken=this->insertIOV(pooldb, myrecord,payloadToken,firstTillTime);
	  pooldb.commit();
	  cond::CoralTransaction& coraldb=m_connection->coralTransaction();
	  cond::MetaData metadata(coraldb);
	  coraldb.start(false);
	  metadata.addMapping(myrecord.m_tag,iovToken);
	  coraldb.commit();
	  myrecord.m_isNewTag=false;
	  myrecord.m_iovtoken=iovToken;
	  m_newtags.push_back( std::make_pair<std::string,std::string>(myrecord.m_tag,iovToken) );
	}catch(const std::exception& er){
	  throw cond::Exception("PoolDBOutputService::createNewIOV "+std::string(er.what()));
	}
      }
      void createNewIOV( const std::string& firstPayloadToken, 
			cond::Time_t firstTillTime,
			const std::string& EventSetupRecordName );
      template<typename T>
	void appendTillTime( T* payloadObj, 
			     cond::Time_t tillTime,
			     const std::string& EventSetupRecordName
			     ){
	cond::service::serviceCallbackRecord& myrecord=this->lookUpRecord(EventSetupRecordName);
	if (!m_dbstarted) this->initDB();
	cond::PoolTransaction& pooldb=m_connection->poolTransaction();
	try{
	  pooldb.start(false);
	  cond::TypedRef<T> myPayload(pooldb,payloadObj);
	  myPayload.markWrite(EventSetupRecordName);
	  std::string payloadToken=myPayload.token();
	  std::string iovToken=this->insertIOV(pooldb,myrecord,payloadToken,tillTime);
	  pooldb.commit();
	}catch(const std::exception& er){
	  throw cond::Exception("PoolDBOutputService::appendTillTime "+std::string(er.what()));
	}
      }
      void appendTillTime( const std::string& payloadToken, 
			   cond::Time_t tillTime,
			   const std::string& EventSetupRecordName
			    );
      
      template<typename T>
	void appendSinceTime( T* payloadObj, 
			      cond::Time_t sinceTime,
			      const std::string& EventSetupRecordName ){
	cond::service::serviceCallbackRecord& myrecord=this->lookUpRecord(EventSetupRecordName);
	if (!m_dbstarted) {
	  this->initDB();
	}
	cond::PoolTransaction& pooldb=m_connection->poolTransaction();
	try{	  
	  pooldb.start(false);
	  cond::TypedRef<T> myPayload(pooldb,payloadObj);
	  myPayload.markWrite(EventSetupRecordName);
	  std::string payloadToken=myPayload.token();
	  this->appendIOV(pooldb,myrecord,payloadToken,sinceTime);
	  pooldb.commit();
	}catch(const std::exception& er){
	  throw cond::Exception("PoolDBOutputService::appendSinceTime "+std::string(er.what()));
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
			   const std::string& EventSetupRecordName );

      //
      // Service time utility callback method 
      // return the infinity value according to the given timetype
      // It is the IOV closing boundary
      //
      cond::Time_t endOfTime() const;
      //
      // Service time utility callback method 
      // return the current conditions time value according to the 
      // given timetype
      //
      cond::Time_t currentTime() const;
      virtual ~PoolDBOutputService();
    private:
      void connect();    
      void disconnect();
      void initDB();
      size_t callbackToken(const std::string& EventSetupRecordName ) const ;
      void appendIOV(cond::PoolTransaction&,
		     cond::service::serviceCallbackRecord& record,
		     const std::string& payloadToken, 
		     cond::Time_t sinceTime);
      std::string insertIOV(cond::PoolTransaction& pooldb,
			    cond::service::serviceCallbackRecord& record,
			    const std::string& payloadToken, 
			    cond::Time_t tillTime);
      //			    const std::string& EventSetupRecordName);
      serviceCallbackRecord& lookUpRecord(const std::string& EventSetupRecordName);
    private:
      std::string m_timetype; 
      cond::Time_t m_currentTime;
      cond::DBSession* m_session;
      cond::Connection* m_connection;
      std::map<size_t, cond::service::serviceCallbackRecord> m_callbacks;
      std::vector< std::pair<std::string,std::string> > m_newtags;
      bool m_dbstarted;
      //cond::IOVService* m_iovservice;
      //edm::ParameterSet m_connectionPset;
    };//PoolDBOutputService
  }//ns service
}//ns cond
#endif
