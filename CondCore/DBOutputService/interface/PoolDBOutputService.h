#ifndef CondCore_PoolDBOutputService_h
#define CondCore_PoolDBOutputService_h
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "CondCore/DBCommon/interface/DbConnection.h"
#include "CondCore/DBCommon/interface/Time.h"
#include "CondCore/DBCommon/interface/ClassInfoLoader.h"
#include "CondCore/MetaDataService/interface/MetaData.h"
#include "CondCore/DBCommon/interface/Logger.h"
#include "CondCore/DBCommon/interface/UserLogInfo.h"
#include "CondCore/DBCommon/interface/TagInfo.h"
#include <string>
#include <map>
#include "CondFormats/Common/interface/PayloadWrapper.h"

// many many clients do not include explicitely!
#ifndef COND_EXCEPTION_H
#include "CondCore/DBCommon/interface/Exception.h"
// #warning please include  "CondCore/DBCommon/interface/Exception.h" explicitely
// #define COND_EXP_WARNING
#endif

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

    /** transaction and data consistency
	for create new tag, 
	start write metadata transaction only if the first pool commit 
	successful;
	for append,start readonly metadata transaction. start pool transaction only if metadata transaction successful.
	
	
    */
    
    
    struct GetToken {
      virtual std::string operator()(cond::DbSession&, bool) const =0;
      static unsigned int sizeDSW();
    };
    
    struct GetTrivialToken : public GetToken {
      
      GetTrivialToken(std::string token) :
	m_token(token){}
      virtual ~GetTrivialToken(){}
      virtual std::string operator()(cond::DbSession&, bool) const {
	return m_token;
      }
      std::string m_token;
    };
    
    template<typename T>
    struct GetTokenFromPointer : public GetToken {
      typedef cond::DataWrapper<T> Wrapper;
      
      GetTokenFromPointer(T * p, Summary * s) :
	m_p(p), m_s(s){}
      //m_w(new Wrapper(p,s)){}
      
      virtual std::string operator()(cond::DbSession& pooldb, bool withWrapper) const {
	std::string ret("");
	if (withWrapper) {
	  Wrapper* wrapper = new Wrapper(m_p,m_s);
	  std::string className = cond::classNameForPointer( wrapper );
	  pool::Ref<Wrapper> myPayload = pooldb.storeObject(wrapper,className.replace(0,sizeDSW(),"DSW"));
	  ret = myPayload.toString();
	} else {
	  std::string className = cond::classNameForPointer( m_p );
	  pool::Ref<T> myPayload = pooldb.storeObject(m_p,className);
	  ret = myPayload.toString();
	}
	return ret;
      }
      
      T * m_p;
      cond::Summary * m_s;
      // Wrapper * m_w;
      //const std::string& m_recordName;
    };
    
    
    
    
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
      void preBeginLumi(const edm::LuminosityBlockID&, 
			const edm::Timestamp& );
      //
      // return the database session in use
      //
      cond::DbSession session() const;
      //
      std::string tag( const std::string& recordName );
      bool isNewTagRequest( const std::string& recordName );
      const cond::Logger& queryLog() const;
      
      // BW-compatible signature
      template<typename T>
      void writeOne(T * payload,
		    Time_t time, const std::string& recordName, 
                    bool withlogging=false) {
        this->writeOne<T>(payload, 0, time, recordName, withlogging);
      }

      /* write one (either create or append)
       * The ONE and ONLY interface supported in future!
       */
      template<typename T>
      void writeOne(T * payload, Summary * summary, 
		    Time_t time, const std::string& recordName, 
                    bool withlogging=false) {
	if (isNewTagRequest(recordName) ){
	  createNewIOV<T>(payload, summary,
                          time, endOfTime(), recordName, withlogging);
        }else{
	  appendSinceTime<T>(payload, summary, time, recordName, withlogging);
        }	
      }
      
      // close the IOVSequence setting lastTill
      void closeIOV(Time_t lastTill, const std::string& recordName, 
                    bool withlogging=false) {


      }

      // BW-compatible signature
      template<typename T>
      void createNewIOV( T* firstPayloadObj,
			 cond::Time_t firstSinceTime,
			 cond::Time_t firstTillTime,
			 const std::string& recordName,
                         bool withlogging=false){
        this->createNewIOV(firstPayloadObj, 0, firstSinceTime, firstTillTime, recordName,withlogging);
      }
      
      //
      // insert the payload and its valid since time into the database
      // Note: user looses the ownership of the pointer to the payloadObj
      // The payload object will be stored as well
      // 
      template<typename T> void createNewIOV( T* firstPayloadObj,  Summary * summary,
                                              cond::Time_t firstSinceTime,
                                              cond::Time_t firstTillTime,
                                              const std::string& recordName,
                                              bool withlogging=false){
        createNewIOV( GetTokenFromPointer<T>(firstPayloadObj, summary),
                      firstSinceTime,
                      firstTillTime,
                      recordName,
                      withlogging);	
      }
      
      void createNewIOV( const std::string& firstPayloadToken,
                         cond::Time_t firstSinceTime,
                         cond::Time_t firstTillTime,
                         const std::string& recordName,
                         bool withlogging=false) {
        createNewIOV( GetTrivialToken(firstPayloadToken),
                      firstSinceTime,
                      firstTillTime,
                      recordName,
                      withlogging);
      }

      
      // BW-compatible signature
      template<typename T> void appendSinceTime( T* payloadObj,
                                                 cond::Time_t sinceTime,
                                                 const std::string& recordName,
                                                 bool withlogging=false){
        this->appendSinceTime<T>(payloadObj, 0, sinceTime, recordName, withlogging);
      }
      
      template<typename T>
      void appendSinceTime( T* payloadObj, Summary * summary,
			    cond::Time_t sinceTime,
                              const std::string& recordName,
			    bool withlogging=false){
        add( GetTokenFromPointer<T>(payloadObj,summary),
	     sinceTime,
	     recordName,
	     withlogging);
      }
      
      // Append the payload and its valid sinceTime into the database
      // Note: user looses the ownership of the pointer to the payloadObj
      // Note: the iov index appended to MUST pre-existing and the existing 
      // conditions data are retrieved from the DB
      // 
      void appendSinceTime( const std::string& payloadToken,
                            cond::Time_t sinceTime,
                            const std::string& recordName,
                            bool withlogging=false) {
        add(GetTrivialToken(payloadToken),
            sinceTime,
            recordName,
            withlogging);
      }
     
      // set last till so that the iov sequence is "closed"
      // void closeSequence(cond::Time_t tillTime,
      //                 const std::string& recordName,
      //                 bool withlogging=false);


      //
      // Service time utility method 
      // return the infinity value according to the given timetype
      //
      cond::Time_t endOfTime() const;
      //
      // Service time utility method 
      // return beginning of time value according to the given timetype
      //
      cond::Time_t beginOfTime() const;
      //
      // Service time utility method 
      // return the time value of the current edm::Event according to the 
      // given timetype
      //
      cond::Time_t currentTime() const;

      // optional. User can inject additional information into the log associated with a given record
      void setLogHeaderForRecord(const std::string& recordName,
				 const std::string& provenance,
				 const std::string& usertext);
      // 
      // Retrieve tag information of the data
      // 
      void tagInfo(const std::string& recordName,
		   cond::TagInfo& result );
      
      virtual ~PoolDBOutputService();  
      
    private:

      struct Record{
	Record(): m_tag(),m_isNewTag(false),
		  m_idName(),
		  m_iovtoken(),
                  m_closeIOV(false),
		  m_freeInsert(false),
		  m_withWrapper(false)
	{}

	std::string timetypestr() const { return cond::timeTypeSpecs[m_timetype].name;}
	std::string m_tag;
	bool m_isNewTag;
	std::string m_idName;
	std::string m_iovtoken;
	cond::TimeType m_timetype;
        bool m_closeIOV;
	bool m_freeInsert;
	bool m_withWrapper;
    };      



      void fillRecord( edm::ParameterSet & pset);
      
      void createNewIOV( GetToken const & token, 
			 cond::Time_t firstSinceTime, 
			 cond::Time_t firstTillTime,
			 const std::string& recordName,
			 bool withlogging=false);
      
      void add( GetToken const & token,  
		cond::Time_t time,
		const std::string& recordName,
		bool withlogging=false);
      
      
      void connect();    
      void disconnect();
      void initDB();
      unsigned int appendIOV(cond::DbSession&,
                             Record& record,
                             const std::string& payloadToken,
                             cond::Time_t sinceTime);
      
      /// Returns payload location index 
      unsigned int 
      insertIOV(cond::DbSession& pooldb,
		Record& record,
		const std::string& payloadToken,
		cond::Time_t tillTime);
      //			    const std::string& recordName);
      
      Record & lookUpRecord(const std::string& recordName);
      cond::UserLogInfo& lookUpUserLogInfo(const std::string& recordName);
      
    private:
      cond::TimeType m_timetype; 
      std::string m_timetypestr;
      cond::Time_t m_currentTime;
      cond::DbConnection m_connection;
      cond::DbSession m_session;
      cond::DbSession m_logSession;
      std::map<std::string, Record> m_callbacks;
      std::vector< std::pair<std::string,std::string> > m_newtags;
      bool m_dbstarted;
      cond::Logger* m_logdb;
      bool m_logdbOn;

      bool m_closeIOV;
      
      bool m_freeInsert;
      
      bool m_withWrapper;
      
      std::map<std::string, cond::UserLogInfo> m_logheaders;

    };//PoolDBOutputService
  }//ns service
}//ns cond
#endif
