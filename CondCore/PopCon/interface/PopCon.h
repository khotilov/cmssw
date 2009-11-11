
#ifndef POPCON_POPCON_H
#define POPCON_POPCON_H
//
// Author: Vincenzo Innocente
// Original Author:  Marcin BOGUSZ
// 


#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CondCore/DBCommon/interface/TagInfo.h"
#include "CondCore/DBCommon/interface/LogDBEntry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"

#include "CondCore/DBCommon/interface/Time.h"

#include "CondCore/DBCommon/interface/IOVInfo.h"


#include <boost/bind.hpp>
#include <algorithm>
#include <vector>
#include <string>


#include<iostream>

namespace cond {
  class Summary;
}


namespace popcon {


  /* Populator of the Condition DB
   *
   */
  class PopCon {
  public:
    typedef cond::Time_t Time_t;
    typedef cond::Summary Summary;

    PopCon(const edm::ParameterSet& pset);
     
     virtual ~PopCon();

     template<typename Source>
       void write(Source const & source);

     template<typename T>
     void writeOne(T * payload, Summary * summary, Time_t time);

   
    
     
  private:
     void initialize();
     void finalize();


  private:

    edm::Service<cond::service::PoolDBOutputService> m_dbService;
    
    std::string  m_record;
    
    std::string m_payload_name;
    
    bool m_since;
    
    bool m_LoggingOn;
    
    bool m_IsDestDbCheckedInQueryLog;

    std::string m_tag;
    
    cond::TagInfo m_tagInfo;
    
    cond::LogDBEntry m_logDBEntry;
    
  };


  template<typename T>
  void PopCon::writeOne(T * payload, Summary * summary, Time_t time) {
    m_dbService->writeOne(payload, summary, time, m_record, m_LoggingOn, m_since);
  }

  
  template<typename Container>
  void displayHelper(Container const & payloads, bool sinceAppend) {
    typename Container::const_iterator it;
    for (it = payloads.begin(); it != payloads.end(); it++)
      edm::LogInfo ("PopCon")<< (sinceAppend ? "Since " :" Till ") << (*it).time << std::endl;
  }     
  
  
  template<typename Container>
  const std::string displayIovHelper(Container const & payloads, bool sinceAppend) {
    if (payloads.empty()) return "Nothing to transfer;";
    std::ostringstream s;    
    // when only 1 payload is transferred; 
    if ( payloads.size()==1)  
      s <<(sinceAppend ? "Since " :" Till ") << (*payloads.begin()).time <<  "; " ;
    else{
      // when more than one payload are transferred;  
      s <<   "first payload " << (sinceAppend ? "Since " :" Till ") <<  (*payloads.begin()).time <<  ","
	<< "last payload " << (sinceAppend ? "Since " :" Till ") << (*payloads.rbegin()).time <<  ";" ;  
    }  
    return s.str();
  }
  
  
  
  
  
  template<typename Source>
  void PopCon::write(Source const & source) {
    typedef typename Source::value_type value_type;
    typedef typename Source::Container Container;
    
    initialize();
    std::pair<Container const *, std::string const> ret = source(m_dbService->session(),
								 m_tagInfo,m_logDBEntry); 
    Container const & payloads = *ret.first;
    
    
    
    
    m_dbService->setLogHeaderForRecord(m_record,source.id(),"PopCon v2.1; " + cond::userInfo() + displayIovHelper(payloads,m_since ) +  ret.second);
    
    
    displayHelper(payloads,m_since);
    
    std::for_each(payloads.begin(),payloads.end(),
		  boost::bind(&popcon::PopCon::writeOne<value_type>,this,
			      boost::bind(&Container::value_type::payload,_1),
			      boost::bind(&Container::value_type::summary,_1),
			      boost::bind(&Container::value_type::time,_1)
			      )
		   );
    
    
    finalize();
  }
  
  
  
}

#endif //  POPCON_POPCON_H


