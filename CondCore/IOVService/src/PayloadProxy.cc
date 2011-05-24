#include "CondCore/IOVService/interface/PayloadProxy.h"

#include "CondCore/ORA/interface/Exception.h"
#include "CondCore/DBCommon/interface/Exception.h"
#include "CondCore/ORA/interface/OId.h"

namespace {

  void fillIt(cond::BasePayloadProxy::ObjId & obj, std::string const & token, cond::Time_t since) {
    obj.since = since;
    ora::OId oid;
    oid.fromString( token );
    obj.oid1 = oid.containerId();
    obj.oid2 = oid.itemId();
  }

}

namespace cond {

  BasePayloadProxy::Stats BasePayloadProxy::gstats = {0,0,0};


  BasePayloadProxy::BasePayloadProxy(cond::DbSession& session,
                                     const std::string & token,
                                     bool errorPolicy) :
    m_doThrow(errorPolicy), m_iov(session,token),m_session(session) {
    ++gstats.nProxy;
    BasePayloadProxy::Stats s = {0,0,0,ObjIds()};
    stats = s;
  }


  BasePayloadProxy::~BasePayloadProxy(){}

  cond::ValidityInterval BasePayloadProxy::loadFor(cond::Time_t time) {
    m_element = *m_iov.find(time);
    make();
    return cond::ValidityInterval(m_element.since(),m_element.till());
  }

  cond::ValidityInterval BasePayloadProxy::loadFor(size_t n) {
    m_element.set(m_iov.iov(),n);
    make();
    return cond::ValidityInterval(m_element.since(),m_element.till());
  }


  void  BasePayloadProxy::make() {
    ++gstats.nMake;
    ++stats.nMake;
    bool ok = false;
    if ( isValid()) {
      // check if (afterall) the payload is still the same...
      if (m_element.token()==token()) return;
      try {
        ok = load( m_session ,m_element.token());
	if (ok) m_token = m_element.token();
      }	catch( const ora::Exception& e) {
        if (m_doThrow) throw cond::Exception(std::string("Condition Payload loader: ")+ e.what());
        ok = false;
      }
    }

    if (!ok) {
      m_element.set(cond::invalidTime,cond::invalidTime,"");
      m_token.clear();
      if (m_doThrow)
        throw cond::Exception("Condition Payload loader: invalid data");
    }
    if (ok) { 
      ++gstats.nLoad; ++stats.nLoad;
      stats.ids.push_back(ObjId());
      fillIt(stats.ids.back(),m_token, m_element.since());
    }
  }


  cond::ValidityInterval BasePayloadProxy::setIntervalFor(cond::Time_t time) {
    //FIXME: shall handle truncation...
    if ( (!(time<m_element.till())) || time<m_element.since() )
      m_element = *m_iov.find(time);
    return cond::ValidityInterval(m_element.since(),m_element.till());
  }
    
  bool BasePayloadProxy::isValid() const {
    return m_element.since()!=cond::invalidTime && m_element.till()!=cond::invalidTime
      &&  !m_element.token().empty();
  }


  bool  BasePayloadProxy::refresh() {
    bool anew = m_iov.refresh();
    if (anew)  m_element = IOVElementProxy();
    return anew;
  }





}
