#ifndef CondCore_IOVService_PayloadProxy_h
#define CondCore_IOVService_PayloadProxy_h

#include "CondCore/IOVService/interface/IOVProxy.h"
#include "CondCore/DBCommon/interface/PayloadRef.h"
#include "CondFormats/Common/interface/PayloadWrapper.h"

namespace pool{
  class IDataSvc;
}


namespace cond {

  /* get iov by name (key, tag, ...)

   */
  class CondGetter {
  public:
    virtual ~CondGetter(){}
    virtual IOVProxy get(std::string name) const=0;

  };

  /* implementation detail: 
    implements the not templated part...

   */
  class BasePayloadProxy {
  public:

    struct Stats {
      int nProxy;
      int nMake;
      int nLoad;
    };

    static Stats stats;

    // errorPolicy=true will throw on load, false will set interval and token to invalid
    BasePayloadProxy(cond::DbSession& session,
                     const std::string & token,
                     bool errorPolicy);
    
    virtual ~BasePayloadProxy();

    virtual void invalidateCache()=0;

    // current cached object token
    std::string const & token() const { return m_token;}
    
    // load Element valid at time
    void loadFor(cond::Time_t time);

    // find ad return interval (does not load)
    cond::ValidityInterval setIntervalFor(cond::Time_t time);
    
    // load element if interval is valid
    void make();

    bool isValid() const;

    TimeType timetype() const { return m_iov.timetype();}

    IOVProxy const & iov() const { return m_iov;}

    virtual void loadMore(CondGetter const &){}

    // reload the iov return true if size has changed
    bool refresh();


  private:
    virtual bool load(pool::IDataSvc * svc, std::string const & token) =0;   


  protected:
    bool m_doThrow;
    IOVProxy m_iov;
    IOVElementProxy m_element;

  private:
    // current loaded payload
    std::string  m_token;

  };


  /* proxy to the payload valid at a given time...

   */
  template<typename DataT>
  class PayloadProxy : public BasePayloadProxy {
  public:
 
    PayloadProxy(cond::DbSession& session,
                 const std::string & token,
                 bool errorPolicy,
                 const char * source=0) :
      BasePayloadProxy(session, token, errorPolicy) {}
    
    virtual ~PayloadProxy(){}

    // dereference (does not load)
    const DataT & operator()() const {
      return (*m_data); 
    }
        
    virtual void invalidateCache() {
      // don't, preserve data for future access
      // m_data.clear();
    }


  protected:
    virtual bool load(pool::IDataSvc * svc, std::string const & itoken) {
      return m_data.load(svc,itoken);
    }

  private:
     cond::PayloadRef<DataT> m_data;
  };

}
#endif // CondCore_IOVService_PayloadProxy_h
