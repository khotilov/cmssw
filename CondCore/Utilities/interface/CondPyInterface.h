#ifndef CondCore_Utilities_CondPyInterface_h
#define CondCore_Utilities_CondPyInterface_h

/*  common utilities of the CondCore Python buiding
 *
 */

#include "CondCore/IOVService/interface/IOVProxy.h"
#include "CondCore/DBCommon/interface/LogDBEntry.h"


#include<boost/shared_ptr.hpp>
#include<string>

namespace cond {

  class DBSession;
  class Connection;
  class Logger;

  namespace impl {
    struct FWMagic;
  }

  // initialize framework
  class FWIncantation {
  public:
    FWIncantation();
    // FWIncantation(FWIncantation const & other );
    ~FWIncantation();
    
  private:
    boost::shared_ptr<impl::FWMagic> magic;
  };

  // a readonly CondDB and its transaction
 class CondDB {
  public:
   CondDB();
   CondDB(const CondDB & other);
   CondDB & operator=(const CondDB & other);
   CondDB(Connection * conn, boost::shared_ptr<cond::Logger> ilog );
   ~CondDB();
   std::string allTags() const;

   IOVProxy iov(std::string const & tag) const;
   IOVProxy iovWithLib(std::string const & tag) const;

   IOVElementProxy payLoad(std::string const & token) const;

   std::string iovToken(std::string const & tag) const;
   
   cond::LogDBEntry lastLogEntry(std::string const & tag) const;
   cond::LogDBEntry lastLogEntryOK(std::string const & tag) const;

 private:
   mutable Connection * me;
   boost::shared_ptr<cond::Logger> logger;
 };

  // initializ cond, coral etc
  class RDBMS {
  public:
    RDBMS();
    ~RDBMS();
    explicit RDBMS(std::string const & authPath);
    RDBMS(std::string const & user,std::string const & pass);
    void setLogger(std::string const & connstr);

    CondDB getDB(std::string const & db);

  private:
    boost::shared_ptr<DBSession> session;
    boost::shared_ptr<cond::Logger> logger;

  };


}


#endif //  CondCore_Utilities_CondPyInterface_h
