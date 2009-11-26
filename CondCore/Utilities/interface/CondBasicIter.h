#ifndef CondIter_CondBasicIter_h
#define CondIter_CondBasicIter_h
#include "CondCore/Utilities/interface/CondPyInterface.h"
#include <string>


class CondBasicIter{

protected:
  cond::RDBMS rdmbs;
  cond::CondDB db;
  cond::IOVProxy iov;
  bool initialized;
  cond::IOVProxy::const_iterator iter;
 
 public:
  
  CondBasicIter();    
  ~CondBasicIter();    
  
  
  /**
     tell Iter to point to a database. After this call Iter can be used.
     Direct Access to database through frontier
     It needs:
     \li \c NameDB -> name of the database (connection string)
     \li \c Tag -> Tag human-readable of the content of the database
     \li \c User -> name of the User (if you don't need to authenticate don't write anything here)
     \li \c Pass -> Password to access database (if you don't need to authenticate don't write anything here)
     \li \c nameBlob -> to handle blob type of data (if it is not needed this field has to be left empty)
  */
  
  CondBasicIter(const std::string & NameDB,
		const std::string & Tag,
		const std::string & User,
		const std::string & Pass,
		const std::string & nameBlob = ""
		);

  CondBasicIter(const std::string & NameDB,
		const std::string & Tag,
		const std::string & auth = ""
		);

  void create(const std::string & NameDB,
	      const std::string & Tag,
	      const std::string & User,
	      const std::string & Pass,
	      const std::string & nameBlob = ""
	      );

  void create(const std::string & NameDB,
	      const std::string & Tag,
	      const std::string & auth = ""
	      );
  
  
  /**
     Set the range of interest of the Iterator of the IOVs.
  */ 
  void setRange(unsigned int min,unsigned int max);
  
  
  /**
     Set the minimum of the range of interest of the Iterator of the IOVs.
  */  
  
  void setMin(unsigned int min);
  
  /**
     Set the maximum of the range of interest of the Iterator of the IOVs.
  */  
  
  void setMax(unsigned int max);
 
  /**
     Get the mean time of the Iterval of Validity.
  */  
  unsigned int getTime();
  
  /**
     Get the SINCE TIME of the Interval of Validity.
  */
  unsigned int getStartTime();
  
  /**
     Get the TILL TIME of the Interval of Validity.
  */
  unsigned int getStopTime();
  


  bool init();
  bool forward();
  boll make();
  virtual bool load(pool::IDataSvc * svc, std::string const & token) =0;

};


#endif
