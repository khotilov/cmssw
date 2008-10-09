#include "CondFormats/Common/interface/IOVSequence.h"
#include <algorithm>
#include <boost/bind.hpp>


namespace cond {
  
  
  IOVSequence::IOVSequence(){}
  
  IOVSequence::IOVSequence(int type, cond::Time_t since, 
			   std::string const& imetadata) :
    timetype(type), firstsince(since),
    m_metadata(imetadata){}
    
  IOVSequence::~IOVSequence(){}
  
  
  size_t IOVSequence::add(cond::Time_t time, 
			  std::string const & payloadToken,
			  std::string const & metadataToken) {
    iovs.push_back(Item(time,  payloadToken, metadataToken));
    return iovs().size()-1;
  }
  
  IOVSequence::iterator IOVSequence::find(cond::Time_t time) {
    return std::lower_bound(iovs().begin(),iovs().end(),Item(time),
			    boost::bind(std::less<cond::Time_t>(),
					boost::bind(&Item::tillTime,_1),
					  boost::bind(&Item::tillTime,_2)
					)
			    );
  }
  
  IOVSequence::const_iterator IOVSequence::find(cond::Time_t time) const {
    return std::lower_bound(iovs().begin(),iovs().end(),Item(time),
			    boost::bind(std::less<cond::Time_t>(),
					boost::bind(&Item::tillTime,_1),
					boost::bind(&Item::tillTime,_2)
					)
			    );
  }
  

}









}
