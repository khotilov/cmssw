//$Id: Utils.cc,v 1.8 2009/08/26 15:19:58 mommsen Exp $
/// @file: Utils.cc

#include "EventFilter/StorageManager/interface/Exception.h"
#include "EventFilter/StorageManager/interface/Utils.h"

#include <iomanip>
#include <sstream>

#include <sys/time.h>
#include <sys/stat.h>


namespace stor
{
  namespace utils
  {
    
    namespace
    {
      /**
	 Convert a (POSIX) timeval into a time_point_t.
      */
      inline void timeval_to_timepoint(timeval const& in, 
				       time_point_t& out)
      {
	// First set the seconds.
	out = static_cast<time_point_t>(in.tv_sec);
	
	// Then set the microseconds.
	out += static_cast<time_point_t>(in.tv_usec)/(1000*1000);
      }
    } // anonymous namespace
      
      
    time_point_t getCurrentTime()
    {
      time_point_t result = -1.0;
      timeval now;
      if (gettimeofday(&now, 0) == 0) timeval_to_timepoint(now, result);
      return result;
    }

    int sleep(duration_t interval)
    {
      if (interval < 0) return -1;
      timespec rqtp;
      rqtp.tv_sec = static_cast<time_t>(interval); // truncate
      rqtp.tv_nsec = static_cast<long>((interval-rqtp.tv_sec)*1000*1000*1000);
      return nanosleep(&rqtp, 0);
    }

    int sleepUntil(time_point_t theTime)
    {
      time_point_t now = getCurrentTime();
      duration_t interval = theTime - now;
      return sleep(interval);
    }
    
    std::string timeStamp(time_point_t theTime)
    {
      time_t rawtime = (time_t)theTime;
      tm * ptm;
      ptm = localtime(&rawtime);
      std::ostringstream timeStampStr;
      std::string colon(":");
      std::string slash("/");
      timeStampStr << std::setfill('0') << std::setw(2) << ptm->tm_mday      << slash 
                   << std::setfill('0') << std::setw(2) << ptm->tm_mon+1     << slash
                   << std::setfill('0') << std::setw(4) << ptm->tm_year+1900 << colon
                   << std::setfill('0') << std::setw(2) << ptm->tm_hour      << slash
                   << std::setfill('0') << std::setw(2) << ptm->tm_min       << slash
                   << std::setfill('0') << std::setw(2) << ptm->tm_sec;
      return timeStampStr.str();
    }

    
    std::string getIdentifier(xdaq::ApplicationDescriptor *appDesc)
    {
      std::ostringstream identifier;
      identifier << appDesc->getClassName() << appDesc->getInstance() << "/";
      return identifier.str();
    }
    
    
    void checkDirectory(const std::string& path)
    {
      struct stat64 buf;
      
      int retVal = stat64(path.c_str(), &buf);
      if( retVal !=0 )
      {
        std::ostringstream msg;
        msg << "Directory " << path << " does not exist. Error=" << errno;
        XCEPT_RAISE(stor::exception::NoSuchDirectory, msg.str());
      }
    }


    void getStdVector(xdata::Vector<xdata::String>& x, std::vector<std::string>& s)
    {
      s.clear();
      s.reserve(x.elements());
      for(xdata::Vector<xdata::String>::iterator it = x.begin(),
            itEnd = x.end();
          it != itEnd;
          ++it)
      {
        s.push_back( it->toString() );
      }
    }


    void getXdataVector(const std::vector<std::string>& v, xdata::Vector<xdata::String>& x)
    {
      x.clear();
      x.reserve(v.size());
      for(std::vector<std::string>::const_iterator it = v.begin(),
            itEnd = v.end();
          it != itEnd;
          ++it)
      {
        x.push_back( static_cast<xdata::String>(*it) );
      }
    }


  } // namespace utils

} // namespace stor

/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
