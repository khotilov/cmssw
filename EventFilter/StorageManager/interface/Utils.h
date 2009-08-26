// $Id: Utils.h,v 1.6 2009/08/20 13:40:03 mommsen Exp $
/// @file: Utils.h 

#ifndef StorageManager_Utils_h
#define StorageManager_Utils_h

#include <string>
#include <vector>

#include "xcept/Exception.h"
#include "xdaq/ApplicationDescriptor.h"
#include "xdata/String.h"
#include "xdata/Vector.h"


namespace stor {

  namespace utils {

    /**
     * Collection of utility functions used in the storage manager
     *
     * $Author: mommsen $
     * $Revision: 1.6 $
     * $Date: 2009/08/20 13:40:03 $
     */

    /**
       time_point_t is used to represent a specific point in time,
       measured by some specific clock. We rely on the "system
       clock". The value is represented by the number of seconds
       (including a fractional part that depends on the resolution of
       the system clock) since the beginning of the "epoch" (as defined
       by the system clock).
    */
    typedef double time_point_t;

    /**
       durtion_t is used to represent a duration (the "distance" between
       two points in time). The value is represented as a number of
       seconds (including a fractional part).
    */
    typedef double duration_t;

    /**
       Returns the current point in time. A negative value indicates
       that an error occurred when fetching the time from the operating
       system.
    */
    time_point_t getCurrentTime();

    /**
       Sleep for at least the given duration. Note that the underlying
       system will round the interval up to an integer multiple of the
       system's sleep resolution. (The underlying implementation
       relies upon nanosleep, so see documentation for nanosleep for
       details). Negative intervals are not allowed, and result in an
       error (returning -1, and no sleeping).
     */
    int sleep(duration_t interval);

    /**
       Sleep until at least the given time_point_t. Uses internally
       sleep(duration_t). See notes to this method.
    */
    int sleepUntil(time_point_t);

    /**
       Converts a time_point_t into a string.
       Note: the string formatting is used by the file summary catalog and
       may or may not depend on the actual formatting
    */
    std::string timeStamp(time_point_t);

    /**
       Returns an identifier string composed of class name and instance
    */
    std::string getIdentifier(xdaq::ApplicationDescriptor*);
    
    /**
       Throws a stor::exception::NoSuchDirectory when the directory does not exist
     */
    void checkDirectory(const std::string&);
 
    /**
       Conversions between std::vector<std::string> and xdata::Vector<xdata::String>
    */
    void getStdVector(xdata::Vector<xdata::String>&, std::vector<std::string>&);
    void getXdataVector(const std::vector<std::string>&, xdata::Vector<xdata::String>&);

  } // namespace utils
  
} // namespace stor

#endif // StorageManager_Utils_h 


/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
