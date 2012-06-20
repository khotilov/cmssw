#ifndef DDException_h
#define DDException_h

#include "FWCore/Utilities/interface/Exception.h"

//! An exception for DDD errors
/** @class DDException DDException.h
 *
 *  @author:  Martin Liendl               Initial Version
 *  @version: 0.0
 *  @date:    
 * 
 *  Description:
 *  
 *  Provides an exception for DDD errors.
 *
 *  Modifications:
 *  MEC:   8 June 2005 Michael Case: changed to inherit from seal::Error
 *  MEC:   25 April 2007 Michael Case: changed to inherit from cms:Exception
 */

class DDException : public cms::Exception 
{
 public: 
  //! constructor takes simply an error message via a std::string
  explicit DDException(const std::string & s);
  DDException();
  DDException(const DDException& dde);
  void swap(DDException& other) throw();
  
  DDException& operator=(DDException const& other);
  
  virtual ~DDException() throw();
};    

#endif
