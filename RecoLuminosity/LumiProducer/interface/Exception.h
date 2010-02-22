#ifndef RecoLuminosity_LumiProducer_Exception_H
#define RecoLuminosity_LumiProducer_Exception_H
#include <string>
#include <exception>
namespace lumi{
  class Exception : public std::exception{
  public:
    Exception( const std::string& message,
	       const std::string& methodname,
	       const std::string& moduleName);
    virtual ~Exception() throw(){}
    virtual char const* what() const throw(){
      return m_message.c_str();
    }
  private:
    std::string m_message;
  };

}//ns lumi
#endif
