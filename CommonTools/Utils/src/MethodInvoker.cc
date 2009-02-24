#include "CommonTools/Utils/src/MethodInvoker.h"
#include "CommonTools/Utils/src/findMethod.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include <algorithm>
using namespace reco::parser;
using namespace Reflex;
using namespace std;

MethodInvoker::MethodInvoker(const Member & method, const vector<AnyMethodArgument> & ints) :
  method_(method), ints_(ints), isFunction_(method.IsFunctionMember())
{ 
  setArgs();
  /*std::cout << "Booking " << method_.Name() 
            << " from " << method_.DeclaringType().Name() 
            << " with " << args_.size() << " arguments"
            << " (were " << ints.size() << ")"
            << std::endl; */
}

MethodInvoker::MethodInvoker(const MethodInvoker & other) :
  method_(other.method_), ints_(other.ints_), isFunction_(other.isFunction_) {
  setArgs();
}

MethodInvoker & MethodInvoker::operator=(const MethodInvoker & other) {
  method_ = other.method_;
  ints_ = other.ints_;
  isFunction_ = other.isFunction_;
  setArgs();
  return *this;
}

void MethodInvoker::setArgs() {
  for(size_t i = 0; i < ints_.size(); ++i) {
      args_.push_back( boost::apply_visitor( AnyMethodArgument2VoidPtr(), ints_[i] ) );
  }
}

std::pair<Object,bool> MethodInvoker::value(const Object & o) const {
  Object ret;
  bool   owned = false;
  /* no need to check the type at run time
  if(method_.IsVirtual()) {
    Type dynType = o.DynamicType();
    Member met = reco::findMethod(dynType, method_.Name(), args_.size()).first;
    if(!met)
      throw edm::Exception(edm::errors::InvalidReference)
	<< "method \"" << method_.Name() << "\" not found in dynamic type \"" 
	<< dynType.Name() << "\"\n";
    ret = met.Invoke(Object(dynType, o.Address()), args_);
    } else */
  /*std::cout << "Invoking " << method_.Name() 
            << " from " << method_.DeclaringType().Name() 
            << " on an instance of " << o.DynamicType().Name() 
            << " with " << args_.size() << " arguments"
            << std::endl;*/
  if(isFunction_) {
     method_.Invoke(o, ret, args_);
  } else {
     ret = method_.Get(o);
  }
  void * addr = ret.Address(); 
  if(addr==0)
    throw edm::Exception(edm::errors::InvalidReference)
      << "method \"" << method_.Name() << "\" called with " << args_.size() 
      << " arguments returned a null pointer ";   
  Type retType = ret.TypeOf();
  if (retType.IsClass() && !retType.IsPointer() && !retType.IsReference()) {
    //std::cout << "Object, and not pointer, at " << addr << ", type " <<  retType.Name() 
    //          << ", from " << method_.Name() << ": I need to delete it" << std::endl;
    owned = true;
  }
  bool stripped = false;
  while(retType.IsTypedef()) { 
    retType = retType.ToType(); stripped = true; 
  }
  bool isPtr = retType.IsPointer() /* isRef = retType.IsReference() */;
  if(isPtr) {
    if(!stripped) {
      stripped = true;
      retType = retType.ToType();
      while(retType.IsTypedef()) {
	retType = retType.ToType();
      }
    }
  }
  if(stripped)
    ret = Object(retType, addr);
  if(!ret) 
     throw edm::Exception(edm::errors::Configuration)
      << "method \"" << method_.Name() 
      << "\" returned void invoked on object of type \"" 
      << o.TypeOf().Name() << "\"\n";
  return std::make_pair(ret,owned);
}
