#ifndef Tracker_ProxyBase_H
#define Tracker_ProxyBase_H

#include "TrackingTools/TrajectoryParametrization/interface/TrajectoryStateExceptions.h"
#include "TrackingTools/TrajectoryState/interface/CopyUsingNew.h"

/** A base class for reference counting proxies.
 *  The class to which this one is proxy must inherit from
 *  ReferenceCounted.
 *  A ProxyBase has value semantics, in contrast to ReferenceCountingPointer,
 *  which has pointer semantics.
 *  The Proxy class inheriting from ProxyBase must duplicate the
 *  part of the interface of the reference counted class that it whiches to expose.
 */

template <class T, class Cloner > 
class ProxyBase {
public:

  ProxyBase() : theData(0) {}

  ProxyBase( T* p) : theData(p) {if (theData) theData->addReference();}

  ProxyBase( const ProxyBase& other) {
    theData = other.theData;
    if (theData) theData->addReference();
  }

  virtual ~ProxyBase() { 
    destroy();
  }

  ProxyBase& operator=( const ProxyBase& other) {
    if ( theData != other.theData) { 
      destroy();
      theData = other.theData;
      if (theData) theData->addReference();
    }
    return *this;
  }

  const T& data() const { check(); return *theData;}

  T& unsharedData() {
    check(); 
    if ( theData->references() > 1) {
      theData->removeReference();
      theData = Cloner().clone( *theData);
      if (theData) theData->addReference();
    }
    return *theData;
  }

  T& sharedData() { check(); return *theData;}

  bool isValid() const { return theData != 0;}

  void check() const {
    if (theData == 0) {
      throw TrajectoryStateException("Error: uninitialized ProxyBase used");
    }
  }

  void destroy() { if (isValid()) theData->removeReference();}

  int  references() const {return theData->references();}  

private:
  T* theData;
};


#endif // Tracker_ProxyBase_H
