#ifndef CondIter_CondIter_h
#define CondIter_CondIter_h


#include "CondCore/Utilities/interface/CondBasicIter.h"

#include "CondCore/DBCommon/interface/PayloadRef.h"

#include <vector>


template <class DataT>
class CondIter : public  CondBasicIter{
  
protected:
  virtual bool load(pool::IDataSvc * svc, std::string const & itoken) {
    if (useCache)
      if (n>=cache.size()) {
	cache.resize(n+1); 
	return cache.back().load(svc,itoken);
      } else return true;
    else return data.load(svc,itoken);
  }
  
private:
  bool initialized;
  bool useCache;
  cond::PayloadRef<DataT> data;
  std::vector<cond::PayloadRef<DataT> > cache;
  size_t n;

public:
  
  
  CondIter(bool cacheIt=false) : initialized(false), useCache(cacheIt),n(0){}
  ~CondIter(){}
  
  
  void reset() { initialized=false;}

  virtual void clear() {
    reset();
    cache.clear();
  }
  
 
  /**
     Obtain the pointer to an object T. If it is the last T the method returns a null pointer.
  */ 
  DataT const * next() {
    bool ok=false;
    if (!initialized) {
      n=0;
      ok =init();
    }
    else {
      ++n;
      ok = forward();
    }
    if (!ok) return 0;
    ok = make();
    if (!ok) return 0;
    return  useCache ?  &(*cache[n]) : &(*data()); 

  }

};




#endif

