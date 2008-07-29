#ifndef  PayLoadInspector_H
#define  PayLoadInspector_H
#include "CondCore/DBCommon/interface/TypedRef.h"
#include "CondCore/IOVService/interface/IOVProxy.h"
#include <string>
#include <vector>

// to be moved in src
#include "CondCore/DBCommon/interface/PoolTransaction.h"


namespace cond {
  // to be moved elsewhere
  class PoolTransactionSentry {
  public:
    PoolTransactionSentry(){}
    PoolTransactionSentry(cond::PoolTransaction & db) : 
      elem(new Elem(db)){}      
  private:
    struct Elem {
      Elem(cond::PoolTransaction & db) : pooldb(db){
	pooldb.start(true);
      }
      ~Elem() { pooldb.commit();}
      cond::PoolTransaction & pooldb;
    };
    boost::shared_ptr<Elem> elem;
      
  };

  template<typename T>
  class BaseValueExtractor {
  public:
    typedef T Class;

    BaseValueExtractor(){}
    
    virtual ~BaseValueExtractor(){}
    void computeW(Class const &o){
      reset();
      compute(o);
    }
    std::vector<float> const & values() const { return m_values;}
  protected:
    void add(float v) { m_values.push_back(v); }
  private:
    void reset() { m_values.clear();}
    virtual void compute(Class const &){}
    
    
  private:
    std::vector<float> m_values;
  };


  // a collection of enumerators...
  template<typename T>
  struct ExtractWhat {

  };



  template<typename T>
  class ValueExtractor : public  BaseValueExtractor<T> {
  public:
    typedef T Class;
    typedef ExtractWhat<Class> What;
    ValueExtractor(){}
    ValueExtractor(What const &, std::vector<int> const&){}
    static What what() { return What();}
  private:
    void compute(Class const &){}
  private:
  };

  template<typename T>
  class PayLoadInspector : PoolTransactionSentry {
  public:
    typedef T Class;
    typedef ValueExtractor<T> Extractor;
        
    PayLoadInspector() {}
    PayLoadInspector(const cond::IOVElement & elem) : 
      PoolTransactionSentry(*elem.db()),
      object(*elem.db(),elem.payloadToken()){}

    std::string dump() const { return ""; }

    std::string summary() const {return ""; }

    // return the real name of the file including extension...
    std::string plot(std::string const & /* filename */,
		     std::string const &, 
		     std::vector<int> const&, std::vector<float> const& ) const {return "";}

    void extract(Extractor & extractor) const {extractor.computeW(*object); }

  private:
    cond::TypedRef<Class> object;    

  };

}

#endif //   PayLoadInspector_H
