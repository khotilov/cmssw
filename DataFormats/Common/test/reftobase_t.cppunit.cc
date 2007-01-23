// $Id: reftobase_t.cppunit.cc,v 1.4 2006/08/03 13:07:01 llista Exp $
#include <cppunit/extensions/HelperMacros.h>
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/Ref.h"

#include <vector>

class testRefToBase : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(testRefToBase);
  CPPUNIT_TEST(check);
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp() {}
  void tearDown() {}
  void check();
};

CPPUNIT_TEST_SUITE_REGISTRATION(testRefToBase);
namespace testreftobase {
  struct Base {
    virtual ~Base() {}
    virtual int val() const=0;
  };  

  struct Inherit1 : public Base {
    virtual int val() const { return 1;}
  };
  struct Inherit2 : public Base {
    virtual int val() const {return 2;}
  };
  
  template<class T>
    struct TestHandle {
      TestHandle(const edm::ProductID& iId, const T* iProd) : id_(iId), prod_(iProd) {}
      const edm::ProductID& id() const { return id_;}
      const T* product() const { return prod_;}
    private:
      edm::ProductID id_;
      const T* prod_;
    };
}

using namespace testreftobase;

void
testRefToBase::check()
{
  using namespace edm;

  std::vector<Inherit1> v1(2,Inherit1());
  std::vector<Inherit2> v2(2,Inherit2());
  
  TestHandle<std::vector<Inherit1> > h1(ProductID(1), &v1);
  Ref<std::vector<Inherit1> > r1(h1, 1);
  RefToBase<Base> b1(r1);
  CPPUNIT_ASSERT(&(*b1) == static_cast<Base*>(&(v1[1])));
  CPPUNIT_ASSERT(b1.operator->() == b1.get());
  CPPUNIT_ASSERT(b1.get() == static_cast<Base*>(&(v1[1])));
  CPPUNIT_ASSERT(b1.id() == ProductID(1));
  
  //copy constructor
  RefToBase<Base> b2(b1);
  CPPUNIT_ASSERT(&(*b2) == static_cast<Base*>(&(v1[1])));
  CPPUNIT_ASSERT(b2.id() == b1.id());

  //operator=
  RefToBase<Base> b3;
  CPPUNIT_ASSERT(b3.isNull());
  CPPUNIT_ASSERT(!(b3.isNonnull()));
  CPPUNIT_ASSERT(!b3);
  b3 = b1;
  CPPUNIT_ASSERT(&(*b3) == static_cast<Base*>(&(v1[1])));
  CPPUNIT_ASSERT(b3.id() == b1.id());
  CPPUNIT_ASSERT(!(b3.isNull()));
  CPPUNIT_ASSERT(b3.isNonnull());
  CPPUNIT_ASSERT(!(!b3));

  CPPUNIT_ASSERT(b1.castTo<Ref<std::vector<Inherit1> > >() == r1);
  bool throwed = false;
  try { b1.castTo<Ref<std::vector<Inherit2> > >(); } 
  catch (edm::Exception e) { throwed = true; }
  CPPUNIT_ASSERT(throwed);
}
