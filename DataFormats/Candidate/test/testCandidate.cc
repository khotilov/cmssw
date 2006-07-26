// $Id: testCandidate.cpp,v 1.1 2006/06/20 14:54:25 llista Exp $
#include <cppunit/extensions/HelperMacros.h>
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include <memory>

class testCandidate : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(testCandidate);
  CPPUNIT_TEST(checkAll);
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp() {}
  void tearDown() {}
  void checkAll(); 
};

CPPUNIT_TEST_SUITE_REGISTRATION(testCandidate);

namespace test {
  
  struct DummyComponent {
    int x;
  };

  struct DummyComponent2 {
    int x;
  };

  class DummyCandidate1 : public reco::LeafCandidate {
  public:
    DummyCandidate1( const LorentzVector & p, Charge q, int x, 
		     int y1, int y2) : reco::LeafCandidate( q, p ) { 
      c.x = x;
      cc[0].x = y1;
      cc[1].x = y2;
    }
    virtual DummyCandidate1 * clone() const { return new DummyCandidate1( * this ); }
    DummyComponent cmp() const { return c; }
    DummyComponent2 cmp2( size_t i ) const { return cc[i]; }
    size_t cmpSize2() const { return 2; }
  private:
    DummyComponent c;
    DummyComponent2 cc[2];
  };

}

namespace reco {
  GET_CANDIDATE_COMPONENT( test::DummyCandidate1, test::DummyComponent, DefaultComponentTag, cmp );
  GET_CANDIDATE_MULTIPLECOMPONENTS( test::DummyCandidate1, test::DummyComponent2, DefaultComponentTag, cmp2, cmpSize2 );
}

void testCandidate::checkAll() {
  reco::Particle::LorentzVector p( 1.0, 2.0, 3.0, 4.0 );
  reco::Particle::Charge q( 1 );
  int x = 123, y0 = 111, y1 = 222;
  std::auto_ptr<reco::Candidate> c( new test::DummyCandidate1( p, q, x, y0, y1 ) );
  CPPUNIT_ASSERT( c->charge() == q );
  CPPUNIT_ASSERT( c->p4() == p );
  CPPUNIT_ASSERT( c->numberOfDaughters() == 0 );
  CPPUNIT_ASSERT( c->get<test::DummyComponent>().x == x );
  CPPUNIT_ASSERT( c->numberOf<test::DummyComponent>() == 1 );
  CPPUNIT_ASSERT( c->get<test::DummyComponent2>( 0 ).x == y0 );
  CPPUNIT_ASSERT( c->get<test::DummyComponent2>( 1 ).x == y1 );
  CPPUNIT_ASSERT( c->numberOf<test::DummyComponent2>() == 2 );
}
