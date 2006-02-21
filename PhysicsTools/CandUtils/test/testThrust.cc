// $Id: testThrust.cc,v 1.3 2006/01/31 11:23:56 llista Exp $
#include <cppunit/extensions/HelperMacros.h>
#include "PhysicsTools/Candidate/interface/LeafCandidate.h"
#include "PhysicsTools/CandUtils/interface/Thrust.h"
using namespace reco;

class testParticle : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(testParticle);
  CPPUNIT_TEST(checkAll);
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp() {}
  void tearDown() {}
  void checkAll(); 
};

CPPUNIT_TEST_SUITE_REGISTRATION(testParticle);

void testParticle::checkAll() {
  {
    CandidateCollection cands;
    cands.push_back( new LeafCandidate( +1, Particle::LorentzVector( 1, 1, 0, 1 ) ) );
    cands.push_back( new LeafCandidate( +1, Particle::LorentzVector( -1, -1, 0, 1 ) ) );
    Thrust t( cands.begin(), cands.end() );
    CPPUNIT_ASSERT( fabs( t.thrust() - 1.0 ) < 1.e-6 );
  }
  {
    CandidateCollection cands;
    cands.push_back( new LeafCandidate( +1, Particle::LorentzVector( 0.7, 0.5, 0, 1 ) ) );
    cands.push_back( new LeafCandidate( +1, Particle::LorentzVector( -0.7, -0.5, 0, 1 ) ) );
    Thrust t( cands.begin(), cands.end() );
    CPPUNIT_ASSERT( fabs( t.thrust() - 1.0 ) < 1.e-6 );
  }
  {
    CandidateCollection cands;
    cands.push_back( new LeafCandidate( +1, Particle::LorentzVector( 1, 0, 0, 1 ) ) );
    cands.push_back( new LeafCandidate( +1, Particle::LorentzVector( 0, 1, 0, 1 ) ) );
    Thrust t( cands.begin(), cands.end() );
    CPPUNIT_ASSERT( fabs( t.thrust() - sqrt( 2.0 )/2 ) < 1.e-6 );
  }
  {
    CandidateCollection cands;
    cands.push_back( new LeafCandidate( +1, Particle::LorentzVector( 1, 0, 0, 1 ) ) );
    cands.push_back( new LeafCandidate( +1, Particle::LorentzVector( 0, 1, 0, 1 ) ) );
    cands.push_back( new LeafCandidate( +1, Particle::LorentzVector( 1, 1, 0, 1 ) ) );
    Thrust t( cands.begin(), cands.end() );
    CPPUNIT_ASSERT( t.thrust() > 0.5 );
    CPPUNIT_ASSERT( fabs( t.axis().z() ) < 1.e-6 );
  }
}
