#include <cppunit/extensions/HelperMacros.h>
#include "CommonTools/Utils/interface/cutParser.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include <iostream>
#include <Reflex/Object.h>
#include <Reflex/Type.h>
#include <typeinfo>

class testCutParser : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(testCutParser);
  CPPUNIT_TEST(checkAll);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}
  void checkAll(); 
  void check(const std::string &, bool);
  void checkHit(const std::string &, bool, const SiStripRecHit2D &);
  void checkMuon(const std::string &, bool, const reco::Muon &);
  reco::Track trk;
  SiStripRecHit2D hitOk, hitThrow;
  Reflex::Object o;
  reco::parser::SelectorPtr sel;
};

CPPUNIT_TEST_SUITE_REGISTRATION(testCutParser);

void testCutParser::check(const std::string & cut, bool res) {
  for (int lazy = 0; lazy <= 1; ++lazy) {
  std::cerr << "parsing " << (lazy ? "lazy " : "") << "cut: \"" << cut << "\"" << std::endl;
  sel.reset();
  CPPUNIT_ASSERT(reco::parser::cutParser<reco::Track>(cut, sel, lazy));
  CPPUNIT_ASSERT((*sel)(o) == res);
  StringCutObjectSelector<reco::Track> select(cut, lazy);
  CPPUNIT_ASSERT(select(trk) == res);
  }
}

void testCutParser::checkHit(const std::string & cut, bool res, const SiStripRecHit2D &hit) {
  Reflex::Type t = Reflex::Type::ByTypeInfo(typeid(SiStripRecHit2D));
  o = Reflex::Object(t, const_cast<void *>(static_cast<const void *>(&hit)));
  for (int lazy = 0; lazy <= 1; ++lazy) {
  std::cerr << "parsing " << (lazy ? "lazy " : "") << "cut: \"" << cut << "\"" << std::endl;
  sel.reset();
  CPPUNIT_ASSERT(reco::parser::cutParser<SiStripRecHit2D>(cut, sel, lazy));
  CPPUNIT_ASSERT((*sel)(o) == res);
  StringCutObjectSelector<SiStripRecHit2D> select(cut, lazy);
  CPPUNIT_ASSERT(select(hit) == res);
  }
}

void testCutParser::checkMuon(const std::string & cut, bool res, const reco::Muon &mu) {
  Reflex::Type t = Reflex::Type::ByTypeInfo(typeid(reco::Muon));
  o = Reflex::Object(t, const_cast<void *>(static_cast<const void *>(&mu)));
  sel.reset();
  for (int lazy = 0; lazy <= 1; ++lazy) {
  std::cerr << "parsing " << (lazy ? "lazy " : "") << "cut: \"" << cut << "\"" << std::endl;
  CPPUNIT_ASSERT(reco::parser::cutParser<reco::Muon>(cut, sel, lazy));
  CPPUNIT_ASSERT((*sel)(o) == res);
  StringCutObjectSelector<reco::Muon> select(cut, lazy);
  CPPUNIT_ASSERT(select(mu) == res);
  }
}



void testCutParser::checkAll() {
  using namespace reco;
  using namespace Reflex;
  const double chi2 = 20.0;
  const int ndof = 10;
  reco::Track::Point v(1, 2, 3);
  reco::Track::Vector p(0, 3, 10);
  double e[] = { 1.1,
                 1.2, 2.2,
                 1.3, 2.3, 3.3,
                 1.4, 2.4, 3.4, 4.4,
                 1.5, 2.5, 3.5, 4.5, 5.5 };
  reco::TrackBase::CovarianceMatrix cov(e, e + 15);
  trk = reco::Track(chi2, ndof, v, p, -1, cov);

  hitOk = SiStripRecHit2D(LocalPoint(1,1), LocalError(1,1,1), 0, SiStripRecHit2D::ClusterRef());

  Reflex::Type t = Reflex::Type::ByTypeInfo(typeid(reco::Track));
  o = Reflex::Object(t, & trk);

  std::cerr << "Track pt: " << trk.pt() << std::endl;
  // note: pt = 3, charge = -1
  check( "", true );
  check( "  ", true );
  check("pt", true);
  check("px", false);
  check( "pt > 2", true );
  check( "charge < 0", true );
  check( "pt < 2", false );
  check( "pt >= 2", true );
  check( "pt <= 2", false );
  check( "pt = 3", true );
  check( "pt == 3", true );
  check( "pt != 3", false );
  check( "! pt == 3", false );
  check( "2.9 < pt < 3.1", true );
  check( "pt > 2 & charge < 0", true );
  check( "pt > 2 && charge < 0", true );
  check( "pt < 2 & charge < 0", false );
  check( "pt > 2 & charge > 0", false );
  check( "pt < 2 & charge > 0", false );
  check( "pt > 2 || charge > 0", true );
  check( "pt > 2 | charge > 0", true );
  check( "pt > 2 | charge < 0", true );
  check( "pt < 2 | charge < 0", true );
  check( "pt < 2 | charge > 0", false );
  check( "pt > 2 | charge > 0 | pt < 2", true );
  check( "pt > 2 | charge > 0 | (pt < 2 && charge < 0)", true );
  check( "pt > 2 | charge < 0 | (pt < 2 && charge < 0)", true );
  check( "pt > 2 | charge > 0 | (pt < 2 && charge > 0)", true );
  check( "(pt) > 2", true );
  check( "-pt < -2", true );
  check( "3.9 < pt + 1 < 4.1", true );
  check( "1.9 < pt - 1 < 2.1", true );
  check( "5.9 < 2 * pt < 6.1", true );
  check( "0.9 < pt / 3 < 1.1", true );
  check( "8.9 < pt ^ 2 < 9.1", true );
  check( "26.9 < 3 * pt ^ 2 < 27.1", true );
  check( "27.9 < 3 * pt ^ 2 + 1 < 28.1", true );
  check( " 0.99 < sin( phi ) < 1.01", true );
  check( " -0.01 < cos( phi ) < 0.01", true );
  check( " 8.9 < pow( pt, 2 ) < 9.1", true );
  check( "( 0.99 < sin( phi ) < 1.01 ) & ( -0.01 < cos( phi ) < 0.01 )", true );
  check( "( 3.9 < pt + 1 < 4.1 ) | ( pt < 2 )", true );
  check( " pt = 3 &  pt > 2 | pt < 2", true );
  check( "( ( pt = 3 &  pt > 2 ) | pt < 2 ) & 26.9 < 3 * pt ^ 2 < 27.1", true );
  check( "! pt > 2", false );
  check( "! pt < 2", true );
  check( "! (( 0.99 < sin( phi ) < 1.01 ) & ( -0.01 < cos( phi ) < 0.01 ))", false );
  check( "pt && pt > 1",true);
  // check trailing space
  check( "pt > 2 ", true );

  // check handling of errors 
  //   first those who are the same in lazy and non lazy parsing
  for (int lazy = 0; lazy <= 1; ++lazy) {
  sel.reset();
  CPPUNIT_ASSERT(!reco::parser::cutParser<reco::Track>("1abc",sel, lazy));
  sel.reset();
  CPPUNIT_ASSERT_THROW(reco::parser::cutParser<reco::Track>("(pt < 1",sel, lazy), edm::Exception);
  sel.reset();
  CPPUNIT_ASSERT_THROW(reco::parser::cutParser<reco::Track>("pt < #",sel, lazy), edm::Exception);
  sel.reset();
  CPPUNIT_ASSERT_THROW(reco::parser::cutParser<reco::Track>("pt <> 5",sel, lazy), edm::Exception);
  sel.reset();

  // FIXME: these three pass instead of failing - except the first one which throws 'unknown method "cos"' in non-lazy parsing
  //std::cerr << "============ BEGIN " << (lazy ? "LAZY ":"") << "===================" << std::endl;
  //CPPUNIT_ASSERT_THROW(reco::parser::cutParser<reco::Track>("cos pt < .5",sel, lazy), edm::Exception);
  //CPPUNIT_ASSERT_THROW(reco::parser::cutParser<reco::Track>("pt() pt < .5",sel, lazy), edm::Exception);
  //CPPUNIT_ASSERT_THROW(reco::parser::cutParser<reco::Track>("pt pt < .5",sel, lazy), edm::Exception);
  //std::cerr << "============ END ===================" << std::endl;

  sel.reset();
  CPPUNIT_ASSERT_THROW(reco::parser::cutParser<reco::Track>("cos( pt < .5",sel, lazy), edm::Exception);
  sel.reset();
  CPPUNIT_ASSERT_THROW(reco::parser::cutParser<reco::Track>(" 2 * (pt + 1 < .5",sel, lazy), edm::Exception);
  }
  // then those which are specific to non lazy parsing
  sel.reset();
  CPPUNIT_ASSERT_THROW(reco::parser::cutParser<reco::Track>("doesNotExist < 1",sel), edm::Exception);
  sel.reset();
  CPPUNIT_ASSERT_THROW(reco::parser::cutParser<reco::Track>("pt.pt < 1",sel), edm::Exception);
  // and afterwards check that in lazy parsing they throw at runtime
  sel.reset();
  CPPUNIT_ASSERT(reco::parser::cutParser<reco::Track>("doesNotExist < 1",sel,true)); // with lazy parsing this doesn't throw
  CPPUNIT_ASSERT_THROW((*sel)(o), reco::parser::Exception);                          // but it throws here!
  sel.reset();
  CPPUNIT_ASSERT(reco::parser::cutParser<reco::Track>("pt.pt < 1",sel, true));  // same for this
  CPPUNIT_ASSERT_THROW((*sel)(o), reco::parser::Exception);                     // it throws wen called


  // check hits (for re-implemented virtual functions and exception handling)
  CPPUNIT_ASSERT(hitOk.hasPositionAndError());
  checkHit( "hasPositionAndError" , true, hitOk );
  CPPUNIT_ASSERT(!hitThrow.hasPositionAndError());
  checkHit( "hasPositionAndError" , false, hitThrow );
  CPPUNIT_ASSERT(hitOk.localPosition().x() == 1);
  checkHit( ".99 < localPosition.x < 1.01", true, hitOk);
  CPPUNIT_ASSERT_THROW(hitThrow.localPosition().x(), cms::Exception);
  CPPUNIT_ASSERT_THROW( checkHit(".99 < localPosition.x < 1.01", true, hitThrow) , cms::Exception);

  // check underscores
  checkHit("cluster_regional.isNull()",true,hitOk);
  checkHit("cluster.isNull()",true,hitOk);
  checkHit("cluster_regional.isNonnull()",false,hitOk);
  checkHit("cluster.isNonnull()",false,hitOk);  

  // check short cirtcuit logics
  CPPUNIT_ASSERT( hitOk.hasPositionAndError() && (hitOk.localPosition().x() == 1) );
  CPPUNIT_ASSERT( !( hitThrow.hasPositionAndError() && (hitThrow.localPosition().x() == 1) ) );
  checkHit( "hasPositionAndError && (localPosition.x = 1)", true,  hitOk    );
  checkHit( "hasPositionAndError && (localPosition.x = 1)", false, hitThrow );
  CPPUNIT_ASSERT( (!hitOk.hasPositionAndError()   ) || (hitOk.localPosition().x()    == 1) );
  CPPUNIT_ASSERT( (!hitThrow.hasPositionAndError()) || (hitThrow.localPosition().x() == 1) );
  checkHit( "!hasPositionAndError || (localPosition.x = 1)", true,  hitOk    );
  checkHit( "!hasPositionAndError || (localPosition.x = 1)", true, hitThrow );

}
