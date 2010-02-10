#include <cppunit/extensions/HelperMacros.h>
#include "CommonTools/Utils/interface/expressionParser.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <iostream>
#include <Reflex/Object.h>
#include <Reflex/Type.h>
#include <typeinfo>
#include "DataFormats/Common/interface/TestHandle.h"

class testExpressionParser : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(testExpressionParser);
  CPPUNIT_TEST(checkAll);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}
  void checkAll(); 
  void checkTrack(const std::string &, double);
  void checkCandidate(const std::string &, double, bool lazy=false);
  void checkJet(const std::string &, double);
  void checkMuon(const std::string &, double);
  reco::Track trk;
  reco::CompositeCandidate cand;
  ROOT::Reflex::Object o;
  reco::parser::ExpressionPtr expr;
  pat::Jet jet;
  pat::Muon muon;
};

CPPUNIT_TEST_SUITE_REGISTRATION(testExpressionParser);

void testExpressionParser::checkTrack(const std::string & expression, double x) {
  for (int lazyMode = 0; lazyMode <= 1; ++lazyMode) {
  std::cerr << "checking " << (lazyMode ? "lazy " : "") << "expression: \"" << expression << "\"" << std::flush; 
  expr.reset();
  CPPUNIT_ASSERT(reco::parser::expressionParser<reco::Track>(expression, expr, lazyMode));
  CPPUNIT_ASSERT(expr.get() != 0);
  double res = expr->value(o);
  StringObjectFunction<reco::Track> f(expression, lazyMode);
  CPPUNIT_ASSERT(fabs(f(trk) - res) < 1.e-6);
  CPPUNIT_ASSERT(fabs(f(trk) - x) < 1.e-6);
  std::cerr << " = " << res << std::endl;
  }
}

void testExpressionParser::checkCandidate(const std::string & expression, double x, bool lazyMode) {
  std::cerr << "checking " << (lazyMode ? "lazy " : "") << "expression: \"" << expression << "\"" << std::flush; 
  expr.reset();
  CPPUNIT_ASSERT(reco::parser::expressionParser<reco::Candidate>(expression, expr, lazyMode));
  CPPUNIT_ASSERT(expr.get() != 0);
  double res = 0;
  CPPUNIT_ASSERT_NO_THROW(res = expr->value(o));
  StringObjectFunction<reco::Candidate> f(expression, lazyMode);
  CPPUNIT_ASSERT(fabs(f(cand) - res) < 1.e-6);
  CPPUNIT_ASSERT(fabs(f(cand) - x) < 1.e-6);
  std::cerr << " = " << res << std::endl;
}

void testExpressionParser::checkJet(const std::string & expression, double x) {
  for (int lazyMode = 0; lazyMode <= 1; ++lazyMode) {
  std::cerr << "checking " << (lazyMode ? "lazy " : "") << "expression: \"" << expression << "\"" << std::flush; 
  expr.reset();
  CPPUNIT_ASSERT(reco::parser::expressionParser<pat::Jet>(expression, expr, lazyMode));
  CPPUNIT_ASSERT(expr.get() != 0);
  double res = expr->value(o);
  StringObjectFunction<pat::Jet> f(expression, lazyMode);
  CPPUNIT_ASSERT(fabs(f(jet) - res) < 1.e-6);
  CPPUNIT_ASSERT(fabs(f(jet) - x) < 1.e-6);
  std::cerr << " = " << res << std::endl;
  }
}

void testExpressionParser::checkMuon(const std::string & expression, double x) {
  for (int lazyMode = 0; lazyMode <= 1; ++lazyMode) {
  std::cerr << "checking " << (lazyMode ? "lazy " : "") << "expression: \"" << expression << "\"" << std::flush; 
  expr.reset();
  CPPUNIT_ASSERT(reco::parser::expressionParser<pat::Muon>(expression, expr, lazyMode));
  CPPUNIT_ASSERT(expr.get() != 0);
  double res = expr->value(o);
  StringObjectFunction<pat::Muon> f(expression, lazyMode);
  std::cerr << " = " << x << " (reference), " << res << " (bare), " << f(muon) << " (full)" << std::endl;
  CPPUNIT_ASSERT(fabs(f(muon) - res) < 1.e-6);
  CPPUNIT_ASSERT(fabs(f(muon) - x) < 1.e-6);
  }
}


void testExpressionParser::checkAll() {
  using namespace reco;
  using namespace ROOT::Reflex;
  const double chi2 = 20.0;
  const int ndof = 10;
  reco::Track::Point v(1, 2, 3);
  reco::Track::Vector p(5, 3, 10);
  double e[] = { 1.1,
                 1.2, 2.2,
                 1.3, 2.3, 3.3,
                 1.4, 2.4, 3.4, 4.4,
                 1.5, 2.5, 3.5, 4.5, 5.5 };
  reco::TrackBase::CovarianceMatrix cov(e, e + 15);
  cov(0, 0) = 1.2345;
  cov(1, 0) = 123.456;
  cov(1, 1) = 5.4321;
  trk = reco::Track(chi2, ndof, v, p, -1, cov);

  edm::ProductID const pid(1);
  reco::TrackExtraCollection trkExtras;
  reco::Track::Point outerV(100, 200, 300), innerV(v);
  reco::Track::Vector outerP(0.5, 3.5, 10.5), innerP(p);
  reco::Track::CovarianceMatrix  outerC, innerC;
  unsigned int outerId = 123, innerId = 456;
  reco::TrackExtra trkExtra(outerV, outerP, true, innerV, innerP, true,
			    outerC, outerId, innerC, innerId,
			    anyDirection);
  trkExtras.push_back(trkExtra);
  edm::TestHandle<reco::TrackExtraCollection> h(&trkExtras, pid);
  reco::TrackExtraRef trkExtraRef(h, 0);
  trk.setExtra(trkExtraRef);
  trk.setAlgorithm(reco::Track::iter2);
  {
    ROOT::Reflex::Type t = ROOT::Reflex::Type::ByTypeInfo(typeid(reco::Track));
    o = ROOT::Reflex::Object(t, & trk);
    checkTrack("pt", trk.pt());
    checkTrack("charge", trk.charge());
    checkTrack("pt/3", trk.pt()/3);
    checkTrack("covariance(0, 0)", trk.covariance(0, 0));
    checkTrack("covariance(1, 0)", trk.covariance(1, 0));
    checkTrack("covariance(1, 1)", trk.covariance(1, 1));
    checkTrack("momentum.x", trk.momentum().x());
    checkTrack("hitPattern.numberOfValidHits", trk.hitPattern().numberOfValidHits());
    checkTrack("extra.outerPhi", trk.extra()->outerPhi());
    checkTrack("referencePoint.R", trk.referencePoint().R());
    checkTrack("algo", reco::Track::iter2);
  }
  reco::Candidate::LorentzVector p1(1, 2, 3, 4);
  reco::Candidate::LorentzVector p2(1.1, 2.2, 3.3, 4.4);
  reco::LeafCandidate c1(+1, p1);
  reco::LeafCandidate c2(-1, p2);
  cand.addDaughter(c1);
  cand.addDaughter(c2);
  CPPUNIT_ASSERT(cand.numberOfDaughters()==2);
  CPPUNIT_ASSERT(cand.daughter(0)!=0);
  CPPUNIT_ASSERT(cand.daughter(1)!=0);
  {
    ROOT::Reflex::Type t = ROOT::Reflex::Type::ByTypeInfo(typeid(reco::Candidate));
    o = ROOT::Reflex::Object(t, & cand);  
    // these can be checked in both modes
    for (int lazyMode = 0; lazyMode <= 1; ++lazyMode) {
    checkCandidate("numberOfDaughters", cand.numberOfDaughters(), lazyMode);
    checkCandidate("daughter(0).isStandAloneMuon", cand.daughter(0)->isStandAloneMuon(), lazyMode);  
    checkCandidate("daughter(1).isStandAloneMuon", cand.daughter(1)->isStandAloneMuon(), lazyMode);  
    checkCandidate("daughter(0).pt", cand.daughter(0)->pt(), lazyMode);
    checkCandidate("daughter(1).pt", cand.daughter(1)->pt(), lazyMode);
    checkCandidate("min(daughter(0).pt, daughter(1).pt)", std::min(cand.daughter(0)->pt(), cand.daughter(1)->pt()), lazyMode);
    checkCandidate("max(daughter(0).pt, daughter(1).pt)", std::max(cand.daughter(0)->pt(), cand.daughter(1)->pt()), lazyMode);
    }
    // these can be checked only in lazy mode
    checkCandidate("name.empty()", true, true);
    checkCandidate("roles.size()", 0, true);

  }

  std::vector<reco::LeafCandidate> cands;
  cands.push_back(c1);  cands.push_back(c2); 
  edm::TestHandle<std::vector<reco::LeafCandidate> > constituentsHandle(&cands, edm::ProductID(42));
  reco::Jet::Constituents constituents;
  constituents.push_back( reco::Jet::Constituent(constituentsHandle, 0) );
  constituents.push_back( reco::Jet::Constituent(constituentsHandle, 1) );
  reco::CaloJet::Specific caloSpecific; caloSpecific.mMaxEInEmTowers = 0.5;
  jet = pat::Jet(reco::CaloJet(p1+p2, reco::Jet::Point(), caloSpecific, constituents));
  CPPUNIT_ASSERT(jet.nConstituents() == 2);
  CPPUNIT_ASSERT(jet.nCarrying(1.0)  == 2);
  CPPUNIT_ASSERT(jet.nCarrying(0.1)  == 1);
  {
    ROOT::Reflex::Type t = ROOT::Reflex::Type::ByTypeInfo(typeid(pat::Jet));
    o = ROOT::Reflex::Object(t, & jet);
    checkJet("nCarrying(1.0)", jet.nCarrying(1.0));
    checkJet("nCarrying(0.1)", jet.nCarrying(0.1));
  }
  {
    // Check a method that starts with the same string as a predefined function
    CPPUNIT_ASSERT(jet.maxEInEmTowers()== 0.5);
    checkJet("maxEInEmTowers", jet.maxEInEmTowers());
  }

  std::pair<std::string,float> bp;
  bp = std::pair<std::string,float>("aaa", 1.0); jet.addBDiscriminatorPair(bp); 
  bp = std::pair<std::string,float>("b c", 2.0); jet.addBDiscriminatorPair(bp); 
  bp = std::pair<std::string,float>("d " , 3.0); jet.addBDiscriminatorPair(bp); 
  pat::JetCorrFactors::FlavourCorrections dummyfc;
  pat::JetCorrFactors jcfs("a",2.f,2.f,2.f,1.f,dummyfc,dummyfc,dummyfc,std::vector<float>());
  jet.setCorrFactors(jcfs);
  jet.setCorrStep(pat::JetCorrFactors::L3);
  CPPUNIT_ASSERT(jet.bDiscriminator("aaa") == 1.0);
  CPPUNIT_ASSERT(jet.bDiscriminator("b c") == 2.0);
  CPPUNIT_ASSERT(jet.bDiscriminator("d ")  == 3.0);
  CPPUNIT_ASSERT(jet.getPairDiscri().size() == 3 );
  {
    ROOT::Reflex::Type t = ROOT::Reflex::Type::ByTypeInfo(typeid(pat::Jet));
    o = ROOT::Reflex::Object(t, & jet);
    checkJet("bDiscriminator(\"aaa\")", jet.bDiscriminator("aaa"));
    checkJet("bDiscriminator('aaa')"  , jet.bDiscriminator("aaa"));
    checkJet("bDiscriminator(\"b c\")", jet.bDiscriminator("b c"));
    checkJet("bDiscriminator(\"d \")" , jet.bDiscriminator("d " ));
    checkJet("correctedJet('RAW').pt",  jet.correctedJet("RAW").pt());
    checkJet("correctedJet('Raw').pt",  jet.correctedJet(pat::JetCorrFactors::Raw).pt());
  }

  {
     ROOT::Reflex::Type t = ROOT::Reflex::Type::ByTypeInfo(typeid(pat::Jet));
     std::vector<reco::SecondaryVertexTagInfo::IndexedTrackData> trackData;
     std::vector<reco::SecondaryVertexTagInfo::VertexData> vertexData(1);
     
     reco::SoftLeptonTagInfo dummyInfo;
     reco::SoftLeptonProperties props;
     props.quality = 10;
     dummyInfo.insert(edm::RefToBase<reco::Track>(), props);
     edm::Ptr<reco::BaseTagInfo> ptrDummyInfo(edm::ProductID(1),&dummyInfo,0);
     jet.addTagInfo("dummy", ptrDummyInfo);
     o = ROOT::Reflex::Object(t, & jet);
     checkJet("tagInfoSoftLepton.properties(0).quality",jet.tagInfoSoftLepton()->properties(0).quality);
  }
  muon = pat::Muon(reco::Muon(+1, p1+p2));
  muon.setUserIso(2.0);
  muon.setUserIso(42.0, 1);
  CPPUNIT_ASSERT( muon.userIso()  == 2.0 );
  CPPUNIT_ASSERT( muon.userIso(0) == 2.0 );
  CPPUNIT_ASSERT( muon.userIso(1) == 42.0 );
  reco::IsoDeposit dep; dep.addCandEnergy(2.0);
  CPPUNIT_ASSERT(dep.candEnergy() == 2.0);
  muon.setIsoDeposit(pat::TrackIso, dep);
  CPPUNIT_ASSERT(muon.trackIsoDeposit() != 0);
  CPPUNIT_ASSERT(muon.trackIsoDeposit()->candEnergy() == 2.0);
  {
    ROOT::Reflex::Type t = ROOT::Reflex::Type::ByTypeInfo(typeid(pat::Muon));
    o = ROOT::Reflex::Object(t, & muon);
    checkMuon("userIso"    , muon.userIso() );
    checkMuon("userIso()"  , muon.userIso() );
    checkMuon("userIso(0)" , muon.userIso(0));
    checkMuon("userIso(1)" , muon.userIso(1));
    checkMuon("trackIsoDeposit.candEnergy", muon.trackIsoDeposit()->candEnergy());
  }
}
