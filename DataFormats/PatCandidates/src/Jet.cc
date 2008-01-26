//
// $Id: Jet.cc,v 1.4 2008/01/22 21:58:15 lowette Exp $
//

#include "DataFormats/PatCandidates/interface/Jet.h"


using namespace pat;


/// default constructor
Jet::Jet() :
  PATObject<JetType>(JetType(reco::Particle::LorentzVector(0, 0, 0, 0), reco::Particle::Point(0,0,0), reco::CaloJet::Specific(), reco::Jet::Constituents())),
  partonFlavour_(0), lrPhysicsJetLRval_(-999.), lrPhysicsJetProb_(-1),
  jetCharge_(0.0), associatedTracks_() {
}


/// constructor from a JetType
Jet::Jet(const JetType & aJet) :
  PATObject<JetType>(aJet),
  partonFlavour_(0), lrPhysicsJetLRval_(-999.), lrPhysicsJetProb_(-1) {
}


/// constructor from ref to JetType
Jet::Jet(const edm::RefToBase<JetType> & aJetRef) : PATObject<JetType>(aJetRef) {
}


/// destructor
Jet::~Jet() {
}


/// return the matched generated parton
const reco::Particle * Jet::genParton() const {
  return (genParton_.size() > 0 ? &genParton_.front() : 0);
}


/// return the matched generated jet
const reco::GenJet * Jet::genJet() const {
  return (genJet_.size() > 0 ? &genJet_.front() : 0);
}


/// return the flavour of the parton underlying the jet
int Jet::partonFlavour() const {
  return partonFlavour_;
}


/// return the correction factor to go to a non-calibrated jet
float Jet::noCorrF() const {
  return noCorrF_;
}


/// return the correction factor to go to a uds-calibrated jet
float Jet::udsCorrF() const {
  return udsCorrF_;
}


/// return the correction factor to go to a gluon-calibrated jet
float Jet::gluCorrF() const {
  return gluCorrF_;
}


/// return the correction factor to go to a c-calibrated jet
float Jet::cCorrF() const {
  return cCorrF_;
}


/// return the correction factor to go to a b-calibrated jet
float Jet::bCorrF() const {
  return bCorrF_;
}


/// return the associated non-calibrated jet
JetType Jet::recJet() const {
  JetType recJet(*this);
  recJet.setP4(noCorrF_*this->p4());
  return recJet;
}


/// return the associated non-calibrated jet
Jet Jet::noCorrJet() const {
  Jet jet(*this);
  jet.setP4(noCorrF_*this->p4());
  // fix the factor to uncalibrate for the fact that we change the scale of the actual jet
  jet.setScaleCalibFactors(1., udsCorrF_, gluCorrF_, cCorrF_, bCorrF_);
  return jet;
}


/// return the associated uds-calibrated jet
Jet Jet::udsCorrJet() const {
  Jet jet(*this);
  jet.setP4(udsCorrF_*noCorrF_*this->p4());
  // fix the factor to uncalibrate for the fact that we change the scale of the actual jet
  jet.setScaleCalibFactors(1./udsCorrF_, udsCorrF_, gluCorrF_, cCorrF_, bCorrF_);
  return jet;
}


/// return the associated gluon-calibrated jet
Jet Jet::gluCorrJet() const {
  Jet jet(*this);
  jet.setP4(gluCorrF_*noCorrF_*this->p4());
  // fix the factor to uncalibrate for the fact that we change the scale of the actual jet
  jet.setScaleCalibFactors(1./gluCorrF_, udsCorrF_, gluCorrF_, cCorrF_, bCorrF_);
  return jet;
}


/// return the associated c-calibrated jet
Jet Jet::cCorrJet() const {
  Jet jet(*this);
  jet.setP4(cCorrF_*noCorrF_*this->p4());
  // fix the factor to uncalibrate for the fact that we change the scale of the actual jet
  jet.setScaleCalibFactors(1./cCorrF_, udsCorrF_, gluCorrF_, cCorrF_, bCorrF_);
  return jet;
}


/// return the associated b-calibrated jet
Jet Jet::bCorrJet() const {
  Jet jet(*this);
  // set the corrected 4-vector
  jet.setP4(bCorrF_*noCorrF_*this->p4());
  // fix the factor to uncalibrate for the fact that we change the scale of the actual jet
  jet.setScaleCalibFactors(1./bCorrF_, udsCorrF_, gluCorrF_, cCorrF_, bCorrF_);
  // set the resolutions assuming this jet to be a b-jet
  jet.setResolutionA(bResA_);
  jet.setResolutionB(bResB_);
  jet.setResolutionC(bResC_);
  jet.setResolutionD(bResD_);
  jet.setResolutionET(bResET_);
  jet.setResolutionEta(bResEta_);
  jet.setResolutionPhi(bResPhi_);
  jet.setResolutionTheta(bResTheta_);
  jet.setCovMatrix(bCovM_);
  return jet;
}


/// return the jet calibrated according to the MC flavour truth
Jet Jet::mcFlavCorrJet() const {
  // determine the correction factor to use depending on MC flavour truth
  float corrF = gluCorrF_; // default, also for unidentified flavour
  if (abs(partonFlavour_) == 1 || abs(partonFlavour_) == 2 || abs(partonFlavour_) == 3) corrF = udsCorrF_;
  if (abs(partonFlavour_) == 4) corrF = cCorrF_;
  if (abs(partonFlavour_) == 5) corrF = bCorrF_;
  Jet jet(*this);
  jet.setP4(corrF*noCorrF_*this->p4());
  // fix the factor to uncalibrate for the fact that we change the scale of the actual jet
  jet.setScaleCalibFactors(1./corrF, udsCorrF_, gluCorrF_, cCorrF_, bCorrF_);
  return jet;
}


/// return the jet calibrated with weights assuming W decay
Jet Jet::wCorrJet() const {
  Jet jet(*this);
  // set the corrected 4-vector weighting for the c-content in W decays
  jet.setP4((3*udsCorrF_+cCorrF_)/4*noCorrF_*this->p4());
  // fix the factor to uncalibrate for the fact that we change the scale of the actual jet
  jet.setScaleCalibFactors(4./(3*udsCorrF_+cCorrF_), udsCorrF_, gluCorrF_, cCorrF_, bCorrF_);
  return jet;
}


/// get b discriminant from label name
float Jet::bDiscriminator(std::string theLabel) const {
  float discriminator = -10.;
  if (theLabel == "" || theLabel == "default") theLabel = "trackCountingJetTags";
  for(unsigned int i=0; i!=pairDiscriVector_.size(); i++){
    if(pairDiscriVector_[i].first == theLabel){
      discriminator = pairDiscriVector_[i].second;
    }
  }
  return discriminator;
}


/// get JetTagRef from labael name
reco::JetTagRef Jet::bJetTagRef(std::string theLabel) const {
  reco::JetTagRef theJetTagRef ;
  //if(pairDiscriJetTagRef.size() == 0){
  //  cout << "no JetTagRef found" << endl;
  //}
  for(unsigned int i=0; i!=pairJetTagRefVector_.size(); i++){
    if(pairJetTagRefVector_[i].first == theLabel){
       theJetTagRef= pairJetTagRefVector_[i].second;
    }
  } 
  return theJetTagRef;
}


/// get the value of the i'th jet cleaning variable
float Jet::lrPhysicsJetVar(unsigned int i) const {
  return (i < lrPhysicsJetVarVal_.size() ? lrPhysicsJetVarVal_[i].first  : 0);
}


/// get the likelihood ratio corresponding to the i'th jet cleaning variable
float Jet::lrPhysicsJetVal(unsigned int i) const {
  return (i < lrPhysicsJetVarVal_.size() ? lrPhysicsJetVarVal_[i].second : 1);
}


/// get the overall jet cleaning likelihood ratio
float Jet::lrPhysicsJetLRval() const {
  return lrPhysicsJetLRval_;
}


/// get the overall jet cleaning probability
float Jet::lrPhysicsJetProb() const {
  return lrPhysicsJetProb_;
}


/// method to return the JetCharge computed when creating the Jet
float Jet::jetCharge() const {
  return jetCharge_;
}


/// method to return a vector of refs to the tracks associated to this jet
const reco::TrackRefVector & Jet::associatedTracks() const {
  return associatedTracks_;
}


/// method to set the matched parton
void Jet::setGenParton(const reco::Particle & gp) {
  genParton_.clear();
  genParton_.push_back(gp);
}


/// method to set the matched generated jet
void Jet::setGenJet(const reco::GenJet & gj) {
  genJet_.clear();
  genJet_.push_back(gj);
}


/// method to set the flavour of the parton underlying the jet
void Jet::setPartonFlavour(int partonFl) {
  partonFlavour_ = partonFl;
}


/// method to set the energy scale correction factors
void Jet::setScaleCalibFactors(float noCorrF, float udsCorrF, float gluCorrF, float cCorrF, float bCorrF) {
  noCorrF_ = noCorrF;
  udsCorrF_ = udsCorrF;
  gluCorrF_ = gluCorrF;
  cCorrF_ = cCorrF;
  bCorrF_ = bCorrF;
}


/// method to set the resolutions under the assumption this is a b-jet
void Jet::setBResolutions(float bResET, float bResEta, float bResPhi, float bResA, float bResB, float bResC, float bResD, float bResTheta) {
  bResET_ = bResET;
  bResEta_ = bResEta;
  bResPhi_ = bResPhi;
  bResA_ = bResA;
  bResB_ = bResB;
  bResC_ = bResC;
  bResD_ = bResD;
  bResTheta_ = bResTheta;
}


/// method to add a algolabel-discriminator pair
void Jet::addBDiscriminatorPair(std::pair<std::string, float> & thePair) {
  pairDiscriVector_.push_back(thePair);
}


/// method to add a algolabel-jettagref pair
void Jet::addBJetTagRefPair(std::pair<std::string, reco::JetTagRef> & thePair) {
  pairJetTagRefVector_.push_back(thePair);
}


/// method to set all jet cleaning variable + LR pairs
void Jet::setLRPhysicsJetVarVal(const std::vector<std::pair<float, float> > & varValVec) {
  for (size_t i = 0; i<varValVec.size(); i++) lrPhysicsJetVarVal_.push_back(varValVec[i]);
}


/// method to set the combined jet cleaning likelihood ratio value
void Jet::setLRPhysicsJetLRval(float clr) {
  lrPhysicsJetLRval_ = clr;
}


/// method to set the jet cleaning probability
void Jet::setLRPhysicsJetProb(float plr) {
  lrPhysicsJetProb_ = plr;
}


/// method to set the jet charge
void Jet::setJetCharge(float jetCharge) {
  jetCharge_ = jetCharge;
}
