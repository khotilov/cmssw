//
// $Id: MET.cc,v 1.11 2008/10/07 18:22:18 lowette Exp $
//

#include "DataFormats/PatCandidates/interface/MET.h"


using namespace pat;


/// default constructor
MET::MET(): uncorInfo_(0) {
}


/// constructor from METType
MET::MET(const METType & aMET) : PATObject<METType>(aMET), uncorInfo_(0) {
    const reco::CaloMET * calo = dynamic_cast<const reco::CaloMET *>(&aMET);
    if (calo != 0) caloMET_.push_back(calo->getSpecific());
}


/// constructor from ref to METType
MET::MET(const edm::RefToBase<METType> & aMETRef) : PATObject<METType>(aMETRef), uncorInfo_(0) {
    const reco::CaloMET * calo = dynamic_cast<const reco::CaloMET *>(aMETRef.get());
    if (calo != 0) caloMET_.push_back(calo->getSpecific());
}

/// constructor from ref to METType
MET::MET(const edm::Ptr<METType> & aMETRef) : PATObject<METType>(aMETRef), uncorInfo_(0) {
    const reco::CaloMET * calo = dynamic_cast<const reco::CaloMET *>(aMETRef.get());
    if (calo != 0) caloMET_.push_back(calo->getSpecific());
}


/// destructor
MET::~MET() {
}


/// return the generated MET from neutrinos
const reco::GenMET * MET::genMET() const {
  return (genMET_.size() > 0 ? &genMET_.front() : 0 );
}

/// method to set the generated MET
void MET::setGenMET(const reco::GenMET & gm) {
  genMET_.clear();
  genMET_.push_back(gm);
}

//! return uncorrrection related stuff
unsigned int MET::nCorrections() const { checkUncor_(); return nCorrections_; }

float MET::corEx(UncorrectionType ix) const { checkUncor_(); return uncorInfo_[ix].corEx; }
float MET::corEy(UncorrectionType ix) const { checkUncor_(); return uncorInfo_[ix].corEy; }
float MET::corSumEt(UncorrectionType ix) const { checkUncor_(); return uncorInfo_[ix].corSumEt; }
float MET::uncorrectedPt(UncorrectionType ix) const { checkUncor_(); return uncorInfo_[ix].pt; }
float MET::uncorrectedPhi(UncorrectionType ix) const { checkUncor_(); return uncorInfo_[ix].phi; }


//! check and set transients
void MET::checkUncor_() const {
  if (uncorInfo_.size() == uncorrMAXN && oldPt_ == pt() ) return;

  oldPt_ = pt();
  std::vector<CorrMETData> corrs(mEtCorr());
  nCorrections_ = corrs.size();

  uncorInfo_.resize(uncorrMAXN);
  UncorrectionType ix;

  //! ugly
  //! ALL
  ix = uncorrALL;
  uncorInfo_[ix] = UncorInfo();
  for (unsigned int iC=0; iC < nCorrections_; ++iC){
    uncorInfo_[ix].corEx +=    corrs[iC].mex;
    uncorInfo_[ix].corEy +=    corrs[iC].mey;
    uncorInfo_[ix].corSumEt += corrs[iC].sumet;
  }
  setPtPhi_(uncorInfo_[ix]);

  //! JES
  ix = uncorrJES;
  uncorInfo_[ix] = UncorInfo();
  if (nCorrections_ >=1 ){
    unsigned int iC = 0;
    uncorInfo_[ix].corEx +=    corrs[iC].mex;
    uncorInfo_[ix].corEy +=    corrs[iC].mey;
    uncorInfo_[ix].corSumEt += corrs[iC].sumet;
  }
  setPtPhi_(uncorInfo_[ix]);

  //! MUON
  ix = uncorrMUON;
  uncorInfo_[ix] = UncorInfo();
  if (nCorrections_ >=2 ){
    unsigned int iC = 1;
    uncorInfo_[ix].corEx +=    corrs[iC].mex;
    uncorInfo_[ix].corEy +=    corrs[iC].mey;
    uncorInfo_[ix].corSumEt += corrs[iC].sumet;
  }
  setPtPhi_(uncorInfo_[ix]);

}

void MET::setPtPhi_(UncorInfo& uci) const {
  float lpx = px() - uci.corEx;
  float lpy = py() - uci.corEy;
  uci.pt = sqrt(lpx*lpx + lpy*lpy);
  uci.phi = atan2(lpy, lpx);  
}
