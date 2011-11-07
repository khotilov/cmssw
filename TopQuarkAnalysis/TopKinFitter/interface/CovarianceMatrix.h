#ifndef CovarianceMatrix_h
#define CovarianceMatrix_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TopQuarkAnalysis/TopKinFitter/interface/TopKinFitter.h"
#include "TopQuarkAnalysis/TopObjectResolutions/interface/MET.h"
#include "TopQuarkAnalysis/TopObjectResolutions/interface/Jet.h"
#include "TopQuarkAnalysis/TopObjectResolutions/interface/Muon.h"
#include "TopQuarkAnalysis/TopObjectResolutions/interface/Electron.h"

/*
  \class   CovarianceMatrix CovarianceMatrix.h "TopQuarkAnalysis/TopKinFitter/interface/CovarianceMatrix.h"
  
  \brief   Helper class used to setup covariance matrices for given objects and known resolutions

  More details to be added here...
  
**/

class CovarianceMatrix {

 public:

  enum ObjectType{kUdscJet, kBJet, kMuon, kElectron, kMet};
  
  /// default constructor
  CovarianceMatrix(){};
  /// constructor for the fully-hadronic channel
  CovarianceMatrix(const std::vector<edm::ParameterSet> udscResolutions, const std::vector<edm::ParameterSet> bResolutions);
  /// constructor for the lepton+jets channel
  CovarianceMatrix(const std::vector<edm::ParameterSet> udscResolutions, const std::vector<edm::ParameterSet> bResolutions,
		   const std::vector<edm::ParameterSet> lepResolutions, const std::vector<edm::ParameterSet> metResolutions);
  // destructor
  ~CovarianceMatrix(){};

  /// return covariance matrix for a PAT object
  template <class T>
    TMatrixD setupMatrix(const pat::PATObject<T>& object, const TopKinFitter::Param param, const std::string resolutionProvider = "");
  /// return covariance matrix for a plain 4-vector
  TMatrixD setupMatrix(const TLorentzVector& object, const ObjectType objType, const TopKinFitter::Param param);
  /// get resolution for a given component of an object
  double getResolution(const TLorentzVector& object, const ObjectType objType, const std::string whichResolution = "");

 private:

  /// vector of strings for the binning of the resolutions
  std::vector<std::string> binsUdsc_, binsB_, binsLep_, binsMet_;
  /// vectors for the resolution functions
  std::vector<std::string> funcEtUdsc_ , funcEtB_ , funcEtLep_ , funcEtMet_;
  std::vector<std::string> funcEtaUdsc_, funcEtaB_, funcEtaLep_, funcEtaMet_;
  std::vector<std::string> funcPhiUdsc_, funcPhiB_, funcPhiLep_, funcPhiMet_;

};

template <class T>
TMatrixD CovarianceMatrix::setupMatrix(const pat::PATObject<T>& object, TopKinFitter::Param param, std::string resolutionProvider)
{
  // This part is for pat objects with resolutions embedded
  if(object.hasKinResolution()) {
    TMatrixD CovM3 (3,3); CovM3.Zero();
    TMatrixD CovM4 (4,4); CovM4.Zero();
    TMatrixD* CovM = &CovM3;
    switch(param){
    case TopKinFitter::kEtEtaPhi :
      CovM3(0,0) = pow(object.resolEt(resolutionProvider) , 2);
      if( dynamic_cast<const reco::MET*>(&object) )
	CovM3(1,1) = pow(9999., 2);
      else
	CovM3(1,1) = pow(object.resolEta(resolutionProvider), 2);
      CovM3(2,2) = pow(object.resolPhi(resolutionProvider), 2);
      CovM = &CovM3;
      break;
    case TopKinFitter::kEtThetaPhi :
      CovM3(0,0) = pow(object.resolEt(resolutionProvider)   , 2);
      CovM3(1,1) = pow(object.resolTheta(resolutionProvider), 2);
      CovM3(2,2) = pow(object.resolPhi(resolutionProvider)  , 2);
      CovM = &CovM3;
      break;
    case TopKinFitter::kEMom :
      CovM4(0,0) = pow(1, 2);
      CovM4(1,1) = pow(1, 2);
      CovM4(2,2) = pow(1, 2);
      CovM4(3,3) = pow(1, 2);
      CovM = &CovM4;
      break;
    }
    return *CovM;
  }
  // This part is for objects without resolutions embedded
  else {
    ObjectType objType;
    // jets
    if( dynamic_cast<const reco::Jet*>(&object) ) {
      if(resolutionProvider == "bjets")
	objType = kBJet;
      else
	objType = kUdscJet;
    }
    // muons
    else if( dynamic_cast<const reco::Muon*>(&object) )
      objType = kMuon;
    // electrons
    else if( dynamic_cast<const reco::GsfElectron*>(&object) ) 
      objType = kElectron;
     // MET
    else if( dynamic_cast<const reco::MET*>(&object) )
      objType = kMet;
    const TLorentzVector p4(object.px(), object.py(), object.pz(), object.energy());
    return setupMatrix(p4, objType, param);
  }
}

#endif
