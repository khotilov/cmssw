#include "TopQuarkAnalysis/TopKinFitter/interface/CovarianceMatrix.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

CovarianceMatrix::CovarianceMatrix(const std::vector<edm::ParameterSet> udscResolutions, const std::vector<edm::ParameterSet> bResolutions){
  for(std::vector<edm::ParameterSet>::const_iterator iSet = udscResolutions.begin(); iSet != udscResolutions.end(); ++iSet){
    if(iSet->exists("bin")) binsUdsc_.push_back(iSet->getParameter<std::string>("bin"));
    else if(udscResolutions.size()==1) binsUdsc_.push_back("");
    else throw cms::Exception("WrongConfig") << "Parameter 'bin' is needed if more than one PSet is specified!\n";
    
    funcEtUdsc_.push_back(iSet->getParameter<std::string>("et"));
    funcEtaUdsc_.push_back(iSet->getParameter<std::string>("eta"));
    funcPhiUdsc_.push_back(iSet->getParameter<std::string>("phi"));
  }
  for(std::vector<edm::ParameterSet>::const_iterator iSet = bResolutions.begin(); iSet != bResolutions.end(); ++iSet){
    if(iSet->exists("bin")) binsB_.push_back(iSet->getParameter<std::string>("bin"));
    else if(bResolutions.size()==1) binsB_.push_back("");
    else throw cms::Exception("WrongConfig") << "Parameter 'bin' is needed if more than one PSet is specified!\n";
    
    funcEtB_.push_back(iSet->getParameter<std::string>("et"));
    funcEtaB_.push_back(iSet->getParameter<std::string>("eta"));
    funcPhiB_.push_back(iSet->getParameter<std::string>("phi"));
  }
}

CovarianceMatrix::CovarianceMatrix(const std::vector<edm::ParameterSet> udscResolutions, const std::vector<edm::ParameterSet> bResolutions,
				   const std::vector<edm::ParameterSet> lepResolutions, const std::vector<edm::ParameterSet> metResolutions){
  for(std::vector<edm::ParameterSet>::const_iterator iSet = udscResolutions.begin(); iSet != udscResolutions.end(); ++iSet){
    if(iSet->exists("bin")) binsUdsc_.push_back(iSet->getParameter<std::string>("bin"));
    else if(udscResolutions.size()==1) binsUdsc_.push_back("");
    else throw cms::Exception("WrongConfig") << "Parameter 'bin' is needed if more than one PSet is specified!\n";
    
    funcEtUdsc_.push_back(iSet->getParameter<std::string>("et"));
    funcEtaUdsc_.push_back(iSet->getParameter<std::string>("eta"));
    funcPhiUdsc_.push_back(iSet->getParameter<std::string>("phi"));
  }
  for(std::vector<edm::ParameterSet>::const_iterator iSet = bResolutions.begin(); iSet != bResolutions.end(); ++iSet){
    if(iSet->exists("bin")) binsB_.push_back(iSet->getParameter<std::string>("bin"));
    else if(bResolutions.size()==1) binsB_.push_back("");
    else throw cms::Exception("WrongConfig") << "Parameter 'bin' is needed if more than one PSet is specified!\n";
    
    funcEtB_.push_back(iSet->getParameter<std::string>("et"));
    funcEtaB_.push_back(iSet->getParameter<std::string>("eta"));
    funcPhiB_.push_back(iSet->getParameter<std::string>("phi"));
  }
  for(std::vector<edm::ParameterSet>::const_iterator iSet = lepResolutions.begin(); iSet != lepResolutions.end(); ++iSet){
    if(iSet->exists("bin")) binsLep_.push_back(iSet->getParameter<std::string>("bin"));
    else if(lepResolutions.size()==1) binsLep_.push_back("");
    else throw cms::Exception("WrongConfig") << "Parameter 'bin' is needed if more than one PSet is specified!\n";
    
    funcEtLep_.push_back(iSet->getParameter<std::string>("et"));
    funcEtaLep_.push_back(iSet->getParameter<std::string>("eta"));
    funcPhiLep_.push_back(iSet->getParameter<std::string>("phi"));
  }
  for(std::vector<edm::ParameterSet>::const_iterator iSet = metResolutions.begin(); iSet != metResolutions.end(); ++iSet){
    if(iSet->exists("bin")) binsMet_.push_back(iSet->getParameter<std::string>("bin"));
    else if(metResolutions.size()==1) binsMet_.push_back("");
    else throw cms::Exception("WrongConfig") << "Parameter 'bin' is needed if more than one PSet is specified!\n";
    
    funcEtMet_.push_back(iSet->getParameter<std::string>("et"));
    funcEtaMet_.push_back(iSet->getParameter<std::string>("eta"));
    funcPhiMet_.push_back(iSet->getParameter<std::string>("phi"));
  }
}

double CovarianceMatrix::getResolution(const TLorentzVector& object, const ObjectType objType, const std::string& whichResolution)
{
  std::vector<std::string> * bins_, * funcEt_, * funcEta_, * funcPhi_;

  switch(objType) {
  case kUdscJet:
    bins_    = &binsUdsc_;
    funcEt_  = &funcEtUdsc_;
    funcEta_ = &funcEtaUdsc_;
    funcPhi_ = &funcPhiUdsc_;
    break;
  case kBJet:
    bins_    = &binsB_;
    funcEt_  = &funcEtB_;
    funcEta_ = &funcEtaB_;
    funcPhi_ = &funcPhiB_;
    break;
  case kMuon:
    bins_    = &binsLep_;
    funcEt_  = &funcEtLep_;
    funcEta_ = &funcEtaLep_;
    funcPhi_ = &funcPhiLep_;
    break;
  case kElectron:
    bins_    = &binsLep_;
    funcEt_  = &funcEtLep_;
    funcEta_ = &funcEtaLep_;
    funcPhi_ = &funcPhiLep_;
    break;
  case kMet:
    bins_    = &binsMet_;
    funcEt_  = &funcEtMet_;
    funcEta_ = &funcEtaMet_;
    funcPhi_ = &funcPhiMet_;
    break;
  default:
    throw cms::Exception("UnsupportedObject") << "The object given is not supported!\n";
  }

  int selectedBin=-1;
  reco::LeafCandidate candidate;
  for(unsigned int i=0; i<bins_->size(); ++i){
    StringCutObjectSelector<reco::LeafCandidate> select_(bins_->at(i));
    candidate = reco::LeafCandidate( 0, reco::LeafCandidate::LorentzVector(object.Px(), object.Py(), object.Pz(), object.Energy()) );
    if(select_(candidate)){
      selectedBin = i;
      break;
    }
  }
  if(selectedBin>=0){
    if(whichResolution == "et")       return StringObjectFunction<reco::LeafCandidate>(funcEt_ ->at(selectedBin)).operator()(candidate);
    else if(whichResolution == "eta") return StringObjectFunction<reco::LeafCandidate>(funcEta_->at(selectedBin)).operator()(candidate);
    else if(whichResolution == "phi") return StringObjectFunction<reco::LeafCandidate>(funcPhi_->at(selectedBin)).operator()(candidate);
    else throw cms::Exception("ProgrammingError") << "Only 'et', 'eta' and 'phi' resolutions supported!\n";
  }
  return 0;
}

TMatrixD CovarianceMatrix::setupMatrix(const TLorentzVector& object, const ObjectType objType, const TopKinFitter::Param param)
{
  TMatrixD CovM3 (3,3); CovM3.Zero();
  TMatrixD CovM4 (4,4); CovM4.Zero();
  TMatrixD* CovM = &CovM3;
  const double pt  = object.Pt();
  const double eta = object.Eta();
  switch(objType) {
  case kUdscJet:
    // if object is a non-b jet
    {
      res::HelperJet jetRes;
      switch(param) {
      case TopKinFitter::kEMom :
	CovM4(0,0) = pow(jetRes.a (pt, eta, res::HelperJet::kUds), 2);
	CovM4(1,1) = pow(jetRes.b (pt, eta, res::HelperJet::kUds), 2);
	CovM4(2,2) = pow(jetRes.c (pt, eta, res::HelperJet::kUds), 2);
	CovM4(3,3) = pow(jetRes.d (pt, eta, res::HelperJet::kUds), 2);
	CovM = &CovM4;
	break;
      case TopKinFitter::kEtEtaPhi : 
	if(!binsUdsc_.size()){
	  CovM3(0,0) = pow(jetRes.et (pt, eta, res::HelperJet::kUds), 2);
	  CovM3(1,1) = pow(jetRes.eta(pt, eta, res::HelperJet::kUds), 2);
	  CovM3(2,2) = pow(jetRes.phi(pt, eta, res::HelperJet::kUds), 2);
	}
	else{
	  CovM3(0,0) = pow(getResolution(object, objType, "et") , 2); 
	  CovM3(1,1) = pow(getResolution(object, objType, "eta"), 2); 
	  CovM3(2,2) = pow(getResolution(object, objType, "phi"), 2);
	}   
	CovM = &CovM3;
	break;
      case TopKinFitter::kEtThetaPhi :
	CovM3(0,0) = pow(jetRes.et   (pt, eta, res::HelperJet::kUds), 2);
	CovM3(1,1) = pow(jetRes.theta(pt, eta, res::HelperJet::kUds), 2);
	CovM3(2,2) = pow(jetRes.phi  (pt, eta, res::HelperJet::kUds), 2);
	CovM = &CovM3;
	break; 
      }
    }
    break;
  case kBJet:
    // if object is a b jet
    {
      res::HelperJet jetRes;
      switch(param) {
      case TopKinFitter::kEMom :
	CovM4(0,0) = pow(jetRes.a (pt, eta, res::HelperJet::kB), 2);
	CovM4(1,1) = pow(jetRes.b (pt, eta, res::HelperJet::kB), 2);
	CovM4(2,2) = pow(jetRes.c (pt, eta, res::HelperJet::kB), 2);
	CovM4(3,3) = pow(jetRes.d (pt, eta, res::HelperJet::kB), 2);
	CovM = &CovM4;
	break;
      case TopKinFitter::kEtEtaPhi : 
	if(!binsUdsc_.size()){
	  CovM3(0,0) = pow(jetRes.et (pt, eta, res::HelperJet::kB), 2);
	  CovM3(1,1) = pow(jetRes.eta(pt, eta, res::HelperJet::kB), 2);
	  CovM3(2,2) = pow(jetRes.phi(pt, eta, res::HelperJet::kB), 2);
	}
	else{
	  CovM3(0,0) = pow(getResolution(object, objType, "et") , 2); 
	  CovM3(1,1) = pow(getResolution(object, objType, "eta"), 2); 
	  CovM3(2,2) = pow(getResolution(object, objType, "phi"), 2);
	}   
	CovM = &CovM3;
	break;
      case TopKinFitter::kEtThetaPhi :
	CovM3(0,0) = pow(jetRes.et   (pt, eta, res::HelperJet::kB), 2);
	CovM3(1,1) = pow(jetRes.theta(pt, eta, res::HelperJet::kB), 2);
	CovM3(2,2) = pow(jetRes.phi  (pt, eta, res::HelperJet::kB), 2);
	CovM = &CovM3;
	break; 
      }
    }
    break;
  case kMuon:
    // if object is a muon
    {
      res::HelperMuon muonRes;
      switch(param){
      case TopKinFitter::kEMom :
	CovM3(0,0) = pow(muonRes.a (pt, eta), 2);
	CovM3(1,1) = pow(muonRes.b (pt, eta), 2); 
	CovM3(2,2) = pow(muonRes.c (pt, eta), 2);
	CovM = &CovM3;
	break;
      case TopKinFitter::kEtEtaPhi :
	if(!binsLep_.size()){
	  CovM3(0,0) = pow(muonRes.et (pt, eta), 2);
	  CovM3(1,1) = pow(muonRes.eta(pt, eta), 2); 
	  CovM3(2,2) = pow(muonRes.phi(pt, eta), 2);
	}
	else{
	  CovM3(0,0) = pow(getResolution(object, objType, "et") , 2);
	  CovM3(1,1) = pow(getResolution(object, objType, "eta"), 2);
	  CovM3(2,2) = pow(getResolution(object, objType, "phi"), 2);
	}
	CovM = &CovM3;
	break;
      case TopKinFitter::kEtThetaPhi :
	CovM3(0,0) = pow(muonRes.et   (pt, eta), 2);
	CovM3(1,1) = pow(muonRes.theta(pt, eta), 2); 
	CovM3(2,2) = pow(muonRes.phi  (pt, eta), 2);
	CovM = &CovM3;
	break;
      }
    }
    break;
  case kElectron:
    {
      // if object is an electron
      res::HelperElectron elecRes;
      switch(param){
      case TopKinFitter::kEMom :
	CovM3(0,0) = pow(elecRes.a (pt, eta), 2);
	CovM3(1,1) = pow(elecRes.b (pt, eta), 2); 
	CovM3(2,2) = pow(elecRes.c (pt, eta), 2);
	CovM = &CovM3;
	break;
      case TopKinFitter::kEtEtaPhi :
	if(!binsLep_.size()){
	  CovM3(0,0) = pow(elecRes.et (pt, eta), 2);
	  CovM3(1,1) = pow(elecRes.eta(pt, eta), 2); 
	  CovM3(2,2) = pow(elecRes.phi(pt, eta), 2);
	}
	else{
	  CovM3(0,0) = pow(getResolution(object, objType, "et") , 2);
	  CovM3(1,1) = pow(getResolution(object, objType, "eta"), 2);
	  CovM3(2,2) = pow(getResolution(object, objType, "phi"), 2);
	}
	CovM = &CovM3;
	break;
      case TopKinFitter::kEtThetaPhi :
	CovM3(0,0) = pow(elecRes.et   (pt, eta), 2);
	CovM3(1,1) = pow(elecRes.theta(pt, eta), 2); 
	CovM3(2,2) = pow(elecRes.phi  (pt, eta), 2);
	CovM = &CovM3;
	break;
      }
    }
    break;
  case kMet:
    // if object is met
    {
      res::HelperMET metRes;
      switch(param){
      case TopKinFitter::kEMom :
	CovM3(0,0) = pow(metRes.a(pt), 2);
	CovM3(1,1) = pow(metRes.b(pt), 2);
	CovM3(2,2) = pow(metRes.c(pt), 2);
	CovM = &CovM3;
	break;
      case TopKinFitter::kEtEtaPhi :
	if(!binsMet_.size()){
	  CovM3(0,0) = pow(metRes.et(pt) , 2);
	  CovM3(1,1) = pow(        9999. , 2);
	  CovM3(2,2) = pow(metRes.phi(pt), 2);
	}
	else{
	  CovM3(0,0) = pow(getResolution(object, objType, "et") , 2);
	  CovM3(1,1) = pow(getResolution(object, objType, "eta"), 2);
	  CovM3(2,2) = pow(getResolution(object, objType, "phi"), 2);
	}
	CovM = &CovM3;
	break;
      case TopKinFitter::kEtThetaPhi :
	CovM3(0,0) = pow(metRes.et(pt) , 2);
	CovM3(1,1) = pow(        9999. , 2);
	CovM3(2,2) = pow(metRes.phi(pt), 2);
	CovM = &CovM3;
	break;
      }
    }
    break;
  }
  return *CovM;
}
