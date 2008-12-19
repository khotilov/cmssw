// Adapted from D0 Experiment jetcorr/jetcorr/RjetCorr.hpp

#ifndef __L3CORR__
#define __L3CORR__

//
// Purpose: low level routines for L3 absolute jet response correction
//



#include "JetDefs.hpp"

#include <iostream>

class L3Corr {
public:
  enum JetAlg {IC5_DATA, SIS5_DATA, SIS7_DATA, KT5_DATA, KT7_DATA,
	       IC5_MC, SIS5_MC, SIS7_MC, KT5_MC, KT7_MC};  

  enum PhotonID {kLoose, kMedium, kTight, kNN, kMedium005, kMedium020};

  L3Corr(const JetAlg& jetalg = L3Corr::IC5_DATA, 
	 const jec::ErrorTypes& errType = jec::kAll,
	 const PhotonID& id = L3Corr::kMedium);
  ~L3Corr(){};
  
  double Rjet(const double pTprime, double& err);

  void SetJetAlg(const JetAlg& jetalg = L3Corr::IC5_DATA);
  void SetErrType(const jec::ErrorTypes& errType = jec::kAll);
  void SetPhoID(const PhotonID& id = L3Corr::kMedium);

  jec::ErrorTypes GetErrType();

  // Half-internal calculations
  // Keep accessible for plotting scripts
  double _Rjet(const double pTprime) const;

  double _SystErr(const double pTprime) const;
  double _StatErr(const double pTprime) const;

  double _purity(const double pTprime, const PhotonID& id) const;
  //const PhotonID& id = L3Corr::kMedium) const;
  double _SystPurity(const double pTprime, const PhotonID& id) const;
  double _SystPurityID(const double pTprime, const PhotonID& id) const;
  double _SystPurityXsec(const double pTprime, const PhotonID& id) const;
  //double _SystPurity2ndJet(const double pTprime, const PhotonID& id) const;
  double _StatPurity(const double pTprime, const PhotonID& id) const;

  double _partonFrag(const double pTprime) const;
  double _partonUE(const double pTprime) const;

  double _flavorMap(const double pTprime) const;

  double _deltaC(const double pTprime,
		 double& syserr,
		 const PhotonID& id, // = L3Corr::kMedium,
		 double (*syserrs)[6] = 0) const;

  double _deltaRjet(const double pTprime, const PhotonID& id) const;

  //double _SystParton(const double pTprime) const;


private:
  JetAlg _jetalg;
  jec::ErrorTypes _errType;
  PhotonID _phoID;
  bool _isMC;

  static const int npar = 3;
  double _E0;
  double _par[npar];

  double _chi2;
  double _pstat[npar];
  double _eig[npar][npar];

  double _powerlaw(const double& pTprime, const double* par) const;
  double _eigpowerlaw(const double& pTprime, const double* par,
		      const double (*eig)[npar], const int& ieig) const;

};

#endif /* __L3CORR__ */
