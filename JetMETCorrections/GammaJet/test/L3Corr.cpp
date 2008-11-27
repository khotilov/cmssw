// Adapted from D0 Experiment jetcorr/src/RjetCorr.cpp

#include "L3Corr.hpp"
#include <math.h>

using namespace std;

L3Corr::L3Corr(const JetAlg& jetalg, 
	       const jec::ErrorTypes& errType) :
  _jetalg(jetalg), _errType(errType)
{

  _E0 = 100;

  this->SetJetAlg(jetalg);
}

// Uncertainty is returned as _absolute_ uncertainty
double L3Corr::Rjet(const double pTprime, double& err)  {

  double rjet = this->_Rjet(pTprime);

  double err2 = 0;

  if (_errType & jec::kL3Sys) {
    // Systematics are a relative uncertainty so need to multiply by rjet
    double errSyst = this->_SystErr(pTprime);
    err2 += rjet * rjet * errSyst * errSyst;
  }

  if (_errType & jec::kL3Stat) {
    // Statistics are reported as relative so need to multiply by rjet
    double errStat = this->_StatErr(pTprime);
    err2 += rjet * rjet * errStat * errStat;
  }
  
  err = sqrt(err2);

  return rjet;
}

void L3Corr::SetJetAlg(const JetAlg& jetalg) {

  _jetalg = jetalg;

  switch (_jetalg) {

  // This uses fit to reconstructed MC
  case L3Corr::IC5_DATA:
    _isMC = false;
    // These are from signal reco MC (medium+deltaeta)
    _par[0] = 0.275; _par[1] = 0.4576; _par[2] = 0.9269;
    // These are from staterr.C toyMC study
    // TRandom3 seed 771156000
    // Signal sample (toyMC 10 pb-1, 1.5-sigma peak):
    _pstat[0] = 0.3473; _pstat[1] = 0.6239; _pstat[2] = 1.004; // [2]+/-0.079
    _chi2 = 0.7393; // 8.9/12
    _eig[0][0] = -0.1097; _eig[1][0] = -0.1064; _eig[2][0]  = -0.1106;
    _eig[0][1] = -0.001698; _eig[1][1] = 0.008222; _eig[2][1]  = -0.006228;
    _eig[0][2] = 0.001132; _eig[1][2] = -0.000357; _eig[2][2]  = -0.00078;
    /*
    // These are from signal reco MC (medium+deltaeta)
    _par[0] = 0.275; _par[1] = 0.4576; _par[2] = 0.9269;
    _pstat[0] = 0.275; _pstat[1] = 0.4576; _pstat[2] = 0.9269;
    _chi2 = 0.5349;
    _eig[0][0] = 0.01217; _eig[1][0] = 0.03932; _eig[2][0]  = 0.01236;
    _eig[0][1] = 0.003087; _eig[1][1] = -0.001454; _eig[2][1]  = 0.001587;
    _eig[0][2] = -0.0003632; _eig[1][2] = -8.516e-05; _eig[2][2]  = 0.0006284;
    */
    /*
    // These are from mixed reco MC (medium+deltaeta)
    _par[0] = 0.2795; _par[1] = 0.4658; _par[2] = 0.931;
    _pstat[0] = 0.2795; _pstat[1] = 0.4658; _pstat[2] = 0.931;
    _chi2 = 0.8428;
    _eig[0][0] = 0.01693; _eig[1][0] = 0.05511; _eig[2][0]  = 0.01752;
    _eig[0][1] = 0.004658; _eig[1][1] = -0.002093; _eig[2][1]  = 0.002083;
    _eig[0][2] = -0.0004272; _eig[1][2] = -0.0001307; _eig[2][2]  = 0.0008241;
    */
    /*
    // This was also supposed to be mixed reco MC (medium) with same settings?
    _par[0] = 0.2848; _par[1] = 0.4885; _par[2] = 0.933;
    _pstat[0] = 0.2848; _pstat[1] = 0.4885; _pstat[2] = 0.933;
    _chi2 = 0.8486;
    _eig[0][0] = -0.02129; _eig[1][0] = -0.06299; _eig[2][0]  = -0.02191;
    _eig[0][1] = -0.00469; _eig[1][1] = 0.002374; _eig[2][1]  = -0.002268;
    _eig[0][2] = -0.0004994; _eig[1][2] = -0.0001395; _eig[2][2]  = 0.0008864;
    */
    break;
  // All other cases use IC5 MC truth until better available
  case L3Corr::SIS5_DATA:
  case L3Corr::SIS7_DATA:
  case L3Corr::KT5_DATA:
  case L3Corr::KT7_DATA:
  //
  case L3Corr::SIS5_MC:
  case L3Corr::SIS7_MC:
  case L3Corr::KT5_MC:
  case L3Corr::KT7_MC:
  //
  case L3Corr::IC5_MC:
    _isMC = true;
    // currently MC truth
    // "1 - 2.316*pow(x,0.6005-1)"
    // => "1 - 0.3679*pow(pT/100., 0.6005-1)
    _par[0] = 0.3679; _par[1] = 0.6005; _par[2] = 1.;
    // These are from signal reco MC (medium+deltaeta)
    _pstat[0] = 0.275; _pstat[1] = 0.4576; _pstat[2] = 0.9269;
    _chi2 = 0.5349;
    _eig[0][0] = 0.01217; _eig[1][0] = 0.03932; _eig[2][0]  = 0.01236;
    _eig[0][1] = 0.003087; _eig[1][1] = -0.001454; _eig[2][1]  = 0.001587;
    _eig[0][2] = -0.0003632; _eig[1][2] = -8.516e-05; _eig[2][2]  = 0.0006284;
    /*
    _pstat[0] = 0.2869; _pstat[1] = 0.4872; _pstat[2] = 0.9356;
    _chi2 = 0.8321;
    _eig[0][0] = -0.01227; _eig[1][0] = -0.03484; _eig[2][0]  = -0.01241;
    _eig[0][1] = -0.002727; _eig[1][1] = 0.001456; _eig[2][1]  = -0.001392;
    _eig[0][2] = -0.0003259; _eig[1][2] = -8.209e-05; _eig[2][2]  = 0.0005524;
    */
    break;
  default:
    cout << "Unknown jet algorithm!" << endl << flush;
    exit(1);
  };
}

void L3Corr::SetErrType(const jec::ErrorTypes& errType) {
  _errType = errType;
}

jec::ErrorTypes L3Corr::GetErrType() {
  return _errType;
}


double L3Corr::_Rjet(const double pTprime) const {

  return _powerlaw(pTprime, _par);
}


// NB: Systematic errors are reported as _relative_ uncertainties
double L3Corr::_SystErr(const double pTprime) const {

  if (pTprime); // supress warning
  double err2 = 0;
 
  if (_errType & jec::kL3PhotonES) {

    // Photon energy scale error from the PhotonID group
    // The 2% estimate on EM scale mentioned at the photon ID workshop
    // in autumn 2007. Relative photon / electron scale is less than about 0.5%
    // Consult photon ID group again once we get data 
    double errEMScale = 0.02;
    double errPhotonRel = 0.005;

    err2 += errEMScale * errEMScale + errPhotonRel * errPhotonRel;
  }

  if (_errType & jec::kL3QCDBackground) {

    // Fit purity from output_mixed_all.root purity;2
    // - Take difference between medium, tight and loose as a conservative
    //   photon ID systematic for P and dC
    // - Assign 10% uncertainty on QCD and photon+jet relative
    //   NLO cross sections (better estimate from Tevatron papers) for P
    double P = _purity(pTprime, kMedium);//0.9;
    double errP = _SystPurity(pTprime, kMedium);//0.1
    double errdC; // 0.1
    double dC = _deltaC(pTprime, errdC, kMedium); // 0.1
    errdC *= (1+dC) * errdC; // turn into absolute uncertainty

    double errQCD_P = errP * fabs(dC);
    double errQCD_dC = (1-P) * errdC;
    double errQCD_cross = errP * errdC;

    // Quadratic addition assuming errP and dC fully unforrelated
    err2 += errQCD_P * errQCD_P + errQCD_dC * errQCD_dC
      + errQCD_cross * errQCD_cross;
    // Linear addition assuming full anticorrelation dC=-errP
    //err2 += pow(-errQCD_P + errQCD_dC - errQCD_cross, 2);
  }

  if (_errType & jec::kL3Parton) {

    // photon+jet balanced at parton level: photon already at parton level,
    // but jet not, so account for uncertainty in jet parton correction
    // Best numbers probably from Attilio's study looking into difference
    // in parton correction between Pythia and Herwig, and between Pythia tunes
    // One piece is effectively coming from fragmentation modeling (showering),
    // the other from underlying event modeling
    // - temporarily use 100% of jet showering in photon+jets
    // NB: disagreement with Attilio's results, should be 2-3% bigger effect
    // [- estimate UE at ~0.5 GeV per 0.7 cone; take 100% as uncertainty
    //    (D0 estimated UE at ~0.2 GeV per 0.7 cone at reco level,
    //     which is around 2.5 times more at generator level)]
    // - QCD group's CR2008/034 gives dN/dphideta~1.7, dpT/dphideta~2.2-3
    //   for charged tracks (I guess charged+neutral is 1.5 times larger)
    //   in Fig. 3; however, Fig. 4 has significantly lower UE of
    //   dpT/dphideta~0.5-1.5, dN~0.5 for drell-Yan? See CMS Note 2006/067
    //   => on average one charged track with 0.75 GeV in 0.5 cone for dN~1.5
    //      guess dpT/dphideta~1.5, take 50% uncertainty
    // NB: UE apparently only affects the mean, not the peak!!
    if (_errType & jec::kL3PartonFrag) {

      double kjets = -4.123 + 5.135 * pow(pTprime, -0.0006598);
      double dkjet = 1.00 * (kjets - 1.);

      err2 += dkjet * dkjet;
    }
    if (_errType & jec::kL3PartonUE) {

      //double due = 0.5 * (0.5/0.7)*(0.5/0.7) / pTprime;
      // 50% * [(charged+neutral)/charged=1.5] * [dpT/dphideta=1.5] * pi * R2
      double due = 0.50 * 1.5 * 1.5 * 3.14 * 0.5 *0.5 / pTprime;

      err2 += due * due;
    }
  }

  if (_errType & jec::kL3HardGluon) {

    // This is the uncertainty on the cut on second jet pT
    // Due to different response in data and MC this cut may be
    // effectively different
    // Test the sensitivity by varying the cut by a large factor
    // (x2 for ratio, x1.5 for pT threshold)
    // These numbers are for mixed sample
  //double rjet  = 0.9310 - 0.2795*pow(pTprime, 0.4658-1.); //pT2<0.10*pT,10GeV
  //double rjet2 = 0.8977 - 0.2783*pow(pTprime, 0.3471-1.); //pT2<0.20*pT,15GeV
  //double rjet3 = 0.9477 - 0.2877*pow(pTprime, 0.4931-1.); //pT2<0.05*pT,6GeV
    // These are for pure signal sample (more stable for rjet2)
  //double rjet  = 0.9269 - 0.2750*pow(pTprime, 0.4576-1.); //pT2<0.10*pT,10GeV
  //double rjet2 = 0.9868 - 0.3509*pow(pTprime, 0.6204-1.); //pT2<0.20*pT,20GeV
  //double rjet3 = 0.9459 - 0.2859*pow(pTprime, 0.4928-1.); //pT2<0.05*pT,5GeV
    // These are for pure signal with deltaeta cut
    double rjet05 = 0.9391 - 0.2758*pow(pTprime,0.4703-1.); //pT2<0.05*pT,5GeV
    double rjet10 = 0.9269 - 0.2750*pow(pTprime,0.4576-1.); //pT2<0.10*pT,10GeV
    double rjet20 = 0.9853 - 0.3465*pow(pTprime,0.6151-1.); //pT2<0.20*pT,20GeV
    double eg = 0.5 * (fabs(rjet05/rjet10 - 1.) + fabs(rjet20/rjet10 - 1.));

    err2 += eg * eg;
  }

  double err = sqrt(err2);

  return err;
}

// NB: Statistical uncertainty is reported as relative uncertainty
double L3Corr::_StatErr(const double pTprime) const {

  if (pTprime); // supress warning
  double err2 = 0;

  if (_errType & jec::kL3Stat) {

    /*
    // Statistical uncertainty per 25 GeV bin
    const double lumi = 0.010; // 10 pb-1 = 0.010 fb-1
    const double N = 1., S = 2., C = 0.06;
    double x = pTprime;
    double res = sqrt(N*N/(x*x) + S*S/x + C*C);
    double npfb = 7.729e+10*pow(x,-3.104)*exp(-0.003991*x); // events per fb-1
    double n = lumi * npfb;
    //double stat = res / sqrt(n);

    err2 += res * res / n;
    */

    // Calculate stat uncertainty from fit (provided as rel. eigenfunctions)
    double e0 = _eigpowerlaw(pTprime, _pstat, &_eig[0], 0);
    double e1 = _eigpowerlaw(pTprime, _pstat, &_eig[0], 1);
    double e2 = _eigpowerlaw(pTprime, _pstat, &_eig[0], 2);
    
    err2 += max(_chi2, 1.) * (e0 * e0 + e1 * e1 + e2 * e2); 

    /*
    // staterr directly
    double f
		 "sqrt(pow(0.01*x,2*[1]-2)*[2]"
		 "+2*pow(0.01*x,2*[1]-2)*[0]*log(0.01*x)*[3]"
		 "-2*pow(0.01*x,[1]-1)*[4]"
		 "+pow(0.01*x,2*[1]-2)*[0]*[0]*log(0.01*x)*log(0.01*x)*[5]"
		 "-2*pow(0.01*x,[1]-1)*[0]*log(0.01*x)*[6]"
		 "+[7])",
    */
  }

  double err = sqrt(err2);

  return err;
}

double L3Corr::_purity(const double pTprime, const PhotonID& id) const {

  double P = 1.;
  double x = log(0.01*pTprime);
  // The 'no deltaeta' were fitted by hand from output_mixed_all.root
  // from plots in the range 30-700 GeV
  if (id == kLoose) {
    P = 0.6339 + x*(0.213 + x*-0.06296); // deltaeta, no rebin
    //P = 0.5514 + x*(0.2773 + x*-0.06944); // deltaeta, automatic
    //P = 6.11287e-01 + x*(1.80336e-01 + x*-1.69342e-02); // no deltaeta?
  }
  if (id == kMedium) {
    P = 0.9673 + x*(-0.01242 + x*0.01187); // deltaeta, no rebin
    //P = 9.66655e-01 + x*(-1.28295e-02 + x*1.29304e-02); // no deltaeta?
  }
  if (id == kTight) {
    P = 0.9797 + x*(-0.006598 + x*0.004506); // deltaeta, no rebin
    //P = 9.78925e-01 + x*(-7.17726e-03 + x*5.88710e-03); // no deltaeta?
  }
  if (id == kMedium005) {
    P = 0.9921 + x*(-0.0128 + x*0.005335); // deltaeta, no rebin
    //P = 0.9876 + x*(-0.03873 + x*0.01767);
  }
  if (id == kMedium020) {
    P = 0.9581 + x*(-0.01365 + x*0.007251); // deltaeta, no rebin
    //P = 0.6982 + x*(0.2757 + x*-0.08513);
  }

  if (P<0) P = 0;
  if (P>1) P = 1;

  return P;
}

// NB: all purity uncertainties are _absolute_ values
double L3Corr::_SystPurity(const double pTprime, const PhotonID& id) const {

  double err2 = 0;
  
  double errstat = _StatPurity(pTprime, id);
  if (_errType & jec::kL3PurityStat)
    err2 += errstat * errstat;

  if (!_isMC && (_errType & jec::kL3PurityID))
    err2 += pow(_SystPurityID(pTprime, id),2);

  if (!_isMC && (_errType & jec::kL3PurityXsec))
    err2 += pow(_SystPurityXsec(pTprime, id), 2);

  if (!_isMC && (_errType & jec::kL3Purity2ndJet)) {
    err2 += pow(_SystPurity2ndJet(pTprime, id), 2);
  }

  double err = sqrt(err2);
		    
  return err;
}

// Systematic on photon+jet purity from photon ID
// (average difference between medium, loose and tight, but could use
//  more detailed estimates for each cut also)
// => Change the impact to 30% of the difference, now it seems too large
//    for medium ID
double L3Corr::_SystPurityID(const double pTprime, const PhotonID& id) const {

  //double dP = 0.5*(fabs(_purity(pTprime, id) - _purity(pTprime, kLoose))
  //	   + fabs(_purity(pTprime, id) - _purity(pTprime, kTight)));

  double dP = 0.30*(fabs(_purity(pTprime, id) - _purity(pTprime, kLoose))
		    + fabs(_purity(pTprime, id) - _purity(pTprime, kTight)));
  
  return dP;
}

// Systematic on photon+jet purity from NLO cross sections
// Assume ratio of QCD and photon+jet is known to 10%
// See e.g. arXiv:/0804.1107 [hep-ex] (photon+jet)
//          arXiv:/0802.2400 [hep-ex] (inclusive jet~dijet)
// => 20% uncertainty on B/S (theory has 10%, data/theory discrepancy 20%)
// P' = P / (P + x*(1-P)), where x = B'/B
// dP'/dx(at x=1) = -(1-P)*P
double L3Corr::_SystPurityXsec(const double pTprime, const PhotonID& id) const{

  double P = _purity(pTprime, id);
  double dxsec = 0.2;

  return (dxsec * (1-P) * P);
}

double L3Corr::_SystPurity2ndJet(const double pTprime,
				 const PhotonID& id) const {

  if (id); // Only for medium now, but suppress error
  double P = _purity(pTprime, kMedium);
  double P2 = _purity(pTprime, kMedium005);
  double P3 = _purity(pTprime, kMedium020);
  double err = 0.5 * (fabs(P2/P-1.) + fabs(P3/P-1.));

  return err;
}

// Statistical uncertainty from MC fit correlation matrix
// (quadratic logarithmic fit)
double L3Corr::_StatPurity(const double pTprime, const PhotonID& id) const {

  double x = log(0.01*pTprime);
  double stat = 0;
  if (id == kLoose)
    stat = sqrt(0.0009192 + 2*-0.0002421*x + 2*-0.0002416*x*x
		+ 0.002093*x*x + 2*-0.001193*x*x*x
		+ 0.0009278*x*x*x*x);
  //stat = sqrt(0.000866 + 2*-0.0003189*x + 2*-0.0001334*x*x
  //	+ 0.001982*x*x + 2*-0.001037*x*x*x
  //	+ 0.0007077*x*x*x*x);
  if (id == kMedium)
    stat = 0.08077; // top of err.bars
  //stat = sqrt(0.0003825 + 2*2.183e-05*x + 2*-0.0001547*x*x
  //	+ 0.0005671*x*x + 2*-0.0003245*x*x*x
  //	+ 0.0002734*x*x*x*x); // deltaeta, no rebin
    //stat = sqrt(0.0003744 + 2*1.692e-05*x + 2*-0.000142*x*x
    //	+ 0.0005641*x*x + 2*-0.0003168*x*x*x
    //	+ 0.0002534*x*x*x*x);
  if (id == kTight)
    stat = sqrt(0.0004265 + 2*-1.484e-05*x + 2*-0.0001517*x*x
		+ 0.0008268*x*x + 2*-0.0004603*x*x*x
		+ 0.0003639*x*x*x*x); // deltaeta, no rebin
    //stat = sqrt(0.0004191 + 2*-2.073e-05*x + 2*-0.0001376*x*x
    //	+ 0.0008222*x*x + 2*-0.0004491*x*x*x
    //	+ 0.0003374*x*x*x*x);
  if (id == kMedium005)
    stat = sqrt(0.0005035 + 2*-2.595e-05*x + 2*-0.0001639*x*x
		+ 0.001006*x*x + 2*-0.0005597*x*x*x
		+ 0.0004097*x*x*x*x); // deltaeta, no rebin
    //stat = sqrt(0.000621 + 2*-0.0001296*x + 2*-0.000138*x*x
    //	+ 0.002332*x*x + 2*-0.001283*x*x*x
    //	+ 0.0008167*x*x*x*x);
  if (id == kMedium020)
    stat = sqrt(0.0003208 + 2*2.855e-05*x + 2*-0.0001461*x*x
		+ 0.000454*x*x + 2*-0.0002643*x*x*x
		+ 0.0002498*x*x*x*x); // deltaeta, no rebin
    //stat = sqrt(0.004524 + 2*-0.006549*x + 2*0.002254*x*x
    //	+ 0.01395*x*x + 2*-0.00594*x*x*x
    //	+ 0.002773*x*x*x*x);

  return stat;
}

// Returns deltaC and _relative_ uncertainties
double L3Corr::_deltaC(const double pTprime, double& syserr,
		       const PhotonID& id,
		       double (*syserrs)[6]) const {
  
  double x = log(0.01*pTprime);
  double y = pTprime;
  double rjets = 1. - 1.881 * pow(y, -0.3802);
  double rjetb = 1. - 2.332 * pow(y, -0.4005);
  double kjets = -4.123 + 5.135 * pow(y, -0.0006598);
  double kjetb = 0.9239 + x * (0.0313 + x * -0.004696);
  double kjetb_p = 0.9726 + x * (0.01841 + x * -0.003336);
  double ktopos = 0.9973 + x * (-0.0005675 + x * 0);
  double ktopob = 1.002 + x * (-0.005414 + x * 0.001552);
  
  double rphos=1., rphob=1., kphos=1., kphob=1.;
  if (id == kLoose) {
    rphos = 1 - 0.007997 * pow(y, 0.07654);
    rphob = 1 - 0.8911 * pow(y, -0.4005);
    kphos = 1 - 1.877e-12 * pow(y, -0.4549);
    kphob = 0.9236 + x * (0.0313 + x * -0.004696);
  }
  if (id == kMedium) {
    rphos = 1 - 0.003445 * pow(y, 0.1838);
    rphob = 1 - 0.226 * pow(y, -0.4005);
    kphos = 1 - 1.397e-12 * pow(y, -0.4546);
    kphob = 0.8957 + x * (0.0313 + x * -0.004696);
  }
  if (id == kTight) {
    rphos = 1 - 0.07228 * pow(y, -0.4365);
    rphob = 1 - 0.6667 * pow(y, -0.4005);
    kphos = 1 - 3.707e-12 * pow(y, -0.4573);
    kphob = 0.9313 + x * (0.0313 + x * -0.004696);
  }
  
  double dC = ((rjetb * kjetb * ktopob) / (rphob * kphob))
    / ((rjets * kjets * ktopos) / (rphos * kphos)) - 1.;
  
  double err2 = 0;

  // Statistical uncertainty of the DeltaC fit in MC
  // Get the shape from factorized corrections which has good
  // effective statistics and fit to the total deltaC
  // Allow some scaling in both energy dependence and magnitude:
  // dC = a*dC0 + b
  // SOURCE: MC statistics
  // NB: Currently just placeholder
  if (_errType & jec::kL3DeltaCStat) {
    double stat = 0.02;//5;
    err2 += stat * stat;
    if (syserrs) (*syserrs)[0] = stat;
  }

  // EM-jet/photon response
  // => 50% of the difference between EM-jet and a real photon in MC
  //    (CMS EM shower MC simulation not yet certified on data)
  // SOURCE: detector simulation of EM shower (+photonID cuts)
  if (_errType & jec::kL3DeltaCRphot) {
    double drpho = 0.50 * (rphob / rphos - 1.);
    err2 += drpho * drpho;
    if (syserrs) (*syserrs)[1] = drpho;
  }

  // EM-jet/photon showering compared to jet showering
  // - normally they'd cancel for QCD, but not with very tight ID
  // => 100% of the difference observed in dijet MC (not well understood)
  // SOURCE: photon ID cuts and fragmentation simulation?
  // NB: Could study this by only applying isolation cuts in PhotonID,
  //     but not H/E (would need hollow cones...) and cluster shape
  if (_errType & jec::kL3DeltaCKphot) {
    double dkb = 1.00 * (kphob / kjetb - 1.);
    err2 += dkb * dkb;
    if (syserrs) (*syserrs)[2] = dkb;
  }
  /*
  //    + 50% of the jet showering observed in photon+jet MC (not corrected)
  // NB: photon showering = 0 by definition and jet absolute showering goes
  //     to another signal uncertainty source => avoid double counting
  if (_errType & jec::kL3DeltaCKphotSig) {
    double dks = 0.50 * (kphos / kjets - 1.);
    if (syserrs) (*syserrs)[3] = dks;
    err2 += dks * dks;
  }
  */

  // Jet response (quark and gluon jets different)
  // => 50% of the difference between photon+jet and dijet MC
  // (at D0 data showed twice as big an effect as MC, but CMS MC
  //  should work better)
  // SOURCE: detector simulation of single pion response
  if (_errType & jec::kL3DeltaCRjet) {
    double drjet = 0.50 * (rjetb / rjets - 1.);
    err2 += drjet * drjet;
    if (syserrs) (*syserrs)[3] = drjet;
  }

  // Jet showering (quark and gluon jets different)
  // => 50% of the difference between photon+jet and dijet MC
  // (replace with a Pythia/Herwig difference from Attilio's
  //  L7 parton corrections)
  // NB: jet showering for QCD cancels with EM-jet, except for
  //     uncertainty already counted in kL3DeltaCKphot;
  //     signal showering is a separate uncertainty
  //     => set this to 0. to avoid double-counting
  // SOURCE: fragmentation simulation
  if (_errType & jec::kL3DeltaCKjet) {
    double dkjet = 0.;//0.50 * (kjetb / kjets - 1.);
    err2 += dkjet * dkjet;
    if (syserrs) (*syserrs)[4] = dkjet;
  }

  // Difference between mean and peak (MPV) for dijet showering
  // While the difference should survive for the probe jet, it is
  // not obvious that it is fully preserved for the EM-jet due
  // to bin-to-bin migration effects => increases showering difference
  // => 100% of the difference between mean and peak
  // SOURCE: fragmentation simulation + accounting of migration effects?
  if (_errType & jec::kL3DeltaCKjetPeak) {
    double dpeak = kjetb_p / kjetb - 1.;
    err2 += dpeak * dpeak;
    if (syserrs) (*syserrs)[5] = dpeak;
  }
  
  syserr = sqrt(err2);

  return dC;
}

double L3Corr::_powerlaw(const double& pTprime, const double* par) const {  

  return (_par[2] - par[0] * pow(0.01*pTprime, par[1]-1));
}

double L3Corr::_eigpowerlaw(const double& pTprime, const double* par,
			    const double (*eig)[npar], const int& ieig) const {

  // Partial derivatives of 'p2 - p0*pow(0.01*pTprime, p1-1)'
  // w.r.t. the parameters times eigenvectors give the eigenfunctions
  double f = _powerlaw(pTprime, par);
  double df0 = -pow(0.01*pTprime, par[1]-1);
  double df1 = -par[0] * log(0.01*pTprime) * pow(0.01*pTprime, par[1]-1);
  double df2 = 1.;

  //return (eig[ieig][0]*df0 + eig[ieig][1]*df1 + eig[ieig][2]*df2) / f;
  return (eig[0][ieig]*df0 + eig[1][ieig]*df1 + eig[2][ieig]*df2) / f;
  
}
