#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyCalibration.h"
#include <TMath.h>
#include <math.h>
#include <vector>
#include <TF1.h>

using namespace std;

PFEnergyCalibration::PFEnergyCalibration() {

//--- initialize calibration parameters
//    for energy correction applied to energy deposits of electrons 
//    and photons in ECAL
//  paramECAL_slope_ = 1.;
//  paramECAL_offset_ = 0.;
//  
////--- initialize calibration parameters
////    for energy correction applied to energy deposits of hadrons in HCAL
//  paramHCAL_slope_ = 2.17;
//  paramHCAL_offset_ = 1.73;
//  paramHCAL_damping_ = 2.49;
//
////--- initialize calibration parameters
////    for energy correction applied to combined energy deposits of hadrons in HCAL and ECAL
//  paramECALplusHCAL_slopeECAL_ = 1.05;
//  paramECALplusHCAL_slopeHCAL_ = 1.06;
//  paramECALplusHCAL_offset_ = 6.11;
	
  paramECAL_slope_ = 1.;
  paramECAL_offset_ = 0.;
  
//--- initialize calibration parameters
//    for energy correction applied to energy deposits of hadrons in HCAL
  paramHCAL_slope_ = 1.0;
  paramHCAL_offset_ = 0.0;
  paramHCAL_damping_ = 1.0;

//--- initialize calibration parameters
//    for energy correction applied to combined energy deposits of hadrons in HCAL and ECAL
  paramECALplusHCAL_slopeECAL_ = 1.0;
  paramECALplusHCAL_slopeHCAL_ = 1.0;
  paramECALplusHCAL_offset_ = 0.0;
}


PFEnergyCalibration::PFEnergyCalibration( double e_slope  , 
					  double e_offset , 
					  double eh_eslope,
					  double eh_hslope,
					  double eh_offset,
					  double h_slope  ,
					  double h_offset ,
					  double h_damping,
					  unsigned newCalib):
  
  paramECAL_slope_(e_slope),
  paramECAL_offset_(e_offset),
  paramECALplusHCAL_slopeECAL_(eh_eslope),
  paramECALplusHCAL_slopeHCAL_(eh_hslope),
  paramECALplusHCAL_offset_(eh_offset),
  paramHCAL_slope_(h_slope),
  paramHCAL_offset_(h_offset),
  paramHCAL_damping_(h_damping) 
{
  if ( newCalib == 2 ) initializeCalibrationFunctions();
}

void
PFEnergyCalibration::initializeCalibrationFunctions() {

  // New calibration with RespCorr factors
  // Thresholds
  threshE = 3.1;
  threshH = 2;

  // Barrel
  faBarrel = new TF1("faBarrel","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
  fbBarrel = new TF1("fbBarrel","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
  fcBarrel = new TF1("fcBarrel","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
  faEtaBarrel = new TF1("faEtaBarrel","[0]+[1]*x+[2]*exp(-x/[3])+[4]*[4]*exp(-x*x/([5]*[5]))",1.,1000.);
  fbEtaBarrel = new TF1("fbEtaBarrel","[0]+[1]*x+[2]*exp(-x/[3])+[4]*[4]*exp(-x*x/([5]*[5]))",1.,1000.);
  faBarrel->SetParameter(0,1.10772);
  fbBarrel->SetParameter(0,1.06012);
  fcBarrel->SetParameter(0,0.979137);
  faEtaBarrel->SetParameter(0,0.02);
  fbEtaBarrel->SetParameter(0,-0.02);
  faBarrel->SetParameter(1,0.186273);
  fbBarrel->SetParameter(1,0.273149);
  fcBarrel->SetParameter(1,0.3488);
  faBarrel->SetParameter(2,-0.47812);
  fbBarrel->SetParameter(2,-0.41739);
  fcBarrel->SetParameter(2,-1.04486);
  faEtaBarrel->SetParameter(2,-0.102858);
  fbEtaBarrel->SetParameter(2,0.14503);
  faBarrel->SetParameter(3,62.5754);
  fbBarrel->SetParameter(3,46.0841);
  fcBarrel->SetParameter(3,18.6968);
  faEtaBarrel->SetParameter(3,337.536);
  fbEtaBarrel->SetParameter(3,273.008);
  faBarrel->SetParameter(4,1.31965);
  fbBarrel->SetParameter(4,1.40679);
  fcBarrel->SetParameter(4,0.68429);
  faEtaBarrel->SetParameter(4,0.653127);
  fbEtaBarrel->SetParameter(4,0.0967542);
  faBarrel->SetParameter(5,35.2559);
  fbBarrel->SetParameter(5,30.2377);
  fcBarrel->SetParameter(5,22.5488);
  faEtaBarrel->SetParameter(5,4.73068);
  fbEtaBarrel->SetParameter(5,36.6991);

  // Endcap
  faEndcap = new TF1("faEndcap","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
  fbEndcap = new TF1("fbEndcap","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
  fcEndcap = new TF1("fcEndcap","[0]+([1]+[2]/sqrt(x))*exp(-x/[3])-[4]*exp(-x*x/[5])",1.,1000.);
  faEtaEndcap = new TF1("faEtaEndcap","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))",1.,1000.);
  fbEtaEndcap = new TF1("fbEtaEndcap","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))",1.,1000.);
  faEndcap->SetParameter(0,1.0877);
  fbEndcap->SetParameter(0,0.984949);
  fcEndcap->SetParameter(0,0.932945);
  faEtaEndcap->SetParameter(0,-0.0109133);
  fbEtaEndcap->SetParameter(0,0.025622);
  faEndcap->SetParameter(1,0.28939);
  fbEndcap->SetParameter(1,0.378368);
  fcEndcap->SetParameter(1,0.100645);
  faEtaEndcap->SetParameter(1,-0.103459);
  fbEtaEndcap->SetParameter(1,-2.58821);
  faEndcap->SetParameter(2,-0.57635);
  fbEndcap->SetParameter(2,-1.4482);
  fcEndcap->SetParameter(2,-0.0718337);
  faEtaEndcap->SetParameter(2,180.264);
  fbEtaEndcap->SetParameter(2,5.87477);
  faEndcap->SetParameter(3,86.5501);
  fbEndcap->SetParameter(3,68.4813);
  fcEndcap->SetParameter(3,46.8387);
  faEtaEndcap->SetParameter(3,0.642761);
  fbEtaEndcap->SetParameter(3,0.365089);
  faEndcap->SetParameter(4,1.02296);
  fbEndcap->SetParameter(4,0.596373);
  fcEndcap->SetParameter(4,0.809857);
  faEtaEndcap->SetParameter(4,14.9234);
  fbEtaEndcap->SetParameter(4,236.042);
  faEndcap->SetParameter(5,64.0116);
  fbEndcap->SetParameter(5,117.031);
  fcEndcap->SetParameter(5,57.1931);

}

void 
PFEnergyCalibration::energyEmHad(double t, double& e, double&h, double eta, double phi) const { 
 
  
  // Use calorimetric energy as true energy for neutral particles
  double tt = t;
  double ee = e;
  double hh = h;
  double a = 1.;
  double b = 1.;
  double etaCorr = 1.;
  t = min(1000.,max(tt,e+h));

  // Barrel calibration
  if ( fabs(eta) < 1.48 ) { 

    // Fudge factors for fast sim 
    // Fudge factors in ECAL crack (correction#0, removed)
    /*
    if ( eta>1.4) h *=0.90;
    else if ( eta>1.3 ) h *=0.86;
    else if ( eta>1.2 ) h *=0.95;
    */

    // The energy correction
    a = e>0. ? faBarrel->Eval(t) : 1.;
    b = e>0. ? fbBarrel->Eval(t) : fcBarrel->Eval(t);
    double thresh = e > 0. ? threshE : threshH;

    // Protection against negative calibration - to be tuned
    if ( a < 0. || b < 0. ) { 
      a = 1.;
      b = 1.;
      thresh = 0.;
    }

    // The new estimate of the true energy
    t = min(1000.,max(tt, thresh+a*e+b*h));

    // The angular correction for ECAL hadronic deposits
    etaCorr = 1. + faEtaBarrel->Eval(t) + fbEtaBarrel->Eval(t)*fabs(eta);
    // etaCorr = 1.;
    t = max(tt, thresh+etaCorr*a*e+b*h);

    if ( e > 0. && thresh > 0. ) 
      e = h > 0. ? threshE-threshH + etaCorr * a * e : threshE + etaCorr * a * e;
    if ( h > 0. && thresh > 0. ) 
      h = threshH + b * h;

    if ( etaCorr > 2. || etaCorr < 0.5 ) 
      std::cout << "Warning : Angular correction ! " << std::endl
		<< "etaCorr,eta,t = " << etaCorr << " " << eta << " " << t << std::endl;

  // Endcap calibration   
  } else {

    // Fudge factor in endcap (correction #2)

    // The energy correction
    a = e>0. ? faEndcap->Eval(t) : 1.;
    b = e>0. ? fbEndcap->Eval(t) : fcEndcap->Eval(t);
    double thresh = e > 0. ? threshE : threshH;

    if ( a < 0. || b < 0. ) { 
      a = 1.;
      b = 1.;
      thresh = 0.;
    }
    // Fudge factor fast sim
    // b *= 1.03;

    // The new estimate of the true energy
    t = min(1000.,max(tt, thresh+a*e+b*h));
    
    // The angular correction
    // That's to compensate for a possible bug in HCALRespCorr's coefficients
    if ( fabs(eta)>2.0 && fabs(eta)<2.65 ) h *= (1+0.2*(fabs(eta)-2.0));
    // And that's the real eta dependence
    if ( fabs(eta) < 2.65 ) 
      etaCorr = 1. + faEtaEndcap->Eval(t) + fbEtaEndcap->Eval(t)*(fabs(eta)-1.48);
    /*
    if ( etaCorr > 2. || etaCorr < 0.5 ) 
      std::cout << "Warning : Angular correction ! " << std::endl
		<< "etaCorr,eta,t = " << etaCorr << " " << eta << " " << tt 
		<< " ee,hh,e,h = " << e << " " << h << " " << a*e << " " << b*h  
		<< std::endl;
    */

    t = min(1000.,max(tt, thresh+etaCorr*a*e+b*h));

    if ( e > 0. && thresh > 0. ) 
      e = h > 0. ? threshE-threshH + etaCorr * a * e : threshE + etaCorr * a * e;
    if ( h > 0. && thresh > 0. ) 
      h = threshH + b * h;


  }

  // Protection
  if ( e < 0. || h < 0. ) { 
    std::cout << "Warning : Energy correction ! " << std::endl
	      << "eta,tt,e,h,a,b = " << eta << " " << tt << " " 
	      << ee << "/" << e << " " << hh << "/" << h << " " << a << " " << b << std::endl;
    // Some protection against crazy calibration
    if ( e < 0. ) e = ee;
    if ( h < 0. ) h = hh;
  }

  // And that's it !

  
}



void PFEnergyCalibration::setCalibrationParametersEm(double paramECAL_slope, 
						     double paramECAL_offset) {

//--- set calibration parameters for energy deposits of electrons 
//    and photons in ECAL;
//    this member function is needed by PFRootEvent

  paramECAL_slope_ = paramECAL_slope;
  paramECAL_offset_ = paramECAL_offset;
}



PFEnergyCalibration::~PFEnergyCalibration()
{
//--- nothing to be done yet  
}


double 
PFEnergyCalibration::energyEm(double uncalibratedEnergyECAL, 
			      double eta, double phi) const {

  //--- apply calibration correction
  //    for energy deposits of electrons and photons in ECAL
  //    (eta and phi dependence not implemented yet)
  
  double calibrated = paramECAL_slope_*uncalibratedEnergyECAL;
  calibrated += paramECAL_offset_;

  return calibrated;
}


double
PFEnergyCalibration::energyEm(const reco::PFCluster& clusterEcal,
			      std::vector<double> &EclustersPS1,
			      std::vector<double> &EclustersPS2,
			      bool crackCorrection ){
  double eEcal = clusterEcal.energy();
  //temporaty ugly fix
  reco::PFCluster myPFCluster=clusterEcal;
  myPFCluster.calculatePositionREP();
  double eta = myPFCluster.positionREP().eta();
  double phi = myPFCluster.positionREP().phi();

  double ePS1 = 0;
  double ePS2 = 0;

  for(unsigned i=0;i<EclustersPS1.size();i++) ePS1 += EclustersPS1[i];
  for(unsigned i=0;i<EclustersPS2.size();i++) ePS2 += EclustersPS2[i];

  double calibrated = Ecorr(eEcal,ePS1,ePS2,eta,phi, crackCorrection);
  if(eEcal!=0 && calibrated==0) std::cout<<"Eecal = "<<eEcal<<"  eta = "<<eta<<"  phi = "<<phi<<std::endl; 
  return calibrated; 
}

double PFEnergyCalibration::energyEm(const reco::PFCluster& clusterEcal,
				     std::vector<double> &EclustersPS1,
				     std::vector<double> &EclustersPS2,
				     double& ps1,double& ps2,
				     bool crackCorrection){
  double eEcal = clusterEcal.energy();
  //temporaty ugly fix
  reco::PFCluster myPFCluster=clusterEcal;
  myPFCluster.calculatePositionREP();
  double eta = myPFCluster.positionREP().eta();
  double phi = myPFCluster.positionREP().phi();

  double ePS1 = 0;
  double ePS2 = 0;

  for(unsigned i=0;i<EclustersPS1.size();i++) ePS1 += EclustersPS1[i];
  for(unsigned i=0;i<EclustersPS2.size();i++) ePS2 += EclustersPS2[i];

  double calibrated = Ecorr(eEcal,ePS1,ePS2,eta,phi,ps1,ps2,crackCorrection);
  if(eEcal!=0 && calibrated==0) std::cout<<"Eecal = "<<eEcal<<"  eta = "<<eta<<"  phi = "<<phi<<std::endl; 
  return calibrated; 
}


double 
PFEnergyCalibration::energyHad(double uncalibratedEnergyHCAL, 
			       double eta, double phi) const {
  
  //--- apply calibration correction
  //    for energy deposits of hadrons in HCAL
  //    (eta and phi dependence not implemented yet)
  
  double numerator = paramHCAL_slope_*uncalibratedEnergyHCAL;
  numerator += paramHCAL_offset_;

  double denominator = 1 + exp(paramHCAL_damping_/uncalibratedEnergyHCAL);


  return numerator/denominator;
}


double 
PFEnergyCalibration::energyEmHad(double uncalibratedEnergyECAL, 
				 double uncalibratedEnergyHCAL, 
				 double eta, double phi) const {
//--- apply calibration correction
//    for energy deposits of hadrons in ECAL and HCAL
//    (eta and phi dependence not implemented yet)

  double calibrated = paramECALplusHCAL_slopeECAL_*uncalibratedEnergyECAL;
  calibrated += paramECALplusHCAL_slopeHCAL_*uncalibratedEnergyHCAL;
  calibrated += paramECALplusHCAL_offset_;

  return calibrated;
}
  
std::ostream& operator<<(std::ostream& out, 
			 const PFEnergyCalibration& calib) {

  if(!out ) return out;

  out<<"PFEnergyCalibration -- "<<endl;
  out<<"ecal      = "<<calib.paramECAL_slope_
     <<" x E + "<< calib.paramECAL_offset_<<endl;
  out<<"hcal only = <add formula>"
     <<calib.paramHCAL_slope_<<","
     <<calib.paramHCAL_offset_<<","
     <<calib.paramHCAL_damping_<<endl;
  out<<"ecal+hcal = "<<calib.paramECALplusHCAL_slopeECAL_<<" x E_e + "
     <<calib.paramECALplusHCAL_slopeHCAL_<<" x E_h + "
     <<calib.paramECALplusHCAL_offset_<<endl;

  return out;
}




///////////////////////////////////////////////////////////////
////                                                       ////  
////             CORRECTION OF PHOTONS' ENERGY             ////
////                                                       ////
////              Material effect: No tracker              ////
////       Tuned on CMSSW_2_1_0_pre4, Full Sim events      ////
////                                                       ////
///////////////////////////////////////////////////////////////
////                                                       ////
////            Jonathan Biteau - June 2008                ////
////                                                       ////
///////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////
////                                                       ////  
////  USEFUL FUNCTIONS FOR THE CORRECTION IN THE BARREL    ////
////                                                       ////
///////////////////////////////////////////////////////////////


//useful to compute the signed distance to the closest crack in the barrel
double
PFEnergyCalibration::minimum(double a,double b){
  if(TMath::Abs(b)<TMath::Abs(a)) a=b;
  return a;
}


//compute the unsigned distance to the closest phi-crack in the barrel
double
PFEnergyCalibration::dCrackPhi(double phi, double eta){

  static double pi= M_PI;// 3.14159265358979323846;
  
  //Location of the 18 phi-cracks
  static std::vector<double> cPhi;
  if(cPhi.size()==0)
    {
      cPhi.resize(18,0);
      cPhi[0]=2.97025;
      for(unsigned i=1;i<=17;++i) cPhi[i]=cPhi[0]-2*i*pi/18;
    }

  //Shift of this location if eta<0
  static double delta_cPhi=0.00638;

  double m; //the result

  //the location is shifted
  if(eta<0) phi +=delta_cPhi;

  if (phi>=-pi && phi<=pi){

    //the problem of the extrema
    if (phi<cPhi[17] || phi>=cPhi[0]){
      if (phi<0) phi+= 2*pi;
      m = minimum(phi -cPhi[0],phi-cPhi[17]-2*pi);        	
    }

    //between these extrema...
    else{
      bool OK = false;
      unsigned i=16;
      while(!OK){
	if (phi<cPhi[i]){
	  m=minimum(phi-cPhi[i+1],phi-cPhi[i]);
	  OK=true;
	}
	else i-=1;
      }
    }
  }
  else{
    m=0.;        //if there is a problem, we assum that we are in a crack
    std::cout<<"Problem in dminphi"<<std::endl;
  }
  if(eta<0) m=-m;   //because of the disymetry
  return m;
}

// corrects the effect of phi-cracks
double
PFEnergyCalibration::CorrPhi(double phi, double eta) {

  // we use 3 gaussians to correct the phi-cracks effect
  static double p1=   5.59379e-01;
  static double p2=   -1.26607e-03;
  static double p3=  9.61133e-04;

  static double p4=   1.81691e-01;
  static double p5=   -4.97535e-03;
  static double p6=   1.31006e-03;

  static double p7=   1.38498e-01;
  static double p8=   1.18599e-04;
  static double p9= 2.01858e-03;
  

  double dminphi = dCrackPhi(phi,eta);
  
  double result = (1+p1*TMath::Gaus(dminphi,p2,p3)+p4*TMath::Gaus(dminphi,p5,p6)+p7*TMath::Gaus(dminphi,p8,p9));

  return result;
}   


// corrects the effect of  |eta|-cracks
double
PFEnergyCalibration::CorrEta(double eta){
  
  // we use a gaussian with a screwness for each of the 5 |eta|-cracks
  static std::vector<double> a;  //amplitude
  static std::vector<double> m;  //mean
  static std::vector<double> s;  //sigma
  static std::vector<double> sa; // screwness amplitude
  static std::vector<double> ss; // screwness sigma

  if(a.size()==0)
    {
      a.push_back(6.13349e-01) ;a.push_back(5.08146e-01)  ;a.push_back(4.44480e-01) ;a.push_back(3.3487e-01)   ;a.push_back(7.65627e-01) ;
      m.push_back(-1.79514e-02);m.push_back(4.44747e-01)  ;m.push_back(7.92824e-01) ;m.push_back(1.14090e+00)  ;m.push_back(1.47464e+00) ;
      s.push_back(7.92382e-03) ;s.push_back(3.06028e-03)  ;s.push_back(3.36139e-03) ;s.push_back(3.94521e-03)  ;s.push_back(8.63950e-04) ;
      sa.push_back(1.27228e+01);sa.push_back(3.81517e-02) ;sa.push_back(1.63507e-01);sa.push_back(-6.56480e-02);sa.push_back(1.87160e-01);
      ss.push_back(5.48753e-02);ss.push_back(-1.00223e-02);ss.push_back(2.22866e-03);ss.push_back(4.26288e-04) ;ss.push_back(2.67937e-03);
    }
 double result = 1;

 for(unsigned i=0;i<=4;i++) result+=a[i]*TMath::Gaus(eta,m[i],s[i])*(1+sa[i]*TMath::Sign(1.,eta-m[i])*TMath::Exp(-TMath::Abs(eta-m[i])/ss[i]));

  return result;
}


//corrects the global behaviour in the barrel
double
PFEnergyCalibration::CorrBarrel(double E, double eta) {

  //Energy dependency
  static double p0=1.00000e+00;
  static double p1=3.27753e+01;
  static double p2=2.28552e-02;
  static double p3=3.06139e+00;
  static double p4=2.25135e-01;
  static double p5=1.47824e+00;
  static double p6=1.09e-02;
  static double p7=4.19343e+01;

  //Eta dependency
  static double p8=2.705593e-03;

  double result = (p0+1/(p1+p2*TMath::Power(E,p3))+p4*TMath::Exp(-E/p5)+p6*TMath::Exp(-E*E/(p7*p7)))*(1+p8*eta*eta);

  return result;
}



///////////////////////////////////////////////////////////////
////                                                       ////  
////  USEFUL FUNCTIONS FOR THE CORRECTION IN THE ENDCAPS   ////
////  Parameters tuned for:                                ////
////          dR(ClustersPS1,ClusterEcal) < 0.08           ////
////          dR(ClustersPS2,ClusterEcal) < 0.13           ////
////                                                       ////
///////////////////////////////////////////////////////////////


//Alpha, Beta, Gamma give the weight of each sub-detector (PS layer1, PS layer2 and Ecal) in the areas of the endcaps where there is a PS
// Etot = Beta*eEcal + Gamma*(ePS1 + Alpha*ePS2) 

double
PFEnergyCalibration::Alpha(double eta) {

  //Energy dependency
  static double p0 = 5.97621e-01;

  //Eta dependency
  static double p1 =-1.86407e-01;
  static double p2 = 3.85197e-01; 

  //so that <feta()> = 1
  static double norm = (p1+p2*(2.6+1.656)/2);

  double result = p0*(p1+p2*eta)/norm;

  return result;
}

double
PFEnergyCalibration::Beta(double E, double eta) {

 //Energy dependency
  static double p0 = 0.032;
  static double p1 = 9.70394e-02;
  static double p2 = 2.23072e+01;
  static double p3 = 100;

  //Eta dependency
  static double p4 = 1.02496e+00 ;
  static double p5 = -4.40176e-03 ;

  //so that <feta()> = 1
  static double norm = (p4+p5*(2.6+1.656)/2);

  double result = (1.0012+p0*TMath::Exp(-E/p3)+p1*TMath::Exp(-E/p2))*(p4+p5*eta)/norm;			  
  return result;
}


double
PFEnergyCalibration::Gamma(double etaEcal) {

 //Energy dependency
  static double p0 = 2.49752e-02;

  //Eta dependency
  static double p1 = 6.48816e-02;
  static double p2 = -1.59517e-02; 
 
  //so that <feta()> = 1
  static double norm = (p1+p2*(2.6+1.656)/2);

  double result = p0*(p1+p2*etaEcal)/norm;					  

  return result;
}



///////////////////////////////////////////////////////////////
////                                                       ////  
////   THE CORRECTIONS IN THE BARREL AND IN THE ENDCAPS    ////
////                                                       ////
///////////////////////////////////////////////////////////////


// returns the corrected energy in the barrel (0,1.48)
// Global Behaviour, phi and eta cracks are taken into account
double
PFEnergyCalibration::EcorrBarrel(double E, double eta, double phi,
				 bool crackCorrection ){

  // double result = E*CorrBarrel(E,eta)*CorrEta(eta)*CorrPhi(phi,eta);
  double correction = crackCorrection ? std::max(CorrEta(eta),CorrPhi(phi,eta)) : 1.;
  double result = E * CorrBarrel(E,eta) * correction;

  return result;
}


// returns the corrected energy in the area between the barrel and the PS (1.48,1.65)
double
PFEnergyCalibration::EcorrZoneBeforePS(double E, double eta){

 //Energy dependency
  static double p0 =1; 
  static double p1 =0.18;
  static double p2 =8.;

  //Eta dependency
  static double p3 =0.3;
  static double p4 =1.11;
  static double p5 =0.025;
  static double p6 =1.49;
  static double p7 =0.6;

  //so that <feta()> = 1
  static double norm = 1.21;

  double result = E*(p0+p1*TMath::Exp(-E/p2))*(p3+p4*TMath::Gaus(eta,p6,p5)+p7*eta)/norm;

  return result;
}


// returns the corrected energy in the PS (1.65,2.6)
// only when (ePS1>0)||(ePS2>0)
double
PFEnergyCalibration::EcorrPS(double eEcal,double ePS1,double ePS2,double etaEcal) {

  // gives the good weights to each subdetector
  double E = Beta(1.0155*eEcal+0.025*(ePS1+0.5976*ePS2)/9e-5,etaEcal)*eEcal+Gamma(etaEcal)*(ePS1+Alpha(etaEcal)*ePS2)/9e-5 ;

  //Correction of the residual energy dependency
  static double p0 = 1.00;
  static double p1 = 2.18;
  static double p2 =1.94;
  static double p3 =4.13;
  static double p4 =1.127;

  double result = E*(p0+p1*TMath::Exp(-E/p2)-p3*TMath::Exp(-E/p4));

  return result;
} 

// returns the corrected energy in the PS (1.65,2.6)
// only when (ePS1>0)||(ePS2>0)
double
PFEnergyCalibration::EcorrPS(double eEcal,double ePS1,double ePS2,double etaEcal,double & outputPS1, double & outputPS2) {

  // gives the good weights to each subdetector
  double gammaprime=Gamma(etaEcal)/9e-5;
  outputPS1=gammaprime*ePS1;
  outputPS2=gammaprime*Alpha(etaEcal)*ePS2;
  double E = Beta(1.0155*eEcal+0.025*(ePS1+0.5976*ePS2)/9e-5,etaEcal)*eEcal+outputPS1+outputPS2;

  //Correction of the residual energy dependency
  static double p0 = 1.00;
  static double p1 = 2.18;
  static double p2 =1.94;
  static double p3 =4.13;
  static double p4 =1.127;
  
  double corrfac=(p0+p1*TMath::Exp(-E/p2)-p3*TMath::Exp(-E/p4));
  outputPS1*=corrfac;
  outputPS2*=corrfac;
  double result = E*corrfac;

  return result;
} 


// returns the corrected energy in the PS (1.65,2.6)
// only when (ePS1=0)&&(ePS2=0)
double 
PFEnergyCalibration::EcorrPS_ePSNil(double eEcal,double eta){

  //Energy dependency
  static double p0= 1.02;
  static double p1= 0.165;
  static double p2= 6.5 ;
  static double p3=  2.1 ;

  //Eta dependency
  static double p4 = 1.02496e+00 ;
  static double p5 = -4.40176e-03 ;

  //so that <feta()> = 1
  static double norm = (p4+p5*(2.6+1.656)/2);

  double result = eEcal*(p0+p1*TMath::Exp(-TMath::Abs(eEcal-p3)/p2))*(p4+p5*eta)/norm;
		  
  return result;
}


// returns the corrected energy in the area between the end of the PS and the end of the endcap (2.6,2.98)
double
PFEnergyCalibration::EcorrZoneAfterPS(double E, double eta){

  //Energy dependency
  static double p0 =1; 
  static double p1 = 0.058;
  static double p2 =12.5;
  static double p3 =-1.05444e+00;
  static double p4 =-5.39557e+00;
  static double p5 =8.38444e+00;
  static double p6 = 6.10998e-01  ;

  //Eta dependency
  static double p7 =1.06161e+00;
  static double p8 = 0.41;
  static double p9 =2.918;
  static double p10 =0.0181;
  static double p11= 2.05;
  static double p12 =2.99;
  static double p13=0.0287;

  //so that <feta()> = 1
  static double norm=1.045;

  double result = E*(p0+p1*TMath::Exp(-(E-p3)/p2)+1/(p4+p5*TMath::Power(E,p6)))*(p7+p8*TMath::Gaus(eta,p9,p10)+p11*TMath::Gaus(eta,p12,p13))/norm;
  return result;
}




// returns the corrected energy everywhere
// this work should be improved between 1.479 and 1.52 (junction barrel-endcap)
double
PFEnergyCalibration::Ecorr(double eEcal,double ePS1,double ePS2,
			   double eta,double phi,
			   bool crackCorrection ) {

  static double endBarrel=1.48;
  static double beginingPS=1.65;
  static double endPS=2.6;
  static double endEndCap=2.98;
 
  double result=0;

  eta=TMath::Abs(eta);

  if(eEcal>0){
    if(eta <= endBarrel)                         result = EcorrBarrel(eEcal,eta,phi,crackCorrection);
    else if(eta <= beginingPS)                   result = EcorrZoneBeforePS(eEcal,eta);
    else if((eta < endPS) && ePS1==0 && ePS2==0) result = EcorrPS_ePSNil(eEcal,eta);
    else if(eta < endPS)                         result = EcorrPS(eEcal,ePS1,ePS2,eta);
    else if(eta < endEndCap)                     result = EcorrZoneAfterPS(eEcal,eta); 
    else result =eEcal;
  }
  else result = eEcal;// useful if eEcal=0 or eta>2.98
  //protection
  if(result<eEcal) result=eEcal;
  return result;
}

// returns the corrected energy everywhere
// this work should be improved between 1.479 and 1.52 (junction barrel-endcap)
double
PFEnergyCalibration::Ecorr(double eEcal,double ePS1,double ePS2,double eta,double phi,double& ps1,double&ps2,bool crackCorrection)  {

  static double endBarrel=1.48;
  static double beginingPS=1.65;
  static double endPS=2.6;
  static double endEndCap=2.98;
 
  double result=0;

  eta=TMath::Abs(eta);

  if(eEcal>0){
    if(eta <= endBarrel)                         result = EcorrBarrel(eEcal,eta,phi,crackCorrection);
    else if(eta <= beginingPS)                   result = EcorrZoneBeforePS(eEcal,eta);
    else if((eta < endPS) && ePS1==0 && ePS2==0) result = EcorrPS_ePSNil(eEcal,eta);
    else if(eta < endPS)                         result = EcorrPS(eEcal,ePS1,ePS2,eta,ps1,ps2);
    else if(eta < endEndCap)                     result = EcorrZoneAfterPS(eEcal,eta); 
    else result =eEcal;
  }
  else result = eEcal;// useful if eEcal=0 or eta>2.98
  // protection
  if(result<eEcal) result=eEcal;
  return result;
}
