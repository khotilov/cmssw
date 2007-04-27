
// $Id: BetafuncEvtVtxGenerator.cc,v 1.4 2007/03/22 02:28:46 yarba Exp $
/*
________________________________________________________________________

 BetafuncEvtVtxGenerator

 Smear vertex according to the Beta function on the transverse plane
 and a Gaussian on the z axis. It allows the beam to have a crossing
 angle (slopes dxdz and dydz).

 Based on GaussEvtVtxGenerator
 implemented by Francisco Yumiceva (yumiceva@fnal.gov)

 FERMILAB
 2006
________________________________________________________________________
*/



#include "IOMC/EventVertexGenerators/interface/BetafuncEvtVtxGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Units/SystemOfUnits.h"
//#include "CLHEP/Vector/ThreeVector.h"
#include "HepMC/SimpleVector.h"

#include <iostream>



BetafuncEvtVtxGenerator::BetafuncEvtVtxGenerator(const edm::ParameterSet & p )
: BaseEvtVtxGenerator(p)
{ 
  
  fRandom = new CLHEP::RandGauss(getEngine());

  fX0 =        p.getParameter<double>("X0")*cm;
  fY0 =        p.getParameter<double>("Y0")*cm;
  fZ0 =        p.getParameter<double>("Z0")*cm;
  fSigmaZ =    p.getParameter<double>("SigmaZ")*cm;
  alpha_ =     p.getParameter<double>("Alpha")*radian;
  phi_ =       p.getParameter<double>("Phi")*radian;
  fbetastar =  p.getParameter<double>("BetaStar")*cm;
  femittance = p.getParameter<double>("Emittance")*cm; // this is not the normalized emittance
  
 
  if (fSigmaZ <= 0) {
	  throw cms::Exception("Configuration")
		  << "Error in BetafuncEvtVtxGenerator: "
		  << "Illegal resolution in Z (SigmaZ is negative)";
  }

  
}

BetafuncEvtVtxGenerator::~BetafuncEvtVtxGenerator() 
{
    delete fRandom; 
}

//Hep3Vector* BetafuncEvtVtxGenerator::newVertex() {
HepMC::FourVector* BetafuncEvtVtxGenerator::newVertex() {	
	double X,Y,Z;
	
	double tmp_sigz = fRandom->fire(0., fSigmaZ);
	Z = tmp_sigz + fZ0;

	double tmp_sigx = BetaFunction(tmp_sigz,fZ0); 
	X = fRandom->fire(0.,tmp_sigx) + fX0; // + Z*fdxdz ;

	double tmp_sigy = BetaFunction(tmp_sigz,fZ0);
	Y = fRandom->fire(0.,tmp_sigy) + fY0; // + Z*fdydz;
	  
	//if (fVertex == 0) fVertex = new CLHEP::Hep3Vector;
	//fVertex->set(X, Y, Z);
	if ( fVertex == 0 ) fVertex = new HepMC::FourVector();
	fVertex->set(X,Y,Z,0.);
		
	return fVertex;
}

double BetafuncEvtVtxGenerator::BetaFunction(double z, double z0)
{
	return sqrt(femittance*(fbetastar+(((z-z0)*(z-z0))/fbetastar)));

}


void BetafuncEvtVtxGenerator::sigmaZ(double s) 
{ 
	if (s>=0 ) {
		fSigmaZ=s; 
	}
	else {
		throw cms::Exception("LogicError")
			<< "Error in BetafuncEvtVtxGenerator::sigmaZ: "
			<< "Illegal resolution in Z (negative)";
	}
}

TMatrixD* BetafuncEvtVtxGenerator::GetInvLorentzBoost() {

	//alpha_ = 0;
	//phi_ = 142.e-6;
	
	if (boost_ != 0 ) return boost_;
	
	//boost_.ResizeTo(4,4);
	//boost_ = new TMatrixD(4,4);
	TMatrixD tmpboost(4,4);
	
	//if ( (alpha_ == 0) && (phi_==0) ) { boost_->Zero(); return boost_; }
	
	// Lorentz boost to frame where the collision is head-on
	// phi is the half crossing angle in the plane ZS
	// alpha is the angle to the S axis from the X axis in the XY plane
	
	tmpboost(0,0) = 1./cos(phi_);
	tmpboost(0,1) = - cos(alpha_)*sin(phi_);
	tmpboost(0,2) = - tan(phi_)*sin(phi_);
	tmpboost(0,3) = - sin(alpha_)*sin(phi_);
	tmpboost(1,0) = - cos(alpha_)*tan(phi_);
	tmpboost(1,1) = 1.;
	tmpboost(1,2) = cos(alpha_)*tan(phi_);
	tmpboost(1,3) = 0.;
	tmpboost(2,0) = 0.;
	tmpboost(2,1) = - cos(alpha_)*sin(phi_);
	tmpboost(2,2) = cos(phi_);
	tmpboost(2,3) = - sin(alpha_)*sin(phi_);
	tmpboost(3,0) = - sin(alpha_)*tan(phi_);
	tmpboost(3,1) = 0.;
	tmpboost(3,2) = sin(alpha_)*tan(phi_);
	tmpboost(3,3) = 1.;

	//std::cout << "matrix initialized: "<< std::endl;

	tmpboost.Invert();
	boost_ = new TMatrixD(tmpboost);
	//boost_->Print();

	
	return boost_;
}
