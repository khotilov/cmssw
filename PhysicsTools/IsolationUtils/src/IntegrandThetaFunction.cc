#include "PhysicsTools/IsolationUtils/interface/IntegrandThetaFunction.h"

// -*- C++ -*-
//
// Package:    IntegrandThetaFunction
// Class:      IntegrandThetaFunction
// 
/**\class IntegrandThetaFunction IntegrandThetaFunction.cc PhysicsTools/IsolationUtils/src/IntegrandThetaFunction.cc

 Description: auxialiary class for fixed area isolation cone computation
              (this class performs the integration over the polar angle)

 Implementation:
     imported into CMSSW on 05/18/2007
*/
//
// Original Author:  Christian Veelken, UC Davis
//         Created:  Thu Nov  2 13:47:40 CST 2006
// $Id: IntegrandThetaFunction.cc,v 1.3 2009/01/14 10:53:15 hegner Exp $
//
//

// system include files
#include <iostream>
#include <iomanip>
#include <vector>

// ROOT include files
#include <TMath.h>

// CMSSW include files
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/normalizedPhi.h"

//
// constructors and destructor
//

IntegrandThetaFunction::IntegrandThetaFunction()
  : fPhi_()
{
  theta0_ = 0.; 
  phi0_ = 0.; 
  alpha_ = 0.; 
}

IntegrandThetaFunction::IntegrandThetaFunction(const IntegrandThetaFunction& bluePrint)
{
  theta0_ = bluePrint.theta0_;
  phi0_ = bluePrint.phi0_;
  alpha_ = bluePrint.alpha_;
  
  fPhi_ = bluePrint.fPhi_;
}

IntegrandThetaFunction::~IntegrandThetaFunction()
{}

//
// assignment operator
//

IntegrandThetaFunction& IntegrandThetaFunction::operator=(const IntegrandThetaFunction& bluePrint)
{
  theta0_ = bluePrint.theta0_;
  phi0_ = bluePrint.phi0_;
  alpha_ = bluePrint.alpha_;
  
  fPhi_ = bluePrint.fPhi_;

  return (*this);
}

//
// member functions
//

void IntegrandThetaFunction::SetParameterTheta0(double theta0)
{
  theta0_ = theta0;
}

void IntegrandThetaFunction::SetParameterPhi0(double phi0)
{
  phi0_ = normalizedPhi(phi0); // map azimuth angle into interval [-pi,+pi]
}

void IntegrandThetaFunction::SetParameterAlpha(double alpha)
{
  alpha_ = alpha;
}

double IntegrandThetaFunction::DoEval(double x) const
{
//--- return zero if theta either close  to zero or close to Pi
//    (numerical expressions might become "NaN"s)
  const double epsilon = 1.e-3;
  if ( x < epsilon || x > (TMath::Pi() - epsilon) ) return 0.;

//--- calculate trigonometric expressions
//    (dependend on angle theta;
//     polar angle of point within cone)
  double sinTheta = TMath::Sin(x);
  double cscTheta = 1./sinTheta;

  double detJacobi = -cscTheta; // partial derrivative dEta/dTheta (for constant particle density in tau id. isolation cone)
  //double detJacobi = 1.; // ordinary solid angle (for TESTING ONLY)

//--- evaluate integral over angle phi
//    (azimuth angle of point within cone)
  fPhi_.SetParameterTheta0(theta0_);
  fPhi_.SetParameterPhi0(phi0_);
  fPhi_.SetParameterAlpha(alpha_);

  double integralOverPhi = fPhi_(x);		

  if ( debugLevel_ > 0 ) {
    edm::LogVerbatim("") << "integralOverPhi = " << integralOverPhi << std::endl
			 << " theta0 = " << theta0_ << std::endl
			 << " phi0 = " << phi0_ << std::endl
			 << " alpha = " << alpha_ << std::endl
			 << " theta = " << x << std::endl
			 << std::endl;
  }
  
//--- integrand for integration over theta
//    equals 
//      |dEta/dTheta| * integral over phi * sin(theta)
//
//    (o the factor dEta/dTheta represents the particle density as function of theta, 
//       assuming that the particle density is flat in eta;
//     o the factor sin(theta) originates from the solid angle surface element 
//       expressed in spherical polar coordinates)
//
  //return TMath::Abs(detJacobi)*integralOverPhi*sinTheta; (in case of TESTING with ordinary solid angle ONLY)
  return TMath::Abs(detJacobi)*integralOverPhi;
}  
