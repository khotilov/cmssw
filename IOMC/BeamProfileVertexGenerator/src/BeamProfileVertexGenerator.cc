#include "IOMC/BeamProfileVertexGenerator/interface/BeamProfileVertexGenerator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <iostream>

BeamProfileVertexGenerator::BeamProfileVertexGenerator(const edm::ParameterSet & p,
                                                       const long& seed) : 
  BaseEventVertexGenerator(p,seed), myVertex(0), myRandom(0) {
  
  edm::ParameterSet vgenParam(p);
  meanX(vgenParam.getUntrackedParameter<double>("BeamMeanX",0.0)*mm);
  meanY(vgenParam.getUntrackedParameter<double>("BeamMeanY",0.0)*mm);
  beamPos(vgenParam.getUntrackedParameter<double>("BeamPosition",0.0)*mm);
  sigmaX(vgenParam.getUntrackedParameter<double>("BeamSigmaX",0.0)*mm);
  sigmaY(vgenParam.getUntrackedParameter<double>("BeamSigmaY",0.0)*mm);
  double fMinEta = vgenParam.getUntrackedParameter<double>("MinEta",-5.5);
  double fMaxEta = vgenParam.getUntrackedParameter<double>("MaxEta",5.5);
  double fMinPhi = vgenParam.getUntrackedParameter<double>("MinPhi",-3.14159265358979323846);
  double fMaxPhi = vgenParam.getUntrackedParameter<double>("MaxPhi", 3.14159265358979323846);
  eta(0.5*(fMinEta+fMaxEta));
  phi(0.5*(fMinPhi+fMaxPhi));
  setType(vgenParam.getUntrackedParameter<bool>("GaussianProfile",true));

  edm::LogInfo("VertexGenerator") << "BeamProfileVertexGenerator: with beam "
				  << "along eta = " << myEta << " (Theta = " 
				  << myTheta/deg << ") phi = " << myPhi/deg 
				  << " centred at (" << myMeanX << ", " 
				  << myMeanY << ", "  << myMeanZ << ") and "
				  << "spread (" << mySigmaX << ", " << mySigmaY
				  << ") of type Gaussian = " << myType;
}

BeamProfileVertexGenerator::~BeamProfileVertexGenerator() {
  if (myVertex) delete myVertex;
  //  if (myRandom) delete myRandom;
}

Hep3Vector * BeamProfileVertexGenerator::newVertex() {

  if (myVertex) delete myVertex;
  double aX, aY;
  if (myType) 
    aX = mySigmaX * (dynamic_cast<RandGauss*>(myRandom))->fire() + myMeanX;
  else
    aX = (dynamic_cast<RandFlat*>(myRandom))->fire(-0.5*mySigmaX,0.5*mySigmaX) + myMeanX ;
  double tX = 90.*deg + myTheta;
  double sX = sin(tX);
  if (fabs(sX)>1.e-12) sX = 1./sX;
  else                 sX = 1.;
  double fX = atan2(sX*cos(myTheta)*sin(myPhi),sX*cos(myTheta)*cos(myPhi));
  if (myType) 
    aY = mySigmaY * (dynamic_cast<RandGauss*>(myRandom))->fire() + myMeanY;
  else
    aY = (dynamic_cast<RandFlat*>(myRandom))->fire(-0.5*mySigmaY,0.5*mySigmaY) + myMeanY;
  double fY = 90.*deg + myPhi;
  double xp = aX*sin(tX)*cos(fX) +aY*cos(fY) +myMeanZ*sin(myTheta)*cos(myPhi);
  double yp = aX*sin(tX)*sin(fX) +aY*cos(fY) +myMeanZ*sin(myTheta)*sin(myPhi);
  double zp = aX*cos(tX)                     +myMeanZ*cos(myTheta);

  myVertex = new Hep3Vector(xp, yp, zp);
  LogDebug("VertexGenerator") << "BeamProfileVertexGenerator: Vertex created "
			      << "at " << *myVertex;
  return myVertex;
}

Hep3Vector * BeamProfileVertexGenerator::lastVertex() {
  return myVertex;
}

void BeamProfileVertexGenerator::sigmaX(double s) { 

  if (s>=0) {
    mySigmaX = s; 
  } else {
    edm::LogWarning("VertexGenerator") << "Warning BeamProfileVertexGenerator:"
				       << " Illegal resolution in X " << s
				       << "- set to default value 0 mm";
    mySigmaX = 0;
  }
}

void BeamProfileVertexGenerator::sigmaY(double s) { 

  if (s>=0) {
    mySigmaY = s; 
  } else {
    edm::LogWarning("VertexGenerator") << "Warning BeamProfileVertexGenerator:"
				       << " Illegal resolution in Y " << s
				       << "- set to default value 0 mm";
    mySigmaY = 0;
  }
}

void BeamProfileVertexGenerator::eta(double s) { 
  myEta   = s; 
  myTheta = 2.0*atan(exp(-myEta));
}

void BeamProfileVertexGenerator::setType(bool s) { 

  myType = s;
  if (myRandom) delete myRandom;
  
  if (myType == true)
    myRandom = new RandGauss(m_Engine);
  else
    myRandom = new RandFlat(m_Engine) ;
}
