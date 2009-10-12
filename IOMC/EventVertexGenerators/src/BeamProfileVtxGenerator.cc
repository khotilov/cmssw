
// $Id: BeamProfileVtxGenerator.cc,v 1.8 2009/04/01 22:43:01 sunanda Exp $

#include "IOMC/EventVertexGenerators/interface/BeamProfileVtxGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "HepMC/SimpleVector.h"

#include<fstream>
#include<string>

BeamProfileVtxGenerator::BeamProfileVtxGenerator(const edm::ParameterSet & p) :
  BaseEvtVtxGenerator(p), fRandom(0) {
  
  meanX(p.getParameter<double>("BeamMeanX")*cm);
  meanY(p.getParameter<double>("BeamMeanY")*cm);
  beamPos(p.getParameter<double>("BeamPosition")*cm);
  sigmaX(p.getParameter<double>("BeamSigmaX")*cm);
  sigmaY(p.getParameter<double>("BeamSigmaY")*cm);
  double fMinEta = p.getParameter<double>("MinEta");
  double fMaxEta = p.getParameter<double>("MaxEta");
  double fMinPhi = p.getParameter<double>("MinPhi");
  double fMaxPhi = p.getParameter<double>("MaxPhi");
  eta(0.5*(fMinEta+fMaxEta));
  phi(0.5*(fMinPhi+fMaxPhi));
  nBinx = p.getParameter<int>("BinX");
  nBiny = p.getParameter<int>("BinY");
  ffile = p.getParameter<bool>("UseFile");
  fTimeOffset = p.getParameter<double>("TimeOffset")*ns*c_light;
  
  if (ffile) {
    std::string file = p.getParameter<std::string>("File");
    ifstream is(file.c_str(), std::ios::in);
    if (is) {
      double elem,sum=0;
      while (!is.eof()) {
	is >> elem;
        fdistn.push_back(elem);
	sum += elem;
      }
      if (((int)(fdistn.size())) >= nBinx*nBiny) {
	double last = 0;
	for (unsigned int i=0; i<fdistn.size(); i++) {
	  fdistn[i] /= sum;
	  fdistn[i] += last;
	  last       = fdistn[i];
	}
	setType(false);
      } else {
	ffile = false;
      }
    } else {
      ffile = false;
    }
  } 
  if (!ffile) {
    setType(p.getParameter<bool>("GaussianProfile"));
  }

  edm::LogInfo("VertexGenerator") << "BeamProfileVtxGenerator: with "
				  << "beam along eta = " << fEta 
				  << " (Theta = " << fTheta/deg 
				  << ") phi = " << fPhi/deg 
				  << " centred at (" << fMeanX << ", " 
				  << fMeanY << ", "  << fMeanZ << ") "
				  << "and spread (" << fSigmaX << ", "
				  << fSigmaY << ") of type Gaussian = "
				  << fType << " use file " << ffile;
  if (ffile) 
    edm::LogInfo("VertexGenerator") << "With " << nBinx << " bins "
				    << " along X and " << nBiny 
				    << " bins along Y";
}

BeamProfileVtxGenerator::~BeamProfileVtxGenerator() {
  delete fRandom;
}


//Hep3Vector * BeamProfileVtxGenerator::newVertex() {
HepMC::FourVector* BeamProfileVtxGenerator::newVertex() {
  double aX, aY;
  if (ffile) {
    double r1 = (dynamic_cast<CLHEP::RandFlat*>(fRandom))->fire();
    int ixy = 0, ix, iy;
    for (unsigned int i=0; i<fdistn.size(); i++) {
      if (r1 > fdistn[i]) ixy = i+1;
    }
    if (ixy >= nBinx*nBiny) {
      ix = nBinx-1; iy = nBiny-1;
    } else {
      ix = ixy%nBinx; iy = (ixy-ix)/nBinx;
    }
    aX = 0.5*(2*ix-nBinx+2*(dynamic_cast<CLHEP::RandFlat*>(fRandom))->fire())*fSigmaX + fMeanX ;
    aY = 0.5*(2*iy-nBiny+2*(dynamic_cast<CLHEP::RandFlat*>(fRandom))->fire())*fSigmaY + fMeanY ;
  } else {
    if (fType) {
      aX = fSigmaX*(dynamic_cast<CLHEP::RandGaussQ*>(fRandom))->fire() +fMeanX;
      aY = fSigmaY*(dynamic_cast<CLHEP::RandGaussQ*>(fRandom))->fire() +fMeanY;
    } else {
      aX = (dynamic_cast<CLHEP::RandFlat*>(fRandom))->fire(-0.5*fSigmaX,0.5*fSigmaX) + fMeanX ;
      aY = (dynamic_cast<CLHEP::RandFlat*>(fRandom))->fire(-0.5*fSigmaY,0.5*fSigmaY) + fMeanY;
    }
  }
  double xp = -aX*cos(fTheta)*cos(fPhi) +aY*sin(fPhi) +fMeanZ*sin(fTheta)*cos(fPhi);
  double yp = -aX*cos(fTheta)*sin(fPhi) -aY*cos(fPhi) +fMeanZ*sin(fTheta)*sin(fPhi);
  double zp =  aX*sin(fTheta)                         +fMeanZ*cos(fTheta);

  //if (fVertex == 0) fVertex = new CLHEP::Hep3Vector;
  //fVertex->set(xp, yp, zp);
  if (fVertex == 0 ) fVertex = new HepMC::FourVector() ;
  fVertex->set(xp, yp, zp, fTimeOffset );

  LogDebug("VertexGenerator") << "BeamProfileVtxGenerator: Vertex created "
			      << "at (" << xp << ", " << yp << ", "
			      << zp << ", " << fTimeOffset << ")";
  return fVertex;
}

void BeamProfileVtxGenerator::sigmaX(double s) { 

  if (s>=0) {
    fSigmaX = s; 
  } else {
    edm::LogWarning("VertexGenerator") << "Warning BeamProfileVtxGenerator:"
				       << " Illegal resolution in X " << s
				       << "- set to default value 0 cm";
    fSigmaX = 0;
  }
}

void BeamProfileVtxGenerator::sigmaY(double s) { 

  if (s>=0) {
    fSigmaY = s; 
  } else {
    edm::LogWarning("VertexGenerator") << "Warning BeamProfileVtxGenerator:"
				       << " Illegal resolution in Y " << s
				       << "- set to default value 0 cm";
    fSigmaY = 0;
  }
}

void BeamProfileVtxGenerator::eta(double s) { 
  fEta   = s; 
  fTheta = 2.0*atan(exp(-fEta));
}

void BeamProfileVtxGenerator::setType(bool s) { 

  fType = s;
  delete fRandom;
  
  if (fType == true)
    fRandom = new CLHEP::RandGaussQ(getEngine());
  else
    fRandom = new CLHEP::RandFlat(getEngine());
}
