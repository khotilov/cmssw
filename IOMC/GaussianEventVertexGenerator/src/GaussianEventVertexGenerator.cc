#include "IOMC/GaussianEventVertexGenerator/interface/GaussianEventVertexGenerator.h"
#include "Utilities/General/interface/CMSexception.h"

//#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <iostream>

using std::cout;
using std::endl;

GaussianEventVertexGenerator::GaussianEventVertexGenerator(const edm::ParameterSet & p,
                                                           const long& seed) 
: BaseEventVertexGenerator(p,seed), 
  m_pGaussianEventVertexGenerator(p), 
  myVertex(0)
{ 
  
  myRandom = new RandGauss(m_Engine);
  
  myMeanX = m_pGaussianEventVertexGenerator.getParameter<double>("MeanX")*mm;
  myMeanY = m_pGaussianEventVertexGenerator.getParameter<double>("MeanY")*mm;
  myMeanZ = m_pGaussianEventVertexGenerator.getParameter<double>("MeanZ")*mm;
  mySigmaX = m_pGaussianEventVertexGenerator.getParameter<double>("SigmaX")*mm;
  mySigmaY = m_pGaussianEventVertexGenerator.getParameter<double>("SigmaY")*mm;
  mySigmaZ = m_pGaussianEventVertexGenerator.getParameter<double>("SigmaZ")*mm;

  if (mySigmaX <= 0) {
    cout << "Error in GaussianEventVertexGenerator: "
	 << "Illegal resolution in X - set to default value 1.0 mm"
	 << endl;
    BaseGenexception ex("GaussianEventVertexGenerator:Illegal resolution in X");
    throw ex;
    mySigmaX = 1; 
  }
  if (mySigmaY <= 0) {
    cout << "Error in GaussianEventVertexGenerator: "
	 << "Illegal resolution in Y - set to default value 1.0 mm"
	 << endl;
    BaseGenexception ex("GaussianEventVertexGenerator:Illegal resolution in Y");
    throw ex;
    mySigmaY = 1; 
  }
  if (mySigmaZ <= 0) {
    cout << "Error in GaussianEventVertexGenerator: "
	 << "Illegal resolution in Z - set to default value 1.0 mm"
	 << endl;
    BaseGenexception ex("GaussianEventVertexGenerator:Illegal resolution in Z");
    throw ex;
    mySigmaZ = 1; 
  }
}

GaussianEventVertexGenerator::~GaussianEventVertexGenerator() 
{
  delete myVertex;
  // I'm not deleting this, since the engine seems to have
  // been delete earlier; thus an attempt tp delete RandGauss
  // results in a core dump... 
  // I need to ask Marc/Jim how to do it right...
  //delete myRandom; 
}

Hep3Vector * GaussianEventVertexGenerator::newVertex() {
  delete myVertex;
  double aX,aY,aZ;
  //aX = mySigmaX * RandGauss::shoot() + myMeanX;
  //aY = mySigmaY * RandGauss::shoot() + myMeanY;
  //aZ = mySigmaZ * RandGauss::shoot() + myMeanZ;
  aX = mySigmaX * myRandom->fire() + myMeanX ;
  aY = mySigmaY * myRandom->fire() + myMeanY ;
  aZ = mySigmaZ * myRandom->fire() + myMeanZ ;
  myVertex = new Hep3Vector(aX, aY, aZ);
  return myVertex;
}

Hep3Vector * GaussianEventVertexGenerator::lastVertex() 
{
  return myVertex;
}

void GaussianEventVertexGenerator::sigmaX(double s) 
{ 
  if (s>=0 ) {
    mySigmaX=s; 
  }
  else {
    cout << "Error in GaussianEventVertexGenerator: "
	 << "Illegal resolution in X - set to default value 1.0 mm"
	 << endl;
    BaseGenexception ex("GaussianEventVertexGenerator:Illegal resolution in X");
    throw ex;
    mySigmaX=1.0;
  }
}

void GaussianEventVertexGenerator::sigmaY(double s) 
{ 
  if (s>=0 ) {
    mySigmaY=s; 
  }
  else {
    cout << "Error in GaussianEventVertexGenerator: "
	 << "Illegal resolution in Y - set to default value 1.0 mm"
	 << endl;
    BaseGenexception ex("GaussianEventVertexGenerator:Illegal resolution in Y");
    throw ex;
    mySigmaY=1.0;
  }
}

void GaussianEventVertexGenerator::sigmaZ(double s) 
{ 
  if (s>=0 ) {
    mySigmaZ=s; 
  }
  else {
    cout << "Error in GaussianEventVertexGenerator: "
	 << "Illegal resolution in Z - set to default value 1.0 mm"
	 << endl;
    BaseGenexception ex("GaussianEventVertexGenerator:Illegal resolution in Z");
    throw ex;
    mySigmaZ=1.0;
  }
}
