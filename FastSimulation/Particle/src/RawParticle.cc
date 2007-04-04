// -----------------------------------------------------------------------------
//  Prototype for a particle class
// -----------------------------------------------------------------------------
//  $Date: 2006/10/03 09:14:48 $
//  $Revision: 1.8 $
// -----------------------------------------------------------------------------
//  Author: Stephan Wynhoff - RWTH-Aachen (Email: Stephan.Wynhoff@cern.ch)
// -----------------------------------------------------------------------------
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "FastSimulation/Particle/interface/ParticleTable.h"
#include "FastSimulation/Particle/interface/RawParticle.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <algorithm>

using namespace HepPDT;

RawParticle::RawParticle() {
  init();
}

RawParticle::RawParticle(const HepLorentzVector& p) 
  : HepLorentzVector(p) {
  init();
}

RawParticle::RawParticle(const int id, const HepLorentzVector& p) 
  : HepLorentzVector(p) {
  this->init();
  this->setID(id);
}

RawParticle::RawParticle(const std::string name, const HepLorentzVector& p) 
  : HepLorentzVector(p) {
  this->init();
  this->setID(name);
}

RawParticle::RawParticle(const HepLorentzVector& p, const HepLorentzVector& xStart) {
  init();
  this->setPx(p.px());
  this->setPy(p.py());
  this->setPz(p.pz());
  this->setE(p.e());
  myVertex = xStart;
}

RawParticle::RawParticle(HepDouble px, HepDouble py, HepDouble pz, HepDouble e) {
  init();
  this->setPx(px);
  this->setPy(py);
  this->setPz(pz);
  this->setE(e);
}

RawParticle::RawParticle(const RawParticle &right) : HepLorentzVector() {
  this->setPx(right.px());
  this->setPy(right.py());
  this->setPz(right.pz());
  this->setE(right.e());
  myId     = right.myId; 
  myStatus = right.myStatus;
  myUsed   = right.myUsed;
  myCharge = right.myCharge;
  myMass   = right.myMass;
  myVertex = (right.myVertex);
  tab = (right.tab);
}

RawParticle::~RawParticle() {
  //  nParticles--;
}

RawParticle&  RawParticle::operator = (const RawParticle & right ) {
  //  cout << "Copy assignment " << endl;
  if (this != &right) { // don't copy into yourself
    this->setPx(right.px());
    this->setPy(right.py());
    this->setPz(right.pz());
    this->setE(right.e());
    myId     = right.myId; 
    myStatus = right.myStatus;
    myUsed   = right.myUsed;
    myCharge = right.myCharge;
    myMass   = right.myMass;
    myVertex = right.myVertex;
    tab      = right.tab;
  }
  return *this;
}

void RawParticle::init() {
  myId=0;  
  myStatus=99;
  myUsed=0;
  myCharge=0.;
  myMass=0.;
  tab = ParticleTable::instance();
}

void RawParticle::setID(const int id) {
  myId = id;
  if ( tab && tab->theTable()->particle(ParticleID(id)) != 0 ) { 
    myCharge = tab->theTable()->particle(ParticleID(id))->charge();
    myMass   = tab->theTable()->particle(ParticleID(id))->mass().value();
  }
}

void RawParticle::setID(const std::string name) {
  if ( tab && tab->theTable()->particle(name) != 0) {
    myId = (tab->theTable()->particle(name))->pid();
    myCharge = tab->theTable()->particle(ParticleID(myId))->charge();
    myMass   = tab->theTable()->particle(ParticleID(myId))->mass().value();
  } else {
    myId = 0;
  }
}

void RawParticle::setStatus(int istat) {
  myStatus = istat;
}

void RawParticle::setMass(float m) {
  myMass = m;
}

void RawParticle::setCharge(float q) {
  myCharge = q;
}

void RawParticle::chargeConjugate() {
  if ( tab && tab->theTable()->particle(ParticleID(-myId)) ) {
    myId = -myId;
    myCharge = -1*myCharge;
  }
}

void RawParticle::setT(const double t) {
    myVertex.setT(t);
}

void RawParticle::setVertex(const HepLorentzVector& vtx) {
  myVertex = vtx;
  //  cout << this->vertex() << endl;
  //    { pidInfo a;}
}

const HepLorentzVector&  RawParticle::vertex() const { 
    return   myVertex ; 
}

void RawParticle::rotate(HepDouble rphi, const Hep3Vector &raxis) {
  this->rotate(rphi, raxis);
  myVertex.rotate(rphi, raxis);
}

void RawParticle::rotateX(HepDouble rphi) {
  this->rotateX(rphi);
  myVertex.rotateX(rphi);
}

void RawParticle::rotateY(HepDouble rphi) {
  this->rotateY(rphi);
  myVertex.rotateY(rphi);
}

void RawParticle::rotateZ(HepDouble rphi) {
  this->rotateZ(rphi);
  myVertex.rotateZ(rphi);
}

std::string RawParticle::PDGname() const {
  std::string MyParticleName;
  //  HepString MyParticleName;
  if ( tab && tab->theTable()->particle(ParticleID(myId)) != 0) {
    MyParticleName = (tab->theTable()->particle(ParticleID(myId)))->name();
  } else {
    MyParticleName = "none";
  }
  return (std::string) MyParticleName;}


void RawParticle::printName() const {
  std::string MyParticleName;
  //  CLHep::HepString MyParticleName;
  if ( tab && tab->theTable()->particle(ParticleID(myId)) != 0) {
    MyParticleName = (tab->theTable()->particle(ParticleID(myId)))->name();
  }  if (MyParticleName.length() != 0) {
    std::cout <<  MyParticleName;
    for(unsigned int k=0;k<9-MyParticleName.length() && k<10; k++) 
      std::cout << " " ;
  }
  else
    std::cout << "unknown  ";
}

void RawParticle::print() const {
  printName();
  std::cout << std::setw(3) << status();
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.setf(std::ios::right, std::ios::adjustfield);
  std::cout << std::setw(8) << std::setprecision(2) << px();
  std::cout << std::setw(8) << std::setprecision(2) << py();
  std::cout << std::setw(8) << std::setprecision(2) << pz();
  std::cout << std::setw(8) << std::setprecision(2) << e();
  std::cout << std::setw(8) << std::setprecision(2) << m();
  std::cout << std::setw(8) << std::setprecision(2) << mass();
  std::cout << std::setw(8) << std::setprecision(2) << charge();
  std::cout << std::setw(8) << std::setprecision(2) << x();
  std::cout << std::setw(8) << std::setprecision(2) << y();
  std::cout << std::setw(8) << std::setprecision(2) << z();
  std::cout << std::setw(8) << std::setprecision(2) << t();
  std::cout << std::setw(0) << std::endl;
}

std::ostream& operator <<(std::ostream& o , const RawParticle& p) {

  o << p.pid() << " (";
  o << p.status() << "): ";
  o << (HepLorentzVector) p << " ";
  o << p.vertex();
  return o;

} 

HepDouble RawParticle::PDGcharge() const { 
  HepDouble q=-99999;
  if ( tab && tab->theTable()->particle(ParticleID(myId))) {
    q=tab->theTable()->particle(ParticleID(myId))->charge();
  }
  return q;
}

HepDouble RawParticle::PDGmass() const  { 
  HepDouble m=-99999;
  if (tab && tab->theTable()->particle(ParticleID(myId))) {
    m=tab->theTable()->particle(ParticleID(myId))->mass().value();
  }
  return m;
}

HepDouble RawParticle::PDGcTau() const {
  HepDouble ct=1E99;
  if (tab && tab->theTable()->particle(ParticleID(myId))) {

    // The lifetime is 0. in the Pythia Particle Data Table !
    //    ct=tab->theTable()->particle(ParticleID(myId))->lifetime().value();

    // Get it from the width (apparently Gamma/c!)
    HepDouble w 
      = tab->theTable()->particle(ParticleID(myId))->totalWidth().value();
    ct = w != 0. ? 6.582119e-25 / w / 10. : 1e99;   // ctau in cm 
  }

  /*
  std::cout << setw(20) << setprecision(18) 
       << "myId/ctau/width = " << myId << " " 
       << ct << " " << w << endl;  
  */

  return ct;
}


HepDouble RawParticle::et() const { 
  double mypp, tmpEt=-1.;
  
  mypp = std::sqrt(this->px()*this->px() + this->py()*this->py() + this->pz()*this->pz());
  if (mypp != 0) {
    tmpEt=this->e() * this->perp()/mypp;
  }
  return tmpEt;
}
