#include "IOMC/CosmicMuonGenerator/interface/SingleParticleEvent.h"

void SingleParticleEvent::create(int id, double px, double py, double pz, double e, double m, double vx, double vy, double vz, double t0){
    ID = id;
    Px = px; Py = py; Pz = pz; E = e; M = m;
    Vx = vx; Vy = vy; Vz = vz; T0 = t0;
    HitTarget = false;
}

void SingleParticleEvent::propagate(double ElossScaleFac){
  // calculated propagation direction
  dX = Px/absmom();
  dY = Py/absmom(); 
  dZ = Pz/absmom();
  // propagate with decreasing step size
  tmpVx = Vx;
  tmpVy = Vy;
  tmpVz = Vz;
  HitTarget = true;
  if (HitTarget == true){
    HitTarget = false;
    double stepSize = MinStepSize*100000.;
    double acceptR = RadiusTarget + stepSize;
    double acceptZ = Z_DistTarget + stepSize;
    bool continuePropagation = true;
    while (continuePropagation){
      if (tmpVy < -acceptR) continuePropagation = false;
      if (absVzTmp() < acceptZ && rVxyTmp() < acceptR){
        HitTarget = true;
        continuePropagation = false;
      }
      if (continuePropagation) updateTmp(stepSize);
    }
  }
  if (HitTarget == true){
    HitTarget = false;
    double stepSize = MinStepSize*10000.;
    double acceptR = RadiusTarget + stepSize;
    double acceptZ = Z_DistTarget + stepSize;
    bool continuePropagation = true;
    while (continuePropagation){
      if (tmpVy < -acceptR) continuePropagation = false;
      if (absVzTmp() < acceptZ && rVxyTmp() < acceptR){
        HitTarget = true;
        continuePropagation = false;
      }
      if (continuePropagation) updateTmp(stepSize);
    }
  }
  if (HitTarget == true){
    HitTarget = false;
    double stepSize = MinStepSize*1000.;
    double acceptR = RadiusTarget + stepSize;
    double acceptZ = Z_DistTarget + stepSize;
    bool continuePropagation = true;
    while (continuePropagation){
      if (tmpVy < -acceptR) continuePropagation = false;
      if (absVzTmp() < acceptZ && rVxyTmp() < acceptR){
        HitTarget = true;
        continuePropagation = false;
      }
      if (continuePropagation) updateTmp(stepSize);
    }
  }
  if (HitTarget == true){
    HitTarget = false;
    double stepSize = MinStepSize*100.;
    double acceptR = RadiusTarget + stepSize;
    double acceptZ = Z_DistTarget + stepSize;
    bool continuePropagation = true;
    while (continuePropagation){
      if (tmpVy < -acceptR) continuePropagation = false;
      if (absVzTmp() < acceptZ && rVxyTmp() < acceptR){
        HitTarget = true;
        continuePropagation = false;
      }
      if (continuePropagation) updateTmp(stepSize);
    }
  }
  if (HitTarget == true){
    HitTarget = false;
    double stepSize = MinStepSize*10.;
    double acceptR = RadiusTarget + stepSize;
    double acceptZ = Z_DistTarget + stepSize;
    bool continuePropagation = true;
    while (continuePropagation){
      if (tmpVy < -acceptR) continuePropagation = false;
      if (absVzTmp() < acceptZ && rVxyTmp() < acceptR){
        HitTarget = true;
        continuePropagation = false;
      }
      if (continuePropagation) updateTmp(stepSize);
    }
  }
  if (HitTarget == true){
    HitTarget = false;
    double stepSize = MinStepSize*1.;
    double acceptR = RadiusTarget + stepSize;
    double acceptZ = Z_DistTarget + stepSize;
    bool continuePropagation = true;
    while (continuePropagation){
      if (tmpVy < -acceptR) continuePropagation = false;
      if (absVzTmp() < acceptZ && rVxyTmp() < acceptR){
        HitTarget = true;
        continuePropagation = false;
      }
      if (continuePropagation) updateTmp(stepSize);
    }
  }
  // actual propagation + energy loss
  if (HitTarget == true){
    HitTarget = false;
    int nAir = 0; int nWall = 0; int nRock = 0;
    double stepSize = MinStepSize*1.; // actual step size
    double acceptR = RadiusCMS + stepSize;
    double acceptZ = Z_DistCMS + stepSize;
    bool continuePropagation = true;
    while (continuePropagation){
      if (Vy < -acceptR) continuePropagation = false;
      if (absVz() < acceptZ && rVxy() < acceptR){
        HitTarget = true;
        continuePropagation = false;
      }
      if (continuePropagation) update(stepSize);
      if (inAir(Vx,Vy,Vz))  ++nAir;
      if (inWall(Vx,Vy,Vz)) ++nWall;
      if (inRock(Vx,Vy,Vz)) ++nRock;
    }
    if (HitTarget){
      double lAir  = double(nAir) *stepSize;
      double lWall = double(nWall)*stepSize;
      double lRock = double(nRock)*stepSize;
      double waterEquivalents = (lAir*RhoAir + lWall*RhoWall + lRock*RhoRock) *ElossScaleFac/10.; // [g cm^-2]
      subtractEloss(waterEquivalents);
      if (E < MuonMass) HitTarget = false; // muon stopped in the material around the target
    }
  }
  // end of propagation part
}

void SingleParticleEvent::update(double stepSize){
  Vx += stepSize*dX;
  Vy += stepSize*dY;
  Vz += stepSize*dZ;
}

void SingleParticleEvent::updateTmp(double stepSize){
  tmpVx += stepSize*dX;
  tmpVy += stepSize*dY;
  tmpVz += stepSize*dZ;
}

void SingleParticleEvent::subtractEloss(double waterEquivalents){
  double L10E = log10(E);
  // parameters for standard rock (PDG 2004, page 230)
  double A = (1.91514 + 0.254957*L10E)/1000.;                         // a [GeV g^-1 cm^2]
  double B = (0.379763 + 1.69516*L10E - 0.175026*L10E*L10E)/1000000.; // b [g^-1 cm^2]
  double EPS = A/B;                                                   // epsilon [GeV]
  E = (E + EPS)*exp(-B*waterEquivalents) - EPS; // updated energy
  double oldAbsMom = absmom();
  double newAbsMom = sqrt(E*E - MuonMass*MuonMass);
  Px = Px*newAbsMom/oldAbsMom;                  // updated px
  Py = Py*newAbsMom/oldAbsMom;                  // updated py
  Pz = Pz*newAbsMom/oldAbsMom;                  // updated pz
}

double SingleParticleEvent::absVzTmp(){
  return abs(tmpVz);
}

double SingleParticleEvent::rVxyTmp(){
  return sqrt(tmpVx*tmpVx + tmpVy*tmpVy);
}

bool SingleParticleEvent::hitTarget(){ return HitTarget; }

int    SingleParticleEvent::id(){ return ID; }

double SingleParticleEvent::px(){ return Px; }

double SingleParticleEvent::py(){ return Py; }

double SingleParticleEvent::pz(){ return Pz; }

double SingleParticleEvent::e(){ return E; }

double SingleParticleEvent::m(){ return M; }

double SingleParticleEvent::vx(){ return Vx; }

double SingleParticleEvent::vy(){ return Vy; }

double SingleParticleEvent::vz(){ return Vz; }

double SingleParticleEvent::t0(){ return T0; }

double SingleParticleEvent::phi(){
  double phiXZ = atan2(Px,Pz);
  if (phiXZ < 0.) phiXZ = phiXZ + TwoPi;
  return  phiXZ;
}

double SingleParticleEvent::theta(){
  return atan2(sqrt(Px*Px+Pz*Pz),-Py);
}

double SingleParticleEvent::absmom(){
  return sqrt(Px*Px + Py*Py + Pz*Pz);
}

double SingleParticleEvent::absVz(){
  return abs(Vz);
}

double SingleParticleEvent::rVxy(){
  return sqrt(Vx*Vx + Vy*Vy);
}
