#ifndef MrParticle_h
#define MrParticle_h

/*  \class MrParticle
*
*  Generic particle class (to be used both for monte carlo and reconstucted objects)
*
*  Authors: Luc Pape & Filip Moortgat      Date: August 2005 
*/

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/TrackReco/interface/Track.h"

using namespace reco;
using namespace std;

class MrParticle {

public:

// constructors
MrParticle() : Px(0.), Py(0.), Pz(0.), E(0.), Ecorfact(1.), Corrected(false), Charge(0.),  
Vx(0.), Vy(0.), Vz(0.), D0Error(0.), DzError(0.),
ParticleType(0), ParticleIso(true), ParticleIsoValue(-1.), Hemi(0),
NumTracks(0), Pt_tracks(0.), Et_em(0.), Et_had(0.), EffID(0.), EffToBeFake(1.),
MCParton(0), PID(0), Status(0), Mother1(0), BtagDiscriminator(-100), TauTagDiscriminator(-100) {};

MrParticle(float px, float py, float pz, float e) : Px(px), Py(py), Pz(pz), E(e), 
Ecorfact(1.), Corrected(false), Charge(0.),  
Vx(0.), Vy(0.), Vz(0.), D0Error(0.), DzError(0.),
ParticleType(0), ParticleIso(true), ParticleIsoValue(-1.), Hemi(0),
NumTracks(0), Pt_tracks(0.), Et_em(0.), Et_had(0.), EffID(0.), EffToBeFake(1.),
MCParton(0), PID(0), Status(0), Mother1(0), BtagDiscriminator(-100), TauTagDiscriminator(-100) {};

MrParticle(MrParticle & p) : Px(p.px()), Py(p.py()), Pz(p.pz()), E(p.energy()), 
Ecorfact(p.ecorfactor()), Corrected(p.isECorrected()), Charge(p.charge()),  
Vx(p.vx()), Vy(p.vy()), Vz(p.vz()), D0Error(p.d0Error()), DzError(p.dzError()),
ParticleType(p.particleType()), ParticleIso(p.particleIso()), ParticleIsoValue(p.particleIsoValue()), Hemi(p.hemisphere()),
NumTracks(p.numTracks()), Pt_tracks(p.pt_tracks()), Et_em(p.et_em()), Et_had(p.et_had()), EffID(p.effID()), EffToBeFake(p.effToBeFake()),
MCParton(p.partonIndex()), PID(p.pid()), Status(p.status()), Mother1(p.motherIndex()), BtagDiscriminator(p.getBtagDiscriminator()), TauTagDiscriminator(p.getTauTagDiscriminator()) {};


// destructor
virtual ~MrParticle(){};

// access methods
// object kinematical quantities
float eta(){if (fabs(Pz) <1.0e-5 ) {return 0;}
            float theta = atan(sqrt(Px*Px+Py*Py)/Pz);
            if (theta < 0.) {theta = theta + 3.141592654;}
            return -log(tan(0.5*theta));}
float phi() {return atan2(Py,Px);}
float p() {return sqrt(Px*Px+Py*Py+Pz*Pz);}
float pt() {return sqrt(Px*Px+Py*Py);}
float px() {return Px;}
float py() {return Py;}
float pz() {return Pz;}
float energy() {return E;}
float ecorfactor() {return Ecorfact;}
float mass() {return sqrt(E*E-Px*Px-Py*Py-Pz*Pz);}
float charge() {return Charge;}
// two-particle invariant mass
float invmass(MrParticle* part) {
   return sqrt( (E+part->energy())*(E+part->energy()) -
                (Px+part->px())*(Px+part->px())       -
                (Py+part->py())*(Py+part->py())       -
                (Pz+part->pz())*(Pz+part->pz()) );}
// coordinates of reference point (closest approach to beam)
float vx() {return Vx;}
float vy() {return Vy;}
float vz() {return Vz;}
// impact parameter errors
float d0Error() {return D0Error;}
float dzError() {return DzError;}
// particleType: particle identification   
//                   10 for good e
//                     11 for e with bad H/E
//                     12 for e with bad shower shape
//                     13 for e with bad matching to track
//                   20 for good mu
//                     21 for bad mu
//                   30 for good tau (hadronic)
//                     31 for bad tau
//                   40 for good photon
//                     41 for photon with bad H/E
//                     42 for photon with bad shower shape
//                   50 for good udscg-jet
//                     51 for electromagnetic fraction in udscg-jet
//                     52 for charged fraction in udscg-jet
//                   60 for good b-jet
//                   70 for good top
//                   80 for an Isolated Unidentified Particle (UFO)
//                     81 for bad UFO
//                   90 for invisible particle (only for MC data)
//
int particleType() {return ParticleType;}
// isolation flag
bool particleIso() {return ParticleIso;}
// isolation value
float particleIsoValue() {return ParticleIsoValue;}
// hemisphere (1 or 2) to which the object is associated
int hemisphere() {return Hemi;} 
// number of tracks (useful for jets) and scalar sum of their pt
int numTracks() {return NumTracks;}
float pt_tracks() {return Pt_tracks;}
// electromagnetic and hadronic transverse energies
float et_em() {return Et_em;}
float et_had() {return Et_had;}
float effID() {return EffID;}
float effToBeFake() {return EffToBeFake;}
// get b-tagging discriminator
double getBtagDiscriminator(){return BtagDiscriminator;}
// get tau-tagging discriminator
double getTauTagDiscriminator(){return TauTagDiscriminator;}
// for Reco: index in MCData of the matched MC particle
int partonIndex() {return MCParton;} 
// for MC: PDG particle identification
int pid() {return PID;} 
// for MC: Pythia status
int status() {return Status;}
// for MC: index in MCData of the mother
int motherIndex() {return Mother1;}
// are jet energy corrections applied or not
bool isECorrected() {return Corrected;}


virtual const GsfElectron* electronCandidate() {
//                   cout << "Pointer to electron candidate not defined." << endl;
                   return NULL;}
virtual const Muon* muonCandidate() {
                   cout << "Pointer to muon candidate not defined." << endl;
                   return NULL;}
virtual const Photon* photonCandidate() {
                   cout << "Pointer to photon candidate not defined." << endl;
                   return NULL;}
virtual const Jet* jetCandidate() {
                   cout << "Pointer to jet candidate not defined." << endl;
                   return NULL;}
virtual const JetTag* jetTag() {
                   cout << "Pointer to jetTag candidate not defined." << endl;
                   return NULL;}
virtual const Track* ufoTrack() {
                   cout << "Pointer to UFO track not defined." << endl;
                   return NULL;}

// set methods
void setPx(float px) {Px = px;}
void setPy(float py) {Py = py;}
void setPz(float pz) {Pz = pz;}
void setEnergy(float e) {E = e;}
void setEcorfactor(float ecor) {Ecorfact = ecor;}
void setCharge(float charge) {Charge = charge;}
void setVx(float vx) {Vx = vx;}
void setVy(float vy) {Vy = vy;}
void setVz(float vz) {Vz = vz;}
void setd0Error(float dd0) {D0Error = dd0;}
void setdzError(float ddz) {DzError = ddz;}
virtual void setParticleType(int ptype) {ParticleType = ptype;}
void setParticleIso(bool piso) {ParticleIso = piso;}
void setParticleIsoValue(float pisoval) {ParticleIsoValue = pisoval;}
void setHemisphere(int hem) {Hemi = hem;}
void setNumTracks(int ntk) {NumTracks = ntk;}
void setPt_tracks(float apt_tracks) {Pt_tracks = apt_tracks;}
void setEt_em(float aet_em){Et_em = aet_em;} 
void setEt_had(float aet_had){Et_had = aet_had;}
void setEffID(float effID) {EffID = effID;}
void setEffToBeFake(float effFake) {EffToBeFake = effFake;}
void setPartonIndex(int mcpart) {MCParton = mcpart;}
void setPID(int pid) {PID = pid; 
                      if (abs(PID) == 11) {ParticleType = 10;} 
                      if (abs(PID) == 13) {ParticleType = 20;} 
                      if (abs(PID) == 15) {ParticleType = 30;} 
                      if (abs(PID) == 22) {ParticleType = 40;} 
                      if (abs(PID) >= 1 && abs(PID) < 5) {ParticleType = 50;} 
		      if (abs(PID) == 5) {ParticleType = 60;} 
		      if (abs(PID) == 6) {ParticleType = 70;} 
} 
void setStatus(int stat) {Status = stat;}
void setMotherIndex(int m1) {Mother1 = m1;} 
void setBtagDiscriminator(double discri) {BtagDiscriminator = discri;}
void setTauTagDiscriminator(double discrim) {TauTagDiscriminator = discrim;}
void setFourVectorToCorrected (bool flag) {
   if (Corrected == false && flag == true) {
     E *= Ecorfact; 
     Px *= Ecorfact; 
     Py *= Ecorfact; 
     Pz *= Ecorfact;
     Corrected = true;
   } else if (Corrected == true && flag == false) {
     E = E / Ecorfact; 
     Px = Px / Ecorfact; 
     Py = Py / Ecorfact; 
     Pz = Pz / Ecorfact;
     Corrected = false;
  }
}


private:

// data members
float Px, Py, Pz, E;
float Ecorfact;
bool Corrected;
float Charge; 
float Vx, Vy, Vz, D0Error, DzError;
int ParticleType, ParticleIso;
float ParticleIsoValue;
int Hemi, NumTracks;
float Pt_tracks, Et_em, Et_had; 
float EffID, EffToBeFake;
int MCParton, PID, Status, Mother1;
double BtagDiscriminator, TauTagDiscriminator;

};




class MrElectron : public MrParticle {

public:

// constructors
MrElectron() : MrParticle(),  PCandidate(0) {MrParticle::setParticleType(10);};

MrElectron(float px, float py, float pz, float e, const GsfElectron* mycand) : 
 MrParticle(px,py,pz,e), PCandidate(mycand) 
 {MrParticle::setParticleType(10);
 MrParticle::setCharge(mycand->charge());}

virtual ~MrElectron() {};

//access methods
virtual const GsfElectron* electronCandidate() {return PCandidate;}

//set methods
void setCandidate (const GsfElectron* mycand) {PCandidate = mycand;}
virtual void setParticleType(int ptype) {
               if (ptype/10 != 1){ cout << "Changing type to non-electron not allowed." << endl;}
                else { MrParticle::setParticleType(ptype); }
                }

private:

// data members
const GsfElectron* PCandidate;


};



class MrMuon : public MrParticle {

public:

// constructors
MrMuon() : MrParticle(),  PCandidate(0) {MrParticle::setParticleType(20);};

MrMuon(float px, float py, float pz, float e, const Muon* mycand) : 
 MrParticle(px,py,pz,e), PCandidate(mycand) 
 {MrParticle::setParticleType(20);
 MrParticle::setCharge(mycand->charge());}

virtual ~MrMuon() {};

//access methods
virtual const Muon* muonCandidate() {return PCandidate;}

//set methods
void setCandidate (const Muon* mycand) {PCandidate = mycand;}
virtual void setParticleType(int ptype) {
               if (ptype/10 != 2){ cout << "Changing type to non-muon not allowed." << endl;}
                else { MrParticle::setParticleType(ptype); }      
                }

private:

// data members
const Muon* PCandidate;

};



class MrPhoton : public MrParticle {

public:

// constructors
MrPhoton() : MrParticle(),  PCandidate(0) {MrParticle::setParticleType(40);};

MrPhoton(float px, float py, float pz, float e, const Photon* mycand) : 
 MrParticle(px,py,pz,e), PCandidate(mycand) 
 {MrParticle::setParticleType(40);
 MrParticle::setCharge(0.);}

virtual ~MrPhoton() {};

//access methods
virtual const Photon* photonCandidate() {return PCandidate;}

//set methods
void setCandidate (const Photon* mycand) {PCandidate = mycand;}
virtual void setParticleType(int ptype) {
               if (ptype/10 != 4){ cout << "Changing type to non-photon not allowed." << endl;}
                else { MrParticle::setParticleType(ptype); }
                }

private:

// data members
const Photon* PCandidate;


};


class MrJet : public MrParticle {

public:

// constructors
MrJet() : MrParticle(),  PCandidate(0) {MrParticle::setParticleType(50);};

MrJet(float px, float py, float pz, float e, const Jet* mycand, const JetTag* myjettag) : MrParticle(px,py,pz,e), 
 PCandidate(mycand), PJetTag(myjettag) {MrParticle::setParticleType(50);};


virtual ~MrJet() {};

//access methods
virtual const Jet* jetCandidate() {return PCandidate;}
virtual const JetTag* jetTag() {return PJetTag;}

//set methods
void setCandidate (const Jet* mycand) {PCandidate = mycand;}
virtual void setParticleType(int ptype) {
               int pptype = ptype/10;
               if (pptype <= 2 || pptype == 4 || pptype > 8){ cout << "Changing type to non-jet not allowed." << endl;}
               else { MrParticle::setParticleType(ptype); }
                }
//void setBtagDiscriminator(double discri) {MrParticle::setBtagDiscriminator(discri);}
//void setTauTagDiscriminator(double discrim) {MrParticle::setTauTagDiscriminator(discrim);}

private:

// data members
const Jet* PCandidate;
const JetTag* PJetTag;

};


class MrUFO : public MrParticle {

public:

// constructors
MrUFO() : MrParticle(),  PTrack(0) {MrParticle::setParticleType(80);};

MrUFO(float px, float py, float pz, float e, const Track* mytrack) : MrParticle(px,py,pz,e), 
 PTrack(mytrack) {MrParticle::setParticleType(80);
 MrParticle::setCharge(mytrack->charge());}


virtual ~MrUFO() {};

//access methods
virtual const Track* ufoTrack() {return PTrack;}

//set methods
void setTrack (const Track* mytrack) {PTrack = mytrack;}
virtual void setParticleType(int ptype) {
               if (ptype/10 != 8){ cout << "Changing type to non-UFO not allowed." << endl;}
                else { MrParticle::setParticleType(ptype); }      
                }

private:

// data members
const Track* PTrack;

};

// The structure Config_t has been introduced 
// to hold the entries to the constants,
// because the entry to the config cannot be transmitted to classes
// invoked from SusyAnalyzer.
// It may be a temporary solution.

struct Config_t {     
     edm::ParameterSet InputMC_params;
     edm::ParameterSet InputReco_params;
     edm::ParameterSet rejectEvent_params;
     edm::ParameterSet acceptance_cuts;
     edm::ParameterSet cleaner_params;
     edm::ParameterSet isolator_params;
     edm::ParameterSet objectmatch_params;
     edm::ParameterSet mcproc_params;
     edm::ParameterSet useranalysis_params;
    
};




#endif
