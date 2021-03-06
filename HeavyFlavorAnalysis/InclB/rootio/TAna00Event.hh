#ifndef TANA00EVENT
#define TANA00EVENT


#include "TObject.h"
#include "TClonesArray.h"
#include "TVector3.h"

#include "TGenCand.hh"
#include "TAnaTrack.hh"
#include "TAnaCand.hh"
#include "TAnaVertex.hh"
#include "TAnaJet.hh"

class TAna00Event : public TObject {

public:

  TAna00Event() { };
  TAna00Event(Int_t Option);
  virtual ~TAna00Event() { };
  virtual void  Clear(Option_t *option ="");
  void          dump();

  // ----------------------------------------------------------------------
  // -- Generator block
  int                 nGenCands() {return fnGenCands;}
  TGenCand*           getGenCand(int n);
  virtual TGenCand*   addGenCand();
  void                dumpGenBlock();
  // For a given 'SimTrack', find the index in the generator block by doing a (p,id) matching
  int                 getGenIndex(double px, double py, double pz, int id, double precision = 0.005);
  // Check whether a RecTrack/TAnaTrack has mother with ID in the generator block
  //  int                 isDescendant(TAnaTrack *pTrk, int ID, int matchCharge = 0);


  // -- RecTracks
  int                 nRecTracks() {return fnRecTracks;}
  TAnaTrack*          getRecTrack(int n);
  virtual TAnaTrack*  addRecTrack();

  // -- Signal Tracks (e.g. leptons)
  int                 nSigTracks() {return fnSigTracks;}
  TAnaTrack*          getSigTrack(int n);
  virtual TAnaTrack*  addSigTrack();

  // -- Simulated Tracks 
  int                 nSimTracks() {return fnSimTracks;}
  TAnaTrack*          getSimTrack(int n);
  virtual TAnaTrack*  addSimTrack();

  // -- Signal Candidates (e.g. muon+jet)
  int                 nCands() {return fnCandidates;}
  TAnaCand*           getCand(int n);
  virtual TAnaCand*   addCand();  

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  // -- D0 Candidates 
  int                 nD0Cands() {return fnD0Candidates;}
  TAnaCand*           getD0Cand(int n);
  virtual TAnaCand*   addD0Cand();
  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 
  // -- CaloJets
  int                 nCaloJets() {return fnCaloJets;}
  TAnaJet*            getCaloJet(int n);
  virtual TAnaJet*    addCaloJet();

  // -- GenJets
  int                 nGenJets() {return fnGenJets;}
  TAnaJet*            getGenJet(int n);
  virtual TAnaJet*    addGenJet();

// -- PFJets
  int                 nPFJets() {return fnPFJets;}
  TAnaJet*            getPFJet(int n);
  virtual TAnaJet*    addPFJet();


  // -- TrackJets
  int                 nTrackJets() {return fnTrackJets;}
  TAnaJet*            getTrackJet(int n);
  virtual TAnaJet*    addTrackJet();

  // -- Primary vertices
  int                 nPV()    {return fnPV;}
  TAnaVertex*         bestPV() {return getPV(fBestPV); }
  TAnaVertex*         getPV(int n);
  virtual TAnaVertex* addPV();

  // ----------------------------------------------------------------------
  int               fRunNumber, fEventNumber;
  int               fEventBits;
  double            fPtHat;
  
  int               fL1Decision, fHLTDecision;

  int               fL1w1, fL1w2, fL1w3, fL1w4;
  int               fHLTw1, fHLTw2, fHLTw3, fHLTw4, fHLTw5, fHLTw6, fHLTw7;
  
  int               fNb2,fNb3,fNc2,fNc3,fNs2,fNs3;


  double            fLumi; 
  int               fLumiSection; 
  int               fOrbit; 
  int               fBx; 
  
  int               fProcessID;
  double            fXsec;
  double            fFilterEff;
  double            fEventWeight;

  int               fEventTag;

  TAnaVertex        fPrimaryVertex;
  int               fnPrimaryVertices;
  TAnaVertex        fPrimaryVertex2;

  double          fGenMET, fMET0, fMET1;  // only x and y component are relevant. z could contain type information. 
  double          fGensumET, fsumET0, fsumET1;  // 0 is calo, 1 is PF

  double           fBeamwx, fBeamwy, fBeamex, fBeamey;  // transverse beamspot size and uncertainty
private:

  int               fnGenCands;
  TClonesArray      *fGenCands;

  int               fnRecTracks;
  TClonesArray      *fRecTracks;

  int               fnSigTracks;
  TClonesArray      *fSigTracks;

  int               fnSimTracks;
  TClonesArray      *fSimTracks;

  int               fnCaloJets;
  TClonesArray      *fCaloJets;

  int               fnGenJets;
  TClonesArray      *fGenJets;

  int               fnPFJets;
  TClonesArray      *fPFJets;

  int               fnTrackJets;
  TClonesArray      *fTrackJets;

  int               fnCandidates;
  TClonesArray      *fCandidates;  

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  int               fnD0Candidates;
  TClonesArray      *fD0Candidates;
  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  int               fBestPV;
  int               fnPV;
  TClonesArray      *fPV;

  ClassDef(TAna00Event,1)

};

#endif
