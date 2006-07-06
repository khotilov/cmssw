#ifndef RecoLocalMuon_CSCRecHitHistograms_H
#define RecoLocalMuon_CSCRecHitHistograms_H

/** \class CSCRecHitHistograms
 *  Collection of histograms for 1D DT RecHit test.
 *
 * Author: D. Fortin  - UC Riverside
 */


#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TString.h"
#include <string>
#include <iostream>


class H2DRecHit {
public:
  /// Constructor from collection name
  H2DRecHit(std::string name_) {
    TString N = name_.c_str();
    name=N;
// Local X position    
    hRecPositionX    = new TH1F(N+"_hRecPositionX", "CSCRecHit x local", 200, -24., 24.);
    hSimPositionX    = new TH1F(N+"_hSimPositionX", "CSCSimHit x local", 200, -24., 24.);
    hResPositionX    = new TH1F(N+"_hResPositionX", "CSCHit local Dx", 101, -0.202, 0.202);
    hResPositionXvsX = new TH2F(N+"_hResPositionXvsX", "CSCHit Dx vs x_{sim}", 200, -24., 24., 101, -2.02, 2.02); 
    hRecPositionYvsX = new TH2F(N+"_hRecPositionYvsX", "CSCRecHit y vs x", 600, -60., 60., 600, -60., 60.);
    hSimPositionYvsX = new TH2F(N+"_hSimPositionYvsX", "CSCSimHit y vs x", 600, -60., 60., 600, -60., 60.);
    hRecXvsSimX      = new TH2F(N+"_hRecXvsSimX", "CSCHit x_{rec} vs x_{sim}", 200, -40., 40., 200, -40., 40.);
// Local Y position
    hRecPositionY    = new TH1F(N+"_hRecPositionY", "CSCRecHit y local", 200, -40., 40.);
    hSimPositionY    = new TH1F(N+"_hSimPositionY", "CSCSimHit y local", 200, -40., 40.);
    hResPositionY    = new TH1F(N+"_hResPositionY", "CSCHit Dy", 101, -2.02, 2.02);  
    hResPositionYvsY = new TH2F(N+"_hResPositionYvsY", "CSCHit Dy vs y", 120, -60., 60., 101, -2.02, 2.02);
    hRecYvsSimY      = new TH2F(N+"_hRecYvsSimY", "CSCRecHit y_{rec} vs y_{sim}", 200, -60., 60., 200, -60., 60.);
// Wire group
    hWireGrp         = new TH1F(N+"_hWireGrp", "CSCRecHit WireGrp", 71, -0.5, 70.5);
    hWireGrpDigi     = new TH1F(N+"_hWireGrpDigi", "CSCDigi WireGrp", 71, -0.5, 70.5);
    hResWireGrp      = new TH1F(N+"_hResWireGrp", "CSCHit DWireGrp", 141, -70.5, 70.5);
    hWireGrpvsX      = new TH2F(N+"_hWireGrpvsX", "CSCRecHit x vs WireGrp", 61, -0.5, 60.5, 200, -24., 24.);
    hWireGrpvsY      = new TH2F(N+"_hWireGrpvsY", "CSCRecHit y vs WireGrp", 61, -0.5, 60.5, 200, -60., 60.);
    hWireGrpvsYres   = new TH2F(N+"_hWireGrpvsYres", "CSCRecHit Dy vs WireGrp", 61, -0.5, 60.5, 101, -2.02, 2.02);
    hAllWireGrpDigi  = new TH1F(N+"_hAllWireGrpDigi", "CSCDigi WireGrp", 71, -0.5, 70.5);
    hSimWireGrp      = new TH1F(N+"_hSimWireGrp", "CSCSimHit WireGrp", 71, -0.5, 70.5);

// Strip position
    hStripNu         = new TH1F(N+"_hStripNu", "CSCRecHit StripNumber", 81, -0.5, 80.5);
    hStripNuDigi     = new TH1F(N+"_hStripNuDigi", "CSCDigi StripNumber vs x", 81, -0.5, 80.5);
    hResStripNu      = new TH1F(N+"_hResStripNu", "CSCHit DStripNumber ", 161, -80.5, 80.5);
    hStripNuvsX      = new TH2F(N+"_hStripNuvsX", "CSCRecHit x vs StripNumber", 71, -0.5, 70.5, 200, -24., 24.);
    hAllStripNuDigi  = new TH1F(N+"_hAllStripNuDigi", "CSCDigi StripNumber vs x", 81, -0.5, 80.5);
    hSimStripNu      = new TH1F(N+"_hSimStripNu", "CSCSimHit StripNumber", 81, -0.5, 80.5);
// Phi
    hRecphi          = new TH1F(N+"_hRecphi", "CSCRecHit phi", 180, -180., 180.);
    hSimphi          = new TH1F(N+"_hSimphi", "CSCSimHit phi", 180, -180., 180.);
    hResphivsSimphi  = new TH2F(N+"_hRecphivsSimphi", "CSCHit Dphi vs phi_{sim}", 180, -180., 180., 101, -5.05, 5.05);
    hrDphi           = new TH1F(N+"_hrDphi", "CSCHit r x Dphi", 101, -0.101, 0.101);
// Eta
    hReceta          = new TH1F(N+"_hReceta", "CSCRecHit eta ", 150, 0.9, 2.4);
    hSimeta          = new TH1F(N+"_hSimeta", "CSCSimHit eta ", 150, 0.9, 2.4);
    hDeta            = new TH1F(N+"_hDeta", "CSC Deta ", 101, -0.0101, 0.0101);
    hDetavseta       = new TH2F(N+"_hDetavseta", "CSC Deta vs eta", 150, 0.9, 2.4, 101, -0.0101, 0.0101);
  }

  /// Constructor from collection name and TFile.
  H2DRecHit(TString name_, TFile* file) {
    name=name_;
// X
    hRecPositionX     = (TH1F *) file->Get(name+"_RecPositionX");
    hSimPositionX     = (TH1F *) file->Get(name+"_SimPositionX");
    hResPositionX     = (TH1F *) file->Get(name+"_ResPositionX");
    hResPositionXvsX  = (TH2F *) file->Get(name+"_ResPositionXvsX");
    hRecPositionYvsX  = (TH2F *) file->Get(name+"_RecPositionYvsX");
    hSimPositionYvsX  = (TH2F *) file->Get(name+"_SimPositionYvsX");
    hRecXvsSimX       = (TH2F *) file->Get(name+"_RecXvsSimX");
// Y
    hRecPositionY     = (TH1F *) file->Get(name+"_RecPositionY");
    hSimPositionY     = (TH1F *) file->Get(name+"_SimPositionY");
    hResPositionY     = (TH1F *) file->Get(name+"_ResPositionY");
    hResPositionYvsY  = (TH2F *) file->Get(name+"_ResPositionYvsY");
    hRecYvsSimY       = (TH2F *) file->Get(name+"_RecYvsSimY");  
// Wire group
    hWireGrp          = (TH1F *) file->Get(name+"_WireGrp");
    hWireGrpDigi      = (TH1F *) file->Get(name+"_WireGrpDigi");
    hResWireGrp       = (TH1F *) file->Get(name+"_ResWireGrp"); 
    hWireGrpvsX       = (TH2F *) file->Get(name+"_WireGrpvsX");
    hWireGrpvsY       = (TH2F *) file->Get(name+"_WireGrpvsY");
    hWireGrpvsYres    = (TH2F *) file->Get(name+"_WireGrpvsYres");
    hAllWireGrpDigi   = (TH1F *) file->Get(name+"_AllWireGrpDigi");
    hSimWireGrp       = (TH1F *) file->Get(name+"_SimWireGrp");
// Strip #
    hStripNu          = (TH1F *) file->Get(name+"_StripNu");
    hStripNuDigi      = (TH1F *) file->Get(name+"_StripNuDigi"); 
    hResStripNu       = (TH1F *) file->Get(name+"_ResStripNu");
    hStripNuvsX       = (TH2F *) file->Get(name+"_StripNuvsX");
    hAllStripNuDigi   = (TH1F *) file->Get(name+"_AllStripNuDigi");
    hSimStripNu       = (TH1F *) file->Get(name+"_SimStripNu");
// Phi
    hRecphi           = (TH1F *) file->Get(name+"_Recphi");
    hSimphi           = (TH1F *) file->Get(name+"_Simphi");
    hResphivsSimphi   = (TH2F *) file->Get(name+"_RecphivsSimphi");
    hrDphi            = (TH1F *) file->Get(name+"_rDphi");
// Eta
    hReceta           = (TH1F *) file->Get(name+"_Receta");
    hSimeta           = (TH1F *) file->Get(name+"_Simeta");  
    hDeta             = (TH1F *) file->Get(name+"_Deta");
    hDetavseta        = (TH2F *) file->Get(name+"_Detavseta");

  }


  /// Destructor
  virtual ~H2DRecHit() {
// X
    delete hRecPositionX;
    delete hSimPositionX;
    delete hResPositionX;
    delete hResPositionXvsX;
    delete hRecPositionYvsX;
    delete hSimPositionYvsX;
    delete hRecXvsSimX;
// Y
    delete hRecPositionY;
    delete hSimPositionY;
    delete hResPositionY;
    delete hResPositionYvsY;
    delete hRecYvsSimY;
// Wire group
    delete hWireGrp;
    delete hWireGrpDigi;
    delete hResWireGrp;
    delete hWireGrpvsX;
    delete hWireGrpvsY;
    delete hWireGrpvsYres;
    delete hAllWireGrpDigi;
    delete hSimWireGrp;
// Strip #
    delete hStripNu;
    delete hStripNuDigi;
    delete hResStripNu;
    delete hStripNuvsX;
    delete hAllStripNuDigi;
    delete hSimStripNu;
// Phi
    delete hRecphi;
    delete hSimphi;
    delete hResphivsSimphi;
    delete hrDphi;    
// Eta
    delete hReceta;
    delete hSimeta;
    delete hDeta;
    delete hDetavseta;
  }

  // Operations
  void FillWDigi(int wiregrp_digi) {hAllWireGrpDigi->Fill(wiregrp_digi);}
  void FillSDigi(int stripnum_digi) {hAllStripNuDigi->Fill(stripnum_digi);}

  /// Fill all the histos
  void Fill(float recx, float recy, float simx, float simy, float recphi, float simphi, float rdphi, float receta, 
            float simeta, float deta, int wiregrp, int stripnum, int wiregrp_digi, int stripnum_digi, 
            int strip_shit, int wiregrp_shit) {
	
    float radtodeg = 180./3.1415927;

// X
    hRecPositionX->Fill(recx);
    hSimPositionX->Fill(simx);
    hResPositionX->Fill(recx-simx);
    hResPositionXvsX->Fill(simx,recx-simx);
    hRecPositionYvsX->Fill(recx,recy);
    hSimPositionYvsX->Fill(simx,simy);
    hRecXvsSimX->Fill(simx,recx);
// Y  
    hRecPositionY->Fill(recy);
    hSimPositionY->Fill(simy);
    hResPositionY->Fill(recy-simy);
    hResPositionYvsY->Fill(simy,recy-simy);
    hRecYvsSimY->Fill(simy,recy);
// Wire group
    hWireGrp->Fill(wiregrp);
    hWireGrpDigi->Fill(wiregrp_digi);
    hResWireGrp->Fill(wiregrp-wiregrp_digi); 
    hWireGrpvsX->Fill(wiregrp,recx);
    hWireGrpvsY->Fill(wiregrp,recy);
    hWireGrpvsYres->Fill(wiregrp,recy-simy);
    hSimWireGrp->Fill(wiregrp_shit);
// Strip #
    hStripNu->Fill(stripnum);
    hStripNuDigi->Fill(stripnum_digi);
    hResStripNu->Fill(stripnum-stripnum_digi);
    hStripNuvsX->Fill(stripnum,recx);
    hSimStripNu->Fill(strip_shit);
// Phi
    hRecphi->Fill(recphi*radtodeg);
    hSimphi->Fill(simphi*radtodeg);
    hResphivsSimphi->Fill(simphi*radtodeg,(recphi-simphi)*radtodeg);
    hrDphi->Fill(rdphi);
// Eta
    hReceta->Fill(receta);
    hSimeta->Fill(simeta); 
    hDeta->Fill(receta-simeta);
    hDetavseta->Fill(simeta,receta-simeta); 
  }

  /// Write all the histos to currently opened file
  void Write() {
// X
    hRecPositionX->Write();
    hSimPositionX->Write();
    hResPositionX->Write();
    hResPositionXvsX->Write();
    hRecPositionYvsX->Write();
    hSimPositionYvsX->Write();
    hRecXvsSimX->Write();
// Y
    hRecPositionY->Write();
    hSimPositionY->Write();
    hResPositionY->Write();
    hResPositionYvsY->Write();
    hRecYvsSimY->Write();
// Wire group
    hWireGrp->Write();
    hWireGrpDigi->Write();
    hResWireGrp->Write();
    hWireGrpvsX->Write();
    hWireGrpvsY->Write();
    hWireGrpvsYres->Write();
    hAllWireGrpDigi->Write();
    hSimWireGrp->Write();
// Strip #
    hStripNu->Write();
    hStripNuDigi->Write();
    hResStripNu->Write();
    hStripNuvsX->Write();
    hAllStripNuDigi->Write();
    hSimStripNu->Write();
// Phi
    hRecphi->Write();
    hSimphi->Write();
    hResphivsSimphi->Write();
    hrDphi->Write();
// Eta
    hReceta->Write();
    hSimeta->Write();
    hDeta->Write();
    hDetavseta->Write();
  }

// X
    TH1F *hRecPositionX;
    TH1F *hSimPositionX;
    TH1F *hResPositionX;
    TH2F *hResPositionXvsX;
    TH2F *hRecPositionYvsX;
    TH2F *hSimPositionYvsX;
    TH2F *hRecXvsSimX;
// Y
    TH1F *hRecPositionY;
    TH1F *hSimPositionY;
    TH1F *hResPositionY;
    TH2F *hResPositionYvsY;
    TH2F *hRecYvsSimY;
// Wire group
    TH1F *hWireGrp;
    TH1F *hWireGrpDigi;
    TH1F *hResWireGrp;
    TH2F *hWireGrpvsX;
    TH2F *hWireGrpvsY;
    TH2F *hWireGrpvsYres;
    TH1F *hAllWireGrpDigi;
    TH1F *hSimWireGrp;
// Strip #
    TH1F *hStripNu;
    TH1F *hStripNuDigi;
    TH1F *hResStripNu;
    TH2F *hStripNuvsX;
    TH1F *hAllStripNuDigi;
    TH1F *hSimStripNu;
// Phi
    TH1F *hRecphi;
    TH1F *hSimphi;
    TH2F *hResphivsSimphi;
    TH1F *hrDphi;
// Eta
    TH1F *hReceta;
    TH1F *hSimeta;
    TH1F *hDeta;
    TH2F *hDetavseta;  

  TString name;
};
#endif


