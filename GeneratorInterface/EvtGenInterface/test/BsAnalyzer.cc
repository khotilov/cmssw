// -*- C++ -*-
//
// Package:    GeneratorInterface
// Class:      BsAnalyzer
// 
//
// Description: Module to analyze Pythia-EvtGen HepMCproducts
//
//
// Original Author:  Roberto Covarelli
//         Created:  April 26, 2007
//

#include <iostream>
#include "fstream.h"
 
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
 
// essentials !!!
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
 
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h" 
 
#include "FWCore/Framework/interface/MakerMacros.h"

#include "GeneratorInterface/EvtGenInterface/test/BsAnalyzer.h"

using namespace edm;
using namespace std;
using namespace HepMC;
 
BsAnalyzer::BsAnalyzer( const ParameterSet& pset )
   : fOutputFileName( pset.getUntrackedParameter<string>("HistOutFile",std::string("TestBs.root")) ),
     theSrc( pset.getUntrackedParameter<string>("theSrc",std::string("source")) ), 
     fOutputFile(0)
{
}

void BsAnalyzer::beginJob( const EventSetup& )
{
 
   nevent = 0;
   nbs = 0;

   fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   // fHist2muMass  = new TH1D(  "Hist2muMass", "2-mu inv. mass", 100,  60., 120. ) ;  
   hGeneralId = new TH1D( "hGeneralId","LundIDs of all particles",  100, -1000., 1000.) ;
   hnB = new TH1D( "hnB", "N(B)", 10,  0., 10. ) ;
   hnBz = new TH1D( "hnBz", "N(B0)", 10,  0., 10. ) ;
   hnBzb = new TH1D( "hnBzb", "N(B0bar)", 10,  0., 10. ) ;
   hPtbs = new TH1D( "hPtbs", "Pt Bs", 100,  0., 50. ) ;
   hPbs  = new TH1D( "hPbs",  "P Bs",  100,  0., 200. ) ;
   hPhibs = new TH1D( "hPhibs","Phi Bs",  100,  -3.14, 3.14) ;
   hEtabs = new TH1D( "hEtabs","Eta Bs",  100,  -7.0, 7.0) ;
   hPtmu = new TH1D( "hPtmu", "Pt Mu", 100,  0., 50. ) ;
   hPmu  = new TH1D( "hPmu",  "P Mu",  100,  0., 200. ) ;
   hPhimu = new TH1D( "hPhimu","Phi Mu",  100,  -3.14, 3.14) ;
   hEtamu = new TH1D( "hEtamu","Eta Mu",  100,  -7.0, 7.0) ;
   htbPlus = new TH1D( "htbPlus", "B+ proper decay time", 50, 0., 12. ) ;
   htbUnmix = new TH1D( "htbUnmix", "B0 proper decay time (unmixed)", 50, 0., 12. ) ;
   htbMix = new TH1D( "htbMix", "B0 proper decay time (mixed)", 50, 0., 12. ) ;
   htbMixPlus = new TH1D( "htbMixPlus", "B0 proper decay time (mixed positive)", 50, 0., 12. ) ;
   htbMixMinus = new TH1D( "htbMixMinus", "B0 proper decay time (mixed negative)", 50, 0., 12. ) ;
   htbsUnmix = new TH1D( "htbsUnmix", "Bs proper decay time (unmixed)", 50, 0., 12. ) ;
   htbsMix = new TH1D( "htbsMix", "Bs proper decay time (mixed)", 50, 0., 12. ) ;
   htbJpsiKs = new TH1D( "htbJpsiKs", "B0 -> J/#psiK_{s} decay time (B0 tags)", 50,  0., 12. ) ;
   htbbarJpsiKs = new TH1D( "htbbarJpsiKs", "B0 -> J/#psiK_{s} decay time (B0bar tags)", 50,  0., 12. ) ;
   hmumuMassSqr = new TH1D( "hmumuMassSqr","#mu^{+}#mu^{-} invariant mass squared",  100, -1.0, 25.0) ;
   hmumuMassSqrPlus = new TH1D( "hmumuMassSqrPlus","#mu^{+}#mu^{-} invariant mass squared (cos#theta > 0)",  100, -1.0, 25.0) ;
   hmumuMassSqrMinus = new TH1D( "hmumuMassSqrMinus","#mu^{+}#mu^{-} invariant mass squared (cos#theta < 0)",  100, -1.0, 25.0) ;
   hIdBsDaugs = new TH1D( "hIdBsDaugs","LundIDs of the Bs's daughters",  100, -1000., 1000.) ;
   hIdPhiDaugs = new TH1D( "hIdPhiDaugs","LundIDs of the phi's daughters",  100, -500., 500.) ;
   hIdBDaugs = new TH1D( "hIdBDaugs","LundIDs of the B's daughters",  100, -1000., 1000.) ;
   hCosHelAngleK = new TH1D( "hCosHelAngleK","cos#theta_{K}",  100, 0., 1.) ;
   hCosHelAngleKbkg = new TH1D( "hCosHelAngleKbkg","cos#theta_{K} background",  100, 0., 1.) ;

   decayed = new ofstream("decayed.txt") ;
   undecayed = new ofstream("undecayed.txt") ;
   return ;
}
 
void BsAnalyzer::analyze( const Event& e, const EventSetup& )
{
      
   Handle< HepMCProduct > EvtHandle ;
   
   // find initial HepMCProduct by its label - source
   // OR
   // find HepMCProduct after evtgenlhc by its label - evtgenproducer, that is
   // 
   e.getByLabel( theSrc , EvtHandle ) ;
   
   const GenEvent* Evt = EvtHandle->GetEvent() ;
   if (Evt) nevent++;

   const float mmcToPs = 3.3355; 
   int nB = 0;
   int nBp = 0;
   int nBz = 0;
   int nBzb = 0;
   int nBzmix = 0;
   int nBzunmix = 0;
   int nBzKmumu = 0;
   int nBJpsiKs = 0;
   int nBJpsiKstar = 0;

   for ( GenEvent::particle_const_iterator p = Evt->particles_begin(); p != Evt->particles_end(); ++p ) {

     // General
     GenVertex* endvert = (*p)->end_vertex(); 
     GenVertex* prodvert = (*p)->production_vertex();
     float gamma = (*p)->momentum().e()/(*p)->generated_mass();
     float dectime = 0.0 ;
     int mixed = -1;  // mixed is: -1 = unmixed
                      //           0 = mixed (before mixing)
                      //           1 = mixed (after mixing)
     if (endvert && prodvert) {
       dectime = (endvert->position().t() - prodvert->position().t())*mmcToPs/gamma;

       // Mixed particle ?
       for ( GenVertex::particles_in_const_iterator p2 = prodvert->particles_in_const_begin(); p2 != prodvert->particles_in_const_end(); ++p2 ) { 
	 if ( (*p)->pdg_id() + (*p2)->pdg_id() == 0) {
	   mixed = 1;
           gamma = (*p2)->momentum().e()/(*p2)->generated_mass();
           GenVertex* mixvert = (*p2)->production_vertex();
           dectime = (prodvert->position().t() - mixvert->position().t())*mmcToPs/gamma;
	 }
       }
       for ( GenVertex::particles_out_const_iterator ap = endvert->particles_out_const_begin(); ap != endvert->particles_out_const_end(); ++ap ) {
	 if ( (*p)->pdg_id() + (*ap)->pdg_id() == 0) mixed = 0;
       }
     }

     hGeneralId->Fill((*p)->pdg_id()); 

     // --------------------------------------------------------------
     if ( abs((*p)->pdg_id()/100) == 5 || abs((*p)->pdg_id()/100) == 4 || abs((*p)->pdg_id()/100) == 3) {
       if (!endvert) {
	 *undecayed << (*p)->pdg_id() << endl;
       } else {
	 *decayed << (*p)->pdg_id() << " --> " ;
         for ( GenVertex::particles_out_const_iterator bp = endvert->particles_out_const_begin(); bp != endvert->particles_out_const_end(); ++bp ) {
	   *decayed << (*bp)->pdg_id() << " " ;
         }
         *decayed << endl;
       }
     }
     // --------------------------------------------------------------
     if ( (*p)->pdg_id() == 531 /* && endvert */ ) {  // B_s 
       nbs++;
       hPtbs->Fill((*p)->momentum().perp());
       hPbs->Fill( sqrt ( pow((*p)->momentum().px(),2)+pow((*p)->momentum().py(),2)+
			  pow((*p)->momentum().pz(),2) )) ;
       hPhibs->Fill((*p)->momentum().phi());
       hEtabs->Fill((*p)->momentum().pseudoRapidity());
       
       for ( GenVertex::particles_out_const_iterator ap = endvert->particles_out_const_begin(); ap != endvert->particles_out_const_end(); ++ap ) {
	 hIdBsDaugs->Fill((*ap)->pdg_id());
       }
       
       if (mixed == 1) {
	 htbsMix->Fill( dectime );
       } else if (mixed == -1) {
         htbsUnmix->Fill( dectime );
       }
       
     }
     // --------------------------------------------------------------
     if ( abs((*p)->pdg_id()) == 511 /* && endvert */ ) {  // B0 
       if (mixed != 0) {
	 nB++;
         if ( (*p)->pdg_id() > 0 ) {
	   nBz++;
	   if ( mixed == 1 ) {nBzmix++;} else {nBzunmix++;}
	 } else {nBzb++;}         
       }
       int isJpsiKs = 0;
       int isKmumu = 0;
       int isSemilept = 0;
       for ( GenVertex::particles_out_const_iterator bp = endvert->particles_out_const_begin(); bp != endvert->particles_out_const_end(); ++bp ) {
	 hIdBDaugs->Fill((*bp)->pdg_id());
         if ( (*bp)->pdg_id() == 443 || (*bp)->pdg_id() == 310 ) isJpsiKs++ ; 
         if ( (*p)->pdg_id() > 0 && ( abs((*bp)->pdg_id()) == 313 || abs((*bp)->pdg_id()) == 13 )) isKmumu++ ; 
         if ( abs((*bp)->pdg_id()) == 11 || abs((*bp)->pdg_id()) == 13 || abs((*bp)->pdg_id()) == 15 ) isSemilept++ ;
       }
        
       if (isSemilept) {
	 if (mixed == 1) {
	   htbMix->Fill( dectime );
           if ( (*p)->pdg_id() > 0 ) { 
	     htbMixPlus->Fill( dectime );
	   } else {
             htbMixMinus->Fill( dectime );
	   }
	 } else if (mixed == -1) {
	   htbUnmix->Fill( dectime );
	 }
       }

       if (isJpsiKs == 2) {
         nBJpsiKs++;
	 if ( (*p)->pdg_id()*mixed < 0 ) { 
	   htbbarJpsiKs->Fill( dectime );
	 } else {
           htbJpsiKs->Fill( dectime );
	 }
       }
      
       if (isKmumu == 3) nBzKmumu++; 
     }
     // --------------------------------------------------------------
     if ( abs((*p)->pdg_id()) == 521 /* && endvert */ ) {  // B+
       nBp++;
       htbPlus->Fill( dectime );
       int isJpsiKstar = 0;
       for ( GenVertex::particles_out_const_iterator bp = endvert->particles_out_const_begin(); bp != endvert->particles_out_const_end(); ++bp ) {
         if ( (*bp)->pdg_id() == 443 || abs( (*bp)->pdg_id() ) == 323 ) isJpsiKstar++ ;          
       }
       if (isJpsiKstar == 2) nBJpsiKstar++;
     }
     // --------------------------------------------------------------
     if ( (*p)->pdg_id() == 333 ) {  // phi 
       if (endvert) {
         for ( GenVertex::particles_out_const_iterator cp = endvert->particles_out_const_begin(); cp != endvert->particles_out_const_end(); ++cp ) {
	   hIdPhiDaugs->Fill((*cp)->pdg_id());
	 }
       }
     }
     // --------------------------------------------------------------
     if ( (*p)->pdg_id() == 13 ) { // mu+
       for ( GenVertex::particles_in_const_iterator p2 = prodvert->particles_in_const_begin(); p2 != prodvert->particles_in_const_end(); ++p2 ) {
	 if ( abs((*p2)->pdg_id()) == 511 ) { // B0
	   hPtmu->Fill((*p)->momentum().perp());
	   hPmu->Fill( sqrt ( pow((*p)->momentum().px(),2)+pow((*p)->momentum().py(),2)+
			      pow((*p)->momentum().pz(),2) )) ;
	   hPhimu->Fill((*p)->momentum().phi());
	   hEtamu->Fill((*p)->momentum().pseudoRapidity());
	   for ( GenVertex::particles_out_const_iterator p3 = prodvert->particles_out_const_begin(); p3 != prodvert->particles_out_const_end(); ++p3 ) {
	     if ( (*p3)->pdg_id() == -13 ) { // mu-
	       TLorentzVector pmu1((*p)->momentum().px(), (*p)->momentum().py(),
				   (*p)->momentum().pz(), (*p)->momentum().e());
	       TLorentzVector pmu2((*p3)->momentum().px(), (*p3)->momentum().py(),
				   (*p3)->momentum().pz(), (*p3)->momentum().e());
               TLorentzVector pb0((*p2)->momentum().px(), (*p2)->momentum().py(),
				  (*p2)->momentum().pz(), (*p2)->momentum().e());
	       TLorentzVector ptot = pmu1 + pmu2;
               TVector3 booster = ptot.BoostVector();
	       TLorentzVector leptdir = ( ((*p2)->pdg_id() > 0) ? pmu1 : pmu2 );

	       leptdir.Boost( -booster );
	       pb0.Boost( -booster );
	       hmumuMassSqr->Fill(ptot.M2());
               if ( cos( leptdir.Vect().Angle(pb0.Vect()) ) > 0) {
		 hmumuMassSqrPlus->Fill(ptot.M2());
	       } else {
		 hmumuMassSqrMinus->Fill(ptot.M2());
	       }
	     }
	   }
	 }
       }
     }
     // --------------------------------------------------------------
     if ( (*p)->pdg_id() == 321 ) {   // K+
       GenVertex* phivert = (*p)->production_vertex();
       if (phivert) {
	 for ( GenVertex::particles_in_const_iterator p2 = phivert->particles_in_const_begin(); p2 != phivert->particles_in_const_end(); ++p2 ) {
	   GenVertex* bsvert = (*p2)->production_vertex();
           if (bsvert) {
	     for ( GenVertex::particles_in_const_iterator p3 = bsvert->particles_in_const_begin(); p3 != bsvert->particles_in_const_end(); ++p3 ) {
	       TLorentzVector pK1((*p)->momentum().px(), (*p)->momentum().py(),
				  (*p)->momentum().pz(), (*p)->momentum().e());
	       TLorentzVector pphi((*p2)->momentum().px(), (*p2)->momentum().py(),
				   (*p2)->momentum().pz(), (*p2)->momentum().e());
	       TLorentzVector pds((*p3)->momentum().px(), (*p3)->momentum().py(),
				  (*p3)->momentum().pz(), (*p3)->momentum().e());
	       TVector3 booster = pphi.BoostVector();
	       pK1.Boost( -booster );
	       pds.Boost( -booster );
	       if ( abs((*p3)->pdg_id()) == 431 && (*p2)->pdg_id() == 333) {   
		 // D_s -> phi
		 hCosHelAngleK->Fill( fabs( cos( pK1.Vect().Angle(pds.Vect())) ));
	       } else {
		 hCosHelAngleKbkg->Fill( fabs( cos( pK1.Vect().Angle(pds.Vect())) ));
	       }
	     }
	   }
	 }
       }
     }
   }
   // ---------------------------------------------------------

   hnB->Fill(nB);
   hnBz->Fill(nBz);
   hnBzb->Fill(nBzb);
   *undecayed << "-------------------------------" << endl;
   *decayed << "-------------------------------" << endl;

   // if (nBz > 0) std::cout << "nBz = " << nBz << " nBz (K*mu+mu-) = " << nBzKmumu << " nMix = " << nBzmix << std::endl;
   // if (nBz > 0 && nBzKmumu == 0) Evt->print();
   // if (nB > 0) std::cout << "nB = " << nB << " nBz (JPsi Ks) = " << nBJpsiKs << std::endl;
   // if (nBp > 0) std::cout << "nB = " << nBp << " nBz (JPsi Kstar) = " << nBJpsiKstar << std::endl;
   // if (nBp > 0 && nBJpsiKstar == 0) Evt->print();

   return ;   
}

void BsAnalyzer::endJob()
{
   
   fOutputFile->Write() ;
   fOutputFile->Close() ;
   cout << "N_events = " << nevent << "\n";
   cout << "N_Bs = " << nbs << "\n"; 
   return ;
}
 
DEFINE_FWK_MODULE(BsAnalyzer);
