 


#include <iostream>

#include "IOMC/GeneratorInterface/test/HZZ4muAnalyzer.h"
 
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
 
// essentials !!!
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
 
#include "TFile.h"
#include "TH1.h"
 
#include "FWCore/Framework/interface/MakerMacros.h"
 
using namespace edm;
using namespace std;

 
HZZ4muAnalyzer::HZZ4muAnalyzer( const ParameterSet& pset )
   : fOutputFileName( pset.getUntrackedParameter<string>("HistOutFile",std::string("TestHiggsMass.root")) ),
     fOutputFile(0), fHist2muMass(0), fHist4muMass(0), fHistZZMass(0)
{
}

void HZZ4muAnalyzer::beginJob( const EventSetup& )
{
 
   fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   fHist2muMass  = new TH1D(  "Hist2muMass", "2-mu inv. mass", 100,  60., 120. ) ;
   fHist4muMass  = new TH1D(  "Hist4muMass", "4-mu inv. mass", 100, 170., 210. ) ;
   fHistZZMass   = new TH1D(  "HistZZMass",  "ZZ inv. mass",   100, 170., 210. ) ;    
 
   return ;
}
 
void HZZ4muAnalyzer::analyze( const Event& e, const EventSetup& )
{
   
   Handle< HepMCProduct > EvtHandle ;
   
   // find initial (unsmeared, unfiltered,...) HepMCProduct
   // by its label - PythiaSource, that is
   // 
   // e.getByLabel( "PythiaSource", EvtHandle ) ;
   e.getByLabel( "source", EvtHandle ) ;
   
   const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
   
   // this a pointer - and not an array/vector/... 
   // because this example explicitely assumes
   // that there one and only Higgs in the record
   //
   HepMC::GenVertex* HiggsDecVtx = 0 ;
	 
   // find the 1st vertex with outgoing Higgs 
   // and get Higgs decay vertex from there;
   //
   // in principal, one can look for the vertex 
   // with incoming Higgs as well...
   //
   for ( HepMC::GenEvent::vertex_const_iterator
         vit=Evt->vertices_begin(); vit!=Evt->vertices_end(); vit++ )
   {
      for ( HepMC::GenVertex::particles_out_const_iterator
            pout=(*vit)->particles_out_const_begin();
	    pout!=(*vit)->particles_out_const_end(); pout++ )
      {
         if ( (*pout)->pdg_id() == 25 ) 
         {
	    if ( (*pout)->end_vertex() != 0 )
	    {
	       HiggsDecVtx = (*pout)->end_vertex() ;
	       break ;
	    }
         }
      }
      if ( HiggsDecVtx != 0 )
      {
         // cout << " Higgs decay found ! " << endl ;
	 break ; // break the initial loop over vertices
      }
   }
   
   if ( HiggsDecVtx == 0 ) 
   {
      cout << " There is NO Higgs in this event ! " << endl ;
      return ;
   }
      
   
   // select and store stable descendants of the Higgs
   //   
   vector<HepMC::GenParticle*> StableHiggsDesc ;
   
   for ( HepMC::GenVertex::particle_iterator
         des=HiggsDecVtx->particles_begin(HepMC::descendants);
	 des!=HiggsDecVtx->particles_end(HepMC::descendants); des++ )
   {
      if ( (*des)->status() == 1 ) StableHiggsDesc.push_back(*des) ;
   }
   
   HepLorentzVector Mom2part ;
   HepLorentzVector Mom4part ;
   HepLorentzVector Mom2pairs ;
   double           XMass2part = 0.;
   double           XMass4part = 0.;
   double           XMass2pairs  = 0.;
   vector< HepLorentzVector > Mom2partCont ;
    
   // browse the array of stable descendants
   // and do 2-mu inv.mass
   //
   for ( unsigned int i=0; i<StableHiggsDesc.size(); i++ )
   {
      // skip other than mu
      //
      if ( abs(StableHiggsDesc[i]->pdg_id()) != 13 ) continue ; 
      
      for ( unsigned int j=i+1; j<StableHiggsDesc.size(); j++ )
      {
         // skip other than mu
	 //
	 if ( abs(StableHiggsDesc[j]->pdg_id()) != 13 ) continue ;
	 //
	 // skip same charge combo's
	 //
	 if ( (StableHiggsDesc[i]->pdg_id()*StableHiggsDesc[j]->pdg_id()) > 0 ) 
	    continue ;
	 //
	 // OK, opposite charges, do the job
	 //
	 Mom2part = StableHiggsDesc[i]->momentum() + StableHiggsDesc[j]->momentum() ;
	 XMass2part = Mom2part.m() ;
	 fHist2muMass->Fill( XMass2part ) ;
	 //cout << " counters : " << StableHiggsDesc[i]->barcode() << " " 
	 //                       << StableHiggsDesc[j]->barcode() 
	 //			<< " -> 2-part mass = " << XMass2part << endl ;
	 //
         // store if 2-part. inv. mass fits into (roughly) Z-mass interval 
	 //
	 if ( XMass2part > 80. && XMass2part < 100. )
	 {
	    Mom2partCont.push_back(Mom2part) ;
	 } 
      }
   }
   
   // make 4-part inv.mass
   //
   if ( StableHiggsDesc.size() == 4 )
   {
      for ( unsigned int i=0; i<StableHiggsDesc.size(); i++ )
      {
         Mom4part += StableHiggsDesc[i]->momentum() ;
      }
      XMass4part = Mom4part.m() ;
      fHist4muMass->Fill( XMass4part ) ;
   }
   //cout << " 4-part inv. mass = " << XMass4part << endl ;
   
   // make 2-pairs (ZZ) inv.mass
   //
   //cout << " selected Z-candidates in this event : " << Mom2partCont.size() << endl ;
   for ( unsigned int i=0; i<Mom2partCont.size(); i++ )
   {
      for ( unsigned int j=i+1; j<Mom2partCont.size(); j++ )
      {
         Mom2pairs = Mom2partCont[i] + Mom2partCont[j] ;
	 XMass2pairs = Mom2pairs.m() ;
	 fHistZZMass->Fill( XMass2pairs ) ;
         //cout << " 2-pairs (ZZ) inv. mass = " << XMass2pairs << endl ;
      }
   }
   
   return ;
   
}

void HZZ4muAnalyzer::endJob()
{
       
   fOutputFile->Write() ;
   fOutputFile->Close() ;
   
   return ;
}
 
DEFINE_FWK_MODULE(HZZ4muAnalyzer)
