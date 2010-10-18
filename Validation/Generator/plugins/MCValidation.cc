#include "Validation/Generator/plugins/MCValidation.h"

using namespace std ;
using namespace edm ;
using namespace reco;


const double pi = acos(-1.) ;
const double twopi = 2*pi ;
const double eMass = 0.00511 ;
const double eMass2 = eMass*eMass ;
const double muMass = 0.106 ;
const double muMass2 = muMass*muMass ;
const double WMass = 80.403 ;
const double WMass2 = WMass*WMass ;
const double ZMass = 91.1876 ;
const double ZMass2 = ZMass*ZMass ;

MCValidation::MCValidation(const edm::ParameterSet& iConfig)

{
   
   s_genEvt = iConfig.getUntrackedParameter<string>("genEvt") ;
//   l_genEvt = iConfig.getParameter<edm::InputTag>("genEvt") ;

   theOutfileName = iConfig.getUntrackedParameter<string>("OutfileName") ;

}

MCValidation::~MCValidation()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void
MCValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   nEvt++ ;

   outfile << " anls fn Evt no. " << iEvent.id().event() << endl ;

//   Handle<reco::GenParticleCollection>genPar ;
//   iEvent.getByLabel(s_genEvt,genPar) ;
//   const GenParticleCollection &genColl = (*(genPar.product())) ;
   Handle<reco::CandidateCollection>genPar ;
   iEvent.getByLabel("genParticleCandidates",genPar) ;
   const CandidateCollection &genColl = (*(genPar.product())) ;
   genparticles(genColl) ;

} // end analyse 

void MCValidation::genparticles(const CandidateCollection &genParticles){ //GenParticleCollection &genParticles){

  for(unsigned i = 0; i < genParticles.size(); ++ i) {
  
//    const GenParticle & part = (genParticles)[i];   
    const Candidate & part = (genParticles)[i];
    
    int pid = part.pdgId(), pst = part.status(), pch = part.charge();
    unsigned nDau = part.numberOfDaughters() ;
    //unsigned nMom = part.numberOfMothers() ;
    double ppt = part.pt(), peta = part.eta(), pphi = part.phi(), pmass = part.mass();
    //double pvx = part.vx(), pvy = part.vy(), pvz = part.vz();

    outfile << " pid " << pid << " pmass " << pmass << " pst " << pst << endl ;
    
    switch(std::abs(pid)){
      
      case     11 : { // e+ OR e-
      
	  gen_e_pT_->Fill(pch*ppt) ;
	  
	  if(pch>0){
	    gen_eP_eta_->Fill(peta); gen_eP_phi_->Fill(pphi);	  
          }else if(pch<0){
	    gen_eM_eta_->Fill(peta); gen_eM_phi_->Fill(pphi);
	  }
	  
	  icode = 1 ;
	  break ;
      }

      case     13 : { // mu+ OR mu-
      
	  gen_e_pT_->Fill(pch*ppt) ;
	  
	  if(pch>0){
	    gen_eP_eta_->Fill(peta); gen_eP_phi_->Fill(pphi);	  
          }else if(pch<0){
	    gen_eM_eta_->Fill(peta); gen_eM_phi_->Fill(pphi);
	  }
	  
	  icode = 2 ;
	  break ;
      }

      case     22 : { // gamma
	
	  gen_gamma_pT_->Fill(pch*ppt) ; gen_gamma_eta_->Fill(peta); gen_gamma_phi_->Fill(pphi);  
	  
	  icode = 3 ;
    	  break ;
	
      } 
      
      case     24 : { // W+ OR W-
	
        if(pst!=3) break ;

	gen_W_pT_->Fill(ppt) ; 
        if(pch>0){
	  gen_Wp_eta_->Fill(peta) ; gen_Wp_phi_->Fill(pphi) ;
        }else if(pch<0){   
	  gen_Wm_eta_->Fill(peta) ; gen_Wm_phi_->Fill(pphi) ;
	}
	
	if(nDau > 1){
	
	  for(unsigned i = 0; i < nDau; ++i){
	  
	    const Candidate *dau = part.daughter(i) ;
	    int did = dau->pdgId(), dst = dau->status(), dch = dau->charge() ;
	    double dpt = dau->pt(), deta = dau->eta(), dphi = dau->phi();// dmass = dau->mass();
	    // double dvx = dau->vx(), dvy = dau->vy(), dvz = dau->vz();
	    
	    if(std::abs(did) == 11 && dst==1){
	      dau_e_pT_->Fill(dpt) ; 
	      if(dch>0){
  	        dau_eP_eta_->Fill(deta); dau_eP_phi_->Fill(dphi) ;
              }else if(dch<0){
	        dau_eM_eta_->Fill(deta); dau_eM_phi_->Fill(dphi) ;
	      } 	
	    }  
	    
	    if(std::abs(did)== 13 && dst==1){ 
              dau_mu_pT_->Fill(dpt); 
	      if(dch>0){
  	        dau_muP_eta_->Fill(deta); dau_muP_phi_->Fill(dphi) ;
              }else if(dch<0){
	        dau_muM_eta_->Fill(deta); dau_muM_phi_->Fill(dphi) ;
	      }
	    }  
	    
	  }
	  
	}

	icode = 4 ; 
	break ;
      
      }
      
      case     23 : { // Z0

	gen_Z_pT_->Fill(ppt) ; 
        gen_Z_eta_->Fill(peta) ; gen_Z_phi_->Fill(pphi) ;
	
	if(nDau > 1){
	
	  for(unsigned i = 0; i < nDau; ++i){
	  
	    const Candidate *dau = part.daughter(i) ;
	    int did = dau->pdgId(), dst = dau->status(), dch = dau->charge() ;
	    double dpt = dau->pt(), deta = dau->eta(), dphi = dau->phi(); //dmass = dau->mass();
	    //double dvx = dau->vx(), dvy = dau->vy(), dvz = dau->vz();
	    
	    if(std::abs(did) == 11 && dst==1){
	      dau_e_pT_->Fill(dpt) ; 
	      if(dch>0){
  	        dau_eP_eta_->Fill(deta); dau_eP_phi_->Fill(dphi) ;
              }else if(dch<0){
	        dau_eM_eta_->Fill(deta); dau_eM_phi_->Fill(dphi) ;
	      } 	
	    }  
	    
	    if(std::abs(did)== 13 && dst==1){ 
              dau_mu_pT_->Fill(dpt); 
	      if(dch>0){
  	        dau_muP_eta_->Fill(deta); dau_muP_phi_->Fill(dphi) ;
              }else if(dch<0){
	        dau_muM_eta_->Fill(deta); dau_muM_phi_->Fill(dphi) ;
	      }
	    }  
	    
	  }
	  
	}

	icode = 5 ;
	break ;

      }	
      
      case     35 : { // H0

	gen_H0_pT_->Fill(ppt) ; 
        gen_H0_eta_->Fill(peta) ; gen_H0_phi_->Fill(pphi) ;
	
	if(nDau > 1){
	
	  for(unsigned i = 0; i < nDau; ++i){
	  
	    const Candidate *dau = part.daughter(i) ;
	    int did = dau->pdgId(), dst = dau->status(), dch = dau->charge() ;
	    double dpt = dau->pt(), deta = dau->eta(), dphi = dau->phi();// dmass = dau->mass();
	    //double dvx = dau->vx(), dvy = dau->vy(), dvz = dau->vz();
	    
	    if(std::abs(did) == 11 && dst==1){
	      dau_e_pT_->Fill(dpt) ; 
	      if(dch>0){
  	        dau_eP_eta_->Fill(deta); dau_eP_phi_->Fill(dphi) ;
              }else if(dch<0){
	        dau_eM_eta_->Fill(deta); dau_eM_phi_->Fill(dphi) ;
	      } 	
	    }  
	    
	    if(std::abs(did)== 13 && dst==1){ 
              dau_mu_pT_->Fill(dpt); 
	      if(dch>0){
  	        dau_muP_eta_->Fill(deta); dau_muP_phi_->Fill(dphi) ;
              }else if(dch<0){
	        dau_muM_eta_->Fill(deta); dau_muM_phi_->Fill(dphi) ;
	      }
	    }  
	    
	    if(did== 22 && dst==1){ 
              dau_gamma_pT_->Fill(dpt); dau_gamma_eta_->Fill(deta); dau_gamma_phi_->Fill(dphi) ;
	    }  
	    
	  }
	  
	}

	icode = 6 ;
	break ;

      }	

      case     32 : { // Z'

	gen_Zprime_pT_->Fill(ppt) ; 
        gen_Zprime_eta_->Fill(peta) ; gen_Zprime_phi_->Fill(pphi) ;
	
	if(nDau > 1){
	
	  for(unsigned i = 0; i < nDau; ++i){
	  
	    const Candidate *dau = part.daughter(i) ;
	    int did = dau->pdgId(), dst = dau->status(), dch = dau->charge() ;
	    double dpt = dau->pt(), deta = dau->eta(), dphi = dau->phi(); //dmass = dau->mass();
	      //	    double dvx = dau->vx(), dvy = dau->vy(), dvz = dau->vz();
	    
	    if(std::abs(did) == 11 && dst==1){
	      dau_e_pT_->Fill(dpt) ; 
	      if(dch>0){
  	        dau_eP_eta_->Fill(deta); dau_eP_phi_->Fill(dphi) ;
              }else if(dch<0){
	        dau_eM_eta_->Fill(deta); dau_eM_phi_->Fill(dphi) ;
	      } 	
	    }  
	    
	    if(std::abs(did)== 13 && dst==1){ 
              dau_mu_pT_->Fill(dpt); 
	      if(dch>0){
  	        dau_muP_eta_->Fill(deta); dau_muP_phi_->Fill(dphi) ;
              }else if(dch<0){
	        dau_muM_eta_->Fill(deta); dau_muM_phi_->Fill(dphi) ;
	      }
	    }  
	    
	  }
	  
	}

	icode = 7 ;
	break ;

      }	
      
      default : icode = 0 ; break ;
      
    } // switch ends

  } // loop over gen particles ends    
    
} // end genparticles

void 
MCValidation::beginJob()
{
  // what is this needed for ?????
  //iSetup.get<IdealMagneticFieldRecord>().get(theMagneticField) ;

  outfile.open(theOutfileName.c_str()) ;

  edm::Service<TFileService> fs;

  Tgen = new TTree("Tgen", "GEN") ; 
  Tgen->Branch("iEvt", &iEvt, "iEvt/I") ;
  Tgen->Branch("iRun", &iRun, "iRun/I") ;

  mom_pT_          = fs->make<TH1F>("mom_pT_",         "p_{T} mother",              100, -50.,  50.)  ;
  mom_eta_         = fs->make<TH1F>("mom_eta_",        "#eta mother",                50,  -3.5,  3.5) ;
  mom_phi_         = fs->make<TH1F>("mom_phi_",        "#phi mother",                50,  -3.5,  3.5) ;
  
  gen_e_pT_        = fs->make<TH1F>("gen_e_pT_",       "Gen e p_{T}",               100, -50.,  50.)  ;
  gen_eP_eta_      = fs->make<TH1F>("gen_eP_eta_",     "Gen e^{+} #eta spectrum",    50,  -3.5,  3.5) ;
  gen_eP_phi_      = fs->make<TH1F>("gen_eP_phi_",     "Gen e^{+} #phi spectrum",    50,  -0.5,  3.5) ;
  gen_eM_eta_      = fs->make<TH1F>("gen_eM_eta_",     "Gen e^{-} #eta spectrum",    50,  -3.5,  3.5) ;
  gen_eM_phi_      = fs->make<TH1F>("gen_eM_phi_",     "Gen e^{-} #phi spectrum",    50,  -0.5,  3.5) ;
 
  gen_mu_pT_       = fs->make<TH1F>("gen_mu_pT_",      "Gen mu p_{T}",              100, -50.,  50.)  ;
  gen_muP_eta_     = fs->make<TH1F>("gen_muP_eta_",    "Gen mu^{+} #eta spectrum",   50,  -3.5,  3.5) ;
  gen_muP_phi_     = fs->make<TH1F>("gen_muP_phi_",    "Gen mu^{+} #phi spectrum",   50,  -0.5,  3.5) ;
  gen_muM_eta_     = fs->make<TH1F>("gen_muM_eta_",    "Gen mu^{-} #eta spectrum",   50,  -3.5,  3.5) ;
  gen_muM_phi_     = fs->make<TH1F>("gen_muM_phi_",    "Gen mu^{-} #phi spectrum",   50,  -0.5,  3.5) ;
 
  gen_gamma_pT_    = fs->make<TH1F>("gen_gamma_pT_",   "Gen #gamma p_{T}",           50,   0.,  50.)  ;
  gen_gamma_eta_   = fs->make<TH1F>("gen_gamma_eta_",  "Gen #gamma #eta spectrum",   50,  -3.5,  3.5) ;
  gen_gamma_phi_   = fs->make<TH1F>("gen_gamma_phi_",  "Gen #gamma #phi spectrum",   50,  -0.5,  3.5) ;
 
  gen_W_pT_        = fs->make<TH1F>("gen_W_pT_",       "Gen W p_{T}",               100, -50.,  50.)  ;
  gen_Wp_eta_      = fs->make<TH1F>("gen_Wp_eta_",     "Gen W^{+} #eta spectrum",    50,  -3.5,  3.5) ;
  gen_Wp_phi_      = fs->make<TH1F>("gen_Wp_phi_",     "Gen W^{+} #phi spectrum",    50,  -0.5,  3.5) ;
  gen_Wm_eta_      = fs->make<TH1F>("gen_Wm_eta_",     "Gen W^{-} #eta spectrum",    50,  -3.5,  3.5) ;
  gen_Wm_phi_      = fs->make<TH1F>("gen_Wm_phi_",     "Gen W^{-} #phi spectrum",    50,  -0.5,  3.5) ;
 
  gen_Z_pT_        = fs->make<TH1F>("gen_Z_pT_",       "Gen Z p_{T}",                50,   0.,  50.)  ;
  gen_Z_eta_       = fs->make<TH1F>("gen_Z_eta_",      "Gen Z #eta spectrum",        50,  -3.5,  3.5) ;
  gen_Z_phi_       = fs->make<TH1F>("gen_Z_phi_",      "Gen Z #phi spectrum",        50,  -0.5,  3.5) ;
  
  gen_H0_pT_       = fs->make<TH1F>("gen_H0_pT_",      "Gen H^{0} p_{T}",            50,   0.,  50.)  ;
  gen_H0_eta_      = fs->make<TH1F>("gen_H0_eta_",     "Gen H^{0} #eta spectrum",    50,  -3.5,  3.5) ;
  gen_H0_phi_      = fs->make<TH1F>("gen_H0_phi_",     "Gen H^{0} #phi spectrum",    50,  -0.5,  3.5) ;

  gen_Zprime_pT_   = fs->make<TH1F>("gen_Zprime_pT_",  "Gen Z^{'} p_{T}",            50,   0.,  50.)  ;
  gen_Zprime_eta_  = fs->make<TH1F>("gen_Zprime_eta_", "Gen Z^{'} #eta spectrum",    50,  -3.5,  3.5) ;
  gen_Zprime_phi_  = fs->make<TH1F>("gen_Zprime_phi_", "Gen Z^{'} #phi spectrum",    50,  -0.5,  3.5) ;
  
  dau_e_pT_        = fs->make<TH1F>("dau_e_pT_",       "Dau e p_{T}",               100, -50.,  50.)  ;
  dau_eP_eta_      = fs->make<TH1F>("dau_eP_eta_",     "Dau e^{+} #eta",             50,  -3.5,  3.5) ;
  dau_eP_phi_      = fs->make<TH1F>("dau_eP_phi_",     "Dau e^{+} #phi",             50,  -0.5,  3.5) ;
  dau_eM_eta_      = fs->make<TH1F>("dau_eM_eta_",     "Dau e^{-} #eta",             50,  -3.5,  3.5) ;
  dau_eM_phi_      = fs->make<TH1F>("dau_eM_phi_",     "Dau e^{-} #phi",             50,  -0.5,  3.5) ;
 
  dau_mu_pT_       = fs->make<TH1F>("dau_mu_pT_",      "Dau mu p_{T}",              100, -50.,  50.)  ;
  dau_muP_eta_     = fs->make<TH1F>("dau_muP_eta_",    "Dau mu^{+} #eta",            50,  -3.5,  3.5) ;
  dau_muP_phi_     = fs->make<TH1F>("dau_muP_phi_",    "Dau mu^{+} #phi",            50,  -0.5,  3.5) ;
  dau_muM_eta_     = fs->make<TH1F>("dau_muM_eta_",    "Dau mu^{-} #eta",            50,  -3.5,  3.5) ;
  dau_muM_phi_     = fs->make<TH1F>("dau_muM_phi_",    "Dau mu^{-} #phi",            50,  -0.5,  3.5) ;
 
  dau_gamma_pT_   = fs->make<TH1F>("dau_gamma_pT_",    "Dau #gamma p_{T}",           50,   0.,  50.)  ;
  dau_gamma_eta_  = fs->make<TH1F>("dau_gamma_eta_",   "Dau #gamma #eta",            50,  -2.5,  2.5) ;
  dau_gamma_phi_  = fs->make<TH1F>("dau_gamma_phi_",   "Dau #gamma #phi",            50,  -0.5,  3.5) ;

   nEvt = iEvt = iRun = 0 ;

}

void 
MCValidation::endJob() {

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MCValidation);

