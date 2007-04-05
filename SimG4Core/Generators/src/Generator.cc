#include "SimG4Core/Generators/interface/Generator.h"
#include "SimG4Core/Generators/interface/HepMCParticle.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HEPEvtParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"

using namespace edm;
using std::cout;
using std::endl;
using std::string;

Generator::Generator(const ParameterSet & p) : 
  fPtCuts(p.getParameter<bool>("ApplyPtCuts")),
  fEtaCuts(p.getParameter<bool>("ApplyEtaCuts")), 
  fPhiCuts(p.getParameter<bool>("ApplyPhiCuts")),
  //theMinPhiCut(p.getParameter<double>("MinPhiCut")*deg),
  //theMaxPhiCut(p.getParameter<double>("MaxPhiCut")*deg),
  theMinPhiCut(p.getParameter<double>("MinPhiCut")),   // now operates in radians (CMS standard)
  theMaxPhiCut(p.getParameter<double>("MaxPhiCut")),
  theMinEtaCut(p.getParameter<double>("MinEtaCut")),
  theMaxEtaCut(p.getParameter<double>("MaxEtaCut")),
  //theMinPtCut(p.getParameter<double>("MinPtCut")*MeV),
  //theMaxPtCut(p.getParameter<double>("MaxPtCut")*MeV),   
  theMinPtCut(p.getParameter<double>("MinPtCut")),    // now operates in GeV (CMS standard)
  theMaxPtCut(p.getParameter<double>("MaxPtCut")),   
  theDecLenCut(p.getUntrackedParameter<double>("DecLenCut",2.9)*cm),
  verbose(p.getUntrackedParameter<int>("Verbosity",0)),
  evt_(0),
  vtx_(0) ,
  weight_(0)
{
  edm::LogInfo("SimG4CoreGenerator") << " Generator constructed " ;
}

Generator::~Generator() 
{ 
}

void Generator::HepMC2G4(const HepMC::GenEvent * evt, G4Event * g4evt)
{

  //M. Vander Donckt : modified to take the generator event weight  
  if ( evt->weights().size() > 0 )
    {
      weight_ = evt->weights()[0] ;
      for ( int iw=1; iw<evt->weights().size(); iw++ )
        {
          // terminate if the versot of weights contains a zero-weight
          if ( evt->weights()[iw] <= 0 ) break;
          weight_ *= evt->weights()[iw] ;
        }     
    }
  // end modification
  
  // in the future, we probably want to skip events of zero-weight
  // but at this point, it's zero in most cases anyway... 
  // just a note for future... (JY)  

  if (vtx_ != 0) delete vtx_;
  vtx_ = new HepLorentzVector((*(evt->vertices_begin()))->position());
  
  if (verbose > 0)
    {
      evt->print();
      cout << " " << endl;
      cout << " Prim.Vtx : " << vtx_->x() << " " 
           << vtx_->y() << " "
           << vtx_->z() << endl;
    }
    
      for(HepMC::GenEvent::vertex_const_iterator vitr= evt->vertices_begin();
      vitr != evt->vertices_end(); ++vitr ) 
    { // loop for vertex ...
        
      // real vertex?
      G4bool qvtx=false;
      for (HepMC::GenVertex::particle_iterator 
             pitr= (*vitr)->particles_begin(HepMC::children);
           pitr != (*vitr)->particles_end(HepMC::children); ++pitr) 
        {
          if (!(*pitr)->end_vertex() && (*pitr)->status()==1) 
            {
              qvtx=true;
              break;
              //bug fix as we need to keep long decaying particles even if the are intermediate
            } else 
            if ( (*pitr)->status()== 2 ) {
              HepLorentzVector xvtx = (*vitr)->position() ;
              HepLorentzVector dvtx=(*pitr)->end_vertex()->position();
              double dd=(xvtx-dvtx).rho();
              if (dd>theDecLenCut){
                qvtx=true;
                break;
              }
            } 
        }
      if (!qvtx) 
        {
          continue;
        }
        
      // check world boundary
      //G4LorentzVector xvtx= (*vitr)-> position();
      HepLorentzVector xvtx = (*vitr)->position() ;
      //fix later
      //if (! CheckVertexInsideWorld(xvtx.vect()*mm)) continue;
        
      // create G4PrimaryVertex and associated G4PrimaryParticles
      G4PrimaryVertex* g4vtx= 
        new G4PrimaryVertex(xvtx.x()*mm, xvtx.y()*mm, xvtx.z()*mm, 
                            xvtx.t()*mm/c_light);
        
      for (HepMC::GenVertex::particle_iterator 
             vpitr= (*vitr)->particles_begin(HepMC::children);
           vpitr != (*vitr)->particles_end(HepMC::children); ++vpitr) 
        {
        
          // M. Vander Donckt modification: to take also decay mother
          // in case decay length is large; decay procuts get setup
          // as daughters of G4Particle in this case, through the method
          // particleAssignDaughters, and they get marked 1000+status 
          // in the generator product (this seem to "violate" the idea 
          // that a product can't be modified one it's in edm::Event... 
          // but it seems to fly...)
          double decay_length=-1;
          if ( (*vpitr)->status() == 2 ) 
            {
              if ( (*vpitr)->end_vertex() != 0 ) // over-protection
                {
                  HepLorentzVector dvtx=(*vpitr)->end_vertex()->position();
                  decay_length=(dvtx-xvtx).rho();
                }
              // end modification
            }           
          if( (*vpitr)->status() == 1 || ((*vpitr)->status() == 2 && decay_length > theDecLenCut ) ) {
            
            //G4LorentzVector p= (*vpitr)->momentum();
            HepLorentzVector p = (*vpitr)->momentum() ;
            
                if ( !particlePassesPrimaryCuts( p ) ) 
                  {
                    continue ;
                  }
            
                G4int pdgcode= (*vpitr)-> pdg_id();
                G4PrimaryParticle* g4prim= 
                  new G4PrimaryParticle(pdgcode, p.x()*GeV, p.y()*GeV, p.z()*GeV);
            
                if ( g4prim->GetG4code() != 0 )
                  { 
                    g4prim->SetMass( g4prim->GetG4code()->GetPDGMass() ) ;
                    g4prim->SetCharge( g4prim->GetG4code()->GetPDGCharge() ) ;  
                  }
            
                g4prim->SetWeight( 10000*(*vpitr)->barcode() ) ;
                setGenId( g4prim, (*vpitr)->barcode() ) ;
                if ( (*vpitr)->status() == 2) particleAssignDaughters(g4prim,(HepMC::GenParticle *) *vpitr, decay_length);
                g4vtx->SetPrimary(g4prim);
              }
            }
          g4evt->AddPrimaryVertex(g4vtx);
        }

      return ;
   
    }
  void Generator::particleAssignDaughters( G4PrimaryParticle* g4p, HepMC::GenParticle* vp, double decaylength)
  {
 
    if ( !(vp->end_vertex())  ) return ;
   
    HepLorentzVector p = vp->momentum() ;
    Hep3Vector cmboost=p.findBoostToCM();
    double proper_time=decaylength/(p.beta()*p.gamma()*c_light);
    g4p->SetProperTime(proper_time*ns); // the particle will decay after the same length if it has not interacted before
    HepLorentzVector xvtx=vp->end_vertex()->position();
    for (HepMC::GenVertex::particle_iterator 
           vpdec= vp->end_vertex()->particles_begin(HepMC::children);
         vpdec != vp->end_vertex()->particles_end(HepMC::children); ++vpdec) {

      //transform decay products such that in the rest frame of mother
      HepLorentzVector pdec = ((*vpdec)->momentum()).boost(cmboost) ;
      G4PrimaryParticle * g4daught= 
        new G4PrimaryParticle((*vpdec)->pdg_id(), pdec.x()*GeV, pdec.y()*GeV, pdec.z()*GeV);
      if ( g4daught->GetG4code() != 0 )
        { 
          g4daught->SetMass( g4daught->GetG4code()->GetPDGMass() ) ;
          g4daught->SetCharge( g4daught->GetG4code()->GetPDGCharge() ) ;  
        }
      g4daught->SetWeight( 10000*(*vpdec)->barcode() ) ;
      setGenId( g4daught, (*vpdec)->barcode() ) ;
      if ( (*vpdec)->status() == 2 && (*vpdec)->end_vertex() != 0 ) 
        {
          HepLorentzVector dvtx=(*vpdec)->end_vertex()->position();
          double dd=(dvtx.x()-xvtx.x())*(dvtx.x()-xvtx.x())
            +(dvtx.y()-xvtx.y())*(dvtx.y()-xvtx.y())
            +(dvtx.z()-xvtx.z())*(dvtx.z()-xvtx.z());
          dd=sqrt(dd);
          particleAssignDaughters(g4daught,*vpdec,dd);
        }
      // children should only be taken into account once
      (*vpdec)->set_status(1000+(*vpdec)->status()); //bug fix position 
      g4p->SetDaughter(g4daught);
    }
    return;
  }

  bool Generator::particlePassesPrimaryCuts( const HepLorentzVector& mom ) const 
  {

    double phi = mom.phi() ;   
    double pt  = sqrt( mom.x()*mom.x() + mom.y()*mom.y() ) ;
    double eta = -log( tan(mom.theta()/2.) ) ;
      
    if ( (fPtCuts)  && (pt  < theMinPtCut  || pt  > theMaxPtCut) )  return false ;
    if ( (fEtaCuts) && (eta < theMinEtaCut || eta > theMaxEtaCut) ) return false ;
    if ( (fPhiCuts) && (phi < theMinPhiCut || phi > theMaxPhiCut) ) return false ;
   
    return true;   
  }
 
  bool Generator::particlePassesPrimaryCuts(const G4PrimaryParticle * p) const
  {
    G4ThreeVector mom = p->GetMomentum();
    double        phi = mom.phi() ;
    double        pt  = sqrt(p->GetPx()*p->GetPx() + p->GetPy()*p->GetPy());
    pt /= GeV ;  // need to convert, since Geant4 operates in MeV
    double        eta = -log(tan(mom.theta()/2));
    if (((fPtCuts)  && (pt  < theMinPtCut  || pt  > theMaxPtCut))           ||
        ((fEtaCuts) && (eta < theMinEtaCut || eta > theMaxEtaCut))          ||
        ((fPhiCuts) && (phi < theMinPhiCut || phi > theMaxPhiCut)))
      return false;
    else return true;
  }

  void Generator::nonBeamEvent2G4(const HepMC::GenEvent * evt, G4Event * g4evt)
  {
    int i = 0; 
    for(HepMC::GenEvent::particle_const_iterator it = evt->particles_begin(); 
        it != evt->particles_end(); ++it )
      {
        i++;
        HepMC::GenParticle * g = (*it); 
        int g_status = g->status();
        // storing only particle with status == 1       
        if (g_status == 1)
          {
            HepLorentzVector mom  = g->momentum();
            int g_id = g->pdg_id();         
            G4PrimaryParticle * g4p = 
              new G4PrimaryParticle(g_id,mom.x()*GeV,mom.y()*GeV,mom.z()*GeV);
            if (g4p->GetG4code() != 0)
              { 
                g4p->SetMass(g4p->GetG4code()->GetPDGMass());
                g4p->SetCharge(g4p->GetG4code()->GetPDGCharge()) ;
              }
            g4p->SetWeight(i*10000);
            setGenId(g4p,i);
            if (particlePassesPrimaryCuts(g4p))
              {
                HepLorentzVector vtx = g->production_vertex()->position();
                G4PrimaryVertex * v = new 
                  G4PrimaryVertex(vtx.x()*mm,vtx.y()*mm,vtx.z()*mm,vtx.t()*mm/c_light);
                v->SetPrimary(g4p);
                g4evt->AddPrimaryVertex(v);
              }
          }
      } // end loop on HepMC particles
  }

