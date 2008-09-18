#include "PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


using namespace reco;

namespace GenParticlesHelper {

  void 
  findParticles(const reco::GenParticleCollection& sourceParticles, reco::GenParticleRefVector& particleRefs, int pdgId, int status ) {

    unsigned index = 0;
    for(IG ig = sourceParticles.begin(); 
	ig!= sourceParticles.end(); ++ig, ++index) {
 
      const GenParticle& gen = *ig;
    
      // status has been specified, and this one does not have the correct
      // status
      if(status && gen.status()!=status ) continue;
    
      if( abs(gen.pdgId()) == pdgId ) {
	GenParticleRef genref( &sourceParticles, index );
	particleRefs.push_back( genref );
      }
    }
  }


  void 
  findDescendents(const reco::GenParticleRef& base, 
		  reco::GenParticleRefVector& descendents, 
		  int status, int pdgId ) {


    const GenParticleRefVector& daughterRefs = base->daughterRefVector();
  
    for(IGR idr = daughterRefs.begin(); 
	idr!= daughterRefs.end(); ++idr ) {
    
      if( (*idr)->status() == status && 
	  (!pdgId || abs((*idr)->pdgId()) == pdgId) ) {
      
	// cout<<"adding "<<(*idr)<<endl;
	descendents.push_back(*idr);
      }
      else 
	findDescendents( *idr, descendents, status, pdgId );
    }
  }



  void 
  findSisters(const reco::GenParticleRef& baseSister, 
	      GenParticleRefVector& sisterRefs) {
  
    assert( baseSister->numberOfMothers() > 0 );
  
    // get first mother 
    const GenParticleRefVector& mothers = baseSister->motherRefVector();
  
    // get sisters 
    const GenParticleRefVector allRefs 
      = mothers[0]->daughterRefVector();
  
    typedef GenParticleRefVector::const_iterator IT;
    for(IT id = allRefs.begin();
	id != allRefs.end();
	++id ) {
    
      if( *id == baseSister ) { 
	continue; // this is myself
      }
      else 
	sisterRefs.push_back( *id );
    }
  }
}
