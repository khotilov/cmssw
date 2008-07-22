// This file was removed but it should not have been.
// This comment is to restore it. 

#include "PhysicsTools/HepMCCandAlgos/interface/FlavorHistoryProducer.h"

// #include "DataFormats/Common/interface/ValueMap.h"
// #include <iterators>

using namespace std;
using namespace reco;
using namespace edm;


ostream & operator<<( ostream & out, Candidate const & cand) 
{
  char buff[1000];
  sprintf(buff, "%5d, status = %5d, nmo = %5d, nda = %5d, pt = %6.2f, eta = %6.2f, phi = %6.2f, m = %6.2f", 
	  cand.pdgId(), cand.status(), 
	  cand.numberOfMothers(),
	  cand.numberOfDaughters(),
	  cand.pt(), cand.eta(), cand.phi(), cand.mass() );
  out << buff;
  return out;
}

ostream & operator<<( ostream & out, FlavorHistory const & cand) 
{
  cout << "Source     = " << cand.flavorSource() << endl;
  if ( cand.hasParton() ) 
    cout << "Parton     = " << cand.parton().key() << " : " << *(cand.parton()) << endl;
  if ( cand.hasProgenitor() ) 
    cout << "Progenitor = " << cand.progenitor().key() << " : " << *(cand.progenitor()) << endl;
  if ( cand.hasSister() ) 
    cout << "Sister     = " << cand.sister().key() << " : " << *(cand.sister()) << endl;
  if ( cand.hasParton() ) {
    cout << "Ancestry: " << endl;
    Candidate const * ipar = cand.parton()->mother();
    while ( ipar->numberOfMothers() > 0 ) {
      cout << *ipar << endl;
      ipar = ipar->mother();
    }
  }
  return out;
}

FlavorHistoryProducer::FlavorHistoryProducer( const ParameterSet & p ) :
  src_( p.getParameter<InputTag>( "src" ) ),
  pdgIdToSelect_( p.getParameter<int> ("pdgIdToSelect") ),
  ptMinParticle_( p.getParameter<double>( "ptMinParticle") ),  
  ptMinShower_( p.getParameter<double>( "ptMinShower") ),  
  etaMaxParticle_( p.getParameter<double>( "etaMaxParticle" )),  
  etaMaxShower_( p.getParameter<double>( "etaMaxShower" )),
  flavorHistoryName_( p.getParameter<string>("flavorHistoryName") ),
  verbose_( p.getUntrackedParameter<bool>( "verbose" ) )
{
  produces<vector<FlavorHistory> >(flavorHistoryName_);
}

FlavorHistoryProducer::~FlavorHistoryProducer() { 
}

void FlavorHistoryProducer::beginJob( const EventSetup & es ) {
  ;
}

void FlavorHistoryProducer::produce( Event& evt, const EventSetup& ) 
{
  
  // Get a handle to the particle collection (OwnVector)
  Handle<CandidateView > particlesViewH;
  evt.getByLabel( src_, particlesViewH );

//   const vector<Candidate const *> & particles = *particlesViewH;

  // Copy the View to an vector for easier iterator manipulation convenience
  vector<const Candidate* > particles;
  for( CandidateView::const_iterator p = particlesViewH->begin();  p != particlesViewH->end(); ++p ) {
    particles.push_back(&*p);
  }
  
  // Make a new flavor history vector
  auto_ptr<vector<FlavorHistory> > flavorHistoryVector ( new vector<FlavorHistory> () ) ;

  // ------------------------------------------------------------
  // Loop over partons
  // ------------------------------------------------------------
  vector<const Candidate* >::size_type j;
  vector<const Candidate* >::size_type j_max = particles.size();
  for( j=0; j<j_max; ++j ) {

    // Get the candidate
    const Candidate *p = particles[j];
    // Set up indices that we'll need for the flavor history
    vector<Candidate const *>::size_type partonIndex=j;
    vector<Candidate const *>::size_type progenitorIndex=0;
    vector<Candidate const *>::size_type sisterIndex=0;
    bool foundProgenitor = false; 
    bool foundSister = false;
    FlavorHistory::FLAVOR_T flavorSource=FlavorHistory::FLAVOR_NULL;


    int idabs = abs( (p)->pdgId() );
    int nDa = (p)->numberOfDaughters();

    // Check if we have a status 2 or 3 particle, which is a parton before the string
    if ( p->status() == 2 ) {
//       if(verbose_) cout << "--------------------------" << endl;
//       if(verbose_) cout << "Processing particle " << j  << " = " << *p << endl;
      // Ensure the parton in question has daughters 
      //       if ( nDa > 0 && ( (p)->daughter(0)->pdgId() == 91 || (p)->daughter(0)->pdgId() == 92 ||
      // 			(p)->daughter(0)->pdgId() == 93) ) {
      if ( nDa > 0 ) {
// 	if(verbose_) cout << "Has daughters" << endl;
	// Ensure the parton passes some minimum kinematic cuts
	if((p)->pt() > ptMinShower_ && fabs((p)->eta())<etaMaxShower_) {
// 	  if(verbose_) cout << "Passes kin cuts" << endl;

 	  if(verbose_) cout << "--------------------------" << endl;
 	  if(verbose_) cout << "Processing particle " << j  << " = " << *p << endl;


	  //Main (but simple) workhorse to get all ancestors
	  vector<Candidate const *> allParents;
	  getAncestors( *p, allParents );
	    
	  if(verbose_) {
	    cout << "Parents : " << endl;
	    vector<Candidate const *>::const_iterator iprint = allParents.begin(),
	      iprintend = allParents.end();
	    for ( ; iprint != iprintend; ++iprint ) 
	      cout << **iprint << endl;
	  }
	  
	  // ------------------------------------------------------------
	  // Now identify the flavor and ancestry of the HF Quark
	  // Mother               Origin
	  // ======               =======
	  // incoming quarks      ISR, likely gluon splitting
	  //   light flavor
	  // incoming partons     ISR, likely flavor excitation
	  //   heavy flavor           
	  // outgoing quark       FSR
	  //   light flavor
	  // outgoing quark       Matrix Element b       
	  //   heavy flavor
	  //     no mother
	  // outgoing quark       Resonance b (e.g. top quark decay)
	  //   heavy flavor
	  //     mother
	  // outgoing resonance   Resonance b (e.g. Higgs decay)
	  // ------------------------------------------------------------
	  vector<Candidate const *>::size_type a_size = allParents.size();
	  int parentIndex=0;

	  // 
	  // Loop over all the ancestors of this parton and find the progenitor.
	  // 
	  for( vector<Candidate const *>::size_type i=0 ; i < a_size && !foundProgenitor; ++i,++parentIndex ) {
	    const Candidate * aParent=allParents[i];
	    vector<Candidate const *>::const_iterator found = find(particles.begin(),particles.end(),aParent);


	    // Get the index of the progenitor candidate
	    progenitorIndex = found - particles.begin();

	    int aParentId = abs(aParent->pdgId());

	    // Here we examine particles that were produced after the collision
 	    if( aParent->numberOfMothers() == 2 && progenitorIndex > 5 ) {
	      // Here is where we have a matrix element
	      if( aParentId == pdgIdToSelect_ ) {
		if(verbose_) cout << "Matrix Element progenitor" << endl;
		flavorSource = FlavorHistory::FLAVOR_ME;
	      } 
	      // Here we have a gluon splitting from final state radiation
	      else if( (aParentId > 0 && aParentId < FlavorHistory::tQuarkId ) || aParentId==FlavorHistory::gluonId ) {
		if(verbose_) cout << "Gluon splitting progenitor" << endl;
		flavorSource = FlavorHistory::FLAVOR_GS;
	      }
	      // Here we have a true decay
	      else if( (aParentId>pdgIdToSelect_ && aParentId<FlavorHistory::gluonId) || aParentId > FlavorHistory::gluonId ) {
		if(verbose_) cout << "Flavor decay progenitor" << endl;
		flavorSource = FlavorHistory::FLAVOR_DECAY;
	      }
	      foundProgenitor = true;
	    }

	    // Here we examine particles that were produced before the collision
	    else if( progenitorIndex==2 || progenitorIndex==3 ) {
 	      // Here is a flavor excitation
	      flavorSource = FlavorHistory::FLAVOR_EXC;

	      // Parent has a quark daughter equal and opposite to this
	      if( aParent->numberOfDaughters() > 0 && 
		  aParent->daughter(0)->pdgId() == -1 * p->pdgId()  ) {
		if(verbose_) cout << "Flavor excitation progenitor" << endl;
		flavorSource = FlavorHistory::FLAVOR_EXC;
	      }
	      // Here is gluon splitting from initial state radiation
	      else {		  
		if(verbose_) cout << "Gluon splitting progenitor" << endl;
		flavorSource = FlavorHistory::FLAVOR_GS;
	      }
	      foundProgenitor = true;
	    }
	  }// End loop over all parents of this parton to find progenitor

	    

	  // 
	  // Now find sister of this particle if there is one
	  // 
	  if ( foundProgenitor ) {
	    // Get the progenitor
	    const Candidate * progenitorCand = particles[progenitorIndex];

	    // Here is the normal case of a sister around
	    if ( progenitorCand->numberOfDaughters() >= 2 ) {
	      const Candidate * sisterCand = 0;
	      
	      for ( unsigned int iida = 0; iida < progenitorCand->numberOfDaughters(); ++iida ) {
		const Candidate * dai = progenitorCand->daughter(iida);

		if ( verbose_ ) cout << "Sister " << *dai << endl;
		
		if ( dai->pdgId() == -1 * p->pdgId() ) {
		  if ( verbose_ ) cout << "Found sister" << endl;
		  sisterCand = dai;
		  foundSister = true;
		}
	      }
		
	      if ( foundSister ) {
		// Find index of daughter in master list
		vector<Candidate const *>::const_iterator found = find(particles.begin(),particles.end(),sisterCand);
		sisterIndex = found - particles.begin();
		if(verbose_) cout << "Sister index = " << sisterIndex << endl;
		if ( found != particles.end() )
		  if(verbose_) cout << "Sister = " << **found << endl;
	      } // end if found sister
	    }
	    // Here is if we have a "transient" decay in the code that isn't
	    // really a decay, so we need to look at the parent of the progenitor
	    else {
	      const Candidate * grandProgenitorCand = progenitorCand->mother(0);
	      const Candidate * sisterCand = 0;

	      if ( verbose_ ) cout << "Looking for sister, progenitor is " << *progenitorCand << endl;
	    
	      // Make sure the progenitor has two daughters
	      if ( grandProgenitorCand->numberOfDaughters() >= 2 ) {

		for ( unsigned int iida = 0; iida < grandProgenitorCand->numberOfDaughters(); ++iida ) {
		  const Candidate * dai = grandProgenitorCand->daughter(iida);

		  if ( verbose_ ) cout << "Looking for sister " << *dai << endl;
		
		  if ( dai->pdgId() == -1 * p->pdgId() ) {
		    if ( verbose_ ) cout << "Found sister" << endl;
		    sisterCand = dai;
		    foundSister = true;
		  }
		}
		
		if ( foundSister ) {
		  // Find index of daughter in master list
		  vector<Candidate const *>::const_iterator found = find(particles.begin(),particles.end(),sisterCand);
		  sisterIndex = found - particles.begin();
		  if(verbose_) cout << "Sister index = " << sisterIndex << endl;
		  if ( found != particles.end() )
		    if(verbose_) cout << "Sister = " << **found << endl;
		} // end if found sister
	      } // End of have at least 2 grand progenitor daughters
	    } // End if we have to look at parents of progenitor to find sister

	  } // end if found progenitor
	  
	}// End if this parton passes some minimal kinematic cuts
      }// End if this particle has strings as daughters
    }// End if this particle was a status==2 parton

    // Make sure we've actually found a sister and a progenitor
    if ( !foundProgenitor ) progenitorIndex = 0;
    if ( !foundSister ) sisterIndex = 0;

    // We've found the particle, add to the list (status 2 only)
    if ( idabs == pdgIdToSelect_ && p->status() == 2 ) 
      flavorHistoryVector->push_back( FlavorHistory( flavorSource, particlesViewH, partonIndex, progenitorIndex, sisterIndex ) ); 
  }


//   ValueMap<FlavorHistory>::Filler filler(*flavorHistory);
//   filler.insert( particlesViewH, flavorHistoryVector.begin(), flavorHistoryVector.end()  );
//   filler.fill();
  // Now add the flavor history to the event record
  if ( true ) {
    cout << "Outputting pdg id = " << pdgIdToSelect_ << " with nelements = " << flavorHistoryVector->size() << endl;
    vector<FlavorHistory>::const_iterator i = flavorHistoryVector->begin(),
      iend = flavorHistoryVector->end();
    for ( ; i !=iend; ++i ) {
      cout << *i << endl;
    }
  }
  evt.put( flavorHistoryVector, flavorHistoryName_ );
}

 
// Helper function to get all ancestors of this candidate
void FlavorHistoryProducer::getAncestors(const Candidate &c,
					 vector<Candidate const *> & moms )
{

  if( c.numberOfMothers() == 1 ) {
    const Candidate * dau = &c;
    const Candidate * mom = c.mother();
    while ( dau->numberOfMothers() != 0) {
      moms.push_back( dau );
      dau = mom ;
      mom = dau->mother();
    } 
  } 
}



#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( FlavorHistoryProducer );
