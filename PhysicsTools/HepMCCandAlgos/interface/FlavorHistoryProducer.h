#ifndef PhysicsTools_HepMCCandAlgos_interface_FlavorHistoryProducer_h
#define PhysicsTools_HepMCCandAlgos_interface_FlavorHistoryProducer_h

/** class 
 *
 * \author Stephen Mrenna, FNAL
 *
 * \version $Id: FlavorHistoryProducer.cc,v 1.0
 *
 */


// -------------------------------------------------------------
// Identify the ancestry of the Quark
// 
// 
// Matrix Element:
//    Status 3 parent with precisely 2 "grandparents" that
//    is outside of the "initial" section (0-5) that has the
//    same ID as the status 2 parton in question. 
//    NOTE: This is not the actual ultimate progenitor,
//    but this is the signature of matrix element decays.
//    The ultimate progenitor is the parent of the status 3
//    parton.
//
// Flavor excitation:
//    Almost the same as the matrix element classification,
//    but has only one outgoing parton product instead of two.
//
// Gluon splitting:
//    Parent is a quark of a different flavor than the parton
//    in question, or a gluon. Can come from either ISR or FSR.
//
// True decay:
//    Decays from a resonance like top, Higgs, etc.
// -------------------------------------------------------------

#include "FWCore/Framework/interface/EDProducer.h"
#include <string>
#include <vector>
#include <set>
#include <utility>
#include <algorithm>

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/ShallowClonePtrCandidate.h"
#include "DataFormats/HepMCCandidate/interface/FlavorHistory.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <fstream>


std::ostream & operator<<( std::ostream & out, reco::Candidate const & cand) ;

std::ostream & operator<<( std::ostream & out, reco::FlavorHistory const & cand);

class FlavorHistoryProducer : public edm::EDProducer {
 public:
  /// constructor
  FlavorHistoryProducer( const edm::ParameterSet & );
  /// destructor
  ~FlavorHistoryProducer();

 private:
  /// module init at begin of job
  void beginJob( const edm::EventSetup & );
  /// process one event
  void produce( edm::Event& e, const edm::EventSetup& );

  void getAncestors(const reco::Candidate &c,
		    std::vector<reco::Candidate const * > & moms );

  reco::Candidate const * getSister(const reco::Candidate &c);

  
  
  edm::InputTag src_;               // source collection name 
  int    pdgIdToSelect_;            // pdg of hf partons to select
  double ptMinParticle_;            // minimum pt of the partons
  double ptMinShower_;              // minimum pt of the shower
  double etaMaxParticle_;           // max eta of the parton
  double etaMaxShower_;             // max eta of the shower
  std::string flavorHistoryName_;   // name to give flavor history
  bool verbose_;                    // verbose flag
};

#endif
