/** \class cand::modules::CandReducer
 *
 * Given a collectin of candidates, produced a
 * collection of LeafCandidas identical to the
 * source collection, but removing all daughters
 * and all components. 
 *
 * This is ment to produce a "light" collection
 * of candiadates just containing kimenatics 
 * information for very fast analysis purpose
 *
 * \author Luca Lista, INFN
 *
 * \version $Revision: 1.1 $
 *
 * $Id: CandReducer.cc,v 1.1 2007/01/30 15:16:56 llista Exp $
 *
 */
#include "FWCore/Framework/interface/EDProducer.h"

class CandReducer : public edm::EDProducer {
public:
  /// constructor from parameter set
  explicit CandReducer( const edm::ParameterSet& );
  /// destructor
  ~CandReducer();
private:
  /// process one evevnt
  void produce( edm::Event& evt, const edm::EventSetup& );
  /// label of source candidate collection
  std::string src_;
};

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "FWCore/Utilities/interface/EDMException.h"

using namespace reco;
using namespace edm;

CandReducer::CandReducer( const edm::ParameterSet& cfg ) :
  src_( cfg.getParameter<std::string>("src") ) {
  produces<CandidateCollection>();
}

CandReducer::~CandReducer() {
}

void CandReducer::produce( Event& evt, const EventSetup& ) {
  Handle<CandidateCollection> cands;
  evt.getByLabel( src_, cands );
  std::auto_ptr<CandidateCollection> comp( new CandidateCollection );
  for( CandidateCollection::const_iterator c = cands->begin(); c != cands->end(); ++c ) {
    std::auto_ptr<Candidate> cand( new LeafCandidate( * c ) );
    comp->push_back( cand.release() );
  }
  evt.put( comp );
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( CandReducer );
