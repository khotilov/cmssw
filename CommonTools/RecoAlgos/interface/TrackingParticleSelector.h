#ifndef RecoSelectors_TrackingParticleSelector_h
#define RecoSelectors_TrackingParticleSelector_h
/* \class TrackingParticleSelector
 *
 * \author Giuseppe Cerati, INFN
 *
 *  $Date: 2008/08/09 14:31:09 $
 *  $Revision: 1.4 $
 *
 */
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

class TrackingParticleSelector {

public:
  TrackingParticleSelector(){}
  TrackingParticleSelector ( double ptMin,double minRapidity,double maxRapidity,
			     double tip,double lip,int minHit, bool signalOnly, bool chargedOnly, 
			     std::vector<int> pdgId = std::vector<int>()) :
    ptMin_( ptMin ), minRapidity_( minRapidity ), maxRapidity_( maxRapidity ),
    tip_( tip ), lip_( lip ), minHit_( minHit ), signalOnly_(signalOnly), chargedOnly_(chargedOnly), pdgId_( pdgId ) { }
  
  /// Operator() performs the selection: e.g. if (tPSelector(tp)) {...}
  bool operator()( const TrackingParticle & tp ) const { 

    //quickly reject if it is from pile-up
    if (signalOnly_ && !(tp.eventId().bunchCrossing()==0 && tp.eventId().event()==0) )return false;

    if (chargedOnly_ && tp.charge()==0) return false;//select only if charge!=0
    bool testId = false;
    unsigned int idSize = pdgId_.size();
    if (idSize==0) testId = true;
    else for (unsigned int it=0;it!=idSize;++it){
      if (tp.pdgId()==pdgId_[it]) testId = true;
    }
    return (
	    tp.matchedHit() >= minHit_ &&
	    sqrt(tp.momentum().perp2()) >= ptMin_ && 
	    tp.momentum().eta() >= minRapidity_ && tp.momentum().eta() <= maxRapidity_ && 
	    sqrt(tp.vertex().perp2()) <= tip_ &&
	    fabs(tp.vertex().z()) <= lip_ &&
	    testId
	    );
  }
  
private:
  double ptMin_;
  double minRapidity_;
  double maxRapidity_;
  double tip_;
  double lip_;
  int    minHit_;
  bool signalOnly_;
  bool chargedOnly_;
  std::vector<int> pdgId_;

};

#include "CommonTools/UtilAlgos/interface/ParameterAdapter.h"

namespace reco {
  namespace modules {
    
    template<>
    struct ParameterAdapter<TrackingParticleSelector> {
      static TrackingParticleSelector make( const edm::ParameterSet & cfg ) {
	return TrackingParticleSelector(    
 	  cfg.getParameter<double>( "ptMin" ),
	  cfg.getParameter<double>( "minRapidity" ),
	  cfg.getParameter<double>( "maxRapidity" ),
	  cfg.getParameter<double>( "tip" ),
	  cfg.getParameter<double>( "lip" ),
	  cfg.getParameter<int>( "minHit" ), 
	  cfg.getParameter<bool>( "signalOnly" ),
	  cfg.getParameter<bool>( "chargedOnly" ),
	cfg.getParameter<std::vector<int> >( "pdgId" )); 
      }
    };
    
  }
}

#endif
