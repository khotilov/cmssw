//
// $Id: TauCaloSpecific.h,v 1.1 2008/06/09 09:03:19 gpetrucc Exp $
//

#ifndef DataFormats_PatCandidates_Tau_CaloSpecific_h
#define DataFormats_PatCandidates_Tau_CaloSpecific_h

/**
  \class    pat::tau::CaloSpecific TauCaloSpecific.h "DataFormats/PatCandidates/interface/TauCaloSpecific.h"
  \brief    Structure to hold information specific to a CaloTau inside a pat::Tau

  \author   Giovanni Petrucciani
  \version  $Id: TauCaloSpecific.h,v 1.1 2008/06/09 09:03:19 gpetrucc Exp $
*/

#include "DataFormats/TauReco/interface/CaloTau.h"

namespace pat { namespace tau {

struct TauCaloSpecific {
// dummy constructor for ROOT I/O
    TauCaloSpecific() {}
// constructor from CaloTau
    TauCaloSpecific(const reco::CaloTau &tau) ;
// datamembers 
    reco::CaloTauTagInfoRef CaloTauTagInfoRef_;
    float leadTracksignedSipt_;
    float leadTrackHCAL3x3hitsEtSum_;
    float leadTrackHCAL3x3hottesthitDEta_;
    float signalTracksInvariantMass_;
    float TracksInvariantMass_; 
    float isolationTracksPtSum_;
    float isolationECALhitsEtSum_;
    float maximumHCALhitEt_;
    float etaetaMoment_;
    float phiphiMoment_;
    float etaphiMoment_;
};

} }

#endif
