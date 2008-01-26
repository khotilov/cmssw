//
// $Id: Tau.h,v 1.6 2008/01/23 16:12:19 lowette Exp $
//

#ifndef DataFormats_PatCandidates_Tau_h
#define DataFormats_PatCandidates_Tau_h

/**
  \class    Tau Tau.h "DataFormats/PatCandidates/interface/Tau.h"
  \brief    Analysis-level tau class

   Tau implements the analysis-level tau class within the 'pat' namespace.

  \author   Steven Lowette
  \version  $Id: Tau.h,v 1.6 2008/01/23 16:12:19 lowette Exp $
*/


#include "DataFormats/TauReco/interface/BaseTau.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"


namespace pat {


  typedef reco::BaseTau TauType;


  class Tau : public Lepton<TauType> {

    public:

      Tau();
      Tau(const TauType & aTau);
      Tau(const edm::RefToBase<TauType> & aTauRef);
      virtual ~Tau();

      float emEnergyFraction() const { return emEnergyFraction_; }
      float eOverP() const { return eOverP_; }

      void setEmEnergyFraction(float fraction) { emEnergyFraction_ = fraction; }
      void setEOverP(float EoP) { eOverP_ = EoP; } 

    private:

      float emEnergyFraction_;
      float eOverP_;

  };


}

#endif
