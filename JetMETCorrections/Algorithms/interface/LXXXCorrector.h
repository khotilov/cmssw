// Generic LX jet corrector class. Inherits from JetCorrector.h
#ifndef LXXXCorrector_h
#define LXXXCorrector_h

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

//----- classes declaration -----------------------------------
namespace edm 
{
  class ParameterSet;
}
class FactorizedJetCorrector;
//----- LXXXCorrector interface -------------------------------
class LXXXCorrector : public JetCorrector 
{
  public:
    //----- constructors---------------------------------------
    LXXXCorrector(const JetCorrectorParameters& fConfig, const edm::ParameterSet& fParameters);   

    //----- destructor ----------------------------------------
    virtual ~LXXXCorrector();

    //----- apply correction using Jet information only -------
    virtual double correction(const LorentzVector& fJet) const;

    //----- apply correction using Jet information only -------
    virtual double correction(const reco::Jet& fJet) const;

    //----- if correction needs event information -------------
    virtual bool eventRequired() const {return false;} 

  private:
    //----- member data ---------------------------------------
    unsigned mLevel;
    FactorizedJetCorrector* mCorrector;
};

#endif
