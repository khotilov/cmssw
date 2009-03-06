#ifndef TopObjects_TtFullLeptonicEvent_h
#define TopObjects_TtFullLeptonicEvent_h

#include "AnalysisDataFormats/TopObjects/interface/TtEvent.h"

namespace TtFullLepDaughter{
  /// full leptonic daughter names for common
  /// use and use with the hypotheses
  static const std::string Nu   ="Nu"   , LepBar="LepBar", WPlus ="WPlus" , B   ="B"   , Top   ="Top";
  static const std::string NuBar="NuBar", Lep   ="Lep"   , WMinus="WMinus", BBar="BBar", TopBar="TopBar"; 
}

// ----------------------------------------------------------------------
// derived class for: 
//
//  * TtSemiLeptonicEvent
//
//  the structure holds information on the leptonic decay channels, 
//  all event hypotheses of different classes (user defined during
//  production) and a reference to the TtGenEvent (if available) 
//  and provides access and administration; the derived class 
//  contains a few additional getters with respect to its base class
// ----------------------------------------------------------------------

class TtFullLeptonicEvent: public TtEvent {
  
 public:

  /// empty constructor
  TtFullLeptonicEvent(){};
  /// default destructor
  virtual ~TtFullLeptonicEvent(){};

  /// get top of the given hypothesis
  const reco::Candidate* top        (const HypoClassKey& key, const unsigned& cmb=0) const { return !isHypoValid(key,cmb) ? 0 : eventHypo(key,cmb). daughter(TtFullLepDaughter::Top   ); };
  /// get b of the given hypothesis
  const reco::Candidate* b          (const HypoClassKey& key, const unsigned& cmb=0) const { return !isHypoValid(key,cmb) ? 0 : top      (key,cmb)->daughter(TtFullLepDaughter::B     ); };
  /// get Wplus of the given hypothesis
  const reco::Candidate* wPlus      (const HypoClassKey& key, const unsigned& cmb=0) const { return !isHypoValid(key,cmb) ? 0 : top      (key,cmb)->daughter(TtFullLepDaughter::WPlus ); };
  /// get anti-lepton of the given hypothesis
  const reco::Candidate* leptonBar  (const HypoClassKey& key, const unsigned& cmb=0) const { return !isHypoValid(key,cmb) ? 0 : wPlus    (key,cmb)->daughter(TtFullLepDaughter::LepBar); };
  /// get neutrino of the given hypothesis
  const reco::Candidate* neutrino   (const HypoClassKey& key, const unsigned& cmb=0) const { return !isHypoValid(key,cmb) ? 0 : wPlus    (key,cmb)->daughter(TtFullLepDaughter::Nu    ); };
  /// get anti-top of the given hypothesis
  const reco::Candidate* topBar     (const HypoClassKey& key, const unsigned& cmb=0) const { return !isHypoValid(key,cmb) ? 0 : eventHypo(key,cmb). daughter(TtFullLepDaughter::TopBar); };
  /// get anti-b of the given hypothesis
  const reco::Candidate* bBar       (const HypoClassKey& key, const unsigned& cmb=0) const { return !isHypoValid(key,cmb) ? 0 : topBar   (key,cmb)->daughter(TtFullLepDaughter::BBar  ); };
  /// get Wminus of the given hypothesis
  const reco::Candidate* wMinus     (const HypoClassKey& key, const unsigned& cmb=0) const { return !isHypoValid(key,cmb) ? 0 : topBar   (key,cmb)->daughter(TtFullLepDaughter::WMinus); };
  /// get lepton of the given hypothesis
  const reco::Candidate* lepton     (const HypoClassKey& key, const unsigned& cmb=0) const { return !isHypoValid(key,cmb) ? 0 : wMinus   (key,cmb)->daughter(TtFullLepDaughter::Lep   ); };
  /// get anti-neutrino of the given hypothesis
  const reco::Candidate* neutrinoBar(const HypoClassKey& key, const unsigned& cmb=0) const { return !isHypoValid(key,cmb) ? 0 : wMinus   (key,cmb)->daughter(TtFullLepDaughter::NuBar ); };

  /// get top of the TtGenEvent
  const reco::GenParticle* genTop        () const { return (!genEvt_ ? 0 : this->genEvent()->top()        ); };
  /// get b of the TtGenEvent
  const reco::GenParticle* genB          () const { return (!genEvt_ ? 0 : this->genEvent()->b()          ); };
  /// get Wplus of the TtGenEvent
  const reco::GenParticle* genWPlus      () const { return (!genEvt_ ? 0 : this->genEvent()->wPlus()      ); };
  /// get anti-lepton of the TtGenEvent
  const reco::GenParticle* genLeptonBar  () const { return (!genEvt_ ? 0 : this->genEvent()->leptonBar()  ); };
  /// get neutrino of the TtGenEvent
  const reco::GenParticle* genNeutrino   () const { return (!genEvt_ ? 0 : this->genEvent()->neutrino()   ); };
  /// get anti-top of the TtGenEvent
  const reco::GenParticle* genTopBar     () const { return (!genEvt_ ? 0 : this->genEvent()->topBar()     ); };
  /// get anti-b of the TtGenEvent
  const reco::GenParticle* genBBar       () const { return (!genEvt_ ? 0 : this->genEvent()->bBar()       ); };
  /// get Wminus of the TtGenEvent
  const reco::GenParticle* genWMinus     () const { return (!genEvt_ ? 0 : this->genEvent()->wMinus()     ); };
  /// get lepton of the TtGenEvent
  const reco::GenParticle* genLepton     () const { return (!genEvt_ ? 0 : this->genEvent()->lepton()     ); };
  /// get anti-neutrino of the TtGenEvent
  const reco::GenParticle* genNeutrinoBar() const { return (!genEvt_ ? 0 : this->genEvent()->neutrinoBar()); };

  /// print full content of the structure as formated 
  /// LogInfo to the MessageLogger output for debugging  
  void print();
};

#endif
