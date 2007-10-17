//
// $Id: TtDilepEvtSolution.h,v 1.9 2007/10/12 14:18:37 delaer Exp $
//

#ifndef TopObjects_TtDilepEvtSolution_h
#define TopObjects_TtDilepEvtSolution_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TopJet.h"
#include "AnalysisDataFormats/TopObjects/interface/TopElectron.h"
#include "AnalysisDataFormats/TopObjects/interface/TopMuon.h"
#include "AnalysisDataFormats/TopObjects/interface/TopMET.h"
#include "AnalysisDataFormats/TopObjects/interface/TopTau.h"

#include <vector>
#include <string>


class TtDilepEvtSolution {

  friend class TtDilepKinSolver;
  friend class TtDilepEvtSolutionMaker;

  public:

    TtDilepEvtSolution();
    virtual ~TtDilepEvtSolution();

    // methods to et the original TopObjects
    TopJet      getJetB() const;
    TopJet      getJetBbar() const;
    TopElectron getElectronp() const;
    TopElectron getElectronm() const;
    TopMuon     getMuonp() const;
    TopMuon     getMuonm() const;
    TopTau      getTaup() const;
    TopTau      getTaum() const;
    TopMET      getMET() const;
    // methods to get the MC matched particles
    const TtGenEvent &      getGenEvent() const;
    const reco::Candidate * getGenT() const;
    const reco::Candidate * getGenWp() const;
    const reco::Candidate * getGenB() const;
    const reco::Candidate * getGenLepp() const;
    const reco::Candidate * getGenN() const;
    const reco::Candidate * getGenTbar() const;
    const reco::Candidate * getGenWm() const;
    const reco::Candidate * getGenBbar() const;
    const reco::Candidate * getGenLepm() const;
    const reco::Candidate * getGenNbar() const;
    // methods to explicitly get reconstructed and calibrated objects 
    TopJetType  getRecJetB() const;
    TopJet      getCalJetB() const;
    TopJetType  getRecJetBbar() const;
    TopJet      getCalJetBbar() const;
    // method to get info on the W decays
    std::string getWpDecay() const { return wpDecay_; }
    std::string getWmDecay() const { return wmDecay_; }
    // method to get the matching residual
    double getResidual() const;
    // miscellaneous methods
    bool   getBestSol()      const { return bestSol_; }
    double getRecTopMass()   const {return topmass_; }
    double getRecWeightMax() const {return weightmax_; }
    
  protected:

    // method to set the generated event
    void setGenEvt(const edm::Handle<TtGenEvent> & aGenEvt);
    // methods to set the basic TopObjects
    void setB(const edm::Handle<std::vector<TopJet> > & jh, int i);
    void setBbar(const edm::Handle<std::vector<TopJet> > & jh, int i);
    void setMuonp(const edm::Handle<std::vector<TopMuon> > & mh, int i);
    void setMuonm(const edm::Handle<std::vector<TopMuon> > & mh, int i);
    void setTaup(const edm::Handle<std::vector<TopTau> > & mh, int i);
    void setTaum(const edm::Handle<std::vector<TopTau> > & mh, int i);
    void setElectronp(const edm::Handle<std::vector<TopElectron> > & eh, int i);
    void setElectronm(const edm::Handle<std::vector<TopElectron> > & eh, int i);
    void setMET(const edm::Handle<std::vector<TopMET> > & nh, int i);
    // miscellaneous methods
    void setBestSol(bool bs);
    void setRecTopMass(double j);
    void setRecWeightMax(double j);
    inline bool isValid() const { return wpDecay_!="NotDefined" && wmDecay_!="NotDefined"; }

  private:

    // particle content
    edm::RefProd<TtGenEvent>            theGenEvt_;
    edm::Ref<std::vector<TopElectron> > elecp_, elecm_;
    edm::Ref<std::vector<TopMuon> >     muonp_, muonm_;
    edm::Ref<std::vector<TopTau> >      taup_, taum_;
    edm::Ref<std::vector<TopJet> >      jetB_, jetBbar_;
    edm::Ref<std::vector<TopMET> >      met_;
    // miscellaneous
    std::string wpDecay_;
    std::string wmDecay_;      
    bool bestSol_;
    double topmass_;
    double weightmax_;

};


#endif
