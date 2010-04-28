#include "TauAnalysis/RecoTools/plugins/PATLeptonIsoDepositSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/IsolationAlgos/interface/IsoDepositVetoFactory.h"

#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include <sstream>

std::string getString(double value)
{
  std::ostringstream value_string;
  value_string << value;
  return value_string.str();
}

template <class T>
PATLeptonIsoDepositSelector<T>::PATLeptonIsoDepositSelector(const edm::ParameterSet& cfg)
  : sumPtMethod_(kAbsoluteIso),
    cfgError_(0)
{
//--- translate isolation "keys" from string to integer format
//    (NOTE: the equivalence between strinsg and integers is defined in
//           PhysicsTools/PatAlgos/src/MultiIsolator.cc )
  std::string isoDepositType_string = cfg.getParameter<std::string>("type");
  if ( isoDepositType_string == "tracker" ) {
    isoDepositType_ = pat::TrackIso;
  } else if ( isoDepositType_string == "ecal" ) {
    isoDepositType_ = pat::EcalIso;
  } else if ( isoDepositType_string == "hcal" ) {
    isoDepositType_ = pat::HcalIso;
  } else if ( isoDepositType_string == "pfAllParticles" ) {
    isoDepositType_ = pat::PfAllParticleIso;
  } else if ( isoDepositType_string == "pfChargedHadron" ) {
    isoDepositType_ = pat::PfChargedHadronIso;
  } else if ( isoDepositType_string == "pfNeutralHadron" ) {
    isoDepositType_ = pat::PfNeutralHadronIso;
  } else if ( isoDepositType_string == "pfGamma" ) {
    isoDepositType_ = pat::PfGammaIso;
  } else {
    edm::LogError ("PATLeptonIsoDepositSelector") 
      << " Type = " << isoDepositType_string << " not defined for IsoDeposits !!";
    cfgError_ = 1;
  }

  if ( cfg.exists("dRvetoCone") ) {
    double dRvetoCone = cfg.getParameter<double>("dRvetoCone");
    std::string dRvetoCone_string = getString(dRvetoCone);
    isoParam_.push_back(IsoDepositVetoFactory::make(dRvetoCone_string.data()));
  }

  if ( cfg.exists("ptMin") ) {
    double ptMin = cfg.getParameter<double>("ptMin");
    std::string ptMin_string = std::string("Threshold(").append(getString(ptMin)).append(")");
    isoParam_.push_back(IsoDepositVetoFactory::make(ptMin_string.data()));
  }

  if ( cfg.exists("vetos") ) {
    typedef std::vector<std::string> vstring;
    vstring vetos = cfg.getParameter<vstring>("vetos");
    for ( vstring::const_iterator veto = vetos.begin();
	  veto != vetos.end(); ++veto ) {
      isoParam_.push_back(IsoDepositVetoFactory::make(veto->data()));
    }
  }

  dRisoCone_ = cfg.getParameter<double>("dRisoCone");
  
  sumPtMin_ = cfg.exists("sumPtMin") ? cfg.getParameter<double>("sumPtMin") : -1.;
  numMin_ = cfg.exists("numMin") ? cfg.getParameter<int>("numMin") : -1;

  sumPtMax_ = cfg.exists("sumPtMax") ? cfg.getParameter<double>("sumPtMax") : -1.;
  numMax_ = cfg.exists("numMax") ? cfg.getParameter<int>("numMax") : -1;
  if ( !(cfg.exists("sumPtMax") || cfg.exists("numMax")) ) {
    edm::LogError ("PATLeptonIsoDepositSelector") 
      << " Need to specify either 'sumPtMax' or 'numMax' Parameters !!";
    cfgError_ = 1;
  }

  if ( cfg.exists("sumPtMethod") ) {
    std::string sumPtMethod_string = cfg.getParameter<std::string>("sumPtMethod");
    if ( sumPtMethod_string == "absolute" ) {
      sumPtMethod_ = kAbsoluteIso;
    } else if ( sumPtMethod_string == "relative" ) {
      sumPtMethod_ = kRelativeIso;
    } else {
      edm::LogError("PATLeptonIsoDepositSelector") 
	<< " Configuration parameter 'sumPtMethod' = " << sumPtMethod_string << " invalid !!";
      cfgError_ = 1;
    }
  }
}

template <class T>
PATLeptonIsoDepositSelector<T>::~PATLeptonIsoDepositSelector()
{
  for ( reco::isodeposit::AbsVetos::const_iterator it = isoParam_.begin();
	it != isoParam_.end(); ++it ) {
    delete (*it);
  }
}

template <class T>
void PATLeptonIsoDepositSelector<T>::select(const edm::Handle<collection>& patLeptons, 
					    const edm::Event& evt, const edm::EventSetup& es)
{
  selected_.clear();
  
  if ( cfgError_ ) {
    edm::LogError ("select") << " Error in Configuration ParameterSet --> no pat::Leptons will pass Selection !!";
    return;
  }

  for ( typename collection::const_iterator patLepton = patLeptons->begin();
	patLepton != patLeptons->end(); ++patLepton ) {
    const reco::IsoDeposit* isoDeposit = patLepton->isoDeposit(isoDepositType_);

    if ( isoDeposit ) {
      double sumPt = isoDeposit->depositWithin(dRisoCone_, isoParam_);
      //std::cout << "sumPt = " << sumPt << std::endl;
      double num = isoDeposit->countWithin(dRisoCone_, isoParam_);
      //std::cout << "num = " << num << std::endl;

      if ( sumPtMethod_ == kAbsoluteIso ) {
	if ( sumPtMin_ >  0. && sumPt < sumPtMin_ ) continue;
	if ( sumPtMax_ >  0. && sumPt > sumPtMax_ ) continue;
      } else if ( sumPtMethod_ == kRelativeIso ) {
	double relIso = ( patLepton->pt() > 1. ) ? (sumPt/patLepton->pt()) : sumPt;
	if ( sumPtMin_ >  0. && relIso < sumPtMin_ ) continue;
	if ( sumPtMax_ >  0. && relIso > sumPtMax_ ) continue;
      }

      if ( numMin_   >= 0  && num   < numMin_   ) continue;
      if ( numMax_   >= 0  && num   > numMax_   ) continue;

      selected_.push_back(&(*patLepton)); 
    } else {
      edm::LogError ("select") 
	<< " No IsoDeposit of type = " << isoDepositType_ << " associated to pat::Lepton with "
	<< " Pt = " << patLepton->pt() << "," 
	<< " theta = " << patLepton->theta() << "," 
	<< " phi = " << patLepton->phi() << " --> pat::Lepton fails Selection !!";
      continue;
    }
  }
}

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

typedef ObjectSelector<PATLeptonIsoDepositSelector<pat::Electron> > PATElectronIsoDepositSelector;
typedef ObjectSelector<PATLeptonIsoDepositSelector<pat::Muon> > PATMuonIsoDepositSelector;
typedef ObjectSelector<PATLeptonIsoDepositSelector<pat::Tau> > PATTauIsoDepositSelector;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_ANOTHER_FWK_MODULE(PATElectronIsoDepositSelector);
DEFINE_ANOTHER_FWK_MODULE(PATMuonIsoDepositSelector);
DEFINE_ANOTHER_FWK_MODULE(PATTauIsoDepositSelector);
