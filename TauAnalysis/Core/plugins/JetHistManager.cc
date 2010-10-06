#include "TauAnalysis/Core/plugins/JetHistManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"
#include "TauAnalysis/GenSimTools/interface/genParticleAuxFunctions.h"

#include <TMath.h>

#include <assert.h>

bool matchesGenJet(const pat::Jet& patJet)
{
  bool isGenMatched = false;
// not implemented yet...
  return isGenMatched;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

JetHistManager::JetHistManager(const edm::ParameterSet& cfg)
  : HistManagerBase(cfg)
{
  //std::cout << "<JetHistManager::JetHistManager>:" << std::endl;

  jetSrc_ = cfg.getParameter<edm::InputTag>("jetSource");
  //std::cout << " jetSrc = " << jetSrc_ << std::endl;

  genParticleSrc_ = ( cfg.exists("genParticleSource") ) ? cfg.getParameter<edm::InputTag>("genParticleSource") : edm::InputTag();
  //std::cout << " genParticleSrc = " << genParticleSrc_ << std::endl;

  requireGenJetMatch_ = cfg.getParameter<bool>("requireGenJetMatch");
  //std::cout << " requireGenJetMatch = " << requireGenJetMatch_ << std::endl;

  skipPdgIdsGenParticleMatch_ = cfg.getParameter<vint>("skipPdgIdsGenParticleMatch");

  edm::ParameterSet cfgCentralJetsToBeVetoed = cfg.getParameter<edm::ParameterSet>("centralJetsToBeVetoed");  
  centralJetsToBeVetoedEtMin_ = cfgCentralJetsToBeVetoed.getParameter<vdouble>("etMin");
  centralJetsToBeVetoedEtaMax_ = cfgCentralJetsToBeVetoed.getParameter<vdouble>("etaMax");
  centralJetsToBeVetoedAlphaMin_ = cfgCentralJetsToBeVetoed.getParameter<vdouble>("alphaMin");

  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgBTaggingDiscriminators = cfg.getParameter<vParameterSet>("bTaggingDiscriminators");
  for ( vParameterSet::const_iterator cfgBTaggingDiscriminator = cfgBTaggingDiscriminators.begin();
	cfgBTaggingDiscriminator != cfgBTaggingDiscriminators.end(); ++cfgBTaggingDiscriminator ) {
    std::string bTaggingDiscrName = cfgBTaggingDiscriminator->getParameter<std::string>("name");
    double bTaggingDiscrThreshold = cfgBTaggingDiscriminator->getParameter<double>("threshold");

    bTaggingDiscriminators_.push_back(bTaggingDiscrName);
    bTaggingDiscriminatorThresholds_.push_back(bTaggingDiscrThreshold);
  }
}

JetHistManager::~JetHistManager()
{
//--- nothing to be done yet...
}

void JetHistManager::bookHistogramsImp()
{
  //std::cout << "<JetHistManager::bookHistogramsImp>:" << std::endl;
  
  hNumJets_ = book1D("NumJets", "NumJets", 20, -0.5, 19.5);
  hSumEtJets_ = book1D("SumEtJets", "SumEtJets", 100, -0.5, 999.5);
  
  bookJetHistograms(hJetPt_, hJetEta_, hJetPhi_, "Jet");
  hJetPtVsEta_ = book2D("JetPtVsEta", "Jet #eta vs P_{T}", 24, -3., +3., 30, 0., 150.);

  bookWeightHistograms(*dqmStore_, "JetWeight", "Jet Weight", 
		       hJetWeightPosLog_, hJetWeightNegLog_, hJetWeightZero_, 
		       hJetWeightLinear_);

  hJetMatchingGenParticlePdgId_ = book1D("JetMatchingGenParticlePdgId", "matching gen. Particle PdgId", 26, -1.5, 24.5);
  
  hJetAlpha_ = book1D("JetAlpha", "Jet #alpha", 102, -0.01, +1.01);
  hJetNumTracks_ = book1D("JetNumTracks", "Jet Track Multiplicity", 50, -0.5, 49.5);
  hJetTrkPt_ = book1D("JetTrkPt", "Jet All Tracks P_{T}", 100, 0., 50.);
  hJetLeadTrkPt_ = book1D("JetLeadTrkPt", "Jet Lead Track P_{T}", 75, 0., 75.);

  for ( vdouble::const_iterator etMin = centralJetsToBeVetoedEtMin_.begin();
	etMin != centralJetsToBeVetoedEtMin_.end(); ++etMin ) {
    for ( vdouble::const_iterator etaMax = centralJetsToBeVetoedEtaMax_.begin();
	  etaMax != centralJetsToBeVetoedEtaMax_.end(); ++etaMax ) {
      for ( vdouble::const_iterator alphaMin = centralJetsToBeVetoedAlphaMin_.begin();
	    alphaMin != centralJetsToBeVetoedAlphaMin_.end(); ++alphaMin ) {
	std::ostringstream hName;
	hName.setf(std::ios::fixed);
	hName.precision(1);
	hName << "numJetsEtGt" << (*etMin) << "EtaLt" << (*etaMax) << "AlphaGt" << (*alphaMin);
	int errorFlag = 0;
	std::string hName_mod = replace_string(hName.str(), ".", "_", 0, 3, errorFlag);
	//std::cout << "hName_mod = " << hName_mod << std::endl;
	
	hNumCentralJetsToBeVetoed_.push_back(book1D(hName_mod, hName_mod, 20, -0.5, 19.5));
      }
    }
  }
  
  assert(bTaggingDiscriminators_.size()== bTaggingDiscriminatorThresholds_.size());
  for (unsigned int i=0;i<bTaggingDiscriminators_.size();i++) {
    std::ostringstream hName0;
    std::ostringstream hTitle0;
    hName0 << "BtagDisc_" << bTaggingDiscriminators_[i];
    hTitle0 << "Jet B-Tag Discriminator: " << bTaggingDiscriminators_[i];
  
    hBtagDisc_.push_back(book1D(hName0.str(),hTitle0.str(), 120, -10., 50.));
    
    std::ostringstream hName1;
    std::ostringstream hTitle1;
    hName1 << "NumBtags_" << bTaggingDiscriminators_[i];
    hTitle1 << "Jet B-Tag Count: " << bTaggingDiscriminators_[i] << ">" << bTaggingDiscriminatorThresholds_[i];
  
    hNumBtags_.push_back(book1D(hName1.str(),hTitle1.str(), 15, -0.5, 14.5));
    
    std::ostringstream hName2;
    std::ostringstream hTitle2;
    hName2 << "PtBtags_" << bTaggingDiscriminators_[i];
    hTitle2 << "B-Tagged Jets P_{T}: " << bTaggingDiscriminators_[i] << ">" << bTaggingDiscriminatorThresholds_[i];
  
    hPtBtags_.push_back(book1D(hName2.str(),hTitle2.str(), 75, 0., 150.));
  }
  
}

double JetHistManager::getJetWeight(const pat::Jet& patJet)
{
  return 1.;
}

void JetHistManager::fillHistogramsImp(const edm::Event& evt, const edm::EventSetup& es, double evtWeight)
{  
  //std::cout << "<JetHistManager::fillHistogramsImp>:" << std::endl; 

  edm::Handle<pat::JetCollection> patJets;
  getCollection(evt, jetSrc_, patJets);

  edm::Handle<reco::GenParticleCollection> genParticles;
  if ( genParticleSrc_.label() != "" ) evt.getByLabel(genParticleSrc_, genParticles);

  //std::cout << " patJets.size = " << patJets->size() << std::endl;
  hNumJets_->Fill(patJets->size(), evtWeight);
  
  double sumJetEt = 0.;
  for ( std::vector<pat::Jet>::const_iterator patJet = patJets->begin(); 
	patJet != patJets->end(); ++patJet ) {
    sumJetEt += patJet->et();
  }
  hSumEtJets_->Fill(sumJetEt, evtWeight);

  double jetWeightSum = 0.;
  for ( std::vector<pat::Jet>::const_iterator patJet = patJets->begin(); 
	patJet != patJets->end(); ++patJet ) {
    if ( requireGenJetMatch_ && !matchesGenJet(*patJet) ) continue;

    jetWeightSum += getJetWeight(*patJet);
  }
  
  std::vector<unsigned int> nbtags;
  for ( unsigned int i = 0; i < bTaggingDiscriminators_.size(); ++i ) {
    nbtags.push_back(0);
  }

  for ( std::vector<pat::Jet>::const_iterator patJet = patJets->begin(); 
	patJet != patJets->end(); ++patJet ) {
  
    //bool isGenJetMatched = matchesGenJet(*patJet);
    //std::cout << " Pt = " << patJet->pt() << ", eta = " << patJet->eta() << ", phi = " << patJet->phi() << std::endl;
    //std::cout << " isGenJetMatched = " << isGenJetMatched << std::endl;

    if ( requireGenJetMatch_ && !matchesGenJet(*patJet) ) continue;

    double jetWeight = getJetWeight(*patJet);
    double weight = getWeight(evtWeight, jetWeight, jetWeightSum);

    fillJetHistograms(*patJet, hJetPt_, hJetEta_, hJetPhi_, weight);
    hJetPtVsEta_->Fill(patJet->eta(), patJet->pt(), weight);

    fillWeightHistograms(hJetWeightPosLog_, hJetWeightNegLog_, hJetWeightZero_, 
			 hJetWeightLinear_, jetWeight);
    
    int matchingGenParticlePdgId = -1;
    if( genParticles.isValid() )
		matchingGenParticlePdgId = getMatchingGenParticlePdgId(patJet->p4(), *genParticles, &skipPdgIdsGenParticleMatch_);
    if ( matchingGenParticlePdgId == -1 ) {
      hJetMatchingGenParticlePdgId_->Fill(-1, weight);
    } else if ( abs(matchingGenParticlePdgId) > 22 ) {
      hJetMatchingGenParticlePdgId_->Fill(24, weight);
    } else {
      hJetMatchingGenParticlePdgId_->Fill(abs(matchingGenParticlePdgId), weight);
    }

    hJetAlpha_->Fill(jetAlphaExtractor_(*patJet), weight);
    unsigned numTracks = 0;
    double maxPt = 0.;
    for ( reco::TrackRefVector::const_iterator track = patJet->associatedTracks().begin();
	  track != patJet->associatedTracks().end(); ++track ) {
      ++numTracks;
      hJetTrkPt_->Fill((*track)->pt(), weight);
      if ( (*track)->pt() > maxPt ) maxPt = (*track)->pt();
    }
    hJetNumTracks_->Fill(numTracks, weight);
    if ( numTracks > 0 ) hJetLeadTrkPt_->Fill(maxPt, weight);

    for ( unsigned int i = 0; i < bTaggingDiscriminators_.size(); ++i ) {
      hBtagDisc_[i]->Fill(patJet->bDiscriminator(bTaggingDiscriminators_[i]), weight);

      if (patJet->bDiscriminator(bTaggingDiscriminators_[i])> bTaggingDiscriminatorThresholds_[i]) {
        hPtBtags_[i]->Fill(patJet->pt(), weight);
	++nbtags[i];
      }
    }
  }
  
  for (unsigned int i=0;i<bTaggingDiscriminators_.size();i++) {
    hNumBtags_[i]->Fill(nbtags[i], evtWeight);
  }

  int index = 0;
  for ( vdouble::const_iterator etMin = centralJetsToBeVetoedEtMin_.begin();
	etMin != centralJetsToBeVetoedEtMin_.end(); ++etMin ) {
    for ( vdouble::const_iterator etaMax = centralJetsToBeVetoedEtaMax_.begin();
	  etaMax != centralJetsToBeVetoedEtaMax_.end(); ++etaMax ) {
      for ( vdouble::const_iterator alphaMin = centralJetsToBeVetoedAlphaMin_.begin();
	    alphaMin != centralJetsToBeVetoedAlphaMin_.end(); ++alphaMin ) {
	fillNumCentralJetsToBeVetoesHistograms(*patJets, hNumCentralJetsToBeVetoed_[index], *etMin, *etaMax, *alphaMin, evtWeight);
	++index;
      }
    }
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void JetHistManager::bookJetHistograms(MonitorElement*& hJetPt, MonitorElement*& hJetEta, MonitorElement*& hJetPhi, const char* histoSetName)
{
  std::string hJetPtName = std::string(histoSetName).append("Pt");
  hJetPt = book1D(hJetPtName, hJetPtName, 75, 0., 150.);
  
  std::string hJetEtaName = std::string(histoSetName).append("Eta");
  hJetEta = book1D(hJetEtaName, hJetEtaName, 60, -3., +3.);
  
  std::string hJetPhiName = std::string(histoSetName).append("Phi");
  hJetPhi = book1D(hJetPhiName, hJetPhiName, 36, -TMath::Pi(), +TMath::Pi());
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void JetHistManager::fillJetHistograms(const pat::Jet& patJet, 
				       MonitorElement* hJetPt, MonitorElement* hJetEta, MonitorElement* hJetPhi,
				       double weight)
{
  //std::cout << "<JetHistManager::fillJetHistograms>:" << std::endl;

  hJetPt->Fill(patJet.pt(), weight);
  hJetEta->Fill(patJet.eta(), weight);
  hJetPhi->Fill(patJet.phi(), weight);
}

void JetHistManager::fillNumCentralJetsToBeVetoesHistograms(const std::vector<pat::Jet>& patJets, MonitorElement* hNumJets,
							    double etMin, double etaMax, double alphaMin, double weight)
{
  unsigned numJets = 0;
  for ( std::vector<pat::Jet>::const_iterator patJet = patJets.begin();
	patJet != patJets.end(); ++patJet ) {
    if ( patJet->et() >= etMin && patJet->eta() <= etaMax && jetAlphaExtractor_(*patJet) >= alphaMin ) ++numJets;
  }
  
  hNumJets->Fill(numJets, weight);
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, JetHistManager, "JetHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, JetHistManager, "JetHistManager");

#include "TauAnalysis/Core/interface/HistManagerAdapter.h"

typedef HistManagerAdapter<JetHistManager> JetAnalyzer;

DEFINE_FWK_MODULE(JetAnalyzer);
