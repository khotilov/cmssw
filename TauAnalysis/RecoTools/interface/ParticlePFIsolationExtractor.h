#ifndef TauAnalysis_RecoTools_ParticlePFIsolationExtractor_h
#define TauAnalysis_RecoTools_ParticlePFIsolationExtractor_h

/** \class ParticlePFIsolationExtractor
 *
 * Auxiliary class to compute isolation of pat::Leptons based on charged hadrons, 
 * neutral hadrons and photons/pi0s reconstructed by particle-flow algorithm
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.14 $
 *
 * $Id: ParticlePFIsolationExtractor.h,v 1.14 2011/11/17 14:28:05 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TauAnalysis/RecoTools/interface/PATLeptonTrackVectorExtractor.h"
#include "TauAnalysis/RecoTools/interface/pfCandAuxFunctions.h"

#include <TMath.h>

#include <vector>

template <class T>
class ParticlePFIsolationExtractor
{
 public:
  enum { kNone, kBeta, kDeltaBeta, kRho };
  enum { kDirP4, kDirTrack };
  explicit ParticlePFIsolationExtractor(const edm::ParameterSet& cfg):
    pfChargedParticleIso_(0),
    addChargedParticleIso_(false),
    pfChargedParticleIsoConeSize_(0.),
    pfNeutralHadronIso_(0),
    addNeutralHadronIso_(false),
    pfNeutralHadronIsoConeSize_(0.),
    pfPhotonIso_(0),
    addPhotonIso_(false),
    pfPhotonIsoConeSize_(0.),
    methodPUcorr_(kNone),
    trackExtractor_(0),
    pfNeutralHadronIsoPUcorr_(0),
    pfPhotonIsoPUcorr_(0)
  {    
    if( cfg.exists("chargedParticleIso") ) {
      edm::ParameterSet cfgChargedParticleIso = cfg.getParameter<edm::ParameterSet>("chargedParticleIso");
      pfChargedParticleIso_ = new pfIsoConfigType(reco::PFCandidate::h, cfgChargedParticleIso);
      addChargedParticleIso_ = cfgChargedParticleIso.exists("add") ? cfgChargedParticleIso.getParameter<bool>("add") : true;
      pfChargedParticleIsoConeSize_ = cfgChargedParticleIso.getParameter<double>("dRisoCone");
    }

    if( cfg.exists("neutralHadronIso") ) {
      edm::ParameterSet cfgNeutralHadronIso = cfg.getParameter<edm::ParameterSet>("neutralHadronIso");
      pfNeutralHadronIso_ = new pfIsoConfigType(reco::PFCandidate::h0, cfgNeutralHadronIso);
      addNeutralHadronIso_ = cfgNeutralHadronIso.exists("add") ? cfgNeutralHadronIso.getParameter<bool>("add") : true;
      pfNeutralHadronIsoConeSize_ = cfgNeutralHadronIso.getParameter<double>("dRisoCone");
    }

    if( cfg.exists("photonIso") ) {
      edm::ParameterSet cfgPhotonIso = cfg.getParameter<edm::ParameterSet>("photonIso");
      pfPhotonIso_ = new pfIsoConfigType(reco::PFCandidate::gamma, cfgPhotonIso);
      addPhotonIso_ = cfgPhotonIso.exists("add") ? cfgPhotonIso.getParameter<bool>("add") : true;
      pfPhotonIsoConeSize_ = cfgPhotonIso.getParameter<double>("dRisoCone");
    }

    if ( cfg.exists("pileUpCorr") ) {
      edm::ParameterSet cfgPUcorr = cfg.getParameter<edm::ParameterSet>("pileUpCorr");

      std::string method_string = cfgPUcorr.getParameter<std::string>("method");
      if      ( method_string == "beta"      ) methodPUcorr_ = kBeta;
      else if ( method_string == "deltaBeta" ) methodPUcorr_ = kDeltaBeta;
      else if ( method_string == "rho"       ) methodPUcorr_ = kRho;
      else if ( method_string == "none"      ) methodPUcorr_ = kNone;
      else 
	throw cms::Exception("ParticlePFIsolationExtractor")
	  << "Invalid Configuration parameter method = " << method_string << "!!\n";
      
      if ( methodPUcorr_ != kNone ) {
	trackExtractor_ = new PATLeptonTrackVectorExtractor<T>();
	
	if ( methodPUcorr_ == kDeltaBeta ) chargedToNeutralFactor_ = cfgPUcorr.getParameter<double>("chargedToNeutralFactor");
	
	if ( methodPUcorr_ == kBeta || methodPUcorr_ == kDeltaBeta ) {
	  if ( pfNeutralHadronIso_ ) {
	    edm::ParameterSet cfgNeutralHadronIso = cfg.getParameter<edm::ParameterSet>("neutralHadronIso");
	    pfNeutralHadronIsoPUcorr_ = new pfIsoConfigType(reco::PFCandidate::h, cfgNeutralHadronIso);
	  }
	  
	  if ( pfPhotonIso_ ) {
	    edm::ParameterSet cfgPhotonIso = cfg.getParameter<edm::ParameterSet>("photonIso");
	    pfPhotonIsoPUcorr_ = new pfIsoConfigType(reco::PFCandidate::h, cfgPhotonIso);
	  }
	  
	  if ( !pfChargedParticleIso_ ) 
	    throw cms::Exception("ParticlePFIsolationExtractor")
	      << "Pile-up correction requires 'chargedParticleIso' Configuration parameter !!\n";
	}

	if ( methodPUcorr_ == kRho ) {	  
	  ueRhoOffset_ = cfgPUcorr.getParameter<double>("ueRhoOffset");

	  if ( !(pfChargedParticleIsoConeSize_ == pfNeutralHadronIsoConeSize_ && 	
		 pfChargedParticleIsoConeSize_ == pfPhotonIsoConeSize_        ) ) 
	    throw cms::Exception("ParticlePFIsolationExtractor")
	      << "Rho FastJet corrections require equal isolation Cone sizes to be used" 
	      << " for PFChargedParticles, PFNeutralHadrons and PFGammas !!\n";								 
	}
      }
    }
  }
  ~ParticlePFIsolationExtractor()
  {
    delete pfChargedParticleIso_;
    delete pfNeutralHadronIso_;
    delete pfPhotonIso_;

    delete trackExtractor_;
    delete pfNeutralHadronIsoPUcorr_;
    delete pfPhotonIsoPUcorr_;
  }

  
  double operator()(const T& lepton, int direction,
		    const reco::PFCandidateCollection& pfCandidates, const reco::PFCandidateCollection& pfNoPileUpCandidates,
		    double rhoFastJetCorrection = 0.)
  {
    reco::Particle::Vector coneAxis;
    if      ( direction == kDirP4    ) coneAxis = lepton.momentum();
    else if ( direction == kDirTrack ) {
      const reco::Track* leadingTrack = 0;
      std::vector<const reco::Track*> signalTracks = (*trackExtractor_)(lepton);
      for ( std::vector<const reco::Track*>::const_iterator signalTrack = signalTracks.begin();
	    signalTrack != signalTracks.end(); ++signalTrack ) {
	if ( leadingTrack == 0 || (*signalTrack)->pt() > leadingTrack->pt() ) leadingTrack = (*signalTrack);
      }
      if ( leadingTrack ) coneAxis = leadingTrack->momentum();
      else                coneAxis = lepton.momentum();
    } else 
      throw cms::Exception("ParticlePFIsolationExtractor")
	<< "Invalid function argument 'direction' = " << direction << " !!\n";

    return this->operator()(lepton, coneAxis, pfCandidates, pfNoPileUpCandidates, rhoFastJetCorrection);
  }

  double operator()(const T& lepton, const reco::Particle::Vector& coneAxis,
		    const reco::PFCandidateCollection& pfCandidates, const reco::PFCandidateCollection& pfNoPileUpCandidates,
		    double rhoFastJetCorrection = 0.)
  {
    std::vector<const reco::PFCandidate*> pfChargedParticles, pfChargedHadrons, pfNeutralHadrons, pfPhotons;
    if ( addChargedParticleIso_ || methodPUcorr_ != kNone ) {
      for ( reco::PFCandidateCollection::const_iterator pfCandidate = pfCandidates.begin();
	    pfCandidate != pfCandidates.end(); ++pfCandidate ) {
	if ( TMath::Abs(pfCandidate->charge()) > 0.5 ) pfChargedParticles.push_back(&(*pfCandidate));
      }
      pfChargedHadrons = getPFCandidatesOfType(pfCandidates, reco::PFCandidate::h);
    }
    if ( addNeutralHadronIso_   ) pfNeutralHadrons = getPFCandidatesOfType(pfCandidates, reco::PFCandidate::h0);
    if ( addPhotonIso_          ) pfPhotons        = getPFCandidatesOfType(pfCandidates, reco::PFCandidate::gamma);

    double sumPt = 0.;
    
    if ( methodPUcorr_ == kNone ) {
      if ( addChargedParticleIso_ ) sumPt += pfChargedParticleIso_->compSumPt(pfChargedParticles, coneAxis);
      if ( addNeutralHadronIso_   ) sumPt += pfNeutralHadronIso_->compSumPt(pfNeutralHadrons, coneAxis);
      if ( addPhotonIso_          ) sumPt += pfPhotonIso_->compSumPt(pfPhotons, coneAxis);
    } else {
      std::vector<const reco::Track*> signalTracks = (*trackExtractor_)(lepton);

      std::vector<const reco::PFCandidate*> pfNoPileUpChargedParticles, pfPileUpChargedParticles;
      getPileUpPFCandidates(pfChargedParticles, pfNoPileUpCandidates, 1.e-3,
			    pfNoPileUpChargedParticles, pfPileUpChargedParticles);
      std::vector<const reco::PFCandidate*> pfNoPileUpChargedHadrons, pfPileUpChargedHadrons;
      getPileUpPFCandidates(pfChargedHadrons, pfNoPileUpCandidates, 1.e-3,
			    pfNoPileUpChargedHadrons, pfPileUpChargedHadrons);
      
      if ( methodPUcorr_ == kBeta ) {
	double pfNoPileUpChargedHadronIsoSumPt = pfChargedParticleIso_->compSumPt(pfNoPileUpChargedHadrons, coneAxis);
	double pfPileUpChargedHadronIsoSumPt  = pfChargedParticleIso_->compSumPt(pfPileUpChargedHadrons, coneAxis);
	sumPt = pfNoPileUpChargedHadronIsoSumPt;
	
	double pfNeutralIsoCorrFactor = ( pfPileUpChargedHadronIsoSumPt > 0. ) ? 
	  (pfNoPileUpChargedHadronIsoSumPt/(pfNoPileUpChargedHadronIsoSumPt + pfPileUpChargedHadronIsoSumPt)) : 1.0;
	
	if ( pfNeutralHadronIso_ ) sumPt += pfNeutralHadronIso_->compSumPt(pfNeutralHadrons, coneAxis)*pfNeutralIsoCorrFactor;
	if ( pfPhotonIso_        ) sumPt += pfPhotonIso_->compSumPt(pfPhotons, coneAxis)*pfNeutralIsoCorrFactor;
      } else if ( methodPUcorr_ == kDeltaBeta ) {
	double pfNoPileUpChargedParticleIsoSumPt = pfChargedParticleIso_->compSumPt(pfNoPileUpChargedParticles, coneAxis);
	sumPt = pfNoPileUpChargedParticleIsoSumPt;

	//std::cout << "pfNoPileUpChargedParticles = " << pfNoPileUpChargedParticleIsoSumPt << std::endl;

	double sumPtNeutralIsoSumPt  = 0.;
	double sumPtNeutralIsoPUcorr = 0.;
	
	if ( pfNeutralHadronIso_ ) {
	  sumPtNeutralIsoSumPt += pfNeutralHadronIso_->compSumPt(pfNeutralHadrons, coneAxis);
	  sumPtNeutralIsoPUcorr += pfNeutralHadronIsoPUcorr_->compSumPt(pfPileUpChargedParticles, coneAxis);
	}
	
	if ( pfPhotonIso_ ) {
	  sumPtNeutralIsoSumPt += pfPhotonIso_->compSumPt(pfPhotons, coneAxis);
	  sumPtNeutralIsoPUcorr += pfPhotonIsoPUcorr_->compSumPt(pfPileUpChargedParticles, coneAxis);
	}
	//std::cout << "sumPtNeutralIsoSumPt = " << sumPtNeutralIsoSumPt << std::endl;
	//std::cout << "sumPtNeutralIsoPUcorr = " << sumPtNeutralIsoPUcorr << std::endl;
	sumPt += TMath::Max(sumPtNeutralIsoSumPt - 0.5*chargedToNeutralFactor_*sumPtNeutralIsoPUcorr, 0.); 
      } else if ( methodPUcorr_ == kRho ) {
	if ( rhoFastJetCorrection == -1. ) 
	  throw cms::Exception("ParticlePFIsolationExtractor")
	    << "Pile-up correction Method = 'rho' requires rhoFastJetCorrection !!\n";

	if ( addChargedParticleIso_ ) sumPt += pfChargedParticleIso_->compSumPt(pfChargedParticles, coneAxis);
	if ( addNeutralHadronIso_   ) sumPt += pfNeutralHadronIso_->compSumPt(pfNeutralHadrons, coneAxis);
	if ( addPhotonIso_          ) sumPt += pfPhotonIso_->compSumPt(pfPhotons, coneAxis);
	
	//std::cout << "before rho FastJet correction: sumPt = " << sumPt << std::endl;

	sumPt -= TMath::Pi()*pfChargedParticleIso_->dRisoCone_*pfChargedParticleIso_->dRisoCone_*(rhoFastJetCorrection - ueRhoOffset_);
	if ( sumPt < 0. ) sumPt = 0.;

	//std::cout << "after rho FastJet correction: sumPt = " << sumPt << std::endl;
      }
      else assert(0);
    }

    return sumPt;
  }

 private:

  struct pfIsoConfigType
  {
    pfIsoConfigType(reco::PFCandidate::ParticleType pfParticleType, const edm::ParameterSet& cfg)
      : pfParticleType_(pfParticleType),
	ptMin_(cfg.getParameter<double>("ptMin")),
	dRvetoCone_(cfg.getParameter<double>("dRvetoCone")),
	dRisoCone_(cfg.getParameter<double>("dRisoCone"))
    {
      dEtaVeto_ = cfg.exists("dEtaVeto") ? 
	cfg.getParameter<double>("dEtaVeto") : -1.;
      dPhiVeto_ = cfg.exists("dPhiVeto") ? 
	cfg.getParameter<double>("dPhiVeto") : -1.;

      vetoNumHighestPtObjects_ = cfg.exists("vetoNumHighestPtObjects") ? 
	cfg.getParameter<unsigned>("vetoNumHighestPtObjects") : 0;
    }
    ~pfIsoConfigType() {}

    bool passesVeto(const reco::PFCandidate& pfCandidate, const reco::Particle::Vector& isoParticleCandidateDirection)
    {
      if ( pfCandidate.particleId() != pfParticleType_ ) return false;

      if ( TMath::IsNaN(pfCandidate.pt()) || pfCandidate.pt() < ptMin_ ) return false;

      double dR = deltaR(pfCandidate.p4(), isoParticleCandidateDirection);
      if ( dR < dRvetoCone_ || dR > dRisoCone_ ) return false;

      if ( TMath::Abs(pfCandidate.eta() - isoParticleCandidateDirection.eta()) < dEtaVeto_ ) return false;
      if ( TMath::Abs(pfCandidate.phi() - isoParticleCandidateDirection.phi()) < dPhiVeto_ ) return false;

      //std::cout << "<ParticlePFIsolationExtractor::passesVeto>:" << std::endl;
      //std::cout << " veto passed: Pt = " << pfCandidate.pt() << "," 
      //	  << " eta = " << pfCandidate.eta() << ", phi = " << pfCandidate.phi() << ", dR = " << dR << std::endl;

      return true;
    }

    double compSumPt(const std::vector<const reco::PFCandidate*>& pfCandidates, 
		     const reco::Particle::Vector& isoParticleCandidateDirection)
    {
      double sumPt = 0.;
      
      if ( vetoNumHighestPtObjects_ > 0 ) {
	std::vector<double> pfCandidatePt;
	
	for ( std::vector<const reco::PFCandidate*>::const_iterator pfCandidate = pfCandidates.begin();
	      pfCandidate != pfCandidates.end(); ++pfCandidate ) {
	  if ( !passesVeto(**pfCandidate, isoParticleCandidateDirection) ) continue;
	  
	  pfCandidatePt.push_back((*pfCandidate)->pt());
	}
	
	// sort transverse momenta of particle-flow candidates
	// ( lowest/highest Pt value will be stored in pfCandidatePt[0]/pfCandidatePt[numPFCandidates - 1];
	//  cf. http://www.cplusplus.com/reference/algorithm/sort/ )
	std::sort(pfCandidatePt.begin(), pfCandidatePt.end());
	
	int numPt = pfCandidatePt.size();
	for ( int iPt = 0; iPt < (numPt - (int)vetoNumHighestPtObjects_); ++iPt ) {
	  sumPt += pfCandidatePt[iPt];
	}
      } else {
	for ( std::vector<const reco::PFCandidate*>::const_iterator pfCandidate = pfCandidates.begin();
	      pfCandidate != pfCandidates.end(); ++pfCandidate ) {
	  if ( !passesVeto(**pfCandidate, isoParticleCandidateDirection) ) continue;
	  
	  sumPt += (*pfCandidate)->pt();
	}   
      } 
      
      return sumPt;
    }

    reco::PFCandidate::ParticleType pfParticleType_;

    double ptMin_;

    double dRvetoCone_;
    double dRisoCone_;

    double dEtaVeto_;
    double dPhiVeto_;

    unsigned vetoNumHighestPtObjects_;
  };

  pfIsoConfigType* pfChargedParticleIso_;
  bool addChargedParticleIso_;
  double pfChargedParticleIsoConeSize_;
  pfIsoConfigType* pfNeutralHadronIso_;
  bool addNeutralHadronIso_;
  double pfNeutralHadronIsoConeSize_;
  pfIsoConfigType* pfPhotonIso_;
  bool addPhotonIso_;
  double pfPhotonIsoConeSize_;

  int methodPUcorr_;
  double chargedToNeutralFactor_;
  PATLeptonTrackVectorExtractor<T>* trackExtractor_;			  
  pfIsoConfigType* pfNeutralHadronIsoPUcorr_;
  pfIsoConfigType* pfPhotonIsoPUcorr_;  
  
  double ueRhoOffset_;
};

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

typedef ParticlePFIsolationExtractor<pat::Electron> PATElectronPFIsolationExtractor;
typedef ParticlePFIsolationExtractor<pat::Muon> PATMuonPFIsolationExtractor;
typedef ParticlePFIsolationExtractor<pat::Tau> PATTauPFIsolationExtractor;

#endif


