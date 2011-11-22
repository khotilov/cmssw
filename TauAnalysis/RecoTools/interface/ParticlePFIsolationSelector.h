#ifndef TauAnalysis_RecoTools_ParticlePFIsolationSelector_h
#define TauAnalysis_RecoTools_ParticlePFIsolationSelector_h

/** \class ParticlePFIsolationSelector
 *
 * Flexible selection of pat::Leptons based on charged hadrons,
 * neutral hadrons and photons/pi0s reconstructed by particle-flow algorithm;
 * a pat::Electron, pat::Muon or pat::Tau passed the selection
 * if the sum of deposits of the speficied type and with Pt > ptMin
 * in an annulus within dRvetoCone (inner cone) and dRisoCone (outer cone)
 * is more than sumPtMin and less than sumPtMax
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.15 $
 *
 * $Id: ParticlePFIsolationSelector.h,v 1.15 2011/08/31 18:47:59 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TauAnalysis/RecoTools/interface/ParticlePFIsolationExtractor.h"

#include <TMath.h>

#include <vector>
#include <string>
#include <algorithm>

template <class T>
class ParticlePFIsolationSelector
{
 public:
  typedef std::vector<T> collection;
  
  explicit ParticlePFIsolationSelector(const edm::ParameterSet& cfg)
    : extractor_(cfg),     
      direction_(ParticlePFIsolationExtractor<T>::kDirP4),
      sumPtMethod_(kAbsoluteIso),
      cfgError_(0)
  {
    pfCandidateSrc_ = cfg.getParameter<edm::InputTag>("pfCandidateSource");
    pfNoPileUpCandidateSrc_ = cfg.getParameter<edm::InputTag>("pfNoPileUpCandidateSource");

    if ( cfg.exists("rhoFastJetSource") ) rhoFastJetSrc_ = cfg.getParameter<edm::InputTag>("rhoFastJetSource");

    if ( cfg.exists("direction") ) {
      std::string direction_string = cfg.getParameter<std::string>("direction");
      if      ( direction_string == "p4"    ) direction_ = ParticlePFIsolationExtractor<T>::kDirP4;
      else if ( direction_string == "track" ) direction_ = ParticlePFIsolationExtractor<T>::kDirTrack;
      else throw cms::Exception("ParticlePFIsolationSelector")
	<< "Invalid Configuration parameter direction = " << direction_string << "!!\n";
    }

    sumPtMaxEB_ = 1000;
    sumPtMaxEE_ = 1000;
    sumPtMinEB_ = -1;
    sumPtMinEE_ = -1;

    if( cfg.exists("sumPtMinEB") && 
        cfg.exists("sumPtMinEE") &&
        cfg.exists("sumPtMin") ) {
      edm::LogWarning("ParticlePFIsolationSelector")
	    << " Configuration specifies sumPtMinEB AND sumPtMinEE AND sumPtMin " << std::endl
        << "  --> using sumPtMin, IGNORING sumPtMinEB and sumPtMinEE";
    } 
    
    if( cfg.exists("sumPtMaxEB") && 
        cfg.exists("sumPtMaxEE") &&
        cfg.exists("sumPtMax") ) {
      edm::LogWarning("ParticlePFIsolationSelector")
	    << " Configuration specifies sumPtMaxEB AND sumPtMaxEE AND sumPtMax " << std::endl
        << "  --> using sumPtMax, IGNORING sumPtMaxEB and sumPtMaxEE";
    } 
    
    if ( cfg.exists("sumPtMinEB") ) {
      sumPtMinEB_ = cfg.getParameter<double>("sumPtMinEB");
    }
    if ( cfg.exists("sumPtMinEE") ) {
      sumPtMinEE_ = cfg.getParameter<double>("sumPtMinEE");
    }
    if (cfg.exists("sumPtMin") ) {
      sumPtMinEB_ = cfg.getParameter<double>("sumPtMin");
      sumPtMinEE_ = cfg.getParameter<double>("sumPtMin");
    }

    if ( cfg.exists("sumPtMaxEB") ) {
      sumPtMaxEB_ = cfg.getParameter<double>("sumPtMaxEB");
    }
    if ( cfg.exists("sumPtMaxEE") ) {
      sumPtMaxEE_ = cfg.getParameter<double>("sumPtMaxEE");
    }
    if (cfg.exists("sumPtMax") ) {
      sumPtMaxEB_ = cfg.getParameter<double>("sumPtMax");
      sumPtMaxEE_ = cfg.getParameter<double>("sumPtMax");
    }
    
    if ( cfg.exists("sumPtMethod") ) {
      std::string sumPtMethod_string = cfg.getParameter<std::string>("sumPtMethod");
      if ( sumPtMethod_string == "absolute" ) {
	sumPtMethod_ = kAbsoluteIso;
      } else if ( sumPtMethod_string == "relative" ) {
	sumPtMethod_ = kRelativeIso;
      } else {
	edm::LogError("ParticlePFIsolationSelector")
	  << " Configuration parameter 'sumPtMethod' = " << sumPtMethod_string << " invalid !!";
	cfgError_ = 1;
      }
    }
  }
    ~ParticlePFIsolationSelector() {}
    
    typename std::vector<const T*>::const_iterator begin() const { return selected_.begin(); }
    typename std::vector<const T*>::const_iterator end() const { return selected_.end(); }
    
    void select(const edm::Handle<collection>& isoParticleCandidates, const edm::Event& evt, const edm::EventSetup& es)
    {
      selected_.clear();
      
      if ( cfgError_ ) {
	edm::LogError ("select")
	  << " Error in Configuration ParameterSet --> no iso. particle candidate will pass Selection !!";
	return;
      }
      
      edm::Handle<reco::PFCandidateCollection> pfCandidates;
      evt.getByLabel(pfCandidateSrc_, pfCandidates);

      edm::Handle<reco::PFCandidateCollection> pfNoPileUpCandidates;
      evt.getByLabel(pfNoPileUpCandidateSrc_, pfNoPileUpCandidates);

      double rhoFastJetCorrection = 0.;
      if ( rhoFastJetSrc_.label() != "" ) {
	edm::Handle<double> rhoFastJetHandle;
	evt.getByLabel(rhoFastJetSrc_, rhoFastJetHandle);
	rhoFastJetCorrection = (*rhoFastJetHandle);
      }

      for ( typename collection::const_iterator isoParticleCandidate = isoParticleCandidates->begin();
	    isoParticleCandidate != isoParticleCandidates->end(); ++isoParticleCandidate ) {
	double sumPt = extractor_(*isoParticleCandidate, direction_, *pfCandidates, *pfNoPileUpCandidates, rhoFastJetCorrection);  
	
	// need to fix for correct eta
	double sumPtMax = ( isoParticleCandidate->eta() < 1.479 ) ? sumPtMaxEB_ : sumPtMaxEE_;
	double sumPtMin = ( isoParticleCandidate->eta() < 1.479 ) ? sumPtMinEB_ : sumPtMinEE_;

	if ( sumPtMethod_ == kAbsoluteIso ) {
	  if ( sumPtMin > 0. && sumPt  < sumPtMin ) continue;
	  if ( sumPtMax  > 0. && sumPt  > sumPtMax  ) continue;
	} else if ( sumPtMethod_ == kRelativeIso ) {
	  double relIso = ( isoParticleCandidate->pt() > 1. ) ? (sumPt/isoParticleCandidate->pt()) : sumPt;
	  if ( sumPtMin > 0. && relIso < sumPtMin ) continue;
	  if ( sumPtMax  > 0. && relIso > sumPtMax  ) continue;
	}
	
	selected_.push_back(&(*isoParticleCandidate));
      }
    }

    size_t size() const { return selected_.size(); }
    
 private:
    std::vector<const T*> selected_;
    
    edm::InputTag pfCandidateSrc_;
    edm::InputTag pfNoPileUpCandidateSrc_;
    
    edm::InputTag rhoFastJetSrc_;

    ParticlePFIsolationExtractor<T> extractor_;
    int direction_;

    double sumPtMinEB_;
    double sumPtMinEE_;
    double sumPtMaxEB_;
    double sumPtMaxEE_;    
        
    enum { kAbsoluteIso, kRelativeIso };
    int sumPtMethod_;
    
    int cfgError_;
};

#endif


