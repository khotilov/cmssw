#ifndef TauAnalysis_Core_eventDumpAuxFunctions_h
#define TauAnalysis_Core_eventDumpAuxFunctions_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <vector>
#include <string>

#include <iostream>

std::ostream* getOutputOptions(const edm::ParameterSet&, bool&, int&);

void printEventSelectionInfo(const std::vector<std::pair<std::string, bool> >&, const std::vector<std::pair<std::string, bool> >&, std::ostream*);

void printGenParticleInfo(edm::Handle<edm::View<reco::GenParticle> >&, edm::Handle<edm::View<reco::GenJet> >&, std::ostream*);

void printTrackInfo(const edm::RefToBase<reco::Track>&, const reco::Candidate::Point&, bool, bool, std::ostream*);
void printVertexInfo(const reco::Candidate::Point&, std::ostream*);

void printTrackIsolationInfo(const edm::Handle<reco::TrackCollection>&, const reco::Candidate::Vector&, double, double, double, const reco::Vertex::Point&, std::ostream*);
void printPFCandidateIsolationInfo(const edm::Handle<reco::PFCandidateCollection>&, std::string,
				   const reco::Candidate::Vector&, double, double, double, std::ostream*);

#endif
