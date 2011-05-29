#ifndef AnalysisDataFormats_TauAnalysis_NSVfitSingleParticleHypothesisBase_h
#define AnalysisDataFormats_TauAnalysis_NSVfitSingleParticleHypothesisBase_h

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Ptr.h"

#include <string>

class NSVfitSingleParticleHypothesisBase
{
 public:

  NSVfitSingleParticleHypothesisBase() {}
  NSVfitSingleParticleHypothesisBase(const std::string& name, int barcode)
    : name_(name),
      barcode_(barcode)
  {}
  NSVfitSingleParticleHypothesisBase(const edm::Ptr<reco::Candidate>& particle, const std::string& name, int barcode)
    : name_(name),
      barcode_(barcode),
      particle_(particle)
  {}
  NSVfitSingleParticleHypothesisBase(const NSVfitSingleParticleHypothesisBase&);
  virtual ~NSVfitSingleParticleHypothesisBase() {}

  virtual NSVfitSingleParticleHypothesisBase* clone() const { return new NSVfitSingleParticleHypothesisBase(*this); }

  virtual NSVfitSingleParticleHypothesisBase& operator=(const NSVfitSingleParticleHypothesisBase&);

  const std::string& name() const { return name_; }
  int barcode() const { return barcode_; }

  /// pointer to reco::Candidate from which this hypothesis was made
  virtual const edm::Ptr<reco::Candidate>& particle() const { return particle_; }

  virtual void print(std::ostream& stream) const
  {
    stream << "<NSVfitSingleParticleHypothesisBase::print>:" << std::endl;
    stream << " name = " << name_ << std::endl;
    stream << " barcode = " << barcode_ << std::endl;
  }

 protected:

  /// "human-readable" label for print-out
  std::string name_;

  /// unique identifier to associate hypothesis to likelihood function
  int barcode_;

  /// pointer to reco::Candidate from which this hypothesis was made
  edm::Ptr<reco::Candidate> particle_;
};

// CV: implementation of operator< needed to store NSVfitSingleParticleHypothesisBase objects in edm::OwnVector
bool operator<(const NSVfitSingleParticleHypothesisBase&, const NSVfitSingleParticleHypothesisBase&);

#endif
