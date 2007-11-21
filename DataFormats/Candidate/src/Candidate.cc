// $Id: Candidate.cc,v 1.12 2007/09/21 14:13:05 llista Exp $
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Utilities/interface/EDMException.h"
using namespace reco;

Candidate::~Candidate() { }

const CandidateBaseRef & Candidate::masterClone() const {
  throw cms::Exception("Invalid Reference") 
    << "this Candidate has no master clone reference."
    << "Can't call masterClone() method.\n";
}

bool Candidate::hasMasterClone() const { 
  return false;
}

double Candidate::vertexChi2() const { 
  return 0; 
}

double Candidate::vertexNdof() const { 
  return 0; 
}

double Candidate::vertexNormalizedChi2() const {
  return 0;
}

double Candidate::vertexCovariance(int i, int j) const { 
  throw edm::Exception(edm::errors::UnimplementedFeature) 
    << "reco::Candidate does not implement vertex covariant matrix.\n";
}

void Candidate::fillVertexCovariance(CovarianceMatrix & err) const {
  throw edm::Exception(edm::errors::UnimplementedFeature) 
    << "reco::Candidate does not implement vertex covariant matrix.\n";
}
