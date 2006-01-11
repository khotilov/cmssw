#ifndef Base_PileUp_h
#define Base_PileUp_h

#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/TripleRand.h"

namespace edm {
  class PileUp {
  public:
    explicit PileUp(ParameterSet const& pset);
    ~PileUp() {}

    void readPileUp(std::vector<std::vector<EventPrincipal *> > & result);

  private:
    int const minBunch_;
    int const maxBunch_;
    double const averageNumber_;
    int const intAverage_;
    bool const poisson_;
    long const seed_;
    VectorInputSource * const input_;
    TripleRand eng_;
    RandPoisson distribution_;
  };
}

#endif
