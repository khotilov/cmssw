#include "TauAnalysis/MCEmbeddingTools/interface/ParticleReplacerBase.h"

ParticleReplacerBase::ParticleReplacerBase(const edm::ParameterSet& iConfig):
  tauMass(1.7769)
{}

ParticleReplacerBase::~ParticleReplacerBase() {}

void ParticleReplacerBase::beginJob() {}
void ParticleReplacerBase::beginRun(edm::Run& iRun, const edm::EventSetup& iSetup) {}
void ParticleReplacerBase::endRun() {}
void ParticleReplacerBase::endJob() {}

