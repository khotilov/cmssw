#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "GeneratorInterface/GenFilters/interface/PythiaFilter.h"
#include "GeneratorInterface/GenFilters/interface/PythiaDauFilter.h"
#include "GeneratorInterface/GenFilters/interface/PythiaFilterGammaJet.h"
#include "GeneratorInterface/GenFilters/interface/PythiaFilterGammaGamma.h"
#include "GeneratorInterface/GenFilters/interface/PythiaFilterZJet.h"
#include "GeneratorInterface/GenFilters/interface/PythiaHLTSoupFilter.h"
#include "GeneratorInterface/GenFilters/interface/BsJpsiPhiFilter.h"
#include "GeneratorInterface/GenFilters/interface/JetFlavourFilter.h"
#include "GeneratorInterface/GenFilters/interface/PythiaFilterGammaJetWithBg.h"
#include "GeneratorInterface/GenFilters/interface/PythiaFilterGammaJetWithOutBg.h"
#include "GeneratorInterface/GenFilters/interface/PythiaFilterZJetWithOutBg.h"
#include "GeneratorInterface/GenFilters/interface/MCDijetResonance.h"
#include "GeneratorInterface/GenFilters/interface/MCProcessFilter.h"
#include "GeneratorInterface/GenFilters/interface/MCProcessRangeFilter.h"
#include "GeneratorInterface/GenFilters/interface/MCSingleParticleFilter.h"
#include "GeneratorInterface/GenFilters/interface/MCSmartSingleParticleFilter.h"
#include "GeneratorInterface/GenFilters/interface/MCZll.h"
#include "GeneratorInterface/GenFilters/interface/MinimumBiasFilter.h"
#include "GeneratorInterface/GenFilters/interface/RecoDiMuon.h"
#include "GeneratorInterface/GenFilters/interface/MCLongLivedParticles.h"
#include "GeneratorInterface/GenFilters/interface/MCParticlePairFilter.h"
#include "GeneratorInterface/GenFilters/interface/CosmicGenFilterHelix.h"
#include "GeneratorInterface/GenFilters/interface/CosmicGenFilterLowE.h"
#include "GeneratorInterface/GenFilters/interface/BHFilter.h"
#include "GeneratorInterface/GenFilters/interface/PythiaFilterIsolatedTrack.h"
#include "GeneratorInterface/GenFilters/interface/BCToEFilter.h"
#include "GeneratorInterface/GenFilters/interface/EMEnrichingFilter.h"
#include "GeneratorInterface/GenFilters/interface/doubleEMEnrichingFilter.h"
#include "GeneratorInterface/GenFilters/interface/MCDecayingPionKaonFilter.h"
#include "GeneratorInterface/GenFilters/interface/PythiaFilterEMJetHeep.h"
#include "GeneratorInterface/GenFilters/interface/ComphepSingletopFilter.h"
#include "GeneratorInterface/GenFilters/interface/STFilter.h"


DEFINE_SEAL_MODULE();
  DEFINE_ANOTHER_FWK_MODULE(PythiaFilter);
  DEFINE_ANOTHER_FWK_MODULE(PythiaDauFilter);
  DEFINE_ANOTHER_FWK_MODULE(PythiaFilterGammaJet);
  DEFINE_ANOTHER_FWK_MODULE(PythiaFilterGammaGamma);
  DEFINE_ANOTHER_FWK_MODULE(PythiaFilterZJet);
  DEFINE_ANOTHER_FWK_MODULE(PythiaHLTSoupFilter);
  DEFINE_ANOTHER_FWK_MODULE(BsJpsiPhiFilter);
  DEFINE_ANOTHER_FWK_MODULE(JetFlavourFilter);
  DEFINE_ANOTHER_FWK_MODULE(PythiaFilterGammaJetWithBg);
  DEFINE_ANOTHER_FWK_MODULE(PythiaFilterGammaJetWithOutBg);
  DEFINE_ANOTHER_FWK_MODULE(PythiaFilterZJetWithOutBg);
  DEFINE_ANOTHER_FWK_MODULE(MCDijetResonance);
  DEFINE_ANOTHER_FWK_MODULE(MCProcessFilter);
  DEFINE_ANOTHER_FWK_MODULE(MCProcessRangeFilter);
  DEFINE_ANOTHER_FWK_MODULE(MCSingleParticleFilter);
  DEFINE_ANOTHER_FWK_MODULE(MCSmartSingleParticleFilter);
  DEFINE_ANOTHER_FWK_MODULE(MCZll);
  DEFINE_ANOTHER_FWK_MODULE(MinimumBiasFilter);
  DEFINE_ANOTHER_FWK_MODULE(RecoDiMuon);
  DEFINE_ANOTHER_FWK_MODULE(MCLongLivedParticles);
  DEFINE_ANOTHER_FWK_MODULE(MCParticlePairFilter);
  DEFINE_ANOTHER_FWK_MODULE(CosmicGenFilterHelix);
  DEFINE_ANOTHER_FWK_MODULE(CosmicGenFilterLowE);
  DEFINE_ANOTHER_FWK_MODULE(BHFilter);
  DEFINE_ANOTHER_FWK_MODULE(PythiaFilterIsolatedTrack);
  DEFINE_ANOTHER_FWK_MODULE(BCToEFilter);
  DEFINE_ANOTHER_FWK_MODULE(EMEnrichingFilter);
  DEFINE_ANOTHER_FWK_MODULE(doubleEMEnrichingFilter);
  DEFINE_ANOTHER_FWK_MODULE(MCDecayingPionKaonFilter);
  DEFINE_ANOTHER_FWK_MODULE(PythiaFilterEMJetHeep);
  DEFINE_ANOTHER_FWK_MODULE(ComphepSingletopFilter);
  DEFINE_ANOTHER_FWK_MODULE(STFilter);


