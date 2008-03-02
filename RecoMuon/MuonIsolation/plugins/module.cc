#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_SEAL_MODULE();

#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractor.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractorFactory.h"
#include "TrackExtractor.h"
#include "CaloExtractor.h"
#include "CaloExtractorByAssociator.h"
#include "JetExtractor.h"
#include "ExtractorFromDeposits.h"
DEFINE_EDM_PLUGIN(IsoDepositExtractorFactory, muonisolation::TrackExtractor, "TrackExtractor");
DEFINE_EDM_PLUGIN(IsoDepositExtractorFactory, muonisolation::CaloExtractor, "CaloExtractor");
DEFINE_EDM_PLUGIN(IsoDepositExtractorFactory, muonisolation::CaloExtractorByAssociator, "CaloExtractorByAssociator");
DEFINE_EDM_PLUGIN(IsoDepositExtractorFactory, muonisolation::JetExtractor, "JetExtractor");
DEFINE_EDM_PLUGIN(IsoDepositExtractorFactory, muonisolation::ExtractorFromDeposits, "ExtractorFromDeposits");
