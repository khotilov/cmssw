#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "EgammaAnalysis/EgammaEfficiencyProducers/interface/EmObjectProducer.h"
#include "EgammaAnalysis/EgammaEfficiencyProducers/interface/TagProbeProducer.h"
#include "EgammaAnalysis/EgammaEfficiencyProducers/interface/TagProbeVerification.h"
#include "EgammaAnalysis/EgammaEfficiencyProducers/interface/GeorgiosTagProbeProducer.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(EmObjectProducer);
DEFINE_ANOTHER_FWK_MODULE(TagProbeProducer);
DEFINE_ANOTHER_FWK_MODULE(TagProbeVerification);
DEFINE_ANOTHER_FWK_MODULE(GeorgiosTagProbeProducer);

