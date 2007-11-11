#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoRomanPot/RecoFP420/interface/ClusterizerFP420.h"
#include "RecoRomanPot/RecoFP420/interface/TrackerizerFP420.h"

using cms::ClusterizerFP420;
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(ClusterizerFP420);

using cms::TrackerizerFP420;
DEFINE_ANOTHER_FWK_MODULE(TrackerizerFP420);

