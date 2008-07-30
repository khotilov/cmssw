// $Id: SealModule.cc,v 1.5 2008/07/23 19:43:43 bachtis Exp $

#include "DQMOffline/Trigger/interface/FourVectorHLTOffline.h"
#include "DQMOffline/Trigger/interface/EgammaHLTOffline.h"
#include "DQMOffline/Trigger/interface/EgHLTOfflineClient.h"
#include "DQMOffline/Trigger/interface/HLTTauRefProducer.h"
#include "DQMOffline/Trigger/interface/HLTTauCaloDQMOfflineSource.h"
#include "DQMOffline/Trigger/interface/HLTTauTrkDQMOfflineSource.h"
#include "DQMOffline/Trigger/interface/HLTTauElDQMOfflineSource.h"
#include "DQMOffline/Trigger/interface/HLTTauL1DQMOfflineSource.h"
#include "DQMOffline/Trigger/interface/HLTTauDQMOfflineSource.h"

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(FourVectorHLTOffline);
DEFINE_ANOTHER_FWK_MODULE(EgammaHLTOffline);
DEFINE_ANOTHER_FWK_MODULE(EgHLTOfflineClient);

//Tau Specific
DEFINE_ANOTHER_FWK_MODULE(HLTTauRefProducer);
DEFINE_ANOTHER_FWK_MODULE(HLTTauDQMOfflineSource);
DEFINE_ANOTHER_FWK_MODULE(HLTTauCaloDQMOfflineSource);
DEFINE_ANOTHER_FWK_MODULE(HLTTauTrkDQMOfflineSource);
DEFINE_ANOTHER_FWK_MODULE(HLTTauElDQMOfflineSource);
DEFINE_ANOTHER_FWK_MODULE(HLTTauL1DQMOfflineSource);




