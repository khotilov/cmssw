#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimGeneral/HepPDTESSource/interface/HepPDTESSource.h"


  DEFINE_SEAL_MODULE();
  DEFINE_ANOTHER_FWK_EVENTSETUP_SOURCE( HepPDTESSource );
