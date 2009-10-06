// define the geometry plugin modules for new layers

#include "SLHCUpgradeSimulations/Geometry/interface/DDPixBarStackLinear.h"
#include "SLHCUpgradeSimulations/Geometry/interface/DDPixBarStackLinearGap.h"
#include "SLHCUpgradeSimulations/Geometry/interface/DDPixBarStackLayerAlgo.h"
#include "SLHCUpgradeSimulations/Geometry/interface/DDPixBarStackTrigLayerAlgo.h"
#include "DetectorDescription/Algorithm/interface/DDAlgorithmFactory.h"
#include "SLHCUpgradeSimulations/Geometry/interface/DDPixFwdInnerDisks.h"
#include "SLHCUpgradeSimulations/Geometry/interface/DDPixFwdOuterDisks.h"

#include "SLHCUpgradeSimulations/Geometry/interface/PrintGeomMatInfo.h"
#include "SimG4Core/Watcher/interface/SimWatcherFactory.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#
DEFINE_EDM_PLUGIN (DDAlgorithmFactory, DDPixBarStackLinear,      "track:DDPixBarStackLinear");
DEFINE_EDM_PLUGIN (DDAlgorithmFactory, DDPixBarStackLinearGap,   "track:DDPixBarStackLinearGap");
DEFINE_EDM_PLUGIN (DDAlgorithmFactory, DDPixBarStackLayerAlgo,   "track:DDPixBarStackLayerAlgo");
DEFINE_EDM_PLUGIN (DDAlgorithmFactory, DDPixBarStackTrigLayerAlgo, "track:DDPixBarStackTrigLayerAlgo");
DEFINE_EDM_PLUGIN (DDAlgorithmFactory, DDPixFwdInnerDisks,      "track:DDPixFwdInnerDisks");
DEFINE_EDM_PLUGIN (DDAlgorithmFactory, DDPixFwdOuterDisks,      "track:DDPixFwdOuterDisks");
DEFINE_SEAL_MODULE ();
DEFINE_SIMWATCHER(PrintGeomMatInfo);
