// define the geometry plugin modules for new layers

#include "SLHCUpgradeSimulations/Geometry/interface/DDPixBarStackLayerAlgo.h"
#include "SLHCUpgradeSimulations/Geometry/interface/DDPixBarStackTrigLayerAlgo.h"
#include "DetectorDescription/Algorithm/interface/DDAlgorithmFactory.h"

DEFINE_EDM_PLUGIN (DDAlgorithmFactory, DDPixBarStackLayerAlgo,   "track:DDPixBarStackLayerAlgo");
DEFINE_EDM_PLUGIN (DDAlgorithmFactory, DDPixBarStackTrigLayerAlgo, "track:DDPixBarStackTrigLayerAlgo");
