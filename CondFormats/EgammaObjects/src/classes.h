#include "CondFormats/PhysicsToolsObjects/interface/Histogram.h"
#include "CondFormats/EgammaObjects/interface/ElectronLikelihoodCategoryData.h"
#include "CondFormats/EgammaObjects/interface/ElectronLikelihoodCalibration.h"

namespace {
  struct dictionary {
    ElectronLikelihoodCategoryData a;
 
    ElectronLikelihoodCalibration b;
    ElectronLikelihoodCalibration::Entry c;
    std::vector<ElectronLikelihoodCalibration::Entry> d;
    std::vector<ElectronLikelihoodCalibration::Entry>::iterator d1;
    std::vector<ElectronLikelihoodCalibration::Entry>::const_iterator d2;
  };
}
