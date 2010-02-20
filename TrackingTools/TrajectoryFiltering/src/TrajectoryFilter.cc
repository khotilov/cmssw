#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"

#include "FWCore/Utilities/interface/typelookup.h"
TYPELOOKUP_DATA_REG(TrajectoryFilter);

TrajectoryFilter::~TrajectoryFilter() {}
