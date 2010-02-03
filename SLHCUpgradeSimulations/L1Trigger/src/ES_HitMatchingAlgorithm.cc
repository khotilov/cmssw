
/*********************************/
/*********************************/
/**                             **/
/** Stacked Tracker Simulations **/
/**        Andrew W. Rose       **/
/**             2008            **/
/**                             **/
/*********************************/
/*********************************/

#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
#include "SLHCUpgradeSimulations/L1Trigger/interface/HitMatchingAlgorithm.h"
#include "FWCore/Framework/interface/eventsetupdata_registration_macro.h" 

EVENTSETUP_DATA_REG(cmsUpgrades::HitMatchingAlgorithm<cmsUpgrades::Ref_PSimHit_>);
EVENTSETUP_DATA_REG(cmsUpgrades::HitMatchingAlgorithm<cmsUpgrades::Ref_PixelDigi_>);
EVENTSETUP_DATA_REG(cmsUpgrades::HitMatchingAlgorithm<cmsUpgrades::Ref_TTHit_>);

