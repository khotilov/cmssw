#include "RecoParticleFlow/PFClusterTools/interface/TreeUtility.h"
#include "RecoParticleFlow/PFClusterTools/interface/Calibrator.h"  
#include "RecoParticleFlow/PFClusterTools/interface/DetectorElement.h"      
#include "RecoParticleFlow/PFClusterTools/interface/LinearCalibrator.h"    
#include "RecoParticleFlow/PFClusterTools/interface/Operators.h"         
#include "RecoParticleFlow/PFClusterTools/interface/SpaceManager.h"  
#include "RecoParticleFlow/PFClusterTools/interface/ToString.h"
#include "RecoParticleFlow/PFClusterTools/interface/Deposition.h"  
#include "RecoParticleFlow/PFClusterTools/interface/DetectorElementType.h"  
#include "RecoParticleFlow/PFClusterTools/interface/PFToolsException.h"  
#include "RecoParticleFlow/PFClusterTools/interface/ParticleDeposit.h"  
#include "RecoParticleFlow/PFClusterTools/interface/SpaceVoxel.h"    
#include "RecoParticleFlow/PFClusterTools/interface/TreeUtility.h"
#include "RecoParticleFlow/PFClusterTools/interface/SingleParticleWrapper.h"
//#include "RecoParticleFlow/PFClusterTools/interface/Exercises.h"
#include "RecoParticleFlow/PFClusterTools/interface/Exercises2.h"
#include "RecoParticleFlow/PFClusterTools/interface/CalibrationResultWrapper.h"
#include "RecoParticleFlow/PFClusterTools/interface/CalibrationProvenance.h"
#include "RecoParticleFlow/PFClusterTools/interface/CalibrationTarget.h"
#include "RecoParticleFlow/PFClusterTools/interface/Calibratable.h"
#include "RecoParticleFlow/PFClusterTools/interface/Calibration.h"
#include "RecoParticleFlow/PFClusterTools/interface/LinearCalibration.h"
#include "RecoParticleFlow/PFClusterTools/interface/IO.h"
#include "RecoParticleFlow/PFClusterTools/interface/Utils.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFClusterCalibration.h"

namespace { 
  namespace {
	pftools::SingleParticleWrapper spw;
	pftools::CalibrationResultWrapper crw;
	pftools::Calibratable c;
	pftools::CalibratableElement ce;
	pftools::CandidateWrapper cw;
	std::vector<pftools::CalibratableElement> svce;
	std::vector<pftools::CandidateWrapper> svcw;
	pftools::Calibration calib;
	pftools::LinearCalibration linCalib;
	
  }
}
