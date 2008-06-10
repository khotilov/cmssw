#ifndef RecoECAL_ECALClusters_CosmicClusterAlgo_h
#define RecoECAL_ECALClusters_CosmicClusterAlgo_h
 
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"

#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "RecoCaloTools/Navigation/interface/EcalBarrelNavigator.h"
#include "RecoCaloTools/Navigation/interface/EcalEndcapNavigator.h"
#include "Geometry/CaloTopology/interface/EcalBarrelHardcodedTopology.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalEtaPhiRegion.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalRecHitLess.h"

#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"

// C/C++ headers
#include <string>
#include <vector>
#include <set>

typedef std::map<DetId, EcalRecHit> RecHitsMap;

class CosmicClusterAlgo 
{
 public:
  
  enum EcalPart { barrel = 0, endcap = 1 };
  enum VerbosityLevel { pDEBUG = 0, pWARNING = 1, pINFO = 2, pERROR = 3 }; 

  CosmicClusterAlgo() {
  }

  CosmicClusterAlgo(double ebst, double ebSt , double ebDt, double ecst, double ecSt, double ecDt, const PositionCalc& posCalc, VerbosityLevel the_verbosity = pERROR) : 
    ecalBarrelSeedThreshold(ebst), ecalBarrelSingleThreshold(ebSt), ecalBarrelSecondThreshold(ebDt), ecalEndcapSeedThreshold(ecst), ecalEndcapSingleThreshold(ecSt), ecalEndcapSecondThreshold(ecDt), verbosity(the_verbosity) {
    posCalculator_ = posCalc;
  }

  virtual ~CosmicClusterAlgo()
    {
    }

  void setVerbosity(VerbosityLevel the_verbosity)
    {
      verbosity = the_verbosity;
    }

  // this is the method that will start the clusterisation
  std::vector<reco::BasicCluster> makeClusters(const EcalRecHitCollection* hits,
                                               const CaloSubdetectorGeometry *geometry,
                                               const CaloSubdetectorTopology *topology_p,
                                               const CaloSubdetectorGeometry *geometryES_p,
                                               EcalPart ecalPart,
											   const std::vector<int>& masked = std::vector<int>(),
											   const EcalIntercalibConstantMap& icalMap = 1.0,
					       bool regional = false,
					       const std::vector<EcalEtaPhiRegion>& regions = std::vector<EcalEtaPhiRegion>());

  /// point in the space
  typedef math::XYZPoint Point;

 private: 

  //algo to compute position of clusters
  PositionCalc posCalculator_;


  // Energy required for a seed:
  double ecalBarrelSeedThreshold;
  double ecalBarrelSingleThreshold;
  double ecalBarrelSecondThreshold;
  double ecalEndcapSeedThreshold;
  double ecalEndcapSingleThreshold;
  double ecalEndcapSecondThreshold;
  
  // collection of all rechits
  const EcalRecHitCollection *recHits_;

  // The vector of seeds:
  std::vector<EcalRecHit> seeds;

  // The set of used DetID's
  std::set<DetId> used_s;
  std::set<DetId> canSeed_s; // set of crystals not to be added but which can seed
                                // a new 3x3 (e.g. the outer crystals in a 5x5)


  // The vector of DetId's in the cluster currently reconstructed
  std::vector<DetId> current_v9;
  std::vector<DetId> current_v25;

  // The vector of clusters
  std::vector<reco::BasicCluster> clusters_v;

  // The verbosity level
  VerbosityLevel verbosity;

  void mainSearch(const EcalRecHitCollection* hits,
                  const CaloSubdetectorGeometry *geometry_p,
                  const CaloSubdetectorTopology *topology_p,
		  const CaloSubdetectorGeometry *geometryES_p,
                  EcalPart ecalPart,
				  const EcalIntercalibConstantMap& icalMap);

  // Is the crystal at the navigator position a 
  // local maxiumum in energy?
  bool checkMaxima(CaloNavigator<DetId> &navigator,
		   const EcalRecHitCollection *hits);

  // prepare the 5x5 taking care over which crystals
  // are allowed to seed new clusters and which are not
  // after the preparation is complete
  void prepareCluster(CaloNavigator<DetId> &navigator,
                const EcalRecHitCollection *hits,
                const CaloSubdetectorGeometry *geometry);

  // Add the crystal with DetId det to the current
  // vector of crystals if it meets certain criteria
  void addCrystal(const DetId &det, const bool in9);
  

  
  // take the crystals in the current_v and build 
  // them into a BasicCluster
  void makeCluster(const EcalRecHitCollection* hits,const CaloSubdetectorGeometry *geometry_p,const CaloSubdetectorGeometry *geometryES_p, const EcalIntercalibConstantMap& icalMap);

  //TEMP JHAUPT 4-27
  std::vector<int> maskedChannels_; //TEMP JHAUPT 4-27
 };

#endif
