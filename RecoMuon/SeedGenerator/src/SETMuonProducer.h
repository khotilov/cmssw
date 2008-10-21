#ifndef RecoMuon_SeedGenerator_SETMuonProducer_H
#define RecoMuon_SeedGenerator_SETMuonProducer_H

/** \class SETMuonProducer 
  */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "RecoMuon/TrackingTools/interface/RecoMuonEnumerators.h"


#include <RecoMuon/TrackingTools/interface/MuonServiceProxy.h>

#include "RecoMuon/SeedGenerator/src/SETFilter.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h" 
#include "RecoMuon/SeedGenerator/src/PtExtractor.h"

#include "FWCore/Framework/interface/Event.h"
#include "RecoMuon/TrackingTools/interface/MuonCandidate.h"

typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
typedef std::vector<Trajectory*> TrajectoryContainer;

class TrajectorySeed;
class STAFilter;
class MuonServiceProxy;

//namespace edm {class ParameterSet;}
namespace edm {class ParameterSet; class Event; class EventSetup;}

class SETMuonProducer : public edm::EDProducer {
  
 public:
  /// Constructor with Parameter set 
  // SETMuonProducer (const edm::ParameterSet&, const MuonServiceProxy*);
  SETMuonProducer (const edm::ParameterSet&);
  
  /// Destructor
  virtual ~SETMuonProducer();
  
  // Operations
  virtual void produce(edm::Event&, const edm::EventSetup&);
  
 protected:

 private:
  
  
  // Returns a vector of measurements sets (for later trajectory seed building)
  std::vector < std::pair < TrajectoryStateOnSurface, 
    TransientTrackingRecHit::ConstRecHitContainer > > trajectories(const edm::Event&);

  /// pre-filter
  SETFilter* filter() const {return theFilter;}
  
  /// pT extractor (given two hits)
  PtExtractor * pt_extractor() const {return thePtExtractor;}

  //---- SET 
  /// Build local clusters of segments that are clearly separated from each other in the eta-phi plane 
  std::vector< MuonRecHitContainer > clusterHits( 
						 MuonRecHitContainer muonRecHits,
						 MuonRecHitContainer muonRecHits_DT2D_hasPhi,
						 MuonRecHitContainer muonRecHits_DT2D_hasZed,
						 MuonRecHitContainer muonRecHits_RPC);
  
  std::vector <MuonRecHitContainer> findAllValidSets(std::vector< MuonRecHitContainer > MuonRecHitContainer_perLayer);

  void validSetsPrePruning(std::vector <MuonRecHitContainer>  & allValidSets); 

  std::pair <int, int> checkAngleDeviation(double dPhi_1, double dPhi_2);

  std::vector <seedSet>  fillSeedSets(std::vector <MuonRecHitContainer> & allValidSets);
  //----

  //private:
  
  SETFilter* theFilter;
  void setEvent(const edm::Event&);
 
  //---- SET
  PtExtractor* thePtExtractor;
  bool apply_prePruning;
  bool useSegmentsInTrajectory;
  bool useRPCs;

  edm::InputTag DTRecSegmentLabel;
  edm::InputTag CSCRecSegmentLabel;
  edm::InputTag RPCRecSegmentLabel;

  MuonServiceProxy *theService;
};
#endif
