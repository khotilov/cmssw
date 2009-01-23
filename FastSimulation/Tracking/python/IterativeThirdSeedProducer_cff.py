import FWCore.ParameterSet.Config as cms

import FastSimulation.Tracking.TrajectorySeedProducer_cfi
iterativeThirdSeeds = FastSimulation.Tracking.TrajectorySeedProducer_cfi.trajectorySeedProducer.clone()
iterativeThirdSeeds.firstHitSubDetectorNumber = [2]
##iterativeThirdSeeds.firstHitSubDetectors = [1, 2, 6]
iterativeThirdSeeds.firstHitSubDetectors = [1, 2]
iterativeThirdSeeds.secondHitSubDetectorNumber = [3]
iterativeThirdSeeds.secondHitSubDetectors = [1, 2, 6]
iterativeThirdSeeds.thirdHitSubDetectorNumber = [0]
iterativeThirdSeeds.thirdHitSubDetectors = []
iterativeThirdSeeds.seedingAlgo = ['ThirdMixedPairs']
iterativeThirdSeeds.minRecHits = [3]
iterativeThirdSeeds.pTMin = [0.2]
iterativeThirdSeeds.maxD0 = [1.]
iterativeThirdSeeds.maxZ0 = [30.]
iterativeThirdSeeds.numberOfHits = [2]
iterativeThirdSeeds.originRadius = [0.7]
iterativeThirdSeeds.originHalfLength = [7.00]
iterativeThirdSeeds.originpTMin = [0.6]
iterativeThirdSeeds.zVertexConstraint = [-1.0]
iterativeThirdSeeds.primaryVertices = ['none']

