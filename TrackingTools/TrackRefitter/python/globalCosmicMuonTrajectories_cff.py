import FWCore.ParameterSet.Config as cms

from TrackingTools.TrackRefitter.TracksToTrajectories_cff import *
# IMPORTANT !!! ##
# If you want to revert the fit direction, then
# Case 1 #
# string RefitDirection = "alongMomentum"
# KFTrajectoryFitterESProducer   ---> Fitter = "KFFitterForRefitInsideOut"
# KFTrajectorySmootherESProducer ---> Smoother = "KFSmootherForRefitInsideOut"
# Case 2 #
# string RefitDirection = "oppositeToMomentum"
# KFTrajectoryFitterESProducer   ---> Fitter = "KFFitterForRefitOutsideIn"
# KFTrajectorySmootherESProducer ---> Smoother = "KFSmootherForRefitOutsideIn"
# the propagator must be the same as the one used by the Fitter
#
globalCosmicMuons = cms.EDProducer("TracksToTrajectories",
                                   Tracks = cms.InputTag("globalCosmicMuons"),
                                   TrackTransformer = cms.PSet(DoPredictionsOnly = cms.bool(False),
                                                               Fitter = cms.string('KFFitterForRefitOutsideIn'),
                                                               #        TrackerRecHitBuilder = cms.string('WithTrackAngleAndTemplate'),
                                                               TrackerRecHitBuilder = cms.string('WithTrackAngle'),
                                                               Smoother = cms.string('KFSmootherForRefitOutsideIn'),
                                                               MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
                                                               RefitDirection = cms.string('oppositeToMomentum'),
                                                               RefitRPCHits = cms.bool(True),
                                                               Propagator = cms.string('SmartPropagatorAnyRK')
                                                               )
                                   )



