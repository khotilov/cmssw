import FWCore.ParameterSet.Config as cms

# The services
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
# parametrization for initial pT
from RecoMuon.MuonSeedGenerator.ptSeedParameterization_38T_cfi import *
SETMuonSeed  = cms.EDProducer("SETMuonProducer",
    MuonServiceProxy,
    SETTrajBuilderParameters = cms.PSet(
        ptSeedParameterization, 
        Apply_prePruning = cms.bool(True),
        FilterParameters = cms.PSet(
            DTRecSegmentLabel = cms.InputTag("dt4DSegments"),
            CSCRecSegmentLabel = cms.InputTag("cscSegments"),
            RPCRecSegmentLabel = cms.InputTag("rpcRecHits"),
            Propagator = cms.string('SteppingHelixPropagatorAny'),

            EnableRPCMeasurement = cms.bool(True),
# NOT USED for now
            EnableDTMeasurement = cms.bool(True),
            EnableCSCMeasurement = cms.bool(True)       
        )
    )
)



