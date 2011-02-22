import FWCore.ParameterSet.Config as cms

apd_sim_parameters = cms.PSet(
    apdAddToBarrel  = cms.bool(False),
    apdSeparateDigi = cms.bool(True),
    apdSimToPELow   = cms.double(2.45e6),
    apdSimToPEHigh  = cms.double(88.2e6),
    apdTimeOffset   = cms.double(-12.0),
    apdTimeOffWidth = cms.double(0.8),
    apdDoPEStats    = cms.bool(True),
    apdDigiTag      = cms.string("APD"),
    apdShapeTstart  = cms.double( 74.5 ),
    apdShapeTau     = cms.double( 40.5 ),
    apdNonlParms    = cms.vdouble( 1.64e-3, -3.72e-2, 0.207, 0.812, 13.0, 0.85)
)

