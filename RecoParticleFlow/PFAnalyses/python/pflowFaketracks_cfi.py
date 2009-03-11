import FWCore.ParameterSet.Config as cms

faketracks = cms.EDProducer("PFPretendTrackProducer",
    debug = cms.int32(0),
    runinfo_cuts=cms.string("/afs/cern.ch/user/b/ballin/scratch0/cmssw/src/RecoParticleFlow/PFAnalyses/macros/testbeam_cuts_v1_3X.root"),
    justCreateEmptyCollections = cms.bool(False)
)
