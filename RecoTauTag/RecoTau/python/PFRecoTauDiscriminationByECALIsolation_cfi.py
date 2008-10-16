import FWCore.ParameterSet.Config as cms

pfRecoTauDiscriminationByECALIsolation = cms.EDFilter("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag('pfRecoTauProducer'),
    ManipulateTracks_insteadofChargedHadrCands = cms.bool(False),
    # following parameters are considered when ManipulateTracks_insteadofChargedHadrCands paremeter is set true
    # *BEGIN*
    TrackerIsolAnnulus_Tracksmaxn = cms.int32(0),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    TrackerIsolAnnulus_Candsmaxn = cms.int32(0),
    ECALIsolAnnulus_Candsmaxn = cms.int32(0)
)


