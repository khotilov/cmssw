import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# produce data-formats providing information 
# about distribution of energy deposits in the event 
# with respect to direction of missing Et vector
#--------------------------------------------------------------------------------

metTopologies = cms.EDProducer("MEtTopologyProducer",
    srcMET = cms.InputTag('met'),
    srcEnergyDeposits = cms.InputTag('towerMaker'),                           
    globalThreshold = cms.double(0.5),
    verbosity = cms.untracked.int32(0)
)

selectedMEtTopology04 = cms.EDFilter("MEtTopologySelector",
    src = cms.InputTag('metTopologies'),
    cut = cms.string('Vanti()/Vparallel() < 0.4'),
    filter = cms.bool(False)
)

selectedMEtTopology025 = cms.EDFilter("MEtTopologySelector",
    src = cms.InputTag('metTopologies'),
    cut = cms.string('Vanti()/Vparallel() < 0.25'),
    filter = cms.bool(False)
)


produceMEtTopology = cms.Sequence( metTopologies
                                   *selectedMEtTopology04
                                   *selectedMEtTopology025)
