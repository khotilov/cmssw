import FWCore.ParameterSet.Config as cms

ErsatzMEtParams = cms.PSet(
MCTruthCollection = cms.InputTag("genParticles"),
ElectronCollection = cms.InputTag("gsfElectrons"),
HybridScCollection = cms.InputTag("correctedHybridSuperClusters"),
M5x5ScCollection = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
sigmaElectronicNoise_EB = cms.double(0.03),
hyb_fCorrPSet = cms.PSet(
        brLinearLowThr = cms.double(1.1),
        fBremVec = cms.vdouble(-0.05208, 0.1331, 0.9196, -0.0005735, 1.343),
        brLinearHighThr = cms.double(8.0),
        fEtEtaVec = cms.vdouble(1.0012, -0.5714, 0, 0,
                                0, 0.5549, 12.74, 1.0448,
                                0, 0, 0, 0,
                                8.0, 1.023, -0.00181, 0, 0)
    ),
sigmaElectronicNoise_EE = cms.double(0.15),
m5x5_fCorrPSet = cms.PSet(
        brLinearLowThr = cms.double(0.6),
        fBremVec = cms.vdouble(-0.04163, 0.08552, 0.95048, -0.002308, 1.077),
        brLinearHighThr = cms.double(6.0),
        fEtEtaVec = cms.vdouble(0.9746, -6.512, 0, 0,
                                0.02771, 4.983, 0, 0,
                                -0.007288, -0.9446, 0, 0,
                                0, 0, 0, 1, 1)
),
CElecPtMin = cms.double(1.),
CEB_sigmaIEtaIEta = cms.double(99999.),
CEB_deltaPhiIn = cms.double(99999.),
CEB_deltaEtaIn = cms.double(99999.),
CEB_EcalIso = cms.double(99999.),
CEB_HcalIso = cms.double(99999.),
CEB_TrckIso = cms.double(99999.),
CEE_sigmaIEtaIEta = cms.double(99999.),
CEE_deltaPhiIn = cms.double(99999.),
CEE_deltaEtaIn = cms.double(99999.),
CEE_EcalIso = cms.double(99999.),
CEE_HcalIso = cms.double(99999.),
CEE_TrckIso = cms.double(99999.),
eIsoTrack = cms.InputTag("electronTrackIsolationLcone"),
eIsoEcal = cms.InputTag("electronEcalRecHitIsolationLcone"),
eIsoHcal = cms.InputTag("electronHcalTowerIsolationLcone"),
CaloMEtCollection = cms.InputTag("met"),
GenMEtCollection = cms.InputTag("genMet"),
PfMEtCollection = cms.InputTag("pfMet"),
TcMEtCollection = cms.InputTag("tcMet"),
#T1MEtCollection = cms.InputTag("layer1METs"),
#CaloMEtT1Collection = cms.InputTag("Type1MET"),
#CaloMEtT1Collection = cms.InputTag("corMetType1Icone5"),
#TriggerResults = cms.InputTag("TriggerResults","","HLT8E29"),
#TriggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT8E29"),
#TriggerPath = cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronLWEt15PixelMatchFilter","","HLT8E29"),
#TriggerName = cms.string("HLT_Ele15_LW_L1R"),
TriggerResults = cms.InputTag("TriggerResults","","HLT"),
TriggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
TriggerPath = cms.InputTag("hltL1NonIsoHLTNonIsoSingleElectronEt15TrackIsolFilter", "", "HLT"),
TriggerName = cms.string("HLT_Ele15_SW_LooseTrackIso_L1R"),
mW = cms.double(80.398),
mZ = cms.double(91.188),
HLTPathCheck = cms.bool(False),
mTPmin = cms.double(61.),
mTPmax = cms.double(121.),
BarrelEtaMax = cms.double(1.4442),
EndCapEtaMin = cms.double(1.56),
EndCapEtaMax = cms.double(2.5),
etaWidth = cms.int32(7),
phiWidth = cms.int32(25)
)
