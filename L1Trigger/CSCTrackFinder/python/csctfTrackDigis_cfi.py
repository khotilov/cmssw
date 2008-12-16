import FWCore.ParameterSet.Config as cms

from L1Trigger.CSCCommonTrigger.CSCCommonTrigger_cfi import *
csctfTrackDigis = cms.EDProducer("CSCTFTrackProducer",
    DTproducer = cms.untracked.InputTag("dtTriggerPrimitiveDigis"),
    SectorReceiverInput = cms.untracked.InputTag("cscTriggerPrimitiveDigis","MPCSORTED"),
    SectorProcessor = cms.PSet(
        CSCCommonTrigger,
        SRLUT = cms.PSet(
            Binary = cms.untracked.bool(False),
            ReadLUTs = cms.untracked.bool(False),
            LUTPath = cms.untracked.string('./'),
            UseMiniLUTs = cms.untracked.bool(True)
        ),
        AllowALCTonly = cms.bool(False),
        PTLUT = cms.PSet(
            LowQualityFlag = cms.untracked.uint32(4),
            ReadPtLUT = cms.untracked.bool(False),
            PtMethod = cms.untracked.uint32(1)
        ),
        singlesTrackOutput = cms.uint32(3),
        singlesTrackPt = cms.uint32(255),
        trigger_on_MB1a = cms.bool(False),
        EtaMin = cms.vuint32(22, 22, 14, 14, 14,
            10, 10, 10),
        trigger_on_ME1a = cms.bool(False),
        trigger_on_ME1b = cms.bool(False),
        EtaMax = cms.vuint32(127, 127, 127, 127, 127,
            24, 24, 24),
        CoreLatency = cms.uint32(8),
        PreTrigger = cms.uint32(2),
        trigger_on_MB1d = cms.bool(False),
        run_core = cms.bool(True),
        mindeta_accp = cms.uint32(4),
        EtaWindows = cms.vuint32(4, 4, 4, 4, 4,
            4),
        AllowCLCTonly = cms.bool(False),
        BXAdepth = cms.uint32(2),
        mindphip = cms.uint32(2),
        maxdeta_accp = cms.uint32(16),
        trigger_on_ME4 = cms.bool(False),
        maxdphi_accp = cms.uint32(64),
        trigger_on_ME3 = cms.bool(False),
        trigger_on_ME2 = cms.bool(False),


        kill_fiber         = cms.uint32(0),
        QualityEnableME1a  = cms.uint32(65535),
        QualityEnableME1b  = cms.uint32(65535),
        QualityEnableME1c  = cms.uint32(65535),
        QualityEnableME1d  = cms.uint32(65535),
        QualityEnableME1e  = cms.uint32(65535),
        QualityEnableME1f  = cms.uint32(65535),
        QualityEnableME2a  = cms.uint32(65535),
        QualityEnableME2b  = cms.uint32(65535),
        QualityEnableME2c  = cms.uint32(65535),
        QualityEnableME3a  = cms.uint32(65535),
        QualityEnableME3b  = cms.uint32(65535),
        QualityEnableME3c  = cms.uint32(65535),
        QualityEnableME4a  = cms.uint32(65535),
        QualityEnableME4b  = cms.uint32(65535),
        QualityEnableME4c  = cms.uint32(65535),

        initializeFromPSet = cms.bool(True)
    ),
    isTMB07 = cms.bool(True),
    useDT = cms.bool(True)
)


