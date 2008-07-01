import FWCore.ParameterSet.Config as cms

simSiStripDigis = cms.EDFilter("SiStripDigitizer",
    TOFCutForPeak = cms.double(100.0),
    DeltaProductionCut = cms.double(0.120425),
    Temperature = cms.double(263.0),
    #-------------------------------------
    #-----SiHitDigitizer
    DepletionVoltage = cms.double(140.0),
    SigmaShapePeak = cms.double(52.17),
    TOFCutForDeconvolution = cms.double(50.0),
    SimProducer = cms.string('SimDigitizer'),
    VerbosityLevel = cms.untracked.int32(0),
    #---------------------------------------
    #SiTrivialInduceChargeOnStrips
    # DECOnvolution Mode
    # TIB
    CouplingCostantDecoTIB = cms.vdouble(0.76, 0.12),
    # TID
    CouplingCostantDecoTID = cms.vdouble(0.76, 0.12),
    # TOB
    CouplingCostantDecoTOB = cms.vdouble(0.76, 0.12),
    # TEC
    CouplingCostantDecoTEC = cms.vdouble(0.76, 0.12),
    # PEAK Mode
    # TIB
    CouplingCostantPeakTIB = cms.vdouble(0.94, 0.03),
    # TID
    CouplingCostantPeakTID = cms.vdouble(0.94, 0.03),
    # TOB
    CouplingCostantPeakTOB = cms.vdouble(0.94, 0.03),
    # TEC
    CouplingCostantPeakTEC = cms.vdouble(0.94, 0.03),
    #
    AppliedVoltage = cms.double(150.0),
    ZeroSuppression = cms.bool(True),
    Noise = cms.bool(True), ## NOTE : turning Noise ON/OFF will make a big change

    noDiffusion = cms.bool(False),
    #--------------------------------------
    #SiLinearChargeDivider
    LandauFluctuations = cms.bool(True),
    #--------------------------------
    #---ZeroSuppression
    FedAlgorithm = cms.int32(4),
    AdcFullScale = cms.int32(255),
    DigiModeList = cms.PSet(
        SCDigi = cms.string('ScopeMode'),
        ZSDigi = cms.string('ZeroSuppressed'),
        PRDigi = cms.string('ProcessedRaw'),
        VRDigi = cms.string('VirginRaw')
    ),
    #-----SiStripDigitizer 
    TrackerConfigurationFromDB = cms.bool(False),
    ROUList = cms.vstring("g4SimHitsTrackerHitsTIBLowTof","g4SimHitsTrackerHitsTIBHighTof",
                          "g4SimHitsTrackerHitsTIDLowTof","g4SimHitsTrackerHitsTIDHighTof",
                          "g4SimHitsTrackerHitsTOBLowTof","g4SimHitsTrackerHitsTOBHighTof",
                          "g4SimHitsTrackerHitsTECLowTof","g4SimHitsTrackerHitsTECHighTof"),
    GevPerElectron = cms.double(3.61e-09),
    chargeDivisionsPerStrip = cms.int32(10),
    GeometryType = cms.string('idealForDigi'),
    #-----SiStripDigitizerAlgorithm 
    electronPerAdc = cms.double(250.0),
    APVpeakmode = cms.bool(False),
    SigmaShapeDeco = cms.double(12.06),
    NoiseSigmaThreshold = cms.double(2.0),
    ChargeDistributionRMS = cms.double(6.5e-10),
    CosmicDelayShift = cms.untracked.double(0.0),
    ChargeMobility = cms.double(480.0)
)



