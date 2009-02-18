import FWCore.ParameterSet.Config as cms

EcalTrivialConditionRetriever = cms.ESSource("EcalTrivialConditionRetriever",
    producedChannelStatus = cms.untracked.bool(True),
    #       Values to get correct noise on RecHit amplitude using 3+5 weights
    EBpedRMSX12 = cms.untracked.double(1.26),
    weightsForTB = cms.untracked.bool(False),
    # channel status
    channelStatusFile = cms.untracked.string(''),
    producedEcalPedestals = cms.untracked.bool(True),
    #       If set true reading optimized weights (3+5 weights) from file 
    getWeightsFromFile = cms.untracked.bool(True),
    intercalibErrorsFile = cms.untracked.string(''),
    laserAPDPNMean = cms.untracked.double(1.0),
    laserAPDPNRefMean = cms.untracked.double(1.0),
    #       untracked string amplWeightsFile = "CalibCalorimetry/EcalTrivialCondModules/data/ampWeights_TB.txt"
    # file with intercalib constants - same format used for online and offline DB
    # by default set all inter calib const to 1.0 if no file provided
    intercalibConstantsFile = cms.untracked.string(''),
    producedEcalWeights = cms.untracked.bool(True),
    EEpedRMSX12 = cms.untracked.double(2.87),
    producedEcalIntercalibConstants = cms.untracked.bool(True),
    producedEcalIntercalibErrors = cms.untracked.bool(True),
    producedEcalTimeCalibConstants = cms.untracked.bool(True),
    producedEcalTimeCalibErrors = cms.untracked.bool(True),
    producedEcalLaserCorrection = cms.untracked.bool(True),
    producedEcalGainRatios = cms.untracked.bool(True),
    producedEcalADCToGeVConstant = cms.untracked.bool(True),
    adcToGeVEBConstant = cms.untracked.double(0.035),
    adcToGeVEEConstant = cms.untracked.double(0.06),
    # cluster functions/corrections -- by default no parameters are passed
    producedEcalClusterLocalContCorrParameters = cms.untracked.bool(True),
    localContCorrParameters = cms.untracked.vdouble( ),
    producedEcalClusterCrackCorrParameters = cms.untracked.bool(True),
    crackCorrParameters = cms.untracked.vdouble( ),
    mappingFile = cms.untracked.string('Geometry/EcalMapping/data/EEMap.txt'),
    producedEcalMappingElectronics = cms.untracked.bool(True)                                             
)
