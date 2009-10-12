# Auto generated configuration file
# using: 
# Revision: 1.123 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: alCaRecoSplitting -s ALCA:MuAlCalIsolatedMu+RpcCalHLT+EcalCalPi0Calib+EcalCalEtaCalib+DQM --datatier RECO --eventcontent RECO --conditions FrontierConditions_GlobalTag,MC_31X_V1::All --no_exec --scenario=pp
import FWCore.ParameterSet.Config as cms

process = cms.Process('ALCA')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/AlCaRecoStreams_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('alCaRecoSplitting nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('alCaRecoSplitting_RECO.root')
)

# Additional output definition
process.ALCARECOStreamRpcCalHLT = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECORpcCalHLT')
    ),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_muonDTDigis_*_*', 
        'keep CSCDetIdCSCWireDigiMuonDigiCollection_*_*_*', 
        'keep CSCDetIdCSCStripDigiMuonDigiCollection_*_*_*', 
        'keep DTLayerIdDTDigiMuonDigiCollection_*_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep RPCDetIdRPCDigiMuonDigiCollection_*_*_*', 
        'keep L1MuRegionalCands_*_RPCb_*', 
        'keep L1MuRegionalCands_*_RPCf_*', 
        'keep L1MuGMTCands_*_*_*', 
        'keep L1MuGMTReadoutCollection_*_*_*'),
    fileName = cms.untracked.string('RpcCalHLT.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('StreamRpcCalHLT'),
        dataTier = cms.untracked.string('ALCARECO')
    )
)
process.ALCARECOStreamMuAlCalIsolatedMu = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOMuAlCalIsolatedMu')
    ),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_ALCARECOMuAlCalIsolatedMu_*_*', 
        'keep *_muonCSCDigis_*_*', 
        'keep *_muonDTDigis_*_*', 
        'keep *_muonRPCDigis_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt2DSegments_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*'),
    fileName = cms.untracked.string('MuAlCalIsolatedMu.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('StreamMuAlCalIsolatedMu'),
        dataTier = cms.untracked.string('ALCARECO')
    )
)
process.ALCARECOStreamEcalCalPi0Calib = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOEcalCalPi0Calib')
    ),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_ecalPi0Corrected_pi0EcalRecHitsEB_*', 
        'keep *_ecalPi0Corrected_pi0EcalRecHitsEE_*', 
        'keep *_hltAlCaPi0RegRecHits_pi0EcalRecHitsES_*'),
    fileName = cms.untracked.string('ALCARECOEcalCalPi0Calib.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('StreamALCARECOEcalCalPi0Calib'),
        dataTier = cms.untracked.string('ALCARECO')
    )
)
process.ALCARECOStreamEcalCalEtaCalib = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOEcalCalEtaCalib')
    ),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_ecalEtaCorrected_etaEcalRecHitsEB_*', 
        'keep *_ecalEtaCorrected_etaEcalRecHitsEE_*', 
        'keep *_hltAlCaEtaRegRecHits_etaEcalRecHitsES_*'),
    fileName = cms.untracked.string('ALCARECOEcalCalEtaCalib.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('StreamALCARECOEcalCalEtaCalib'),
        dataTier = cms.untracked.string('ALCARECO')
    )
)

# Other statements
process.GlobalTag.globaltag = 'MC_31X_V1::All'

# Path and EndPath definitions
process.ALCARECOStreamRpcCalHLTOutPath = cms.EndPath(process.ALCARECOStreamRpcCalHLT)
process.ALCARECOStreamMuAlCalIsolatedMuOutPath = cms.EndPath(process.ALCARECOStreamMuAlCalIsolatedMu)
process.ALCARECOStreamEcalCalPi0CalibOutPath = cms.EndPath(process.ALCARECOStreamEcalCalPi0Calib)
process.ALCARECOStreamEcalCalEtaCalibOutPath = cms.EndPath(process.ALCARECOStreamEcalCalEtaCalib)

# Schedule definition
process.schedule = cms.Schedule(process.ALCARECOStreamRpcCalHLTOutPath,process.ALCARECOStreamMuAlCalIsolatedMuOutPath,process.ALCARECOStreamEcalCalPi0CalibOutPath,process.ALCARECOStreamEcalCalEtaCalibOutPath)
