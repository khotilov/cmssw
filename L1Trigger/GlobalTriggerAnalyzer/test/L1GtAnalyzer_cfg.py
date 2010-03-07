#
# cfg file to run the L1 GT test analyzer according to 
#   the options set in "user choices"
#

import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("L1GtAnalyzer")

print '\n'
from L1Trigger.GlobalTriggerAnalyzer.UserOptions_cff import *
if errorUserOptions == True :
    print '\nError returned by UserOptions_cff\n'
    sys.exit()


# source according to data type
if dataType == 'StreamFile' :
    process.source = cms.Source("NewEventStreamFileReader", fileNames=readFiles)
else :        
    process.source = cms.Source ('PoolSource', fileNames=readFiles, secondaryFileNames=secFiles)


# number of events to be processed and source file
process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(maxNumberEvents)
)

#
# load and configure modules via Global Tag
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions

process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = useGlobalTag+'::All'

# processes to be run


process.load("L1Trigger.GlobalTriggerAnalyzer.l1GtAnalyzer_cfi")
#
# input tag for GT readout collection: 
#process.l1GtAnalyzer.L1GtDaqInputTag = 'gtDigis' 
 
# input tags for GT lite record
#process.l1GtAnalyzer.L1GtRecordInputTag = 'l1GtRecord'

# input tag for GT object map collection
#process.l1GtAnalyzer.L1GtObjectMapTag = 'hltL1GtObjectMap'

# physics algorithm name or alias, technical trigger name 
#process.l1GtAnalyzer.AlgorithmName = 'L1_SingleEG5'
process.l1GtAnalyzer.AlgorithmName = 'L1_BscMinBiasOR_BptxPlusORMinus'
#process.l1GtAnalyzer.AlgorithmName = 'L1Tech_BPTX_plus_AND_minus_instance1.v0'
#process.l1GtAnalyzer.AlgorithmName = 'L1Tech_BPTX_quiet.v0'
#process.l1GtAnalyzer.AlgorithmName = 'L1Tech_BPTX_plus_AND_minus.v0'

# condition in the above algorithm to test the object maps
process.l1GtAnalyzer.ConditionName = 'SingleNoIsoEG_0x0A'

# a bit number
process.l1GtAnalyzer.BitNumber = 10




process.load("L1Trigger.GlobalTriggerAnalyzer.l1GtTrigReport_cfi")

# boolean flag to select the input record
#process.l1GtTrigReport.UseL1GlobalTriggerRecord = True

# input tag for the GT record requested: 
#   GT emulator:    gtDigis (DAQ record)
#   GT unpacker:    gtDigis (DAQ record)
#   GT lite record: l1GtRecord 
#process.l1GtTrigReport.L1GtRecordInputTag = "gtDigis"

#process.l1GtTrigReport.PrintVerbosity = 2

# print output: 0 = std::cout; 1 = LogTrace; 2 = LogVerbatim; 3 = LogInfo
#process.l1GtTrigReport.PrintOutput = 1


# for RAW data, run first the RAWTODIGI 
if (dataType == 'RAW') and not (useRelValSample) :
    process.load('Configuration/StandardSequences/RawToDigi_Data_cff')
    process.l1GtTrigReport.L1GtRecordInputTag = "gtDigis"
    process.p = cms.Path(process.RawToDigi+process.l1GtTrigReport+process.l1GtAnalyzer)

elif (dataType == 'RAW') and (useRelValSample) :
    process.load('Configuration/StandardSequences/RawToDigi_cff')
    process.l1GtTrigReport.L1GtRecordInputTag = "gtDigis"
    process.p = cms.Path(process.RawToDigi+process.l1GtTrigReport+process.l1GtAnalyzer)
    
else :        
    # path to be run for RECO
    process.p = cms.Path(process.l1GtTrigReport+process.l1GtAnalyzer)


# Message Logger
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.debugModules = ['l1GtAnalyzer']
process.MessageLogger.categories.append('L1GtAnalyzer')
process.MessageLogger.categories.append('L1GtUtils')

process.MessageLogger.cerr.default.limit = 0
process.MessageLogger.cerr.FwkJob.limit = 0
process.MessageLogger.cerr.FwkReport.limit = 0
process.MessageLogger.cerr.FwkSummary.limit = 0

process.MessageLogger.debugs = cms.untracked.PSet( 
        threshold = cms.untracked.string('DEBUG'),
        DEBUG = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
        INFO = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
        WARNING = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
        ERROR = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
        L1GtAnalyzer = cms.untracked.PSet( limit = cms.untracked.int32(-1) ), 
        L1GtUtils = cms.untracked.PSet( limit = cms.untracked.int32(-1) ) 
        )

process.MessageLogger.warnings = cms.untracked.PSet( 
        threshold = cms.untracked.string('WARNING'),
        WARNING = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
        ERROR = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
        L1GtAnalyzer = cms.untracked.PSet( limit = cms.untracked.int32(-1) ), 
        L1GtUtils = cms.untracked.PSet( limit = cms.untracked.int32(-1) ) 
        )

process.MessageLogger.errors = cms.untracked.PSet( 
        threshold = cms.untracked.string('ERROR'),
        ERROR = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
        L1GtAnalyzer = cms.untracked.PSet( limit = cms.untracked.int32(-1) ), 
        L1GtUtils = cms.untracked.PSet( limit = cms.untracked.int32(-1) ) 
        )