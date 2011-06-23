import FWCore.ParameterSet.Config as cms

process = cms.Process("EcalTimeCalibrationValidator")

# Global Tag -- for original timing calibrations
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_43_V3::All'

# shaping our Message logger to suit our needs
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
    categories = cms.untracked.vstring('*'),
    destinations = cms.untracked.vstring('cout')
)

process.expressValidator = cms.EDAnalyzer("EcalTimeCalibrationValidator",
  InputFileNames = cms.vstring(#"file:/data2/kubota/TimingCalibrationOct384/src/CalibCalorimetry/EcalTiming/test/input_files/ecaltime-run-143953-144114/2ndhalf1.root",
                              #"file:/data2/kubota/TimingCalibrationOct384/src/CalibCalorimetry/EcalTiming/test/input_files/ecaltime-run-143953-144114/EcalTimeTree_999999_170_1_f5s.root"
                       'file:/data/franzoni/data/Ecaltime/skimming-jun17/rhQuinq-SingleEletron-Run2011A-WElectron-PromptSkim-jun16_HADDED.root'
    ),
  OutputFileName = cms.string("file:converted1.root"),
  CalibConstantXMLFileName = cms.string("EcalTimeCalibs.xml"),
  MaxTreeEntriesToProcess = cms.untracked.int32(-1)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("EmptySource",
       numberEventsInRun = cms.untracked.uint32(1),
       firstRun = cms.untracked.uint32(888888) # Use last IOV
)


process.p = cms.Path(process.expressValidator)
