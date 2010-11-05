import FWCore.ParameterSet.Config as cms

process = cms.Process("EcalTimeCalibrationValidator")

# shaping our Message logger to suit our needs
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
    categories = cms.untracked.vstring('*'),
    destinations = cms.untracked.vstring('cout')
)

process.expressValidator = cms.EDAnalyzer("EcalTimeCalibrationValidator",
  InputFileNames = cms.vstring("file:/data2/kubota/TimingCalibrationOct384/src/CalibCalorimetry/EcalTiming/test/input_files/ecaltime-run-143953-144114/2ndhalf1.root",
                              "file:/data2/kubota/TimingCalibrationOct384/src/CalibCalorimetry/EcalTiming/test/input_files/ecaltime-run-143953-144114/EcalTimeTree_999999_170_1_f5s.root"
    ),
  #InputFileNames = cms.vstring("file:EcalTimeTree_147114_9_1_uaG.root"),
  OutputFileName = cms.string("file:converted1.root"),
  CalibConstantXMLFileName = cms.string("myCalibBoth.xml"),
  MaxTreeEntriesToProcess = cms.untracked.int32(100000000)
)

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(1) )
process.source = cms.Source("EmptySource")

process.p = cms.Path(process.expressValidator)
