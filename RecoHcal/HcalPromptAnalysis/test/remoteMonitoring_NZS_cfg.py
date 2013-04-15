import FWCore.ParameterSet.Config as cms

process = cms.Process("testAnalyzer")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     'file:/tmp/zhokin/HcalNZSRAW-199-751.root'
   )
)

process.Analyzer = cms.EDAnalyzer("VeRawAnalyzer",
                                  #
                                  MapCreation = cms.untracked.int32(1),
                                  #
                                  #recordNtuples = cms.untracked.bool(False),
                                  recordNtuples = cms.untracked.bool(True),
                                  maxNeventsInNtuple = cms.int32(100),
                                  #
                                  #recordHistoes = cms.untracked.bool(False),
                                  recordHistoes = cms.untracked.bool(True),
                                  #
                                  ##scripts: zRunRatio34.C, zRunNbadchan.C
                                  studyRunDependenceHist = cms.untracked.bool(True),
                                  #studyRunDependenceHist = cms.untracked.bool(False),
                                  #
                                  ##scripts: zerrors.C
                                  studyCapIDErrorsHist = cms.untracked.bool(True),
                                  #studyCapIDErrorsHist = cms.untracked.bool(False),
                                  #
                                  ##scripts: zrms.C
                                  studyRMSshapeHist = cms.untracked.bool(True),
                                  #studyRMSshapeHist = cms.untracked.bool(False),
                                  #
                                  ##scripts: zratio34.C
                                  studyRatioShapeHist = cms.untracked.bool(True),
                                  #studyRatioShapeHist = cms.untracked.bool(False),
                                  #
                                  ##scripts: zadcamplitude.C
                                  studyADCAmplHist = cms.untracked.bool(True),
                                  #studyADCAmplHist = cms.untracked.bool(False),
                                  #
                                  ##scripts: ztsmean.C
                                  studyTSmeanShapeHist = cms.untracked.bool(True),
                                  #studyTSmeanShapeHist = cms.untracked.bool(False),
                                  #
                                  ##scripts: ztsmaxa.C
                                  studyTSmaxShapeHist = cms.untracked.bool(True),
                                  #studyTSmaxShapeHist = cms.untracked.bool(False),
                                  #
                                  ##scripts: zcalib....C
                                  studyCalibCellsHist = cms.untracked.bool(True),
                                  #studyCalibCellsHist = cms.untracked.bool(False),
                                  #
                                  ##scripts: zdifampl.C
                                  studyDiffAmplHist = cms.untracked.bool(True),
                                  #studyDiffAmplHist = cms.untracked.bool(False),
                                  #
                                  ##DigiCollectionLabel = cms.untracked.InputTag("hcalDigis"),
                                  #Verbosity = cms.untracked.int32(-54),
                                  #Verbosity = cms.untracked.int32(-22),
                                  #Verbosity = cms.untracked.int32(-11),
                                  #Verbosity = cms.untracked.int32(-12),
                                  #Verbosity = cms.untracked.int32(-13),
                                  #Verbosity = cms.untracked.int32(-51),
                                  #Verbosity = cms.untracked.int32(-24),
                                  #Verbosity = cms.untracked.int32(-244),
                                  #Verbosity = cms.untracked.int32(-233),
                                  Verbosity = cms.untracked.int32(0),
                                  #
                                  #
                                  #         Normal channels:
                                  #
                                  # -53 for  BAD HBHEHF channels from study on shape Ratio
                                  #Verbosity = cms.untracked.int32(-53),
                                  ratioHBMin = cms.double(0.50),
                                  ratioHBMax = cms.double(0.88),
                                  ratioHEMin = cms.double(0.50),
                                  ratioHEMax = cms.double(0.93),
                                  ratioHFMin = cms.double(0.30),
                                  ratioHFMax = cms.double(1.02),
                                  ratioHOMin = cms.double(0.54),
                                  ratioHOMax = cms.double(1.04),
                                  # -54 for  BAD HBHEHF channels from study on RMS of shapes
                                  #Verbosity = cms.untracked.int32(-54),
                                  rmsHBMin = cms.double(1.0),
                                  rmsHBMax = cms.double(2.3),
                                  rmsHEMin = cms.double(1.0),
                                  rmsHEMax = cms.double(2.3),
                                  rmsHFMin = cms.double(0.2),
                                  rmsHFMax = cms.double(2.8),
                                  rmsHOMin = cms.double(0.2),
                                  rmsHOMax = cms.double(2.4),
                                  # -55 for  BAD HBHEHF channels from study on TSmean of shapes
                                  #Verbosity = cms.untracked.int32(-55),
                                  TSmeanHBMin = cms.double(1.0),
                                  TSmeanHBMax = cms.double(5.2),
                                  TSmeanHEMin = cms.double(1.0),
                                  TSmeanHEMax = cms.double(5.2),
                                  TSmeanHFMin = cms.double(1.5),
                                  TSmeanHFMax = cms.double(4.5),
                                  TSmeanHOMin = cms.double(2.5),
                                  TSmeanHOMax = cms.double(4.5),
                                  # -55 for  BAD HBHEHF channels from study on TSmax of shapes
                                  #Verbosity = cms.untracked.int32(-55),
                                  TSpeakHBMin = cms.double(1.5),
                                  TSpeakHBMax = cms.double(6.5),
                                  TSpeakHEMin = cms.double(1.5),
                                  TSpeakHEMax = cms.double(6.5),
                                  TSpeakHFMin = cms.double(0.5),
                                  TSpeakHFMax = cms.double(5.5),
                                  TSpeakHOMin = cms.double(0.5),
                                  TSpeakHOMax = cms.double(5.5),
                                  # -56 for  BAD HBHEHOHF channels from study on ADC Amplitude
                                  #Verbosity = cms.untracked.int32(-56),
                                  ADCAmplHBMin = cms.double(100.),
                                  ADCAmplHEMin = cms.double(200.),
                                  ADCAmplHFMin = cms.double(30.),
                                  ADCAmplHOMin = cms.double(300.),
                                  #
                                  #             CALIBRATION channels:
                                  #
                                  # for  BAD HBHEHOHF CALIBRATION channels from study on ADC amplitude
                                  # cuts for Laser runs:
                                  #calibrADCHBMin = cms.double(15.0),
                                  #calibrADCHEMin = cms.double(15.0),
                                  #calibrADCHOMin = cms.double(15.0),
                                  #calibrADCHFMin = cms.double(15.0),
                                  # cuts for LED runs:
                                  calibrADCHBMin = cms.double(200.),
                                  calibrADCHEMin = cms.double(200.),
                                  calibrADCHOMin = cms.double(200.),
                                  calibrADCHFMin = cms.double(100.),
                                  # for  BAD HBHEHOHF CALIBRATION channels from study on shape Ratio
                                  calibrRatioHBMin = cms.double(0.63),
                                  calibrRatioHEMin = cms.double(0.63),
                                  calibrRatioHOMin = cms.double(0.10),
                                  calibrRatioHFMin = cms.double(0.21),
                                  # for  BAD HBHEHOHF CALIBRATION channels from study on TSmax
                                  calibrTSmaxHBMin = cms.double(0.50),
                                  calibrTSmaxHBMax = cms.double(8.50),
                                  calibrTSmaxHEMin = cms.double(0.50),
                                  calibrTSmaxHEMax = cms.double(8.50),
                                  calibrTSmaxHOMin = cms.double(0.50),
                                  calibrTSmaxHOMax = cms.double(9.50),
                                  calibrTSmaxHFMin = cms.double(0.50),
                                  calibrTSmaxHFMax = cms.double(9.50),
                                  # for  BAD HBHEHOHF CALIBRATION channels from study on TSmean
                                  calibrTSmeanHBMin = cms.double(1.00),
                                  calibrTSmeanHBMax = cms.double(4.80),
                                  calibrTSmeanHEMin = cms.double(1.00),
                                  calibrTSmeanHEMax = cms.double(4.20),
                                  calibrTSmeanHOMin = cms.double(1.00),
                                  calibrTSmeanHOMax = cms.double(3.50),
                                  calibrTSmeanHFMin = cms.double(1.00),
                                  calibrTSmeanHFMax = cms.double(5.20),
                                  # for  BAD HBHEHOHF CALIBRATION channels from study on Width
                                  calibrWidthHBMin = cms.double(1.00),
                                  calibrWidthHBMax = cms.double(2.30),
                                  calibrWidthHEMin = cms.double(1.00),
                                  calibrWidthHEMax = cms.double(2.30),
                                  calibrWidthHOMin = cms.double(1.00),
                                  calibrWidthHOMax = cms.double(2.20),
                                  calibrWidthHFMin = cms.double(0.20),
                                  calibrWidthHFMax = cms.double(2.40),
                                  #
                                  # Special task of run or LS quality:
                                  #
                                  # flag for choise of criterion of bad channels:
                                  # =0-CapIdErr, =1-Ratio, =2-Width, =3-TSmax, =4-TSmean, =5-adcAmplitud
                                  flagtodefinebadchannel = cms.int32(0),
                                  #
                                  # to see run dependence of rate of events with Nbadchannels bigger:
                                  nbadchannels1 = cms.int32(5),
                                  nbadchannels2 = cms.int32(30),
                                  nbadchannels3 = cms.int32(60),
                                  #
                                  #
                                  HistOutFile = cms.untracked.string('testNZS.root'),
#                                  HistOutFile = cms.untracked.string('testLaser178165.root'),
#                                  HistOutFile = cms.untracked.string('testLaser141849.root'),
#                                  HistOutFile = cms.untracked.string('testLaser133070.root'),
#                                  HistOutFile = cms.untracked.string('testLaser211875.root'),
                                  MAPOutFile = cms.untracked.string('LogEleMapdb.h')
                                  #
                                  ##OutputFilePath = cms.string('/tmp/zhokin/'),        
                                  ##OutputFileExt = cms.string(''),
                                  #
)		

process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.hcalDigis.FilterDataQuality = cms.bool(False)
				
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_P_V41::All'
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.p = cms.Path(process.hcalDigis*process.Analyzer)

process.MessageLogger = cms.Service("MessageLogger",
     categories   = cms.untracked.vstring(''),
     destinations = cms.untracked.vstring('cout'),
     debugModules = cms.untracked.vstring('*'),
     cout = cms.untracked.PSet(
         threshold = cms.untracked.string('WARNING'),
	 WARNING = cms.untracked.PSet(limit = cms.untracked.int32(0))
     )
 )



