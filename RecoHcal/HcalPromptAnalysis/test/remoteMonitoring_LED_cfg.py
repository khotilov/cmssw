import FWCore.ParameterSet.Config as cms

process = cms.Process("testAnalyzer")

process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(1000)
   input = cms.untracked.int32(-1)
)

#process.source = cms.Source("PoolSource",
process.source = cms.Source("HcalTBSource",
    fileNames = cms.untracked.vstring(
#
## Laser:
#     'file:rfiles/USC_211775.root' ## only HF digis 10K
#     'file:rfiles/USC_211673.root' ## only HBHE digis 10K
#
## LED digis 10K:
     'file:rfiles/USC_194057.root',
#     'file:rfiles/USC_194124.root',
#     'file:rfiles/USC_194165.root',
#     'file:rfiles/USC_194807.root',
#     'file:rfiles/USC_195172.root',
     'file:rfiles/USC_195682.root',
     'file:rfiles/USC_196870.root',
#     'file:rfiles/USC_198084.root',
#     'file:rfiles/USC_199392.root',
#     'file:rfiles/USC_200281.root',
     'file:rfiles/USC_200707.root',
     'file:rfiles/USC_201305.root',
#     'file:rfiles/USC_201571.root',
#     'file:rfiles/USC_202066.root',
     'file:rfiles/USC_202341.root',
     'file:rfiles/USC_203339.root',
#     'file:rfiles/USC_203536.root',
     'file:rfiles/USC_204172.root',
     'file:rfiles/USC_205355.root',
#     'file:rfiles/USC_205928.root',
     'file:rfiles/USC_206456.root',
#     'file:rfiles/USC_207167.root',
     'file:rfiles/USC_207406.root',
#     'file:rfiles/USC_207545.root',
#     'file:rfiles/USC_207914.root',
#     'file:rfiles/USC_208361.root',
     'file:rfiles/USC_208565.root',
#     'file:rfiles/USC_209311.root',
#     'file:rfiles/USC_209997.root',
     'file:rfiles/USC_210378.root',
#     'file:rfiles/USC_211006.root',
     'file:rfiles/USC_211659.root'
#before:     
#     'file:rfiles/USC_195682.root',
#     'file:rfiles/USC_196870.root',
#     'file:rfiles/USC_201305.root',
#     'file:rfiles/USC_202341.root',
#     'file:rfiles/USC_203339.root',
#     'file:rfiles/USC_204172.root',
#     'file:rfiles/USC_205355.root',
#     'file:rfiles/USC_206456.root'
#
   ), 
    streams = cms.untracked.vstring(
		  "HCAL_Trigger",
		  "HCAL_DCC700","HCAL_DCC701","HCAL_DCC702","HCAL_DCC703","HCAL_DCC704","HCAL_DCC705",
		  "HCAL_DCC706","HCAL_DCC707","HCAL_DCC708","HCAL_DCC709","HCAL_DCC710","HCAL_DCC711",
		  "HCAL_DCC712","HCAL_DCC713","HCAL_DCC714","HCAL_DCC715","HCAL_DCC716","HCAL_DCC717",
		  "HCAL_DCC718","HCAL_DCC719","HCAL_DCC720","HCAL_DCC721","HCAL_DCC722","HCAL_DCC723",
		  "HCAL_DCC724","HCAL_DCC725","HCAL_DCC726","HCAL_DCC727","HCAL_DCC728","HCAL_DCC729",
		  "HCAL_DCC730","HCAL_DCC731"
		 )	
)
process.Analyzer = cms.EDAnalyzer("VeRawAnalyzer",
                                  #
                                  MapCreation = cms.untracked.int32(1),
                                  #
                                  recordNtuples = cms.untracked.bool(False),
                                  #recordNtuples = cms.untracked.bool(True),
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
                                  ratioHBMin = cms.double(0.48),
                                  ratioHBMax = cms.double(0.92),
                                  ratioHEMin = cms.double(0.48),
                                  ratioHEMax = cms.double(0.92),
                                  ratioHFMin = cms.double(0.48),
                                  ratioHFMax = cms.double(1.02),
                                  ratioHOMin = cms.double(0.48),
                                  ratioHOMax = cms.double(1.04),
                                  # -54 for  BAD HBHEHF channels from study on RMS of shapes
                                  #Verbosity = cms.untracked.int32(-54),
                                  rmsHBMin = cms.double(0.7),
                                  rmsHBMax = cms.double(2.4),
                                  rmsHEMin = cms.double(0.7),
                                  rmsHEMax = cms.double(2.4),
                                  rmsHFMin = cms.double(0.1),
                                  rmsHFMax = cms.double(2.5),
                                  rmsHOMin = cms.double(0.1),
                                  rmsHOMax = cms.double(2.5),
                                  # -55 for  BAD HBHEHF channels from study on TSmean of shapes
                                  #Verbosity = cms.untracked.int32(-55),
                                  TSmeanHBMin = cms.double(1.0),
                                  TSmeanHBMax = cms.double(5.2),
                                  TSmeanHEMin = cms.double(1.0),
                                  TSmeanHEMax = cms.double(5.2),
                                  TSmeanHFMin = cms.double(1.0),
                                  TSmeanHFMax = cms.double(4.5),
                                  TSmeanHOMin = cms.double(1.0),
                                  TSmeanHOMax = cms.double(5.3),
                                  # -55 for  BAD HBHEHF channels from study on TSmax of shapes
                                  #Verbosity = cms.untracked.int32(-55),
                                  TSpeakHBMin = cms.double(1.5),
                                  TSpeakHBMax = cms.double(6.5),
                                  TSpeakHEMin = cms.double(1.5),
                                  TSpeakHEMax = cms.double(6.5),
                                  TSpeakHFMin = cms.double(0.5),
                                  TSpeakHFMax = cms.double(6.5),
                                  TSpeakHOMin = cms.double(0.5),
                                  TSpeakHOMax = cms.double(6.5),
                                  # -56 for  BAD HBHEHOHF channels from study on ADC Amplitude
                                  #Verbosity = cms.untracked.int32(-56),
                                  ADCAmplHBMin = cms.double(70.),
                                  ADCAmplHEMin = cms.double(70.),
                                  ADCAmplHFMin = cms.double(40.),
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
                                  calibrADCHBMin = cms.double(120.),
                                  calibrADCHEMin = cms.double(120.),
                                  calibrADCHOMin = cms.double(120.),
                                  calibrADCHFMin = cms.double(60.),
                                  # for  BAD HBHEHOHF CALIBRATION channels from study on shape Ratio
                                  calibrRatioHBMin = cms.double(0.60),
                                  calibrRatioHEMin = cms.double(0.64),
                                  calibrRatioHOMin = cms.double(0.25),
                                  calibrRatioHFMin = cms.double(0.25),
                                  # for  BAD HBHEHOHF CALIBRATION channels from study on TSmax
                                  calibrTSmaxHBMin = cms.double(0.50),
                                  calibrTSmaxHBMax = cms.double(8.50),
                                  calibrTSmaxHEMin = cms.double(0.50),
                                  calibrTSmaxHEMax = cms.double(8.50),
                                  calibrTSmaxHOMin = cms.double(0.50),
                                  calibrTSmaxHOMax = cms.double(8.50),
                                  calibrTSmaxHFMin = cms.double(0.50),
                                  calibrTSmaxHFMax = cms.double(9.50),
                                  # for  BAD HBHEHOHF CALIBRATION channels from study on TSmean
                                  calibrTSmeanHBMin = cms.double(1.00),
                                  calibrTSmeanHBMax = cms.double(5.50),
                                  calibrTSmeanHEMin = cms.double(1.00),
                                  calibrTSmeanHEMax = cms.double(5.50),
                                  calibrTSmeanHOMin = cms.double(1.00),
                                  calibrTSmeanHOMax = cms.double(3.50),
                                  calibrTSmeanHFMin = cms.double(1.00),
                                  calibrTSmeanHFMax = cms.double(5.20),
                                  # for  BAD HBHEHOHF CALIBRATION channels from study on Width
                                  calibrWidthHBMin = cms.double(1.00),
                                  calibrWidthHBMax = cms.double(2.30),
                                  calibrWidthHEMin = cms.double(1.00),
                                  calibrWidthHEMax = cms.double(2.30),
                                  calibrWidthHOMin = cms.double(0.50),
                                  calibrWidthHOMax = cms.double(3.00),
                                  calibrWidthHFMin = cms.double(0.20),
                                  calibrWidthHFMax = cms.double(2.30),
                                  #
                                  # Special task of run or LS quality:
                                  #
                                  # flag for choise of criterion of bad channels:
                                  #=0-CapIdErr, =1-Ratio, =2-Width, =3-TSmax, =4-TSmean, =5-adcAmplitud
                                  flagtodefinebadchannel = cms.int32(2),
                                  #how many bins you want on the plots:better to choice (Nruns+1)
                                  howmanybinsonplots = cms.int32(15),
                                  #
                                  # to see run dependence of rate of events with Nbadchannels bigger:
                                  nbadchannels1 = cms.int32(5),
                                  nbadchannels2 = cms.int32(30),
                                  nbadchannels3 = cms.int32(60),
                                  #
                                  #
                                  HistOutFile = cms.untracked.string('test.root'),
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
process.hcalDigis.InputLabel = cms.InputTag("source")
				
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



