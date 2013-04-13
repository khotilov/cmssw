import FWCore.ParameterSet.Config as cms

process = cms.Process("writeGBRForests")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(1) # CV: needs to be set to 1 so that GBRForestWriter::analyze method gets called exactly once         
)

process.source = cms.Source("EmptySource")

process.load('Configuration/StandardSequences/Services_cff')

process.gbrForestWriter = cms.EDAnalyzer("GBRForestWriter",
    jobs = cms.VPSet(
        cms.PSet(
            inputFileName = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_Dec2012.root'),
            gbrForestName = cms.string("U1Correction"),
            outputRecord = cms.string("mvaPFMET_53_Dec2012_U")
        ),
        cms.PSet(
            inputFileName = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_Dec2012.root'),
            gbrForestName = cms.string("PhiCorrection"),
            outputRecord = cms.string("mvaPFMET_53_Dec2012_DPhi")
        ),
        cms.PSet(
            inputFileName = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012.root'),
            gbrForestName = cms.string("CovU1"),
            outputRecord = cms.string("mvaPFMET_53_Dec2012_CovU1")
        ),
        cms.PSet(
            inputFileName = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012.root'),
            gbrForestName = cms.string("CovU2"),
            outputRecord = cms.string("mvaPFMET_53_Dec2012_CovU2")
        )                                             
    )
)

process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.CondDBCommon.connect = 'sqlite_file:../data/mvaPFMEt_53_Dec2012.db'

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDBCommon,
    timetype = cms.untracked.string('runnumber'),
    toPut = cms.VPSet(
        cms.PSet(
            record = cms.string('mvaPFMET_53_Dec2012_U'),
            tag = cms.string('mvaPFMET_53_Dec2012_U')
        ),
        cms.PSet(
            record = cms.string('mvaPFMET_53_Dec2012_DPhi'),
            tag = cms.string('mvaPFMET_53_Dec2012_DPhi')
        ),
        cms.PSet(
            record = cms.string('mvaPFMET_53_Dec2012_CovU1'),
            tag = cms.string('mvaPFMET_53_Dec2012_CovU1')
        ),
        cms.PSet(
            record = cms.string('mvaPFMET_53_Dec2012_CovU2'),
            tag = cms.string('mvaPFMET_53_Dec2012_CovU2')
        )
    )
)

process.p = cms.Path(process.gbrForestWriter)
