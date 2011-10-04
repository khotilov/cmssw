##*****************************************************
##*****************************************************
##******** for MC
##*****************************************************
##*****************************************************
import FWCore.ParameterSet.Config as cms

from SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi import *
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from RecoBTag.ImpactParameter.impactParameter_cfi import *

process = cms.Process("Mistag")

process.source = cms.Source(
    "PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(   
# Summer11 S4 MC
#         'file:FEB0DEC7-DF8D-E011-8D10-00266CF33288.root'
# Summer11 QCD_Pt-30to50_MuPt5Enriched S3 MC
        'file:FEBD6E0A-967F-E011-819A-0024E87699A6.root'

#         '/store/mc/Summer11/QCD_Pt-80to120_TuneZ2_7TeV_pythia6/AODSIM/PU_S3_START42_V11-v2/0002/FA97B0C7-DB7D-E011-990A-485B39800C14.root',
#         '/store/mc/Summer11/QCD_Pt-80to120_TuneZ2_7TeV_pythia6/AODSIM/PU_S3_START42_V11-v2/0002/FA5C01BD-E47D-E011-87CF-485B39800BCA.root',
#         '/store/mc/Summer11/QCD_Pt-80to120_TuneZ2_7TeV_pythia6/AODSIM/PU_S3_START42_V11-v2/0002/FA49A8F9-E17D-E011-B8DD-E0CB4E1A118E.root',
#         '/store/mc/Summer11/QCD_Pt-80to120_TuneZ2_7TeV_pythia6/AODSIM/PU_S3_START42_V11-v2/0002/F8E63805-E27D-E011-A7F2-485B39800C0A.root',
#         '/store/mc/Summer11/QCD_Pt-80to120_TuneZ2_7TeV_pythia6/AODSIM/PU_S3_START42_V11-v2/0002/F86EC6F9-E17D-E011-991B-485B39800BA0.root',
  )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START42_V12::All" 

##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#$$
##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')
##-------------------- Turn-on the FastJet density calculation -----------------------
process.kt6PFJets.doRhoFastjet = True
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
process.ak5PFJets.doAreaFastjet = True
#$$

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("RecoBTau.JetTagComputer.jetTagRecord_cfi")

process.load("SimTracker.TrackHistory.TrackHistory_cff")

# for Impact Parameter based taggers
process.load("RecoBTag.ImpactParameter.negativeOnlyJetBProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetBProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D2ndComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D3rdComputer_cfi")
process.load("RecoBTag.Configuration.RecoBTag_cff")
process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetBProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetBProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCountingHighPur_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCountingHighEffJetTags_cfi")
process.load("RecoBTag.ImpactParameter.jetProbabilityBJetTags_cfi")
process.load("RecoBTag.ImpactParameter.jetBProbabilityBJetTags_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# for Secondary Vertex taggers
process.load("RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi")
process.load("RecoBTag.SecondaryVertex.secondaryVertexNegativeTagInfos_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighEffBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighPurBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexNegativeHighEffBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexNegativeHighPurBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeES_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexPositiveBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexPositiveES_cfi")

# for Soft Muon tagger
process.load("RecoBTag.SoftLepton.negativeSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftMuonBJetTags_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonBJetTags_cfi")

process.load("RecoBTag.SoftLepton.negativeSoftLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftMuonByPtBJetTags_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonByPtBJetTags_cfi")

process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")  

process.load("SimTracker.TrackHistory.TrackClassifier_cff")
process.load("RecoBTag.PerformanceMeasurements.MistagAnalyzer_cff")


process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 
        'keep recoJetTags_*_*_*'),
    fileName = cms.untracked.string('testtt.root')
)


#-------------------------------------
#Produce PFJets from PF2PAT

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

# process.load("CommonTools.ParticleFlow.Sources.source_ZtoMus_DBS_cfi")
runOnMC = True

from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PF2PAT"  # to have only PF2PAT
jetAlgo="AK5"

usePFnoPU = True # before any top projection

# levels to be accessible from the jets
# jets are corrected to L3Absolute (MC), L2L3Residual (data) automatically, if enabled here
# and remain uncorrected, if none of these levels is enabled here
useL1FastJet    = True  # needs useL1Offset being off, error otherwise
useL1Offset     = False # needs useL1FastJet being off, error otherwise
useL2Relative   = True
useL3Absolute   = True
# useL2L3Residual = True  # takes effect only on data; currently disabled for CMSSW_4_2_X GlobalTags!
useL5Flavor     = False
useL7Parton     = False

# JEC set
jecSetBase = jetAlgo
jecSetPF = jecSetBase + 'PF'
if usePFnoPU:
  jecSetPF += 'chs'

# JEC levels
if useL1FastJet and useL1Offset:
  sys.exit( 'ERROR: switch off either "L1FastJet" or "L1Offset"' )
jecLevels = []
if useL1FastJet:
  jecLevels.append( 'L1FastJet' )
if useL1Offset:
  jecLevels.append( 'L1Offset' )
if useL2Relative:
  jecLevels.append( 'L2Relative' )
if useL3Absolute:
  jecLevels.append( 'L3Absolute' )
# if useL2L3Residual and not runOnMC:
#   jecLevelsPF.append( 'L2L3Residual' )
if useL5Flavor:
  jecLevels.append( 'L5Flavor' )
if useL7Parton:
  jecLevels.append( 'L7Parton' )

usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix) 

applyPostfix(process,"patJetCorrFactors",postfix).levels     = cms.vstring('L1FastJet','L2Relative','L3Absolute')
applyPostfix(process,"patJetCorrFactors",postfix).rho        = cms.InputTag("kt6PFJets","rho")
applyPostfix(process,"patJetCorrFactors",postfix).payload    = cms.string('AK5PFchs')
applyPostfix(process,"pfPileUp",postfix).checkClosestZVertex = cms.bool(False) 

from PhysicsTools.PatAlgos.tools.metTools import *

addTcMET(process, 'TC')
process.patMETsTC.addGenMET = False

process.pfPileUpPF2PAT.Enable = True
#process.pfPileUpPF2PAT.checkClosestZVertex = cms.bool(False)
process.pfPileUpPF2PAT.Vertices = cms.InputTag('goodOfflinePrimaryVertices')
process.pfJetsPF2PAT.doAreaFastjet = True
process.pfJetsPF2PAT.doRhoFastjet = False

process.load('RecoJets.Configuration.RecoJets_cff')
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets

process.kt6PFJets               = kt4PFJets.clone()
process.kt6PFJets.rParam        = 0.6     
#process.kt6PFJets.src           = cms.InputTag('pfNoElectron'+postfix)
process.kt6PFJets.Rho_EtaMax    = cms.double( 4.4)
process.kt6PFJets.doRhoFastjet  = True
process.kt6PFJets.doAreaFastjet = True
#process.kt6PFJets.voronoiRfact  = 0.9

#process.patJetCorrFactorsPFlow.rho = cms.InputTag("kt6PFJets", "rho")
getattr(process,"patJetCorrFactors"+postfix).rho = cms.InputTag("kt6PFJets", "rho")

#-------------------------------------
#Redo the primary vertex

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )
    
#-------------------------------------
#JetID
process.load('PhysicsTools.SelectorUtils.pfJetIDSelector_cfi')
process.load('PhysicsTools.SelectorUtils.jetIDSelector_cfi')
process.jetIDSelector.version = cms.string('PURE09')


#-------------------------------------
#Filter for PFJets
process.PFJetsFilter = cms.EDFilter("PFJetSelector",
 src = cms.InputTag("selectedPatJetsPF2PAT"),
 cut = cms.string("pt > 10.0 && abs(eta) < 2.5 && neutralHadronEnergyFraction < 0.99 && neutralEmEnergyFraction < 0.99 && nConstituents > 1 && chargedHadronEnergyFraction > 0.0 && chargedMultiplicity > 0.0 && chargedEmEnergyFraction < 0.99"),
 filter = cms.bool(True)
)
#---------------------------------------
process.load("bTag.CommissioningCommonSetup.caloJetIDFilter_cff")


#Filter for removing scraping events
process.noscraping = cms.EDFilter("FilterOutScraping",
                               applyfilter = cms.untracked.bool(True),
                               debugOn = cms.untracked.bool(False),
                               numtrack = cms.untracked.uint32(10),
                               thresh = cms.untracked.double(0.25)
                               )


#Filter for good primary vertex
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                      minimumNDOF = cms.uint32(4) ,
 		      maxAbsZ = cms.double(24), 
 		      maxd0 = cms.double(2)	
                     )


### from Cristina JP calibration for cmsRun only : 
# from CondCore.DBCommon.CondDBCommon_cfi import *
# process.load("RecoBTag.TrackProbability.trackProbabilityFakeCond_cfi")
# process.trackProbabilityFakeCond.connect =cms.string( "sqlite_fip:RecoBTag/PerformanceMeasurements/test/btagnew_MC_414_2011.db")
# process.es_prefer_trackProbabilityFakeCond = cms.ESPrefer("PoolDBESSource","trackProbabilityFakeCond")

### from Cristina JP calibration for crab only: 
process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
       tag = cms.string("TrackProbabilityCalibration_2D_2011_v1_mc"),
       connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU")),
  cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
       tag = cms.string("TrackProbabilityCalibration_3D_2011_v1_mc"),
       connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU"))
)


#---------------------------------------
# process.TFileService = cms.Service("TFileService", fileName = cms.string("TrackTree.root") )
process.TFileService = cms.Service("TFileService", fileName = cms.string("JetTree.root") )
 
process.mistag.isData              = False
process.mistag.useTrackHistory     = False
process.mistag.produceJetProbaTree = False
process.mistag.triggerTable = 'TriggerResults::HLT' # Data and MC
#---------------------------------------


process.AK5byRef.jets = "selectedPatJetsPF2PAT"

# process.mistag.Jets = 'ak5PFJets'
process.mistag.Jets = 'selectedPatJetsPF2PAT'
#process.mistag.jetCorrector = cms.string('ak5PFL2L3')
process.mistag.jetCorrector = cms.string('ak5PFL1FastL2L3')
#process.mistag.jetCorrector = cms.string('ak5PFL1FastL2L3Residual')

# process.ak5JetTracksAssociatorAtVertex.jets = "ak5PFJets"
process.ak5JetTracksAssociatorAtVertex.jets = "selectedPatJetsPF2PAT"
# process.softMuonTagInfos.jets = "ak5PFJets"
process.softMuonTagInfos.jets = "selectedPatJetsPF2PAT"


process.jetProbabilityMixed = cms.ESProducer("JetProbabilityESProducer",
    impactParameterType = cms.int32(0), ## 0 = 3D, 1 = 2D

    deltaR = cms.double(0.3),
    maximumDistanceToJetAxis = cms.double(0.07),
    trackIpSign = cms.int32(0), ## 0 = use both, 1 = positive only, -1 = negative only

    minimumProbability = cms.double(0.005),
    maximumDecayLength = cms.double(5.0),
    trackQualityClass = cms.string("any")
)


#-------------------------------------
#Jet Probability
process.jetBProbabilityMixed = cms.ESProducer("JetBProbabilityESProducer",
    impactParameterType = cms.int32(0), ## 0 = 3D, 1 = 2D

    deltaR = cms.double(-1.0), ## use cut from JTA

    maximumDistanceToJetAxis = cms.double(0.07),
    trackIpSign = cms.int32(0), ## 0 = use both, 1 = positive only, -1 = negative only

    minimumProbability = cms.double(0.005),
    numberOfBTracks = cms.uint32(4),
    maximumDecayLength = cms.double(5.0),

    trackQualityClass = cms.string("any")
)

process.jetProbabilityBJetTags.jetTagComputer  = 'jetProbabilityMixed'
process.jetBProbabilityBJetTags.jetTagComputer = 'jetBProbabilityMixed'


#-------------------------------------
process.p = cms.Path(
#$$
        process.kt6PFJets
        *process.offlinePrimaryVertices 
	*process.goodOfflinePrimaryVertices
	*getattr(process,"patPF2PATSequence"+postfix)
        #*process.PFJetsFilter
	*process.myPartons*process.AK5Flavour
#$$
        *process.noscraping
        *process.primaryVertexFilter
	*process.ak5JetTracksAssociatorAtVertex
	*process.btagging
	*process.positiveOnlyJetProbabilityJetTags*process.negativeOnlyJetProbabilityJetTags
	*process.negativeTrackCountingHighEffJetTags*process.negativeTrackCountingHighPur
	*process.secondaryVertexTagInfos*process.simpleSecondaryVertexHighEffBJetTags*process.simpleSecondaryVertexHighPurBJetTags
	*process.secondaryVertexNegativeTagInfos*process.simpleSecondaryVertexNegativeHighEffBJetTags*process.simpleSecondaryVertexNegativeHighPurBJetTags
        *process.combinedSecondaryVertexNegativeBJetTags*process.combinedSecondaryVertexPositiveBJetTags
	#*process.negativeSoftMuonBJetTags*process.positiveSoftMuonBJetTags	
	*process.negativeSoftLeptonByPtBJetTags*process.positiveSoftLeptonByPtBJetTags	
	*process.negativeOnlyJetBProbabilityJetTags*process.positiveOnlyJetBProbabilityJetTags
	*process.mistag
	)
	
	## Output Module Configuration (expects a path 'p')




