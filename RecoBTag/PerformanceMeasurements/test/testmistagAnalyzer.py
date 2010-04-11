import FWCore.ParameterSet.Config as cms

from SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi import *
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from JetMETCorrections.Configuration.JetCorrectionsRecord_cfi import *
from RecoBTag.ImpactParameter.impactParameter_cfi import *
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
process = cms.Process("Demo")


process.source = cms.Source(
    "PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(   

    '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/F4C92A98-163C-DF11-9788-0030487C7392.root',
    '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/F427D642-173C-DF11-A909-0030487C60AE.root',
    '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/E27821C3-0C3C-DF11-9BD9-0030487CD718.root',
    '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/D87D5469-2E3C-DF11-A470-000423D99896.root',
    '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/B647CAD9-0E3C-DF11-886F-0030487CD716.root',
    '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/A860D55E-193C-DF11-BE29-0030487C60AE.root',
    '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/9884BC11-0C3C-DF11-8F9C-000423D986C4.root',
    '/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/92684831-233C-DF11-ABA0-0030487CD16E.root'
    
    )
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)


#process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "MC_31X_V5::All" # use a valid global tag here!
process.GlobalTag.globaltag = "GR10_P_V4::All" # use a valid global tag here!
#process.GlobalTag.globaltag = "START3X_V25B::All" # use a valid global tag here!


process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("RecoBTau.JetTagComputer.jetTagRecord_cfi")

process.load("SimTracker.TrackHistory.TrackHistory_cff")


# for Impact Parameter based taggers
process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D2ndComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D3rdComputer_cfi")
process.load("RecoBTag.Configuration.RecoBTag_cff")
process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")


process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCountingHighPur_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCountingHighEffJetTags_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load("JetMETCorrections.Configuration.MCJetCorrections152_cff")

process.load("RecoBTag.ImpactParameter.jetProbabilityBJetTags_cfi")


#for Secondary Vertex taggers
process.load("RecoBTag.SecondaryVertex.secondaryVertexNegativeTagInfos_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexNegativeBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeES_cfi")


#For Soft Muon tagger

process.load("RecoBTag.SoftLepton.negativeSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftMuonBJetTags_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonBJetTags_cfi")

process.load("RecoBTag.SoftLepton.negativeSoftLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftMuonByPtBJetTags_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonByPtBJetTags_cfi")




#process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")

# process.load("SimTracker.TrackHistory.TrackHistory_cff")

# process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")


process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")  

#############   Include the jet corrections ##########
process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer09_7TeV_cff")
# set the record's IOV. Must be defined once. Choose ANY correction service. #
#process.prefer("L2L3JetCorrectorIC5Calo") 
process.prefer("L2L3JetCorrectorAK5PF") 



process.load("SimTracker.TrackHistory.TrackClassifier_cff")
process.load("RecoBTag.PerformanceMeasurements.MistagAnalyzer_cff")


#process.mistag.useTrackHistory = cms.bool(False)
#mistag.jetCorrector    = cms.string('L2L3JetCorrectorIC5Calo')
#process.mcAlgoJetFlavour = cms.Sequence(
#     process.myPartons *
#     process.IC5byRef *
#     process.IC5byValAlgo
# )



process.TFileService = cms.Service("TFileService", fileName = cms.string("analyzerTree.root") )
 

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 
        'keep recoJetTags_*_*_*'),
    fileName = cms.untracked.string('testtt.root')
)

#process.p = cms.Path(process.positiveOnlyJetProbabilityJetTags*process.negativeOnlyJetProbabilityJetTags*process.negativeTrackCountingHigEffJetTags*process.negativeTrackCountingHigPur*process.mistag)

#process.p = cms.Path(process.mcAlgoJetFlavour*process.positiveOnlyJetProbabilityJetTags*process.negativeOnlyJetProbabilityJetTags
#	*process.negativeTrackCountingHighEffJetTags*process.negativeTrackCountingHighPur
#	*process.secondaryVertexNegativeTagInfos*process.simpleSecondaryVertexNegativeBJetTags
#	*process.negativeSoftMuonBJetTags*process.positiveSoftMuonBJetTags	
#	*process.jetProbabilityBJetTags*process.mistag
#	)

process.mistag.isData = True
process.mistag.useTrackHistory = False
process.mistag.Jets = 'ak5PFJets'
process.mistag.jetCorrector = cms.string('L2L3JetCorrectorAK5PF')
process.ak5JetTracksAssociatorAtVertex.jets = "ak5PFJets"

process.softMuonTagInfos.jets = "ak5PFJets"


process.p = cms.Path(
        #*process.myPartons*process.AK5Flavour
	process.ak5JetTracksAssociatorAtVertex
	*process.btagging*process.positiveOnlyJetProbabilityJetTags
	*process.negativeOnlyJetProbabilityJetTags
	*process.negativeTrackCountingHighEffJetTags*process.negativeTrackCountingHighPur
	*process.secondaryVertexNegativeTagInfos*process.simpleSecondaryVertexNegativeBJetTags*process.combinedSecondaryVertexNegativeBJetTags
	#*process.negativeSoftMuonBJetTags*process.positiveSoftMuonBJetTags	
	*process.negativeSoftLeptonByPtBJetTags*process.positiveSoftLeptonByPtBJetTags	
	*process.jetProbabilityBJetTags
	*process.mistag
	)

