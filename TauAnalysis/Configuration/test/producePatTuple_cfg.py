import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('producePatTupleZtoElecTau')

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR10_P_V5::All'


#--------------------------------------------------------------------------------
# import sequence for PAT-tuple production
process.load("TauAnalysis.Configuration.producePatTuple_cff")

# define configuration parameters for submission of jobs to CERN batch system 
# (running over skimmed samples stored on CASTOR)
#from TauAnalysis.Configuration.recoSampleDefinitionsAHtoElecMu_cfi import *
#from TauAnalysis.Configuration.recoSampleDefinitionsWtoTauNu_cfi import *
#from TauAnalysis.Configuration.recoSampleDefinitionsZtoElecMu_cfi import *
from TauAnalysis.Configuration.recoSampleDefinitionsZtoElecTau_cfi import *
from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_10TeV_cfi import *
from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_cfi import *
#from TauAnalysis.BgEstimationTools.recoSampleDefinitionsTauIdEffZtoMuTau_cfi import *

# import event-content definition of products to be stored in patTuple
from TauAnalysis.Configuration.patTupleEventContent_cff import *
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# print event content 
process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")
#--------------------------------------------------------------------------------

process.savePatTuple = cms.OutputModule("PoolOutputModule",
	patTupleEventContent,
	fileName = cms.untracked.string('patTuple.root')
)

process.maxEvents = cms.untracked.PSet(            
	input = cms.untracked.int32(10)    
)

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
		'rfio:/castor/cern.ch/user/j/jkolb/eTauSkims/data/EG_SD/skimElecTau_24_0.root'
		#'rfio:/castor/cern.ch/user/j/jkolb/eTauSkims/Summer09_CMSSW_3_1_4/Ztautau_7TeV/skimElecTau_Ztautau_7TeV_17.root'
	)
	#skipBadFiles = cms.untracked.bool(True)    
)

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system
#
#__process.source.fileNames = #inputFileNames#
#__process.maxEvents.input = cms.untracked.int32(#maxEvents#)
#__process.savePatTuple.fileName = #patTupleOutputFileName#
#
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# add patElectronIsolation, which was removed from standard pat sequence in CMSSW_3_4
from TauAnalysis.Configuration.tools.producePatElectronIsolation import *
producePatElectronIsolation(process)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import function to remove MC-specific modules 
from TauAnalysis.Configuration.tools.switchToData import switchToData
# uncomment to run over data
##switchToData(process)#
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for switching pat::Tau input
# to different reco::Tau collection stored on AOD
from PhysicsTools.PatAlgos.tools.tauTools import * 

# comment-out to take reco::CaloTaus instead of reco::PFTaus
# as input for pat::Tau production
#switchToCaloTau(process)

# comment-out to take shrinking dR = 5.0/Et(PFTau) signal cone
# instead of fixed dR = 0.07 signal cone reco::PFTaus
# as input for pat::Tau production
switchToPFTauShrinkingCone(process)
#switchToPFTauFixedCone(process)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::Jets
from PhysicsTools.PatAlgos.tools.jetTools import *

# switchJetCollection complains if this doesn't exist
process.jetCorrFactors = cms.PSet()

# uncomment to replace caloJets by pfJets
switchJetCollection(process, cms.InputTag("iterativeCone5PFJets") )
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::METs
from TauAnalysis.Configuration.tools.metTools import *

# uncomment to add pfMET
# (set boolean parameter to true/false to enable/disable type-1 MET corrections)
addPFMet(process, False)

# uncomment to replace caloMET by pfMET in all di-tau objects
process.load("TauAnalysis.CandidateTools.diTauPairProductionAllKinds_cff")
replaceMETforDiTaus(process, cms.InputTag('layer1METs'), cms.InputTag('layer1PFMETs'))

from TauAnalysis.Configuration.tools.sysUncertaintyTools import *
# uncomment to disable produceSysErrGenEventReweights sequence from PAT post-production
disableSysUncertainties_patTupleProduction(process)
#--------------------------------------------------------------------------------

process.p = cms.Path(
		process.producePatTuple
		#  + process.printEventContent      
		+ process.savePatTuple
)

#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True)
#)

# print-out all python configuration parameter information
#print process.dumpPython()
