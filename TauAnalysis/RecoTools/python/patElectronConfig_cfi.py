import FWCore.ParameterSet.Config as cms
import copy

from PhysicsTools.PatAlgos.recoLayer0.electronId_cff import *
from PhysicsTools.PatAlgos.recoLayer0.electronIsolation_cff import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerMatcher_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerMatchEmbedder_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi import *
from PhysicsTools.PatAlgos.cleaningLayer1.electronCleaner_cfi import *

#--------------------------------------------------------------------------------  
# PAT layer 0 electron configuration parameters
#--------------------------------------------------------------------------------

# add HLT electron trigger to PAT trigger match sequence
patTriggerElectronMatcher += electronTriggerMatchHLTIsoEle15LWL1I
patTriggerEvent.patTriggerMatches.append("electronTriggerMatchHLTIsoEle15LWL1I")
cleanLayer1ElectronsTriggerMatch.matches = cms.VInputTag("electronTriggerMatchHLTIsoEle15LWL1I")

#--------------------------------------------------------------------------------  
# PAT layer 1 electron configuration parameters
#--------------------------------------------------------------------------------

#
# Track isolation
#
allLayer1Electrons.userIsolation.tracker.src = cms.InputTag("eleIsoDepositTk")
allLayer1Electrons.userIsolation.tracker.deltaR = cms.double(0.6)
allLayer1Electrons.userIsolation.tracker.vetos = cms.vstring(
    '0.015',         # inner radius veto cone 
    'Threshold(0.3)' # threshold on individual track pt
)
allLayer1Electrons.userIsolation.tracker.skipDefaultVeto = cms.bool(True)
#
# ECAL isolation
#
allLayer1Electrons.userIsolation.ecal.src = cms.InputTag("eleIsoDepositEcalFromHits")
allLayer1Electrons.userIsolation.ecal.deltaR = cms.double(0.6)
allLayer1Electrons.userIsolation.ecal.vetos = cms.vstring(
    'EcalBarrel:0.045', 
    'EcalBarrel:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)',
    'EcalEndcaps:0.1',                          # default: 0.07
    'EcalEndcaps:RectangularEtaPhiVeto(-0.05,0.05,-0.5,0.5)',
    'EcalBarrel:ThresholdFromTransverse(0.12)', # default: 0.08
    'EcalEndcaps:ThresholdFromTransverse(0.3)'
)
allLayer1Electrons.userIsolation.ecal.skipDefaultVeto = cms.bool(True)
#
# HCAL isolation
#
allLayer1Electrons.userIsolation.hcal.src = cms.InputTag("eleIsoDepositHcalFromTowers")
allLayer1Electrons.userIsolation.hcal.deltaR = cms.double(0.6)
#
# add IsoDeposit objects for Track, ECAL and HCAL based isolation
#
allLayer1Electrons.isoDeposits = cms.PSet(
   tracker          = allLayer1Electrons.userIsolation.tracker.src,
   ecal             = allLayer1Electrons.userIsolation.ecal.src,
   hcal             = allLayer1Electrons.userIsolation.hcal.src,
   particle         = cms.InputTag("pfeleIsoDepositPFCandidates"),
   pfChargedHadrons = cms.InputTag("pfeleIsoChDepositPFCandidates"),
   pfNeutralHadrons = cms.InputTag("pfeleIsoNeDepositPFCandidates"),
   pfPhotons        = cms.InputTag("pfeleIsoGaDepositPFCandidates")
)
#
# add electron Id. flags
#
allLayer1Electrons.addElectronID = cms.bool(True)
#
# enable matching to HLT trigger information
#
allLayer1Electrons.embedHighLevelSelection = cms.bool(True)
#
# enable matching to generator level information
#
allLayer1Electrons.addGenMatch = cms.bool(True)

allLayer1Electrons.genParticleMatch = cms.InputTag("electronMatch")

# do not remove electrons overlapping with muons
# (instead, leave this removal to the subsequent selector stage)
cleanLayer1Electrons.checkOverlaps = cms.PSet()
