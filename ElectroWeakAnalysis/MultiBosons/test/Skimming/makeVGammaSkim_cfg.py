###############################################################################
## This is the VGamma PAT configuration + Skim for Spring10 MC full sim samples
###############################################################################

## TO-DO List
## * Do we want to use VBTF WLNu candidates instead of ours leptonPlusMETs?
##   They have the acop method.

import copy
import FWCore.ParameterSet.Config as cms
import ElectroWeakAnalysis.MultiBosons.Skimming.VgEventContent as vgEventContent

from PhysicsTools.PatAlgos.patTemplate_cfg import * ## Default PAT process skeleton
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.trigTools import *
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from ElectroWeakAnalysis.MultiBosons.Skimming.egammaUserDataProducts_cff import *
from ElectroWeakAnalysis.MultiBosons.Skimming.jobOptions import *
from ElectroWeakAnalysis.MultiBosons.Skimming.matchHltPaths import matchHltPaths
from ElectroWeakAnalysis.MultiBosons.Skimming.options import options as defaultOptions
from ElectroWeakAnalysis.MultiBosons.tools.skimmingTools import embedTriggerMatches

## See link below for the definition of the selection
## https://twiki.cern.ch/twiki/bin/view/CMS/VGammaFirstPaper#Vgamma_Group_skims
skimVersion = 3  # Do we need this?
basePath = "ElectroWeakAnalysis.MultiBosons.Skimming." # shorthand

options = copy.deepcopy(defaultOptions)

## Define default options specific to this configuration file

#options.jobType = "testSummer10"
#options.jobType = "testPOWHEG"
#options.jobType = "testRealData"
options.jobType = "testJetRealData"
# options.jobType = "testMC"


## Parse (command-line) arguments - this overrides the options given above
# options.parseArguments()
applyJobOptions(options)

## Input
process.source.fileNames = options.inputFiles
if options.isMaxEventsOutput:
  process.maxEvents = cms.untracked.PSet(
    output = cms.untracked.int32(options.maxEvents)
    )
else:
  process.maxEvents.input = options.maxEvents

## Global tag
process.GlobalTag.globaltag = options.globalTag


## Remove MC matching and apply cleaning if we run on data
## (No need to remove pfMET and tcMET explicitly if this is done first
if options.isRealData:
  removeMCMatching(process)

## Add non-default MET flavors
addPfMET(process , 'PF')
addTcMET(process , 'TC')

## Spring10 MC was produced with CMSS_3_5_6 - make sure we can run on it
if options.use35XInput:
  run36xOn35xInput(process, genJets = "ak5GenJets")

## Embed lepton tracks
process.patMuons.embedTrack = True
process.patElectrons.embedTrack = True

## Keep only global and tracker muons
process.selectedPatMuons.cut = cms.string("isGlobalMuon | isTrackerMuon")

## Reject soft jets to reduce event content
process.selectedPatJets.cut = cms.string("pt > 30")

## Overlap of photons and electrons
##+ Set to True to remove overlapping photons from the event, False to embed
##+ embed the overlap information
process.cleanPatPhotons.checkOverlaps.electrons.requireNoOverlaps = False

## Add photon user data
process.load(basePath + "photonUserData_cfi")
process.patDefaultSequence.replace(process.patPhotons,
  process.photonUserData * process.patPhotons
  )
process.patPhotons.userData.userFloats.src = egammaUserDataFloats(
  moduleName = "photonUserData"
  )
process.patPhotons.userData.userInts.src = egammaUserDataInts(
  moduleName = "photonUserData"
  )

process.load("RecoEcal.EgammaClusterProducers.preshowerClusterShape_cfi")
process.load("EgammaAnalysis.PhotonIDProducers.piZeroDiscriminators_cfi")

process.load("ElectroWeakAnalysis.MultiBosons.Skimming.pi0Discriminator_cfi")
process.patDefaultSequence.replace(process.patPhotons,
  process.preshowerClusterShape *
  process.piZeroDiscriminators  *
  process.pi0Discriminator      *
  process.patPhotons
  )
process.patPhotons.userData.userFloats.src.append(
  cms.InputTag("pi0Discriminator") #, "piZeroDiscriminatorsPhotonPi0DiscriminatorAssociationMap")
  )

## Add electron user data
process.load(basePath + "electronUserData_cfi")
process.patDefaultSequence.replace(process.patElectrons,
  process.electronUserData * process.patElectrons
  )
process.patElectrons.userData.userFloats.src = egammaUserDataFloats(
  moduleName = "electronUserData"
  )
process.patElectrons.userData.userInts.src = egammaUserDataInts(
  moduleName = "electronUserData"
  )

## Add electron official electron ID from the Egamma / VBTF
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
process.patDefaultSequence.replace(process.patElectrons,
  process.simpleEleIdSequence * process.patElectrons
  )
process.patElectrons.addElectronID = cms.bool(True)
process.patElectrons.electronIDSources = cms.PSet(
  simpleEleId95relIso = cms.InputTag("simpleEleId95relIso"),
  simpleEleId90relIso = cms.InputTag("simpleEleId90relIso"),
  simpleEleId85relIso = cms.InputTag("simpleEleId85relIso"),
  simpleEleId80relIso = cms.InputTag("simpleEleId80relIso"),
  simpleEleId70relIso = cms.InputTag("simpleEleId70relIso"),
  simpleEleId60relIso = cms.InputTag("simpleEleId60relIso"),
  simpleEleId95cIso   = cms.InputTag("simpleEleId95cIso"),
  simpleEleId90cIso   = cms.InputTag("simpleEleId90cIso"),
  simpleEleId85cIso   = cms.InputTag("simpleEleId85cIso"),
  simpleEleId80cIso   = cms.InputTag("simpleEleId80cIso"),
  simpleEleId70cIso   = cms.InputTag("simpleEleId70cIso"),
  simpleEleId60cIso   = cms.InputTag("simpleEleId60cIso"),
  )

## PAT Trigger
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
switchOnTrigger(process)
process.patTrigger.processName = options.hltProcessName
process.patTriggerEvent.processName = options.hltProcessName
embedTriggerMatches(process, matchHltPaths)
## Drop matched target collections from the event content to only keep the trigger matched versions
for collection in matchHltPaths.keys():
  vgEventContent.extraSkimEventContent.append("drop *_%s_*_*" % collection)

## HLT trigger
process.load(basePath + "hltFilter_cfi")
process.hltFilter.HLTPaths = options.hltPaths
process.hltFilter.TriggerResultsTag = \
  "TriggerResults::" + options.hltProcessName

## Define the path that's used to select events
process.skimFilterSequence = cms.Sequence(process.hltFilter) # Extend below
process.skimFilterPath = cms.Path(process.skimFilterSequence)

if options.skimType == "MuonPhoton":
    process.hltFilter.HLTPaths += ["HLT_Mu9", "HLT_Mu11"]
    process.load(basePath + "muonPhotonSkimFilterSequence_cff")
    process.skimFilterSequence += process.muonPhotonSkimFilterSequence

elif options.skimType == "ElectronPhoton":
    process.hltFilter.HLTPaths += ["HLT_Ele15_LW_L1R", "HLT_Ele15_SW_L1R"]
    process.load(basePath + "electronPhotonSkimFilterSequence_cff")
    process.skimFilterSequence += process.electronPhotonSkimFilterSequence

elif options.skimType == "Dimuon":
    process.hltFilter.HLTPaths += ["HLT_Mu9", "HLT_Mu11", "HLT_DoubleMu3"]
    process.load(basePath + "dimuonSkimFilterSequence_cff")
    process.skimFilterSequence += process.dimuonSkimFilterSequence
elif options.skimType == "Jet":
  process.load('Configuration.StandardSequences.Services_cff')
  process.load('Configuration.StandardSequences.MagneticField_38T_cff')
  process.load('Configuration.StandardSequences.Geometry_cff')
  process.load('Configuration.StandardSequences.Reconstruction_cff')
  #process.load("RecoEgamma.Configuration.RecoEgamma_cff")
  process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
  process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

  process.photonReReco = cms.Sequence(process.conversionSequence*
                                      process.trackerOnlyConversionSequence*
                                      process.photonSequence*
                                      process.photonIDSequence)

  process.hltFilter.HLTPaths += ["HLT_Jet100U","HLT_Jet70U","HLT_Jet50U","HLT_Jet30U","HLT_Jet15U"]
  process.load(basePath + "jetSkimFilterSequence_cff")
  process.skimFilterSequence += process.jetSkimFilterSequence
  # now change the photon reco to much looser settings
  process.photonCore.minSCEt = 10.0
  process.photons.minSCEtBarrel = 10.0
  process.photons.minSCEtEndcap = 10.0
  process.photons.maxHoverEBarrel = 10.0
  process.photons.maxHoverEEndcap = 10.0
  #edit the pat sequence to do the rereco
  process.patDefaultSequence = cms.Sequence(process.photonReReco+process.patDefaultSequence)

else:
    raise RuntimeError, "Illegal skimType option: %s" % options.skimType

process.load(basePath + "VGammaSkimSequences_cff")

## Add cleaning of collision data (no scraping events etc.)
##+ https://twiki.cern.ch/twiki/bin/viewauth/CMS/Collisions2010Recipes
process.load(basePath + "goodCollisionDataSequence_cff")

if options.isRealData:
  ## Remove the hltPhysicsDeclared - it kills some good events, reference?
  process.goodCollisionDataSequence.remove("hltPhysicsDeclared")
  ## Run the hltPhysicsDeclared filter in a separate path to
  ##+ store its result in the triggerEvent product.
  process.hltPhysicsDeclaredPath = cms.Path(process.hltPhysicsDeclared)
  process.skimFilterSequence = cms.Sequence(
    process.goodCollisionDataSequence +
    process.skimFilterSequence
    )
  process.defaultSequence = cms.Sequence(
    process.skimFilterSequence +
    process.patDefaultSequence
  )
else:
  process.primaryVertexFilterPath = cms.Path(process.primaryVertexFilter)
  process.load(basePath + "prunedGenParticles_cfi")
  process.defaultSequence = cms.Sequence(
    process.skimFilterSequence +
    process.prunedGenParticles *
    process.patDefaultSequence
  )
  ## Add parton shower related filters (prevent ISR/FSR double-counting)
  process.isLeadingPhotonPythiaPartonShowerIsr = cms.Path(
    process.defaultSequence *
    process.pythiaPartonShowerIsrSequence
  )
  process.isLeadingPhotonPythiaPartonShowerFsr = cms.Path(
    process.defaultSequence *
    process.pythiaPartonShowerFsrSequence
  )
  process.hasPhotonCandidateNotPythiaPartonShower = cms.Path(
    process.defaultSequence *
    process.pythiaPartonShowerPhotonVeto
  )
# if options.isRealData <-----------------------------------------------------

process.WENuGammaPath  = cms.Path(
  process.defaultSequence * process.WENuGammaSequence
  )
process.WMuNuGammaPath = cms.Path(
  process.defaultSequence * process.WMuNuGammaSequence
  )
process.ZEEGammaPath   = cms.Path(
  process.defaultSequence * process.ZEEGammaSequence
  )
process.ZMuMuGammaPath = cms.Path(
  process.defaultSequence * process.ZMuMuGammaSequence
  )
process.ZInvisibleGammaPath = cms.Path(
  process.defaultSequence * process.ZInvisibleGammaSequence
  )

## Output configuration (add event content, select events, output file name)
process.out.outputCommands += vgEventContent.extraSkimEventContent

if not options.isRealData:
  process.out.outputCommands += ["keep *_prunedGenParticles_*_PAT"]

if options.skimType == "Jet":
  process.out.outputCommands += ["drop *_photonCore_*_RECO",
                                 "drop *_conversions_*_RECO",
                                 "drop *_trackerOnlyConversions_*_RECO",
                                 "drop *_ckf*TracksFromConversions_*_RECO"]

process.out.SelectEvents.SelectEvents = ["skimFilterPath"]
process.out.fileName = options.outputFile

## Logging
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
if not options.isRealData:
  ## Suppress many warnings about missing HLT prescale tables
  process.MessageLogger.categories += ["hltPrescaleTable"]
  process.MessageLogger.cerr.hltPrescaleTable = cms.untracked.PSet(
    limit = cms.untracked.int32(5)
    )

process.options.wantSummary = options.wantSummary

## Add tab completion + history during inspection
if __name__ == "__main__": import user
