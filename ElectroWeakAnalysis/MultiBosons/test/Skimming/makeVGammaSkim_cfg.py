###############################################################################
## This is the VGamma PAT configuration + Skim for Spring10 MC full sim samples
###############################################################################

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
from ElectroWeakAnalysis.MultiBosons.Skimming.matchHltPaths import *
from ElectroWeakAnalysis.MultiBosons.Skimming.options import options as defaultOptions
from ElectroWeakAnalysis.MultiBosons.tools.skimmingTools import *


## See link below for the definition of the selection
## https://twiki.cern.ch/twiki/bin/view/CMS/VGammaFirstPaper#Vgamma_Group_skims
skimVersion = 3  # Do we need this?
basePath = "ElectroWeakAnalysis.MultiBosons.Skimming." # shorthand

options = copy.deepcopy(defaultOptions)

## Define default options specific to this configuration file
#options.jobType = "testSummer10"
#options.jobType = "testPOWHEG"
#options.jobType = "testRealData"
#options.jobType = "testJetRealData"
#options.jobType = "testMC"
#options.jobType = "testMC414"

## Parse (command-line) arguments - this overrides the options given above
applyJobOptions(options)

## Input
process.source.fileNames = options.inputFiles

process.maxEvents.input = -1
if options.outEvents >= 0:
    process.maxEvents.output = cms.untracked.int32(options.outEvents)
    #process.maxEvents = cms.untracked.PSet( output = cms.untracked.int32(options.outEvents) )
if options.maxEvents >= 0:
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
## TODO: rename photon(electron)UserData to photon(electron)ClusterShape
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

if not options.isRealData:
    process.photonGenMatch = cms.EDProducer("PhotonGenMatchUserDataProducer",
        src = cms.InputTag("photons"),
        match = cms.InputTag("photonMatch")
        )
    process.patDefaultSequence.replace(process.patPhotons,
        process.photonGenMatch *
        process.patPhotons
        )
    process.patPhotons.userData.userInts.src.extend([
        cms.InputTag("photonGenMatch", "motherPdgId"),
        cms.InputTag("photonGenMatch", "motherStatus"),
        cms.InputTag("photonGenMatch", "grandMotherPdgId"),
        cms.InputTag("photonGenMatch", "grandMotherStatus"),
        ])

if not options.isAOD:
    process.patDefaultSequence.replace(process.patPhotons,
        process.preshowerClusterShape *
        process.piZeroDiscriminators  *
        process.pi0Discriminator      *
        process.patPhotons
        )
else:
    process.patDefaultSequence.replace(process.patPhotons,
        process.piZeroDiscriminators  *
        process.pi0Discriminator      *
        process.patPhotons
        )
    process.piZeroDiscriminators.preshClusterShapeProducer = cms.string('multi5x5PreshowerClusterShape')
    process.piZeroDiscriminators.preshClusterShapeCollectionX = cms.string('multi5x5PreshowerXClustersShape')
    process.piZeroDiscriminators.preshClusterShapeCollectionY = cms.string('multi5x5PreshowerYClustersShape')
    


process.patPhotons.userData.userFloats.src.append(
    cms.InputTag("pi0Discriminator")
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



## Drop matched target collections from the event content to only keep the trigger matched versions
for collection in matchHltPaths.keys():
    vgEventContent.extraSkimEventContent.append("drop *_%s_*_*" % collection)
## Drop default patTriggerObjectsedmAssociation
vgEventContent.extraSkimEventContent.append(
    "drop patTriggerObjectsedmAssociation_patTriggerEvent_*_*"
    )


## HLT trigger
process.load(basePath + "hltFilter_cfi")
process.hltFilter.TriggerResultsTag = \
    "TriggerResults::" + options.hltProcessName

## Add cleaning of collision data (no scraping events etc.)
##+ https://twiki.cern.ch/twiki/bin/viewauth/CMS/Collisions2010Recipes
process.load(basePath + "goodCollisionDataSequence_cff")
## Remove the hltPhysicsDeclared - it kills some good events, reference?
process.goodCollisionDataSequence.remove(process.hltPhysicsDeclared)
## Run the hltPhysicsDeclared filter in a separate path to
##+ store its result in the triggerEvent product.
process.hltPhysicsDeclaredPath = cms.Path(process.hltPhysicsDeclared)

## Define the path that's used to select events
process.skimFilterSequence = cms.Sequence(
    process.hltFilter +
    process.goodCollisionDataSequence
    ) # Extend below
process.skimFilterPath = cms.Path(process.skimFilterSequence)

if options.skimType == "MuonPhoton":
    removeTriggerPathsForAllBut(matchHltPaths, ["cleanPatMuons"])
    process.hltFilter.HLTPaths = matchHltPaths["cleanPatMuons"]
    if not options.ignoreSkimFilter:
        process.load(basePath + "muonPhotonSkimFilterSequence_cff")
        process.skimFilterSequence += process.muonPhotonSkimFilterSequence

elif options.skimType == "ElectronPhoton":
    removeTriggerPathsForAllBut(matchHltPaths, ["cleanPatElectrons"])
    process.hltFilter.HLTPaths = matchHltPaths["cleanPatElectrons"]
    if not options.ignoreSkimFilter:
        process.load(basePath + "electronPhotonSkimFilterSequence_cff")
        process.skimFilterSequence += process.electronPhotonSkimFilterSequence

elif options.skimType == "Dimuon":
    removeTriggerPathsForAllBut(matchHltPaths, ["cleanPatMuons"])
    ## Don't require any triggers.
    process.skimFilterSequence.remove(process.hltFilter)
    ## Require the Muon PD
    # process.hltFilter.HLTPaths = ["*Mu*"]
    if not options.ignoreSkimFilter:
        process.load(basePath + "dimuonSkimFilterSequence_cff")
        process.skimFilterSequence += process.dimuonSkimFilterSequence
    ## Add the photon re-reco.
    addPhotonReReco(process)
    ## Now change the photon reco to much looser settings.
    process.photonCore.minSCEt = 2.0
    process.photons.minSCEtBarrel = 2.0
    process.photons.minSCEtEndcap = 2.0
    process.photons.maxHoverEBarrel = 10.0
    process.photons.maxHoverEEndcap = 10.0
    ## Remove the pi0 discriminator
    ## (currently doesn't work with extremely loose photons)
    ## FIXME: make the pi0Discriminator work for these weird photons too
    for module in [process.preshowerClusterShape,
                   process.piZeroDiscriminators,
                   process.pi0Discriminator]:
        process.patDefaultSequence.remove( module )
    while cms.InputTag("pi0Discriminator") in process.patPhotons.userData.userFloats.src:
        process.patPhotons.userData.userFloats.src.remove(
            cms.InputTag("pi0Discriminator")
            )

    ## Add more photon-related event content (super clusters, clusters)
    vgEventContent.extraSkimEventContent += \
        vgEventContent.vgExtraPhotonEventContent
    ## Add island basic clusters to the sequence and event content
    process.load("RecoEcal.EgammaClusterProducers.islandBasicClusters_cfi")
    process.patDefaultSequence = cms.Sequence(
       process.islandBasicClusters *
       process.patDefaultSequence
    )
    vgEventContent.extraSkimEventContent.append("keep *_islandBasicClusters_*_*")

elif options.skimType == "Jet":
    removeTriggerPathsForAllBut(matchHltPaths, ["cleanPatJets"])
    process.hltFilter.HLTPaths = matchHltPaths["cleanPatJets"]
    if not options.ignoreSkimFilter:
        process.load(basePath + "jetSkimFilterSequence_cff")
        process.skimFilterSequence += process.jetSkimFilterSequence
    addPhotonReReco(process)
    # now change the photon reco to much looser settings
    process.photonCore.minSCEt = 2.0
    process.photons.minSCEtBarrel = 2.0
    process.photons.minSCEtEndcap = 2.0
    process.photons.maxHoverEBarrel = 10.0
    process.photons.maxHoverEEndcap = 10.0
    #edit the pat sequence to do the rereco
    process.patDefaultSequence = cms.Sequence(process.photonReReco*process.patDefaultSequence)

elif options.skimType == "Inclusive":
    ## Remove all skimming filters
    process.skimFilterSequence = cms.Sequence()

else:
    raise RuntimeError, "Illegal skimType option: %s" % options.skimType

embedTriggerMatches(process, matchHltPaths)

process.load(basePath + "VGammaSkimSequences_cff")

if options.isRealData:
    process.defaultSequence = cms.Sequence(
        process.skimFilterSequence +
        process.patDefaultSequence
    )
else:
    process.primaryVertexFilterPath = cms.Path(process.primaryVertexFilter)
    process.noScrapingPath = cms.Path(process.noScraping)
    process.load(basePath + "prunedGenParticles_cfi")
    if not options.applyCollisionDataCleaningToMC \
        and repr(process.skimFilterSequence) != repr(cms.Sequence()):
        process.skimFilterSequence.remove(process.goodCollisionDataSequence)
    process.defaultSequence = cms.Sequence(
        process.skimFilterSequence +
        process.prunedGenParticles *
        process.patDefaultSequence
    )
    ## Add parton shower related filters (prevent ISR/FSR double-counting)
    ##+ "POWHEG parton shower" means actually Pythia parton shower in Pythia
    ##+ hadronization run on top of POWHEG.
    process.isLeadingPhotonPythiaPartonShowerIsr = cms.Path(
        process.defaultSequence *
        process.pythiaPartonShowerIsrSequence
    )
    process.isLeadingPhotonPythiaPartonShowerFsr = cms.Path(
        process.defaultSequence *
        process.pythiaPartonShowerFsrSequence
    )
    process.isLeadingPhotonPythiaPartonShowerIsrVeto = cms.Path(
        process.defaultSequence *
        process.pythiaPartonShowerIsrSequenceVeto
        )
    process.isLeadingPhotonPythiaPartonShowerFsrVeto = cms.Path(
        process.defaultSequence *
        process.pythiaPartonShowerFsrSequenceVeto
        )
    process.isLeadingPhotonPowhegPartonShowerFsr = cms.Path(
        process.defaultSequence *
        process.powhegPartonShowerFsrSequence
    )
    process.isLeadingPhotonPowhegPartonShowerFsrVeto = cms.Path(
        process.defaultSequence *
        process.powhegPartonShowerFsrSequenceVeto
        )

    process.hasPhotonCandidateNotPythiaPartonShower = cms.Path(
        process.defaultSequence *
        process.pythiaPartonShowerPhotonVeto
        )

    # construct path for select events based on options...
    process.partonShowerFilterPath = cms.Path( process.defaultSequence )

    if options.vetoISR:
        process.partonShowerFilterPath.replace(process.defaultSequence,
                                               process.defaultSequence*
                                               process.pythiaPartonShowerIsrSequenceVeto)

    if options.vetoFSR:
        if options.isPOWHEG:
            process.partonShowerFilterPath.replace(process.defaultSequence,
                                                   process.defaultSequence*
                                                   process.powhegPartonShowerFsrSequenceVeto)
        else:
            process.partonShowerFilterPath.replace(process.defaultSequence,
                                                   process.defaultSequence*
                                                   process.pythiaPartonShowerFsrSequenceVeto)


    if options.skimISR:
        process.partonShowerFilterPath.replace(process.defaultSequence,
                                               process.defaultSequence*
                                               process.pythiaPartonShowerIsrSequence)

    if options.skimFSR:
        if options.isPOWHEG:
            process.partonShowerFilterPath.replace(process.defaultSequence,
                                                   process.defaultSequence*
                                                   process.powhegPartonShowerFsrSequence)
        else:
            process.partonShowerFilterPath.replace(process.defaultSequence,
                                                   process.defaultSequence*
                                                   process.pythiaPartonShowerFsrSequence)
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

if not (options.vetoISR or options.vetoFSR or options.skimISR or options.skimFSR):
    process.out.SelectEvents.SelectEvents = ["skimFilterPath"]
else: #if we want to filter on I/FSR we must process the whole event first... sadface.
    process.out.SelectEvents.SelectEvents = ["partonShowerFilterPath"]

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

## Check for an empty path in the output
if str(process.skimFilterPath) == "None":
    del process.out.SelectEvents

## Test photon re-reco
# process.TFileService = cms.Service("TFileService",
#   fileName = cms.string("histo_withPhotonReReco.root")
# )
# process.load(basePath + "testPhotonReRecoSequence_cff")
# process.testPhotonReRecoPath = cms.Path(process.defaultSequence *
#                                         process.testPhotonReRecoSequence)


## Add tab completion + history during inspection
if __name__ == "__main__": import user
