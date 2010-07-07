import FWCore.ParameterSet.Config as cms

# Primary vertex smearing.
from IOMC.EventVertexGenerators.VtxSmearedGauss_cfi import *

# Conversion to GenParticleCandidates 
from PhysicsTools.HepMCCandAlgos.genParticleCandidatesFast_cfi import *

# Famos PileUp Producer
from FastSimulation.PileUpProducer.PileUpProducer_cff import *

# Famos SimHits producer
from FastSimulation.EventProducer.FamosSimHits_cff import *

# Mixing module 
from FastSimulation.Configuration.mixNoPU_cfi import *

# Gaussian Smearing RecHit producer
from FastSimulation.TrackingRecHitProducer.SiTrackerGaussianSmearingRecHitConverter_cfi import *

# Rec Hit Tranlator to the Full map with DeTId'
from FastSimulation.TrackingRecHitProducer.TrackingRecHitTranslator_cfi import *

# CTF and Iterative tracking (contains pixelTracks and pixelVertices)

# 1) Common algorithms and configuration taken from full reconstruction
# Note: The runge-kutta propagator is not used here 
# (because no magnetic field inhomogeneities are simulated between layers)
from FastSimulation.Tracking.GSTrackFinalFitCommon_cff import *

# 2) Specific cuts - not needed anymore, as a specific KFFittingSmoother deals with that.
# Add a chi**2 cut to retain/reject hits
# KFFittingSmoother.EstimateCut = 15.0
# Request three hits to make a track
# KFFittingSmoother.MinNumberOfHits = 3

# 3) Fast Simulation tracking sequences
# this one is added before 340pre3 to cope with adding SiPixelTemplateDBObjectESProducer and corresponding objects to the ConfDB (MC_3XY_V11, STARTUP3X_V10)
from CalibTracker.SiPixelESProducers.SiPixelTemplateDBObjectESProducer_cfi import *
from FastSimulation.Tracking.GlobalPixelTracking_cff import *
from FastSimulation.Tracking.IterativeTracking_cff import *

# Calo RecHits producer (with no HCAL miscalibration by default)
from FastSimulation.CaloRecHitsProducer.CaloRecHits_cff import *

# ECAL clusters
from RecoEcal.Configuration.RecoEcal_cff import *
reducedRecHitsSequence.remove(seldigis)


# Calo Towers
from RecoJets.Configuration.CaloTowersRec_cff import *

# Particle Flow
from RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff import *
#from RecoParticleFlow.PFTracking.particleFlowTrack_cff import *
from RecoParticleFlow.PFTracking.particleFlowTrackWithDisplacedVertex_cff import *
from RecoParticleFlow.PFProducer.particleFlowSimParticle_cff import *
from RecoParticleFlow.PFProducer.particleFlowBlock_cff import *
from RecoParticleFlow.PFProducer.particleFlow_cff import *
from RecoParticleFlow.PFProducer.pfElectronTranslator_cff import *
from RecoParticleFlow.PFTracking.trackerDrivenElectronSeeds_cff import *

particleFlowSimParticle.sim = 'famosSimHits'

#Deactivate the recovery of dead towers since dead towers are not simulated
particleFlowRecHitHCAL.ECAL_Compensate = cms.bool(False)
#Similarly, deactivate HF cleaning for spikes
particleFlowRecHitHCAL.ShortFibre_Cut = cms.double(1E5)
particleFlowRecHitHCAL.LongFibre_Cut = cms.double(1E5)
particleFlowRecHitHCAL.LongShortFibre_Cut = cms.double(1E5)
particleFlowRecHitHCAL.ApplyLongShortDPG = cms.bool(False)
particleFlowClusterHFEM.thresh_Clean_Barrel = cms.double(1E5)
particleFlowClusterHFEM.thresh_Clean_Endcap = cms.double(1E5)
particleFlowClusterHFHAD.thresh_Clean_Barrel = cms.double(1E5)
particleFlowClusterHFHAD.thresh_Clean_Endcap = cms.double(1E5)

#particleFlowBlock.useNuclear = cms.bool(True)
#particleFlowBlock.useConversions = cms.bool(True)
#particleFlowBlock.useV0 = cms.bool(True)

#particleFlow.rejectTracks_Bad =  cms.bool(False)
#particleFlow.rejectTracks_Step45 = cms.bool(False)

#particleFlow.usePFNuclearInteractions = cms.bool(True)
#particleFlow.usePFConversions = cms.bool(True)
#particleFlow.usePFDecays = cms.bool(True)


famosParticleFlowSequence = cms.Sequence(
    caloTowersRec+
#    pfTrackElec+
    particleFlowTrackWithDisplacedVertex+
    particleFlowBlock+
    particleFlow+
    pfElectronTranslatorSequence    
)

# Reco Jets and MET
from RecoJets.Configuration.RecoJets_cff import *
from RecoJets.Configuration.JetIDProducers_cff import *
from RecoJets.Configuration.RecoPFJets_cff import *
from RecoMET.Configuration.RecoMET_cff import *
from RecoMET.Configuration.RecoPFMET_cff import *


caloJetMet = cms.Sequence(
    recoJets+
    recoJetIds+ 
    metreco
)

PFJetMet = cms.Sequence(
    recoPFJets+
    recoPFMET
)

# Gen Jets
from PhysicsTools.HepMCCandAlgos.genParticles_cfi import *
from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.Configuration.RecoGenJets_cff import *
from RecoMET.Configuration.GenMETParticles_cff import *
from RecoMET.Configuration.RecoGenMET_cff import *
# No longer applicable according to Ronny
#genCandidatesForMET.verbose = False
caloJetMetGen = cms.Sequence(
    genParticles+
    genJetParticles+
    recoGenJets+
    genMETParticles+
    recoGenMET
)

# Muon parametrization
from FastSimulation.ParamL3MuonProducer.ParamL3Muon_cfi import *

# Muon simHit sequence 
from FastSimulation.MuonSimHitProducer.MuonSimHitProducer_cfi import *

# Muon Digi sequence
from SimMuon.Configuration.SimMuon_cff import *
simMuonCSCDigis.strips.doCorrelatedNoise = False ## Saves a little bit of time

simMuonCSCDigis.InputCollection = 'MuonSimHitsMuonCSCHits'
simMuonDTDigis.InputCollection = 'MuonSimHitsMuonDTHits'
simMuonRPCDigis.InputCollection = 'MuonSimHitsMuonRPCHits'

# Muon RecHit sequence
from RecoLocalMuon.Configuration.RecoLocalMuon_cff import *
csc2DRecHits.stripDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigi")
csc2DRecHits.wireDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
rpcRecHits.rpcDigiLabel = 'simMuonRPCDigis'
dt1DRecHits.dtDigiLabel = 'simMuonDTDigis'
dt1DCosmicRecHits.dtDigiLabel = 'simMuonDTDigis'

# Muon reconstruction sequence
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.TrackingTools.MuonTrackLoader_cff import *
KFSmootherForMuonTrackLoader.Propagator = 'SmartPropagatorAny'
from RecoMuon.MuonSeedGenerator.standAloneMuonSeeds_cff import *
from RecoMuon.StandAloneMuonProducer.standAloneMuons_cff import *
from FastSimulation.Configuration.globalMuons_cff import *
globalMuons.GLBTrajBuilderParameters.TrackTransformer.TrackerRecHitBuilder = 'WithoutRefit'
globalMuons.GLBTrajBuilderParameters.TrackerRecHitBuilder = 'WithoutRefit'
globalMuons.GLBTrajBuilderParameters.TransformerOutPropagator = cms.string('SmartPropagatorAny')
globalMuons.GLBTrajBuilderParameters.MatcherOutPropagator = cms.string('SmartPropagator')

from RecoMuon.GlobalMuonProducer.tevMuons_cfi import *
GlobalMuonRefitter.TrackerRecHitBuilder = 'WithoutRefit'
GlobalMuonRefitter.Propagator = 'SmartPropagatorAny'
GlobalTrajectoryBuilderCommon.TrackerRecHitBuilder = 'WithoutRefit'
tevMuons.RefitterParameters.TrackerRecHitBuilder = 'WithoutRefit'
tevMuons.RefitterParameters.Propagator =  'SmartPropagatorAny'
KFSmootherForRefitInsideOut.Propagator = 'SmartPropagatorAny'
KFSmootherForRefitOutsideIn.Propagator = 'SmartPropagator'
KFFitterForRefitInsideOut.Propagator = 'SmartPropagatorAny'
KFFitterForRefitOutsideIn.Propagator = 'SmartPropagatorAny'

famosMuonSequence = cms.Sequence(
    muonDigi+
    muonlocalreco+
    ancientMuonSeed+
    standAloneMuons+
    globalMuons+
    tevMuons
)

#Muon identification sequence
from FastSimulation.Configuration.muonIdentification_cff import *
# Use FastSim tracks and calo hits for muon id
muons.inputCollectionLabels = cms.VInputTag(
    'generalTracks',
    'globalMuons',
    cms.InputTag("standAloneMuons","UpdatedAtVtx")
)
# Use FastSim tracks and calo hits for calo muon id
calomuons.inputTracks = 'generalTracks'

# Muon isolation
from RecoMuon.MuonIsolationProducers.muIsolation_cff import *

famosMuonIdAndIsolationSequence = cms.Sequence(
    ak5CaloJets+
    muonIdProducerSequence+
    muIsolation
)

# Electron reconstruction
from FastSimulation.Tracking.globalCombinedSeeds_cfi import *
from RecoEgamma.EgammaElectronProducers.ecalDrivenElectronSeeds_cfi import *
from FastSimulation.EgammaElectronAlgos.electronGSGsfTrackCandidates_cff import *
from RecoEgamma.EgammaElectronProducers.gsfElectronSequence_cff import *
from TrackingTools.GsfTracking.GsfElectronFit_cff import *
from RecoEgamma.EgammaPhotonProducers.trackerOnlyConversionSequence_cff import *
from TrackingTools.GsfTracking.CkfElectronCandidateMaker_cff import *
from TrackingTools.GsfTracking.FwdElectronPropagator_cfi import *
import TrackingTools.GsfTracking.GsfElectronFit_cfi

electronGsfTracks = TrackingTools.GsfTracking.GsfElectronFit_cfi.GsfGlobalElectronTest.clone()
electronGsfTracks.src = 'electronGSGsfTrackCandidates'
electronGsfTracks.TTRHBuilder = 'WithoutRefit'
electronGsfTracks.TrajectoryInEvent = True

from RecoParticleFlow.PFTracking.mergedElectronSeeds_cfi import *
from RecoEgamma.ElectronIdentification.electronIdSequence_cff import *

famosGsfTrackSequence = cms.Sequence(
    iterativeFirstSeeds+
    newCombinedSeeds+
    particleFlowCluster+ 
    ecalDrivenElectronSeeds+
    trackerDrivenElectronSeeds+
    electronMergedSeeds+
    electronGSGsfTrackCandidates+
    electronGsfTracks
)

# Photon reconstruction
from RecoEgamma.EgammaPhotonProducers.photonSequence_cff import *
photons.hbheInstance = ''
#photons.pixelSeedProducer = 'fastElectronSeeds'
from RecoEgamma.PhotonIdentification.photonId_cff import *

famosPhotonSequence = cms.Sequence(
    photonSequence+
    photonIDSequence
)

# Add pre-calculated isolation sums for electrons (NB for photons they are stored in the Photon. All is done in the
# sequence above
from RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff import *

#Add egamma ecal interesting rec hits
from RecoEgamma.EgammaIsolationAlgos.interestingEleIsoDetIdModule_cff import *
from RecoEgamma.EgammaIsolationAlgos.interestingGamIsoDetIdModule_cff import *

from RecoEgamma.EgammaIsolationAlgos.interestingEgammaIsoDetIdsSequence_cff import *
#import  RecoEgamma.EgammaIsolationAlgos.interestingEgammaIsoDetIdsSequence_cff



# B tagging
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
ak5JetTracksAssociatorAtVertex.tracks = 'generalTracks'
from RecoVertex.Configuration.RecoVertex_cff import *
from RecoVertex.BeamSpotProducer.BeamSpot_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *
offlinePrimaryVerticesWithBS.TrackLabel = 'generalTracks'

famosBTaggingSequence = cms.Sequence(
    btagging
)

#Tau tagging
from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import *
from RecoJets.JetAssociationProducers.ic5PFJetTracksAssociatorAtVertex_cfi import *
ic5JetTracksAssociatorAtVertex.tracks = 'generalTracks'
ic5PFJetTracksAssociatorAtVertex.tracks = 'generalTracks'
from RecoTauTag.Configuration.RecoTauTag_cff import *

famosTauTaggingSequence = cms.Sequence(tautagging)

from RecoTauTag.Configuration.RecoPFTauTag_cff import *

famosPFTauTaggingSequence = cms.Sequence(PFTau)

# The sole simulation sequence
famosSimulationSequence = cms.Sequence(
    offlineBeamSpot+
    famosPileUp+
    famosSimHits+
    MuonSimHits+
    mix
)

# Famos pre-defined sequences (and self-explanatory names)
famosWithTrackerHits = cms.Sequence(
    famosSimulationSequence+
    siTrackerGaussianSmearingRecHits
)

famosWithTrackerAndCaloHits = cms.Sequence(
    famosWithTrackerHits+
    caloRecHits
)

famosWithTracks = cms.Sequence(
    famosWithTrackerHits+
    iterativeTracking
)

famosWithTracksAndMuonHits = cms.Sequence(
    famosSimulationSequence+
    siTrackerGaussianSmearingRecHits+
    iterativeTracking+
    vertexreco+
    famosMuonSequence
)

famosWithTracksAndMuons = cms.Sequence(
    famosSimulationSequence+
    siTrackerGaussianSmearingRecHits+
    iterativeTracking+
    vertexreco+
    famosMuonSequence+
    caloRecHits+
    caloTowersRec+
    famosMuonIdAndIsolationSequence
)

famosWithCaloHits = cms.Sequence(
    famosSimulationSequence+
    caloRecHits
)

famosWithEcalClusters = cms.Sequence(
    famosWithCaloHits+
    ecalClusters+
    particleFlowCluster
)

famosWithTracksAndCaloHits = cms.Sequence(
    famosWithTracks+
    caloRecHits
)

famosWithTracksAndEcalClusters = cms.Sequence(
    famosWithTracksAndCaloHits+
    ecalClusters+
    particleFlowCluster
)

famosWithParticleFlow = cms.Sequence(
    famosWithTracksAndEcalClusters+
    vertexreco+
    trackerOnlyConversionSequence+
    caloTowersRec+ 
    famosParticleFlowSequence+
    PFJetMet
)

famosWithCaloTowers = cms.Sequence(
    famosWithCaloHits+
    caloTowersRec
)

famosElectronSequence = cms.Sequence(
        famosGsfTrackSequence+
        famosWithParticleFlow+
        gsfElectronSequence+
        eIdSequence
        )

famosWithTracksAndCaloTowers = cms.Sequence(
    famosWithTracksAndCaloHits+
    caloTowersRec
)

famosWithTracksAndJets = cms.Sequence(
    famosWithTracksAndCaloTowers+
    vertexreco+
    caloJetMetGen+
    caloJetMet
)

### Standard Jets _cannot_ be done without many other things...
#######################################################################
famosWithJets = cms.Sequence(
    famosWithTracksAndCaloTowers+
    vertexreco+
    ecalClusters+
    particleFlowCluster+
    famosGsfTrackSequence+
    famosMuonSequence+
    famosMuonIdAndIsolationSequence+
    famosParticleFlowSequence+
    gsfElectronSequence+	
    caloJetMetGen+
    caloJetMet
)

##--- simplified IC05 jets only
famosWithSimpleJets = cms.Sequence(
    famosWithTracksAndCaloTowers+
    vertexreco+
    caloJetMetGen+
    iterativeCone5CaloJets+
    ic5JetTracksAssociatorAtVertex
)

famosWithCaloTowersAndParticleFlow = cms.Sequence(
    famosWithParticleFlow+
    caloTowersRec
)

famosWithMuons = cms.Sequence(
    famosWithTracks+
    paramMuons
)

famosWithMuonsAndIsolation = cms.Sequence(
    famosWithTracksAndCaloTowers+
    paramMuons+
    ak5CaloJets+
    muIsolation_ParamGlobalMuons
)

famosWithElectrons = cms.Sequence(
    famosWithTracksAndEcalClusters+
    caloTowersRec+
    famosGsfTrackSequence+
    famosParticleFlowSequence+
    famosElectronSequence+
    interestingEleIsoDetIdEB+
    interestingEleIsoDetIdEE+
    egammaIsolationSequence
)

famosWithPhotons = cms.Sequence(
    famosWithTracks+
    vertexreco+
    caloRecHits+
    ecalClusters+
    famosPhotonSequence+
    interestingGamIsoDetIdEB+
    interestingGamIsoDetIdEE
)

famosWithElectronsAndPhotons = cms.Sequence(
    famosWithTracks+
    vertexreco+
    caloRecHits+
    ecalClusters+
    caloTowersRec+
    famosElectronSequence+
    famosPhotonSequence+
    interestingEgammaIsoDetIds+
    egammaIsolationSequence
)

famosWithBTagging = cms.Sequence(
    famosWithTracksAndCaloTowers+
    vertexreco+
    ak5CaloJets+
    ak5JetTracksAssociatorAtVertex+
    ecalClusters+
    famosMuonSequence+
    reducedRecHitsSequence+ 
    famosBTaggingSequence
    )

famosWithTauTagging = cms.Sequence(
    famosWithTracksAndCaloTowers+
    vertexreco+
    iterativeCone5CaloJets+
    ic5JetTracksAssociatorAtVertex+
    ecalClusters+
    famosTauTaggingSequence
)

famosWithPFTauTagging = cms.Sequence(
    famosWithCaloTowersAndParticleFlow+
    famosPFTauTaggingSequence
)

# The simulation sequence without muon digitization
simulationNoMuonDigiWithFamos = cms.Sequence(
    famosSimulationSequence+
    siTrackerGaussianSmearingRecHits+
    caloRecHits
)

# The simulation and digitization sequence
simulationWithFamos = cms.Sequence(
    famosSimulationSequence+
    muonDigi+
    siTrackerGaussianSmearingRecHits+
    caloRecHits
)


# The reconstruction sequence
reconstructionWithFamos = cms.Sequence(
    iterativeTracking+
    vertexreco+
    trackerOnlyConversionSequence+
    caloTowersRec+
    ecalClusters+
    particleFlowCluster+
    famosGsfTrackSequence+
    famosMuonSequence+
    famosMuonIdAndIsolationSequence+
    famosParticleFlowSequence+
    gsfElectronSequence+
    eIdSequence+
    famosPhotonSequence+
    egammaIsolationSequence+
    interestingEgammaIsoDetIds+
    caloJetMetGen+
    caloJetMet+
    PFJetMet+
#    paramMuons+
#    muIsolation_ParamGlobalMuons+
    ic5JetTracksAssociatorAtVertex+
    ak5JetTracksAssociatorAtVertex+
    famosTauTaggingSequence+
    reducedRecHitsSequence+
    famosBTaggingSequence+
    famosPFTauTaggingSequence
)

# Simulation plus reconstruction
famosWithEverything = cms.Sequence(
    simulationWithFamos+
    reconstructionWithFamos
)

