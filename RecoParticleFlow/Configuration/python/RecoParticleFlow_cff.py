import FWCore.ParameterSet.Config as cms


from RecoParticleFlow.PFTracking.particleFlowTrack_cff import *
#from RecoParticleFlow.PFTracking.particleFlowTrackWithDisplacedVertex_cff import *

from RecoParticleFlow.PFProducer.particleFlowSimParticle_cff import *
from RecoParticleFlow.PFProducer.particleFlowBlock_cff import *

from RecoParticleFlow.PFProducer.particleFlow_cff import *
from RecoParticleFlow.PFProducer.pfElectronTranslator_cff import *
from RecoParticleFlow.PFProducer.pfPhotonTranslator_cff import *

from RecoParticleFlow.PFProducer.pfLinker_cff import * 

from CommonTools.ParticleFlow.ParticleSelectors.pfCandsForIsolation_cff import *

from CommonTools.ParticleFlow.pfPileUpCandidates_cff import *

particleFlowReco = cms.Sequence( particleFlowTrackWithDisplacedVertex*
                                 particleFlowBlock*
                                 particleFlowTmp*
                                 pfElectronTranslatorSequence*
                                 pfPhotonTranslatorSequence*
                                 pfCandsForIsolationSequence*
                                 pfPileUpCandidatesSequence )

particleFlowLinks = cms.Sequence( particleFlow )
