import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.pfMET_cfi  import *
from CommonTools.ParticleFlow.pfParticleSelection_cff import *
from CommonTools.ParticleFlow.pfElectrons_cff import *
from CommonTools.ParticleFlow.pfMuons_cff import *
from CommonTools.ParticleFlow.pfJets_cff import *
from CommonTools.ParticleFlow.pfTaus_cff import *

# sequential top projection cleaning
from CommonTools.ParticleFlow.TopProjectors.pfNoMuon_cfi import * 
from CommonTools.ParticleFlow.TopProjectors.pfNoElectron_cfi import * 
from CommonTools.ParticleFlow.TopProjectors.pfNoJet_cfi import *
from CommonTools.ParticleFlow.TopProjectors.pfNoTau_cfi import *

# generator tools
from CommonTools.ParticleFlow.genForPF2PAT_cff import *


PF2PAT = cms.Sequence(
    pfParticleSelectionSequence + 
    pfMuonSequence + 
    pfNoMuon +
    pfElectronSequence +
    pfNoElectron + 
    pfJetSequence +
    pfNoJet + 
    pfTauSequence +
    pfNoTau +
    pfMET
    )
