import FWCore.ParameterSet.Config as cms

# default JetMET calibration on IC, KT and MC Jets
from JetMETCorrections.Configuration.L2L3Corrections_Summer08_cff import *

# define jet realtive correction services
#from JetMETCorrections.Configuration.L2RelativeCorrection152_cff import *

# define jet absolute correction services
#from JetMETCorrections.Configuration.L3AbsoluteCorrection152_cff import *

# define jet flavour correction services
#from JetMETCorrections.Configuration.L5FlavorCorrections_cff import *

# define jet parton correction services
from JetMETCorrections.Configuration.L7PartonCorrections_cff import *

# produce associated jet correction factors in a valuemap
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *

# Pick the L2+L3 corrections
#es_prefer_L2L3JetCorrectorIcone5 = cms.ESPrefer("JetCorrectionServiceChain","L2L3JetCorrectorIcone5")


# MET corrections from JES
from JetMETCorrections.Type1MET.MetType1Corrections_cff import *

# change corrector to L2+L3
corMetType1Icone5.corrector = cms.string('L2L3JetCorrectorIC5Calo')


# MET corrections from muons
import JetMETCorrections.Type1MET.corMetMuons_cfi
# muon MET correction maker 
corMetType1Icone5Muons = JetMETCorrections.Type1MET.corMetMuons_cfi.corMetMuons.clone()
# JetMET corrections for muons: input jet-corrected MET
corMetType1Icone5Muons.inputUncorMetLabel = cms.InputTag('corMetType1Icone5')

# It would be better to get this config to JetMETCorrections/Type1MET/data/ at some point
corMetType1Icone5Muons.TrackAssociatorParameters.useEcal = False ## RecoHits
corMetType1Icone5Muons.TrackAssociatorParameters.useHcal = False ## RecoHits
corMetType1Icone5Muons.TrackAssociatorParameters.useHO = False ## RecoHits
corMetType1Icone5Muons.TrackAssociatorParameters.useCalo = True ## CaloTowers
corMetType1Icone5Muons.TrackAssociatorParameters.useMuon = False ## RecoHits
corMetType1Icone5Muons.TrackAssociatorParameters.truthMatch = False


# re-key jet energy corrections to layer 0 output
layer0JetCorrFactors = cms.EDFilter("JetCorrFactorsValueMapSkimmer",
    collection  = cms.InputTag("allLayer0Jets"),
    backrefs    = cms.InputTag("allLayer0Jets"),
    association = cms.InputTag("jetCorrFactors"),
)


# default PAT sequence for JetMET corrections before cleaners
patAODJetMETCorrections = cms.Sequence(jetCorrFactors * corMetType1Icone5 * corMetType1Icone5Muons)

# default PAT sequence for JetMET corrections after cleaners
patLayer0JetMETCorrections = cms.Sequence(layer0JetCorrFactors)

