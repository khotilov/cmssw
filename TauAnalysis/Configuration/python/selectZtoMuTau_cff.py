import FWCore.ParameterSet.Config as cms

from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *
from TauAnalysis.CandidateTools.sysErrDefinitions_cfi import *

#--------------------------------------------------------------------------------
# define event selection criteria for Z --> mu + tau-jet channel
#--------------------------------------------------------------------------------

# generator level phase-space selection
# (NOTE: to be used in case of Monte Carlo samples
#        overlapping in simulated phase-space only !!)
cfgGenPhaseSpaceCut = cms.PSet(
    pluginName = cms.string('genPhaseSpaceCut'),
    pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
    src = cms.InputTag('genPhaseSpaceEventInfo'),
    cut = cms.string('')
)

# trigger selection
cfgTrigger = cms.PSet(
    pluginName = cms.string('Trigger'),
    pluginType = cms.string('PATTriggerEventSelector'),
    src = cms.InputTag('patTriggerEvent'),
    hltAcceptPaths = cms.vstring('HLT_Mu9')
)

# muon candidate selection
cfgGlobalMuonCut = cms.PSet(
    pluginName = cms.string('globalMuonCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatMuonsGlobalCumulative'),
    src_individual = cms.InputTag('selectedPatMuonsGlobalIndividual'),
    systematics = cms.vstring(muonSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgMuonEtaCut = cms.PSet(
    pluginName = cms.string('muonEtaCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatMuonsEta21Cumulative'),
    src_individual = cms.InputTag('selectedPatMuonsEta21Individual'),
    systematics = cms.vstring(muonSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgMuonPtCut = cms.PSet(
    pluginName = cms.string('muonPtCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatMuonsPt15Cumulative'),
    src_individual = cms.InputTag('selectedPatMuonsPt15Individual'),
    systematics = cms.vstring(muonSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgMuonVbTfIdCut = cms.PSet(
    pluginName = cms.string('muonVbTfIdCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatMuonsVbTfIdCumulative'),
    src_individual = cms.InputTag('selectedPatMuonsVbTfIdIndividual'),
    systematics = cms.vstring(muonSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgMuonPFRelIsoCut = cms.PSet(
    pluginName = cms.string('muonPFRelIsoCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatMuonsPFRelIsoCumulative'),
    src_individual = cms.InputTag('selectedPatMuonsPFRelIsoIndividual'),
    systematics = cms.vstring(muonSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgMuonTrkIPcut = cms.PSet(
    pluginName = cms.string('muonTrkIPcut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatMuonsTrkIPcumulative'),
    src_individual = cms.InputTag('selectedPatMuonsTrkIPindividual'),
    systematics = cms.vstring(muonSystematics.keys()),
    minNumber = cms.uint32(1)
)

# tau candidate selection
cfgTauAntiOverlapWithMuonsVeto = cms.PSet(
    pluginName = cms.string('tauAntiOverlapWithMuonsVeto'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForMuTauAntiOverlapWithMuonsVetoCumulative'),
    src_individual = cms.InputTag('selectedPatTausForMuTauAntiOverlapWithMuonsVetoIndividual'),
    systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauEtaCut = cms.PSet(
    pluginName = cms.string('tauEtaCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForMuTauEta23Cumulative'),
    src_individual = cms.InputTag('selectedPatTausForMuTauEta23Individual'),
    systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauPtCut = cms.PSet(
    pluginName = cms.string('tauPtCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForMuTauPt20Cumulative'),
    src_individual = cms.InputTag('selectedPatTausForMuTauPt20Individual'),
    systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauLeadTrkCut = cms.PSet(
    pluginName = cms.string('tauLeadTrkCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForMuTauLeadTrkCumulative'),
    src_individual = cms.InputTag('selectedPatTausForMuTauLeadTrkIndividual'),
    systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauLeadTrkPtCut = cms.PSet(
    pluginName = cms.string('tauLeadTrkPtCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForMuTauLeadTrkPtCumulative'),
    src_individual = cms.InputTag('selectedPatTausForMuTauLeadTrkPtIndividual'),
    systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauTaNCdiscrCut = cms.PSet(
    pluginName = cms.string('tauTaNCdiscrCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForMuTauTaNCdiscrCumulative'),
    src_individual = cms.InputTag('selectedPatTausForMuTauTaNCdiscrIndividual'),
    systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauProngCut = cms.PSet(
    pluginName = cms.string('tauProngCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForMuTauProngCumulative'),
    src_individual = cms.InputTag('selectedPatTausForMuTauProngIndividual'),
    systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauChargeCut = cms.PSet(
    pluginName = cms.string('tauChargeCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForMuTauChargeCumulative'),
    src_individual = cms.InputTag('selectedPatTausForMuTauChargeIndividual'),
    systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauMuonVeto = cms.PSet(
    pluginName = cms.string('tauMuonVeto'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForMuTauMuonVetoCumulative'),
    src_individual = cms.InputTag('selectedPatTausForMuTauMuonVetoIndividual'),
    systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauElectronVeto = cms.PSet(
    pluginName = cms.string('tauElectronVeto'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForMuTauElectronVetoCumulative'),
    src_individual = cms.InputTag('selectedPatTausForMuTauElectronVetoIndividual'),
    systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)

# di-tau candidate selection
cfgDiTauCandidateForMuTauAntiOverlapVeto = cms.PSet(
    pluginName = cms.string('diTauCandidateForMuTauAntiOverlapVeto'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedMuTauPairsAntiOverlapVetoCumulative'),
    src_individual = cms.InputTag('selectedMuTauPairsAntiOverlapVetoIndividual'),
    systematics = cms.vstring(muTauPairSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgDiTauCandidateForMuTauMt1METcut = cms.PSet(
    pluginName = cms.string('diTauCandidateForMuTauMt1METcut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedMuTauPairsMt1METcumulative'),
    src_individual = cms.InputTag('selectedMuTauPairsMt1METindividual'),
    systematics = cms.vstring(muTauPairSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgDiTauCandidateForMuTauPzetaDiffCut = cms.PSet(
    pluginName = cms.string('diTauCandidateForMuTauPzetaDiffCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedMuTauPairsPzetaDiffCumulative'),
    src_individual = cms.InputTag('selectedMuTauPairsPzetaDiffIndividual'),
    systematics = cms.vstring(muTauPairSystematics.keys()),
    minNumber = cms.uint32(1)
)

# selection of event vertex associated to muon + tau-jet pair
cfgPrimaryEventVertexForMuTau = cms.PSet(
    pluginName = cms.string('primaryEventVertexForMuTau'),
    pluginType = cms.string('VertexMinEventSelector'),
    src = cms.InputTag('selectedPrimaryVertexForMuTau'),
    minNumber = cms.uint32(1)
)
cfgPrimaryEventVertexQualityForMuTau = cms.PSet(
    pluginName = cms.string('primaryEventVertexQualityForMuTau'),
    pluginType = cms.string('VertexMinEventSelector'),
    src = cms.InputTag('selectedPrimaryVertexQualityForMuTau'),
    minNumber = cms.uint32(1)
)
cfgPrimaryEventVertexPositionForMuTau = cms.PSet(
    pluginName = cms.string('primaryEventVertexPositionForMuTau'),
    pluginType = cms.string('VertexMinEventSelector'),
    src = cms.InputTag('selectedPrimaryVertexPositionForMuTau'),
    minNumber = cms.uint32(1)
)

# "final" selection of di-tau candidates for "OppositeSign" signal region
cfgDiTauCandidateForMuTauZeroChargeCut = cms.PSet(
    pluginName = cms.string('diTauCandidateForMuTauZeroChargeCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedMuTauPairsZeroChargeCumulative'),
    src_individual = cms.InputTag('selectedMuTauPairsZeroChargeIndividual'),
    systematics = cms.vstring(muTauPairSystematics.keys()),
    minNumber = cms.uint32(1)
)

# "final" selection of di-tau candidates for "SameSign" background dominated control region
cfgDiTauCandidateForMuTauNonZeroChargeCut = cms.PSet(
    pluginName = cms.string('diTauCandidateForMuTauNonZeroChargeCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedMuTauPairsNonZeroChargeCumulative'),
    src_individual = cms.InputTag('selectedMuTauPairsNonZeroChargeIndividual'),
    systematics = cms.vstring(muTauPairSystematics.keys()),
    minNumber = cms.uint32(1)
)

# veto events compatible with Z/gamma* --> mu+ mu- hypothesis
# (two paths:two
#   o two (loosely) isolated global muons of Pt > 10 GeV and opposite charge
#   o two global muons of opposite charge and dR(mu+,mu-) < 1.0
cfgDiMuPairZmumuHypothesisVetoByLooseIsolation = cms.PSet(
    pluginName = cms.string('diMuPairZmumuHypothesisVetoByLooseIsolation'),
    pluginType = cms.string('PATCandViewMaxEventSelector'),
    src = cms.InputTag('selectedDiMuPairZmumuHypothesesByLooseIsolation'),
    #systematics = cms.vstring(muonSystematics.keys()),
    maxNumber = cms.uint32(0)
)

cfgDiMuPairDYmumuHypothesisVeto = cms.PSet(
    pluginName = cms.string('diMuPairDYmumuHypothesisVeto'),
    pluginType = cms.string('PATCandViewMaxEventSelector'),
    src = cms.InputTag('selectedDiMuPairDYmumuHypotheses'),
    #systematics = cms.vstring(muonSystematics.keys()),
    maxNumber = cms.uint32(0)
)

zToMuTauEventSelConfiguratorOS = eventSelFlagProdConfigurator(
    [ cfgGenPhaseSpaceCut,
      cfgTrigger,      
      cfgGlobalMuonCut,
      cfgMuonEtaCut,
      cfgMuonPtCut,
      cfgTauAntiOverlapWithMuonsVeto,
      cfgTauEtaCut,
      cfgTauPtCut,
      cfgMuonVbTfIdCut,
      cfgMuonPFRelIsoCut,
      cfgMuonTrkIPcut,
      cfgTauLeadTrkCut,
      cfgTauLeadTrkPtCut,
      cfgTauTaNCdiscrCut,
      cfgTauProngCut,
      cfgTauChargeCut,
      cfgTauMuonVeto,
      cfgTauElectronVeto,
      cfgDiTauCandidateForMuTauAntiOverlapVeto,
      cfgDiTauCandidateForMuTauMt1METcut,
      cfgDiTauCandidateForMuTauPzetaDiffCut,
      cfgDiTauCandidateForMuTauZeroChargeCut,
      cfgPrimaryEventVertexForMuTau,
      cfgPrimaryEventVertexQualityForMuTau,
      cfgPrimaryEventVertexPositionForMuTau,
      cfgDiMuPairZmumuHypothesisVetoByLooseIsolation,
      cfgDiMuPairDYmumuHypothesisVeto ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

produceEventSelFlagsZtoMuTauOS = zToMuTauEventSelConfiguratorOS.configure()

zToMuTauEventSelConfiguratorSS = eventSelFlagProdConfigurator(
    [ cfgDiTauCandidateForMuTauNonZeroChargeCut ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

produceEventSelFlagsZtoMuTauSS = zToMuTauEventSelConfiguratorSS.configure()

produceEventSelFlagsZtoMuTau = cms.Sequence(produceEventSelFlagsZtoMuTauOS * produceEventSelFlagsZtoMuTauSS)

isRecZtoMuTau = cms.EDProducer("BoolEventSelFlagProducer",
    pluginName = cms.string('isRecZtoMuTau'),
    pluginType = cms.string('MultiBoolEventSelFlagSelector'),
    flags = cms.VInputTag(
        cms.InputTag('Trigger'),
        cms.InputTag('muonTrkIPcut', 'cumulative'),
        cms.InputTag('tauElectronVeto', 'cumulative'),
        cms.InputTag('diTauCandidateForMuTauZeroChargeCut', 'cumulative'),
        cms.InputTag('primaryEventVertexPositionForMuTau'),                           
        cms.InputTag('diMuPairZmumuHypothesisVetoByLooseIsolation'),
        cms.InputTag('diMuPairDYmumuHypothesisVeto')                           
    )
)

selectZtoMuTauEvents = cms.Sequence(
    produceEventSelFlagsZtoMuTau
   * isRecZtoMuTau
)

