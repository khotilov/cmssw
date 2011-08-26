import FWCore.ParameterSet.Config as cms

from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *
from TauAnalysis.CandidateTools.sysErrDefinitions_cfi import *

#--------------------------------------------------------------------------------
# define event selection criteria for Z --> e + tau-jet channel
#--------------------------------------------------------------------------------

# 
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
	hltAcceptPaths = cms.vstring(
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2'

    )
)

# electron candidate selection
cfgElectronIdCut = cms.PSet(
    pluginName = cms.string('electronIdCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatElectronsForElecTauIdCumulative'),
    src_individual = cms.InputTag('selectedPatElectronsForElecTauIdIndividual'),
	systematics = cms.vstring(electronSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgElectronAntiCrackCut = cms.PSet(
    pluginName = cms.string('electronAntiCrackCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatElectronsForElecTauAntiCrackCutCumulative'),
    src_individual = cms.InputTag('selectedPatElectronsForElecTauAntiCrackCutIndividual'),
	systematics = cms.vstring(electronSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgElectronEtaCut = cms.PSet(
    pluginName = cms.string('electronEtaCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatElectronsForElecTauEtaCumulative'),
    src_individual = cms.InputTag('selectedPatElectronsForElecTauEtaIndividual'),
	systematics = cms.vstring(electronSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgElectronPtCut = cms.PSet(
    pluginName = cms.string('electronPtCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatElectronsForElecTauPtCumulative'),
    src_individual = cms.InputTag('selectedPatElectronsForElecTauPtIndividual'),
	systematics = cms.vstring(electronSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgElectronIsoCut = cms.PSet(
    pluginName = cms.string('electronIsoCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatElectronsForElecTauIsoCumulative'),
    src_individual = cms.InputTag('selectedPatElectronsForElecTauIsoIndividual'),
    systematics = cms.vstring(electronSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgElectronConversionVeto = cms.PSet(
    pluginName = cms.string('electronConversionVeto'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatElectronsForElecTauConversionVetoCumulative'),
    src_individual = cms.InputTag('selectedPatElectronsForElecTauConversionVetoIndividual'),
	systematics = cms.vstring(electronSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgElectronTrkIPcut = cms.PSet(
    pluginName = cms.string('electronTrkIPcut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatElectronsForElecTauTrkIPcumulative'),
    src_individual = cms.InputTag('selectedPatElectronsForElecTauTrkIPindividual'),
	systematics = cms.vstring(electronSystematics.keys()),
    minNumber = cms.uint32(1)
)

# tau candidate selection
cfgTauAntiOverlapWithElectronsVeto = cms.PSet(
    pluginName = cms.string('tauAntiOverlapWithElectronsVeto'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForElecTauAntiOverlapWithElectronsVetoCumulative'),
    src_individual = cms.InputTag('selectedPatTausForElecTauAntiOverlapWithElectronsVetoIndividual'),
	systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauEtaCut = cms.PSet(
    pluginName = cms.string('tauEtaCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForElecTauEtaCumulative'),
    src_individual = cms.InputTag('selectedPatTausForElecTauEtaIndividual'),
	systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauPtCut = cms.PSet(
    pluginName = cms.string('tauPtCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForElecTauPtCumulative'),
    src_individual = cms.InputTag('selectedPatTausForElecTauPtIndividual'),
	systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauDecayModeFindingCut = cms.PSet(
    pluginName = cms.string('tauDecayModeFindingCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForElecTauDecayModeFindingCumulative'),
    src_individual = cms.InputTag('selectedPatTausForElecTauDecayModeFindingIndividual'),
	systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauLeadTrkPtCut = cms.PSet(
    pluginName = cms.string('tauLeadTrkPtCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForElecTauLeadTrkPtCumulative'),
    src_individual = cms.InputTag('selectedPatTausForElecTauLeadTrkPtIndividual'),
	systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauIsoCut = cms.PSet(
    pluginName = cms.string('tauIsoCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForElecTauIsoCumulative'),
    src_individual = cms.InputTag('selectedPatTausForElecTauIsoIndividual'),
	systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauElectronVeto = cms.PSet(
    pluginName = cms.string('tauElectronVeto'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForElecTauElectronVetoCumulative'),
    src_individual = cms.InputTag('selectedPatTausForElecTauElectronVetoIndividual'),
	systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauEcalCrackVeto = cms.PSet(
    pluginName = cms.string('tauEcalCrackVeto'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForElecTauEcalCrackVetoCumulative'),
    src_individual = cms.InputTag('selectedPatTausForElecTauEcalCrackVetoIndividual'),
	systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgTauMuonVeto = cms.PSet(
    pluginName = cms.string('tauMuonVeto'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedPatTausForElecTauMuonVetoCumulative'),
    src_individual = cms.InputTag('selectedPatTausForElecTauMuonVetoIndividual'),
	systematics = cms.vstring(tauSystematics.keys()),
    minNumber = cms.uint32(1)
)


# di-tau candidate selection
cfgDiTauCandidateForElecTauAntiOverlapVeto = cms.PSet(
    pluginName = cms.string('diTauCandidateForElecTauAntiOverlapVeto'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedElecTauPairsAntiOverlapVetoCumulative'),
    src_individual = cms.InputTag('selectedElecTauPairsAntiOverlapVetoIndividual'),
	systematics = cms.vstring(elecTauPairSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgDiTauCandidateForElecTauMt1METCut = cms.PSet(
    pluginName = cms.string('diTauCandidateForElecTauMt1METCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedElecTauPairsMt1METcumulative'),
    src_individual = cms.InputTag('selectedElecTauPairsMt1METindividual'),
	systematics = cms.vstring(elecTauPairSystematics.keys()),
    minNumber = cms.uint32(1)
)
cfgDiTauCandidateForElecTauPzetaDiffCut = cms.PSet(
    pluginName = cms.string('diTauCandidateForElecTauPzetaDiffCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedElecTauPairsPzetaDiffCumulative'),
    src_individual = cms.InputTag('selectedElecTauPairsPzetaDiffIndividual'),
	systematics = cms.vstring(elecTauPairSystematics.keys()),
    minNumber = cms.uint32(1)
)

## selection of vertex associated with elec + tau
cfgPrimaryEventVertexForElecTau = cms.PSet(
    pluginName = cms.string('primaryEventVertexForElecTau'),
    pluginType = cms.string('VertexMinEventSelector'),
    src = cms.InputTag('selectedPrimaryVertexForElecTau'),
    minNumber = cms.uint32(1)
)
cfgPrimaryEventVertexQualityForElecTau = cms.PSet(
    pluginName = cms.string('primaryEventVertexQualityForElecTau'),
    pluginType = cms.string('VertexMinEventSelector'),
    src = cms.InputTag('selectedPrimaryVertexQualityForElecTau'),
    minNumber = cms.uint32(1)
)
cfgPrimaryEventVertexPositionForElecTau = cms.PSet(
    pluginName = cms.string('primaryEventVertexPositionForElecTau'),
    pluginType = cms.string('VertexMinEventSelector'),
    src = cms.InputTag('selectedPrimaryVertexPositionForElecTau'),
    minNumber = cms.uint32(1)
)
#  defines opposite-sign (OS) workflow
cfgDiTauCandidateForElecTauZeroChargeCut = cms.PSet(
    pluginName = cms.string('diTauCandidateForElecTauZeroChargeCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedElecTauPairsZeroChargeCumulative'),
    src_individual = cms.InputTag('selectedElecTauPairsZeroChargeIndividual'),
	systematics = cms.vstring(elecTauPairSystematics.keys()),
    minNumber = cms.uint32(1)
)
#  defines same-sign (SS) workflow
cfgDiTauCandidateForElecTauNonZeroChargeCut = cms.PSet(
    pluginName = cms.string('diTauCandidateForElecTauNonZeroChargeCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src_cumulative = cms.InputTag('selectedElecTauPairsNonZeroChargeCumulative'),
    src_individual = cms.InputTag('selectedElecTauPairsNonZeroChargeIndividual'),
	systematics = cms.vstring(elecTauPairSystematics.keys()),
    minNumber = cms.uint32(1)
)

# veto events compatible with Z --> e+ e- hypothesis
# (based on reconstructed (visible) invariant mass of e + tau-jet pair)
cfgElecTauPairZeeHypothesisVeto = cms.PSet(
    pluginName = cms.string('elecTauPairZeeHypothesisVeto'),
    pluginType = cms.string('ZllHypothesisElecTauMaxEventSelector'),
    src = cms.InputTag('selectedElecTauPairZeeHypotheses'),
    maxNumber = cms.uint32(0)
)

# veto events compatible with Z -> e+ e- by finding a second loosely-isolated electron
#   with opposite charge w.r.t the primary electron
cfgDiElecPairZeeHypothesisVetoByLooseIsolation = cms.PSet(
    pluginName = cms.string('diElecPairZeeHypothesisVetoByLooseIsolation'),
    pluginType = cms.string('PATCandViewMaxEventSelector'),
    src = cms.InputTag('selectedDiElecPairZeeHypothesesByLooseIsolation'),
    maxNumber = cms.uint32(0)
)

zToElecTauEventSelConfiguratorOS = eventSelFlagProdConfigurator(
    [ cfgGenPhaseSpaceCut,
	  cfgTrigger,
      cfgElectronIdCut,
      cfgElectronAntiCrackCut,
      cfgElectronEtaCut,
      cfgElectronPtCut,
      cfgElectronIsoCut,
      cfgElectronConversionVeto,
      cfgElectronTrkIPcut,
      cfgTauAntiOverlapWithElectronsVeto,
      cfgTauEtaCut,
      cfgTauPtCut,
      cfgTauDecayModeFindingCut,
      cfgTauLeadTrkPtCut,
      cfgTauIsoCut,
      cfgTauElectronVeto,
      cfgTauEcalCrackVeto,
      cfgTauMuonVeto,
      cfgDiTauCandidateForElecTauAntiOverlapVeto,
      cfgDiTauCandidateForElecTauMt1METCut,
      cfgDiTauCandidateForElecTauPzetaDiffCut,
      cfgDiTauCandidateForElecTauZeroChargeCut,
      cfgPrimaryEventVertexForElecTau,
      cfgPrimaryEventVertexQualityForElecTau,
      cfgPrimaryEventVertexPositionForElecTau,
      cfgDiElecPairZeeHypothesisVetoByLooseIsolation ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

produceEventSelFlagsZtoElecTauOS = zToElecTauEventSelConfiguratorOS.configure()

zToElecTauEventSelConfiguratorSS = eventSelFlagProdConfigurator(
      [ cfgDiTauCandidateForElecTauNonZeroChargeCut ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

produceEventSelFlagsZtoElecTauSS = zToElecTauEventSelConfiguratorSS.configure()

produceEventSelFlagsZtoElecTau = cms.Sequence( produceEventSelFlagsZtoElecTauOS * produceEventSelFlagsZtoElecTauSS)

isRecZtoElecTau = cms.EDProducer("BoolEventSelFlagProducer",
    pluginName = cms.string('isRecZtoElecTau'),
    pluginType = cms.string('MultiBoolEventSelFlagSelector'),
    flags = cms.VInputTag(
        cms.InputTag('Trigger'),
        cms.InputTag('genPhaseSpaceCut'),
        cms.InputTag('primaryEventVertexPositionForElecTau'),
        cms.InputTag('electronTrkIPcut', 'cumulative'),
        cms.InputTag('tauMuonVeto', 'cumulative'),
        cms.InputTag('diTauCandidateForElecTauZeroChargeCut', 'cumulative'),
        cms.InputTag('diElecPairZeeHypothesisVetoByLooseIsolation'),
    )
)

selectZtoElecTauEvents = cms.Sequence(
    produceEventSelFlagsZtoElecTau
    + isRecZtoElecTau
)

