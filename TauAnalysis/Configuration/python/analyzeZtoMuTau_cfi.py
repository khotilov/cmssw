import FWCore.ParameterSet.Config as cms
import copy

# import config for histogram manager filling information about phase-space simulated in Monte Carlo sample
from TauAnalysis.Core.genPhaseSpaceEventInfoHistManager_cfi import *

# import config for muon histogram manager
from TauAnalysis.Core.muonHistManager_cfi import *

# import config for tau-jet histogram manager
from TauAnalysis.Core.pftauHistManager_cfi import *

# import config for di-tau histogram manager
from TauAnalysis.Core.diTauCandidateHistManager_cfi import *
diTauCandidateHistManagerForMuTau = copy.deepcopy(diTauCandidateHistManager)
diTauCandidateHistManagerForMuTau.name = cms.string('diTauCandidateHistManagerForMuTau')
diTauCandidateHistManagerForMuTau.type = cms.string('PATMuTauPairHistManager')
diTauCandidateHistManagerForMuTau.diTauCandidateSource = cms.InputTag('allMuTauPairs')

# import config for missing-Et histogram manager
from TauAnalysis.Core.metHistManager_cfi import *

# import config for central jet veto histogram manager
from TauAnalysis.Core.jetHistManager_cfi import *

# import config for primary event vertex histogram manager
from TauAnalysis.Core.vertexHistManager_cfi import *

# import config for L1 & HLT histogram manager
from TauAnalysis.Core.triggerHistManager_cfi import *
triggerHistManager.l1Bits = cms.vstring('L1_SingleMu3', 'L1_SingleMu5', 'L1_SingleMu7', 'L1_SingleMu10', 'L1_SingleMu14')
triggerHistManager.hltPaths = cms.vstring('HLT_Mu15', 'HLT_IsoMu11')

muTauHistManagers = cms.vstring( 'genPhaseSpaceEventInfoHistManager',
                                 'muonHistManager',
                                 'tauHistManager',
                                 'diTauCandidateHistManagerForMuTau',
                                 'metHistManager',
                                 'jetHistManager',
                                 'vertexHistManager',
                                 'triggerHistManager' )

#--------------------------------------------------------------------------------
# define event selection criteria
#--------------------------------------------------------------------------------

# generator level phase-space selection
# (NOTE: to be used in case of Monte Carlo samples
#        overlapping in simulated phase-space only !!)
genPhaseSpaceCut = cms.PSet(
  name = cms.string('genPhaseSpaceCut'),
  type = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)

# generator level selection of Z --> mu + tau-jet events
# passing basic acceptance and kinematic cuts
# (NOTE: to be used for efficiency studies only !!)
#genMuonCut = cms.PSet(
#  name = cms.string('genMuonCut'),
#  type = cms.string('TauGenJetMinEventSelector'),
#  src = cms.InputTag('selectedGenTauDecaysToMuonPt15Cumulative'),
#  minNumber = cms.uint32(1)
#)
#genTauCut = cms.PSet(
#  name = cms.string('genTauCut'),
#  type = cms.string('TauGenJetMinEventSelector'),
#  src = cms.InputTag('selectedGenTauDecaysToHadronsPt20Cumulative'),
#  minNumber = cms.uint32(1)
#)

# trigger selection
Trigger = cms.PSet(
  name = cms.string('Trigger'),
  type = cms.string('TriggerResultEventSelector'),
  src = cms.InputTag('TriggerResults::HLT'),
  triggerPaths = cms.vstring('HLT_Mu15', 'HLT_IsoMu11')
)

# primary event vertex selection
primaryEventVertex = cms.PSet(
  name = cms.string('primaryEventVertex'),
  type = cms.string('VertexMinEventSelector'),
  src = cms.InputTag('selectedPrimaryVertexHighestPtTrackSum'),
  minNumber = cms.uint32(1)
)
primaryEventVertexQuality = cms.PSet(
  name = cms.string('primaryEventVertexQuality'),
  type = cms.string('VertexMinEventSelector'),
  src = cms.InputTag('selectedPrimaryVertexQuality'),
  minNumber = cms.uint32(1)
)
primaryEventVertexPosition = cms.PSet(
  name = cms.string('primaryEventVertexPosition'),
  type = cms.string('VertexMinEventSelector'),
  src = cms.InputTag('selectedPrimaryVertexPosition'),
  minNumber = cms.uint32(1)
)

# muon candidate selection
globalMuonCut = cms.PSet(
  name = cms.string('globalMuonCut'),
  type = cms.string('PATMuonMinEventSelector'),
  src = cms.InputTag('selectedLayer1MuonsGlobal'),
  minNumber = cms.uint32(1)
)
muonEtaCut = cms.PSet(
  name = cms.string('muonEtaCut'),
  type = cms.string('PATMuonMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1MuonsEta21Cumulative'),
  src_individual = cms.InputTag('selectedLayer1MuonsEta21Individual'),
  minNumber = cms.uint32(1)
)
muonPtCut = cms.PSet(
  name = cms.string('muonPtCut'),
  type = cms.string('PATMuonMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1MuonsPt15Cumulative'),
  src_individual = cms.InputTag('selectedLayer1MuonsPt15Individual'),
  minNumber = cms.uint32(1)
)
muonHLTmatchCut = cms.PSet(
  name = cms.string('muonHLTmatchCut'),
  type = cms.string('PATMuonMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1MuonsHLTmatchCumulative'),
  src_individual = cms.InputTag('selectedLayer1MuonsHLTmatchIndividual'),
  minNumber = cms.uint32(1)
)
muonTrkIsoCut = cms.PSet(
  name = cms.string('muonTrkIsoCut'),
  type = cms.string('PATMuonMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1MuonsTrkIsoCumulative'),
  src_individual = cms.InputTag('selectedLayer1MuonsTrkIsoIndividual'),
  minNumber = cms.uint32(1)
)
muonEcalIsoCut = cms.PSet(
  name = cms.string('muonEcalIsoCut'),
  type = cms.string('PATMuonMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1MuonsEcalIsoCumulative'),
  src_individual = cms.InputTag('selectedLayer1MuonsEcalIsoIndividual'),
  minNumber = cms.uint32(1)
)
muonAntiPionCut = cms.PSet(
  name = cms.string('muonAntiPionCut'),
  type = cms.string('PATMuonMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1MuonsPionVetoCumulative'),
  src_individual = cms.InputTag('selectedLayer1MuonsPionVetoIndividual'),
  minNumber = cms.uint32(1)
)
muonTrkIPcut = cms.PSet(
  name = cms.string('muonTrkIPcut'),
  type = cms.string('PATMuonMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1MuonsTrkIPcumulative'),
  src_individual = cms.InputTag('selectedLayer1MuonsTrkIPindividual'),
  minNumber = cms.uint32(1)
)

# tau candidate selection
tauAntiOverlapWithMuonsVeto = cms.PSet(
  name = cms.string('tauAntiOverlapWithMuonsVeto'),
  type = cms.string('PATTauMinEventSelector'),
  src = cms.InputTag('selectedLayer1TausForMuTauAntiOverlapWithMuonsVeto'),
  minNumber = cms.uint32(1)
)
tauEtaCut = cms.PSet(
  name = cms.string('tauEtaCut'),
  type = cms.string('PATTauMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1TausForMuTauEta21Cumulative'),
  src_individual = cms.InputTag('selectedLayer1TausEta21Individual'),
  minNumber = cms.uint32(1)
)
tauPtCut = cms.PSet(
  name = cms.string('tauPtCut'),
  type = cms.string('PATTauMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1TausForMuTauPt20Cumulative'),
  src_individual = cms.InputTag('selectedLayer1TausPt20Individual'),
  minNumber = cms.uint32(1)
)
tauLeadTrkCut = cms.PSet(
  name = cms.string('tauLeadTrkCut'),
  type = cms.string('PATTauMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1TausForMuTauLeadTrkCumulative'),
  src_individual = cms.InputTag('selectedLayer1TausLeadTrkIndividual'),
  minNumber = cms.uint32(1)
)
tauLeadTrkPtCut = cms.PSet(
  name = cms.string('tauLeadTrkPtCut'),
  type = cms.string('PATTauMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1TausForMuTauLeadTrkPtCumulative'),
  src_individual = cms.InputTag('selectedLayer1TausLeadTrkPtIndividual'),
  minNumber = cms.uint32(1)
)
tauTrkIsoCut = cms.PSet(
  name = cms.string('tauTrkIsoCut'),
  type = cms.string('PATTauMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1TausForMuTauTrkIsoCumulative'),
  src_individual = cms.InputTag('selectedLayer1TausTrkIsoIndividual'),
  minNumber = cms.uint32(1)
)
tauEcalIsoCut = cms.PSet(
  name = cms.string('tauEcalIsoCut'),
  type = cms.string('PATTauMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1TausForMuTauEcalIsoCumulative'),
  src_individual = cms.InputTag('selectedLayer1TausEcalIsoIndividual'),
  minNumber = cms.uint32(1)
)
tauProngCut = cms.PSet(
  name = cms.string('tauProngCut'),
  type = cms.string('PATTauMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1TausForMuTauProngCumulative'),
  src_individual = cms.InputTag('selectedLayer1TausProngIndividual'),
  minNumber = cms.uint32(1)
)
tauMuonVeto = cms.PSet(
  name = cms.string('tauMuonVeto'),
  type = cms.string('PATTauMinEventSelector'),
  src_cumulative = cms.InputTag('selectedLayer1TausForMuTauMuonVetoCumulative'),
  src_individual = cms.InputTag('selectedLayer1TausMuonVetoIndividual'),
  minNumber = cms.uint32(1)
)

# di-tau candidate selection
diTauCandidateForMuTauAntiOverlapVeto = cms.PSet(
  name = cms.string('diTauCandidateForMuTauAntiOverlapVeto'),
  type = cms.string('PATMuTauPairMinEventSelector'),
  src = cms.InputTag('selectedMuTauPairsAntiOverlapVeto'),
  minNumber = cms.uint32(1)
)
diTauCandidateForMuTauAcoplanarityCut = cms.PSet(
  name = cms.string('diTauCandidateForMuTauAcoplanarityCut'),
  type = cms.string('PATMuTauPairMinEventSelector'),
  src_cumulative = cms.InputTag('selectedMuTauPairsAcoplanarityCumulative'),
  src_individual = cms.InputTag('selectedMuTauPairsAcoplanarityIndividual'),
  minNumber = cms.uint32(1)
)
diTauCandidateForMuTauZeroChargeCut = cms.PSet(
  name = cms.string('diTauCandidateForMuTauZeroChargeCut'),
  type = cms.string('PATMuTauPairMinEventSelector'),
  src_cumulative = cms.InputTag('selectedMuTauPairsZeroChargeCumulative'),
  src_individual = cms.InputTag('selectedMuTauPairsZeroChargeIndividual'),
  minNumber = cms.uint32(1)
)
diTauCandidateForMuTauMt1METCut = cms.PSet(
  name = cms.string('diTauCandidateForMuTauMt1METCut'),
  type = cms.string('PATMuTauPairMinEventSelector'),
  src_cumulative = cms.InputTag('selectedMuTauPairsMt1METCumulative'),
  src_individual = cms.InputTag('selectedMuTauPairsMt1METIndividual'),
  minNumber = cms.uint32(1)
)


# veto events containing additional central jets with Et > 20 GeV
centralJetVeto = cms.PSet(
  name = cms.string('centralJetVeto'),
  type = cms.string('PATJetMaxEventSelector'),
  src = cms.InputTag('selectedLayer1JetsEt20'),
  maxNumber = cms.uint32(0)
)

#--------------------------------------------------------------------------------
# define event print-out
#--------------------------------------------------------------------------------

muTauEventDump = cms.PSet(
  name = cms.string('muTauEventDump'),
  type = cms.string('MuTauEventDump'),

  l1GtReadoutRecordSource = cms.InputTag('hltGtDigis::HLT'),
  l1GtObjectMapRecordSource = cms.InputTag('hltL1GtObjectMap::HLT'),
  l1BitsToPrint = cms.vstring('L1_SingleMu3', 'L1_SingleMu5', 'L1_SingleMu7', 'L1_SingleMu10', 'L1_SingleMu14'),

  hltResultsSource = cms.InputTag('TriggerResults::HLT'),
  hltPathsToPrint = cms.vstring('HLT_Mu15', 'HLT_IsoMu11'),
      
  genParticleSource = cms.InputTag('genParticles'),
  genTauJetSource = cms.InputTag('tauGenJets'),
  electronSource = cms.InputTag('allLayer1ElectronsSelForTauAnalyses'),
  muonSource = cms.InputTag('allLayer1MuonsSelForTauAnalyses'),
  tauSource = cms.InputTag('allLayer1PFTausSelForTauAnalyses'),
  diTauCandidateSource = cms.InputTag('allMuTauPairs'),
  metSource = cms.InputTag('allLayer1METs'),
  genMEtSource = cms.InputTag('genMETWithMu'),
  jetSource = cms.InputTag('selectedLayer1JetsEt20'),

  #output = cms.string("muTauEventDump.txt"),
  output = cms.string("std::cout"),

  triggerConditions = cms.vstring("diTauCandidateForMuTauZeroChargeCut: passed_cumulative")
  #triggerConditions = cms.vstring("Trigger: rejected_cumulative",
  #                                "globalMuonCut: rejected_cumulative",
  #                                "muonEtaCut: rejected_cumulative",
  #                                "muonPtCut: rejected_cumulative",
  #                                "muonHLTmatchCut: rejected_cumulative",
  #                                "muonTrkIsoCut: rejected_cumulative",
  #                                "muonEcalIsoCut: rejected_cumulative",
  #                                "muonAntiPionCut: rejected_cumulative",
  #                                "muonTrkIPcut: rejected_cumulative",
  #                                "tauEtaCut: rejected_cumulative",
  #                                "tauPtCut: rejected_cumulative",
  #                                "tauLeadTrkCut: rejected_cumulative",
  #                                "tauLeadTrkPtCut: rejected_cumulative",
  #                                "tauTrkIsoCut: rejected_cumulative",
  #                                "tauEcalIsoCut: rejected_cumulative",
  #                                "tauProngCut: rejected_cumulative",
  #                                "tauMuonVeto: rejected_cumulative")
)

#--------------------------------------------------------------------------------
# define analysis sequence
# (ordered list of event selection criteria and histogram filling)
#--------------------------------------------------------------------------------

muTauAnalysisSequence = cms.VPSet(
  # fill histograms for full event sample
  cms.PSet(
    histManagers = muTauHistManagers
  ),

  # generator level phase-space selection
  # (NOTE: to be used in case of Monte Carlo samples
  #        overlapping in simulated phase-space only !!)
  cms.PSet(
    filter = cms.string('genPhaseSpaceCut'),
    title = cms.string('gen. Phase-Space'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers
  ),

  # generator level selection of Z --> mu + tau-jet events
  # passing basic acceptance and kinematic cuts
  # (NOTE: to be used for efficiency studies only !!)
  #cms.PSet(
  #  filter = cms.string('genMuonCut'),
  #  title = cms.string('gen. Muon'),
  #  saveRunEventNumbers = cms.vstring('')
  #),
  #cms.PSet(
  #  filter = cms.string('genTauCut'),
  #  title = cms.string('gen. Tau'),
  #  saveRunEventNumbers = cms.vstring('')
  #),
  #cms.PSet(
  #  histManagers = muTauHistManagers
  #),
    
  # trigger selection
  cms.PSet(
    filter = cms.string('Trigger'),
    title = cms.string('mu15 || isoMu11 Trigger'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers
  ),

  # primary event vertex selection
  cms.PSet(
    filter = cms.string('primaryEventVertex'),
    title = cms.string('Vertex'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    filter = cms.string('primaryEventVertexQuality'),
    title = cms.string('p(chi2Vertex) > 0.01'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    filter = cms.string('primaryEventVertexPosition'),
    title = cms.string('-50 < zVertex < +50 cm'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers
  ),
    
  # selection of muon candidate
  # produced in muonic tau decay
  cms.PSet(
    filter = cms.string('globalMuonCut'),
    title = cms.string('global Muon'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('muonHistManager.muonSource = selectedLayer1MuonsGlobal')
  ),
  cms.PSet(
    filter = cms.string('muonEtaCut'),
    title = cms.string('-2.1 < eta(Muon) < +2.1'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('muonHistManager.muonSource = selectedLayer1MuonsEta21Cumulative')
  ),
  cms.PSet(
    filter = cms.string('muonPtCut'),
    title = cms.string('Pt(Muon) > 15 GeV'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('muonHistManager.muonSource = selectedLayer1MuonsPt15Cumulative')
  ),
  cms.PSet(
    filter = cms.string('muonHLTmatchCut'),
    title = cms.string('Muon Trigger match'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('muonHistManager.muonSource = selectedLayer1MuonsHLTmatchCumulative')
  ),
  cms.PSet(
    filter = cms.string('muonTrkIsoCut'),
    title = cms.string('Muon Track iso.'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('muonHistManager.muonSource = selectedLayer1MuonsTrkIsoCumulative')
  ),
  cms.PSet(
    filter = cms.string('muonEcalIsoCut'),
    title = cms.string('Muon ECAL iso.'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('muonHistManager.muonSource = selectedLayer1MuonsEcalIsoCumulative')
  ),
  cms.PSet(
    filter = cms.string('muonAntiPionCut'),
    title = cms.string('Muon pi-Veto'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('muonHistManager.muonSource = selectedLayer1MuonsPionVetoCumulative')
  ),
  cms.PSet(
    filter = cms.string('muonTrkIPcut'),
    title = cms.string('Muon Track IP'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('muonHistManager.muonSource = selectedLayer1MuonsTrkIPcumulative')
  ),
    
  # selection of tau-jet candidate
  # produced in hadronic tau decay
  cms.PSet(
    filter = cms.string('tauAntiOverlapWithMuonsVeto'),
    title = cms.string('Tau not overlapping with Muon'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('tauHistManager.tauSource = selectedLayer1TausForMuTauAntiOverlapWithMuonsVeto')
  ),
  cms.PSet(
    filter = cms.string('tauEtaCut'),
    title = cms.string('-2.1 < eta(Tau) < +2.1'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('tauHistManager.tauSource = selectedLayer1TausForMuTauEta21Cumulative')
  ),
  cms.PSet(
    filter = cms.string('tauPtCut'),
    title = cms.string('Pt(Tau) > 20 GeV'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('tauHistManager.tauSource = selectedLayer1TausForMuTauPt20Cumulative')
  ),
  cms.PSet(
    filter = cms.string('tauLeadTrkCut'),
    title = cms.string('Tau lead. Track find.'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('tauHistManager.tauSource = selectedLayer1TausForMuTauLeadTrkCumulative')
  ),
  cms.PSet(
    filter = cms.string('tauLeadTrkPtCut'),
    title = cms.string('Tau lead. Track Pt'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('tauHistManager.tauSource = selectedLayer1TausForMuTauLeadTrkPtCumulative')
  ),
  cms.PSet(
    filter = cms.string('tauTrkIsoCut'),
    title = cms.string('Tau Track iso.'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('tauHistManager.tauSource = selectedLayer1TausForMuTauTrkIsoCumulative')
  ),
  cms.PSet(
    filter = cms.string('tauEcalIsoCut'),
    title = cms.string('Tau ECAL iso.'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('tauHistManager.tauSource = selectedLayer1TausForMuTauEcalIsoCumulative')
  ),
  cms.PSet(
    filter = cms.string('tauProngCut'),
    title = cms.string('Tau 1||3-Prong'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('tauHistManager.tauSource = selectedLayer1TausForMuTauProngCumulative')
  ),
  cms.PSet(
    filter = cms.string('tauMuonVeto'),
    title = cms.string('Tau mu-Veto'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('tauHistManager.tauSource = selectedLayer1TausForMuTauMuonVetoCumulative')
  ),

  #selection of muon + tau-jet combinations
  cms.PSet(
    filter = cms.string('diTauCandidateForMuTauAntiOverlapVeto'),
    title = cms.string('dR(Muon-Tau) > 0.7'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('diTauCandidateHistManagerForMuTau.diTauCandidateSource = selectedMuTauPairsAntiOverlapVeto')
  ),
  cms.PSet(
    filter = cms.string('diTauCandidateForMuTauAcoplanarityCut'),
    title = cms.string('dPhi(Muon-MET) < 2.4'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('diTauCandidateHistManagerForMuTau.diTauCandidateSource = selectedMuTauPairsAcoplanarityCumulative')
  ),
  cms.PSet(
    filter = cms.string('diTauCandidateForMuTauZeroChargeCut'),
    title = cms.string('Charge(Muon+Tau) = 0'),
    saveRunEventNumbers = cms.vstring('')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('diTauCandidateHistManagerForMuTau.diTauCandidateSource = selectedMuTauPairsZeroChargeCumulative'),
  ),
  cms.PSet(
    filter = cms.string('diTauCandidateForMuTauMt1METCut'),
    title = cms.string('M_{T}(Muon-MET) < 60 GeV'),
    saveRunEventNumbers = cms.vstring('passed_cumulative')
  ),
  cms.PSet(
    histManagers = muTauHistManagers,
    replace = cms.vstring('diTauCandidateHistManagerForMuTau.diTauCandidateSource = selectedMuTauPairsMt1METCumulative')
  #),

  # veto events containing additional central jets with Et > 20 GeV
  #cms.PSet(
  #  filter = cms.string('centralJetVeto'),
  #  title = cms.string('central Jet Veto'),
  #  saveRunEventNumbers = cms.vstring('passed_cumulative')
  #),
  #cms.PSet(
  #  histManagers = muTauHistManagers,
  )
)
