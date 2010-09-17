import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# import the two configs for event selection, event print-out and analysis sequence
# of Z --> tau-jet + tau-jet events with loose and with tight tau id. criteria applied;
# import config of "regular" Z --> tau-jet + tau-jet analysis module
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoDiTau_factorized_cfi import *
from TauAnalysis.Configuration.analyzeZtoDiTau_cff import *
from TauAnalysis.Configuration.tools.factorizationTools import replaceEventSelections

#--------------------------------------------------------------------------------
# define Z --> tau-jet + tau-jet analysis module
# for the path with "regular" (tight) tau id. criteria applied
#--------------------------------------------------------------------------------

analyzeZtoDiTauEvents_factorizedTight2ndTau = analyzeZtoDiTauEvents.clone(
    name = cms.string('zDiTauAnalyzer_factorizedTight2ndTau'),
    analysisSequence = diTauAnalysisSequence_factorizedTight2ndTau
)    
analyzeZtoDiTauEvents_factorizedTight2ndTau.eventDumps[0] = diTauEventDump_factorizedTight2ndTau

#--------------------------------------------------------------------------------
# define Z --> tau-jet + tau-jet analysis module
# for the path with "loose" tau id. criteria applied
#
# NOTE: modifications to analyzeZtoDiTauEvents_factorizedLoose2ndTau
#       modify the original analyzeZtoDiTauEvents sequence
#
#      --> analyzeZtoDiTauEvents_factorizedTight2ndTau needs to be defined
#          before analyzeZtoDiTauEvents_factorizedLoose2ndTau !!
#
#--------------------------------------------------------------------------------

analyzeZtoDiTauEvents_factorizedLoose2ndTau = analyzeZtoDiTauEvents.clone(
    name = cms.string('zDiTauAnalyzer_factorizedLoose2ndTau')
)
analyzeZtoDiTauEvents_factorizedLoose2ndTau.eventDumps[0] = diTauEventDump_factorizedLoose2ndTau
replaceEventSelections(analyzeZtoDiTauEvents_factorizedLoose2ndTau, 
    [ [ evtSelSecondTauLeadTrkPt, evtSelSecondTauLeadTrkPtLoose ],
      [ evtSelSecondTauTaNCdiscr, evtSelSecondTauTaNCdiscrLoose ],
      [ evtSelSecondTauTrkIso, evtSelSecondTauTrkIsoLoose ],
      [ evtSelSecondTauEcalIso, evtSelSecondTauEcalIsoLoose ],
      [ evtSelSecondTauProng, evtSelSecondTauProngLoose ],
      [ evtSelSecondTauCharge, evtSelSecondTauChargeLoose ],
      [ evtSelSecondTauMuonVeto, evtSelSecondTauMuonVetoLoose ],
      [ evtSelSecondTauElectronVeto, evtSelSecondTauElectronVetoLoose ],
      [ evtSelDiTauCandidateForDiTauAntiOverlapVeto, evtSelDiTauCandidateForDiTauAntiOverlapVetoLoose2ndTau ],
      [ evtSelDiTauCandidateForDiTauZeroCharge, evtSelDiTauCandidateForDiTauZeroChargeLoose2ndTau ],
      [ evtSelDiTauCandidateForDiTauAcoplanarity, evtSelDiTauCandidateForDiTauAcoplanarityLoose2ndTau] ]
      #[ evtSelCentralJetVeto, evtSelCentralJetVetoLoose2ndTau] ]
)                       
analyzeZtoDiTauEvents_factorizedLoose2ndTau.analysisSequence = diTauAnalysisSequence_factorizedLoose2ndTau

