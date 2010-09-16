import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.RecoTools.patPFTauSelection_cfi import *

# require tau candidate to be within geometric acceptance of Pixel + SiTracker detectors
selectedPatTausForDiTauEta21 = copy.deepcopy(selectedPatTausEta21)

# require tau candidate to have transverse energy above threshold
selectedPatTausForDiTauPt20 = copy.deepcopy(selectedPatTausPt20)

# require tau candidate to have a leading track
# (track of Pt > 1. GeV within matching cone of size dR = 0.2 around jet-axis)
selectedPatTausForDiTauLeadTrk = copy.deepcopy(selectedPatTausLeadTrk)

# require leading track of first (second) tau candidate to have Pt > 12. (8.) GeV
selectedPatTausForDiTau1stLeadTrkPt = selectedPatTausLeadTrkPt.clone(
    cut = cms.string('leadPFChargedHadrCand().isNonnull() & leadPFChargedHadrCand().pt() > 12.')
)    

selectedPatTausForDiTau2ndLeadTrkPt = selectedPatTausLeadTrkPt.clone(
    cut = cms.string('leadPFChargedHadrCand().isNonnull() & leadPFChargedHadrCand().pt() > 8.')
) 

# require tau candidate to pass TaNC discriminator
selectedPatTausForDiTau1stTaNCdiscr = copy.deepcopy(selectedPatTausTaNCdiscr)
selectedPatTausForDiTau2ndTaNCdiscr = copy.deepcopy(selectedPatTausTaNCdiscr)

# require tau candidate to have no tracks of Pt > 1. GeV
# in isolation cone of size dR = 0.8, surrounding signal cone of size dR = 5./Et
selectedPatTausForDiTau1stTrkIso = copy.deepcopy(selectedPatTausTrkIso)
selectedPatTausForDiTau2ndTrkIso = copy.deepcopy(selectedPatTausTrkIso)

# require tau candidate to be isolated
# with respect to energy deposits in ECAL
selectedPatTausForDiTau1stEcalIso = copy.deepcopy(selectedPatTausEcalIso)
selectedPatTausForDiTau2ndEcalIso = copy.deepcopy(selectedPatTausEcalIso)

# require tau candidate to have either one or three tracks within signal cone
selectedPatTausForDiTau1stProng = copy.deepcopy(selectedPatTausProng)
selectedPatTausForDiTau2ndProng = copy.deepcopy(selectedPatTausProng)

# require tau candidate to have charge either +1 or -1
# (computed as sum of charges of tracks within signal cone)
selectedPatTausForDiTau1stCharge = copy.deepcopy(selectedPatTausCharge)
selectedPatTausForDiTau2ndCharge = copy.deepcopy(selectedPatTausCharge)

# require tau candidate to pass muon veto
selectedPatTausForDiTau1stMuonVeto = copy.deepcopy(selectedPatTausMuonVeto)
selectedPatTausForDiTau2ndMuonVeto = copy.deepcopy(selectedPatTausMuonVeto)

# require tau candidate to pass electron veto
selectedPatTausForDiTau1stElectronVeto = copy.deepcopy(selectedPatTausElectronVeto)
selectedPatTausForDiTau2ndElectronVeto = copy.deepcopy(selectedPatTausElectronVeto)
