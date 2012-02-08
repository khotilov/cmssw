import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.RecoTools.patPFTauSelection_cfi import *

# require tau candidate not to overlap with selected electrons
# (in order to avoid double-counting one and the same physical particle
#  as electron and as tau candidate)
selectedPatTausForElecTauAntiOverlapWithElectronsVeto = cms.EDFilter("PATTauAntiOverlapSelector",
    srcNotToBeFiltered = cms.VInputTag("selectedPatElectronsForElecTauPtCumulative"),
    dRmin = cms.double(0.3),
    filter = cms.bool(False)                                           
)

# require tau candidate to be within geometric acceptance of Pixel + SiTracker detectors
selectedPatTausForElecTauEta = copy.deepcopy(selectedPatTausEta23)

# require tau candidate to have transverse energy above threshold
selectedPatTausForElecTauPt = copy.deepcopy(selectedPatTausPt20)

# require tau candidate to have a leading track
# (track of Pt > 1. GeV within matching cone of size dR = 0.2 around jet-axis)
selectedPatTausForElecTauLeadTrk = copy.deepcopy(selectedPatTausLeadTrk)

# preform decay mode reconstruction
selectedPatTausForElecTauDecayModeFinding = cms.EDFilter("PATTauSelector",
    cut = cms.string('tauID("decayModeFinding") > 0.5'),
    filter = cms.bool(False)
)

# require leading track of tau candidate to have Pt > 5. GeV
# BUT 'leadingTrackPtCut' discriminant is not defined for HPS taus
# so make the cut manually
selectedPatTausForElecTauLeadTrkPt = copy.deepcopy(selectedPatTausLeadTrkPt)
selectedPatTausForElecTauLeadTrkPt.cut = cms.string('leadPFChargedHadrCand().isNonnull() & leadPFChargedHadrCand().pt() > 5.')

# require tau candidate to pass TaNC discriminator
selectedPatTausForElecTauTaNCdiscr = copy.deepcopy(selectedPatTausTaNCdiscr)

# require tau candidate to pass HPS/TaNC isolation discriminator
selectedPatTausForElecTauIso = cms.EDFilter("PATTauSelector",
    cut = cms.string('tauID("byHPSloose") > 0.5'),
    filter = cms.bool(False)
)

# require tau candidate to have either one or three tracks within signal cone
selectedPatTausForElecTauProng = copy.deepcopy(selectedPatTausProng)

# require tau candidate to have charge either +1 or -1
# (computed as sum of charges of tracks within signal cone)
selectedPatTausForElecTauCharge = copy.deepcopy(selectedPatTausCharge)

# require tau candidate to pass electron veto
selectedPatTausForElecTauElectronVeto = copy.deepcopy(selectedPatTausElectronVeto)
# update cut to MVA algorithm
selectedPatTausForElecTauElectronVeto.cut = cms.string('againstElectronMVA')

# require tau candidate not to be in ECAL barrel/endcap crack
selectedPatTausForElecTauEcalCrackVeto = copy.deepcopy(selectedPatTausEcalCrackVeto)

# require tau candidate to pass muon veto
selectedPatTausForElecTauMuonVeto = copy.deepcopy(selectedPatTausMuonVeto)
