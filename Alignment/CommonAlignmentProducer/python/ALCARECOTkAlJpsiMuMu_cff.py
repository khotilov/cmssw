import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
# AlCaReco for track based alignment using J/Psi->MuMu events
ALCARECOTkAlJpsiMuMuHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
import Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi
ALCARECOTkAlJpsiMuMu = Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi.AlignmentTrackSelector.clone()
seqALCARECOTkAlJpsiMuMu = cms.Sequence(ALCARECOTkAlJpsiMuMuHLT+ALCARECOTkAlJpsiMuMu)
ALCARECOTkAlJpsiMuMuHLT.andOr = True ## choose logical OR between Triggerbits

ALCARECOTkAlJpsiMuMuHLT.HLTPaths = ['HLT_DoubleMu3_JPsi', 'HLT_DoubleMu4_BJPsi']
ALCARECOTkAlJpsiMuMu.filter = True ##do not store empty events

ALCARECOTkAlJpsiMuMu.applyBasicCuts = True
ALCARECOTkAlJpsiMuMu.ptMin = 0.8 ##GeV

ALCARECOTkAlJpsiMuMu.etaMin = -3.5
ALCARECOTkAlJpsiMuMu.etaMax = 3.5
ALCARECOTkAlJpsiMuMu.nHitMin = 0
ALCARECOTkAlJpsiMuMu.GlobalSelector.applyIsolationtest = False
ALCARECOTkAlJpsiMuMu.GlobalSelector.applyGlobalMuonFilter = True
ALCARECOTkAlJpsiMuMu.TwoBodyDecaySelector.applyMassrangeFilter = True
ALCARECOTkAlJpsiMuMu.TwoBodyDecaySelector.minXMass = 3.0 ##GeV

ALCARECOTkAlJpsiMuMu.TwoBodyDecaySelector.maxXMass = 3.2 ##GeV

ALCARECOTkAlJpsiMuMu.TwoBodyDecaySelector.daughterMass = 0.105 ##GeV (Muons)

ALCARECOTkAlJpsiMuMu.TwoBodyDecaySelector.applyChargeFilter = False
ALCARECOTkAlJpsiMuMu.TwoBodyDecaySelector.charge = 0
ALCARECOTkAlJpsiMuMu.TwoBodyDecaySelector.applyAcoplanarityFilter = False
ALCARECOTkAlJpsiMuMu.TwoBodyDecaySelector.acoplanarDistance = 1 ##radian


