import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
# AlCaReco for muon based alignment using ZMuMu events
ALCARECOMuAlZMuMuHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
import Alignment.CommonAlignmentProducer.AlignmentMuonSelector_cfi
ALCARECOMuAlZMuMu = Alignment.CommonAlignmentProducer.AlignmentMuonSelector_cfi.AlignmentMuonSelector.clone()
seqALCARECOMuAlZMuMu = cms.Sequence(ALCARECOMuAlZMuMuHLT+ALCARECOMuAlZMuMu)
ALCARECOMuAlZMuMuHLT.andOr = True ## choose logical OR between Triggerbits

#ALCARECOMuAlZMuMuHLT.HLTPaths = ['HLT2MuonIso', 'HLT2MuonNonIso', 'HLT2MuonZ']
ALCARECOMuAlZMuMuHLT.HLTPaths = ['HLT_DoubleIsoMu3', 'HLT_DoubleMu3', 'HLT_DoubleMu7_Z']

