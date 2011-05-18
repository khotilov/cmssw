import FWCore.ParameterSet.Config as cms

# This attempts to implement exactly the VBTF selection, which is
# documented at
#
# https://twiki.cern.ch/twiki/bin/view/CMS/VbtfZMuMuBaselineSelection
#
# but this is slightly out of date. The most up-to-date info comes
# from AN-10-264, AN-10-345, and e-mail with Michele de Gruttola.
#
# Currently, their golden dimuons are formed from opposite-charged
# pairs in which both muons must pass this selection:
#
# - muon must be a global muon
# - pT > 20
# - |eta| < 2.1
# - |dxy wrt beamspot| < 0.2 cm (abs(dB) < 0.2)
# - relative combined isolation < 15% ((isolationR03.sumPt + isolationR03.emEt + isolationR03.hadEt) / innerTrack.pt < 0.15)
# - number of tracker hits > 10 (innerTrack.hitPattern.numberOfValidTrackerHits > 10)
# - at least one pixel hit (innerTrack.hitPattern.numberOfValidPixelHits >= 1)
# - at least two muon stations in the fit (globalTrack.hitPattern.muonStationsWithValidHits >= 2)
#
# At least one muon must be matched to a trigger object firing the
# single muon HLT path
# (e.g. !triggerObjectMatchesByPath("HLT_Mu9").empty()) The single
# muon HLT path used will change as the trigger menu evolves with
# luminosity.
#
# So we have a LooseTightCandViewShallowCloneCombiner that requires
# both muons to pass the above cuts ("loose" is then a misnomer), and
# at least one must pass the trigger match requirement (the only
# "tight" cut).

loose_cut = 'isGlobalMuon && ' \
            'isTrackerMuon && ' \
            'innerTrack.pt > 20. && ' \
            'abs(innerTrack.eta) < 2.1 && ' \
            'abs(dB) < 0.2 && ' \
            '(isolationR03.sumPt + isolationR03.emEt + isolationR03.hadEt) / innerTrack.pt < 0.15 && ' \
            'globalTrack.hitPattern.numberOfValidTrackerHits > 10 && ' \
            'globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ' \
            'globalTrack.hitPattern.numberOfValidMuonHits > 0 && ' \
            'numberOfMatches >= 2'

# For the trigger match, currently HLT_Mu24_v1/v2 is the lowest-pT
# unprescaled single muon path. MC does not have such a path; emulate
# with Mu15.
trigger_match = '(' \
                '(!triggerObjectMatchesByPath("HLT_Mu15_v1").empty() && triggerObjectMatchesByPath("HLT_Mu15_v1").at(0).pt() > 24) || ' \
                '!triggerObjectMatchesByPath("HLT_Mu24_v1").empty() || ' \
                '!triggerObjectMatchesByPath("HLT_Mu24_v2").empty()' \
                ')'

tight_cut = trigger_match

allDimuons = cms.EDProducer('LooseTightCandViewShallowCloneCombiner',
                            decay = cms.string('leptons:muons@+ leptons:muons@-'),
                            cut = cms.string(''),
                            loose_cut = cms.string(loose_cut),
                            tight_cut = cms.string(tight_cut)
                            )

dimuons = cms.EDProducer('Zprime2muCompositeCandidatePicker',
                         src = cms.InputTag('allDimuons'),
                         cut = cms.string(''),
                         max_candidates = cms.uint32(100),
                         do_remove_overlap = cms.bool(False),
                         )
