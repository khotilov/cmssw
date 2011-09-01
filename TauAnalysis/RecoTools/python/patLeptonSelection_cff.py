import FWCore.ParameterSet.Config as cms

from TauAnalysis.RecoTools.patElectronSelection_cfi import *
from TauAnalysis.RecoTools.patElectronSelectionForElecTau_cfi import *
from TauAnalysis.RecoTools.patElectronSelectionForElecMu_cfi import *
from TauAnalysis.RecoTools.patMuonMomentumCorrection_cfi import *
from TauAnalysis.RecoTools.patMuonSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForElecTau_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForMuTau_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForDiTau_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForWTauNu_cfi import *
from TauAnalysis.RecoTools.patLeptonPFIsolationSelector_cfi import *

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *

#--------------------------------------------------------------------------------
# define selection criteria for pat::Electrons
# (settings made here overwrite values defined in electronPatSelector_cfi)
#--------------------------------------------------------------------------------

selectedPatElectronsTightId.cut = cms.string('(abs(superCluster.eta) < 1.479 & electronID("eidRobustTight") > 0 & eSuperClusterOverP < 1.05 & eSuperClusterOverP > 0.95) | (abs(superCluster.eta) > 1.479 & electronID("eidRobustTight") > 0 & eSuperClusterOverP < 1.12 & eSuperClusterOverP > 0.95)')
selectedPatElectronsLooseId.cut = cms.string('(abs(superCluster.eta) < 1.479 & electronID("eidRobustLoose") > 0 & eSuperClusterOverP < 1.05 & eSuperClusterOverP > 0.95) | (abs(superCluster.eta) > 1.479 & electronID("eidRobustLoose") > 0 & eSuperClusterOverP < 1.12 & eSuperClusterOverP > 0.95)')
selectedPatElectronsAntiCrackCut.cut = cms.string('abs(superCluster.eta) < 1.442 | abs(superCluster.eta) > 1.560')
selectedPatElectronsEta21.cut = cms.string('abs(eta) < 2.1')
selectedPatElectronsPt15.cut = cms.string('pt > 15.')
selectedPatElectronsTrkIso.cut = cms.string('userIsolation("pat::TrackIso") < 1.')
selectedPatElectronsEcalIso.cut = cms.string('(abs(superCluster.eta) < 1.479 & userIsolation("pat::EcalIso") < 2.5) | (abs(superCluster.eta) > 1.479 & userIsolation("pat::EcalIso") < 3.5)')
selectedPatElectronsTrk.cut = cms.string('gsfTrack.isNonnull')
selectedPatElectronsTrkIP.vertexSource = cms.InputTag("selectedPrimaryVertexHighestPtTrackSum")
selectedPatElectronsTrkIP.IpMax = cms.double(0.05)

patElectronSelConfigurator = objSelConfigurator(
    [ selectedPatElectronsTightId,
      selectedPatElectronsAntiCrackCut,
      selectedPatElectronsEta21,
      selectedPatElectronsPt15,
      selectedPatElectronsTrkIso,
      selectedPatElectronsEcalIso,
      selectedPatElectronsTrk,
      selectedPatElectronsTrkIP ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatElectrons = patElectronSelConfigurator.configure(pyNameSpace = locals())

selectedPatElectronsTrkIsoLooseIsolation.cut = cms.string('userIsolation("pat::TrackIso") < 8.')
selectedPatElectronsEcalIsoLooseIsolation.cut = cms.string('userIsolation("pat::EcalIso") < 8.')
selectedPatElectronsTrkLooseIsolation.cut = selectedPatElectronsTrk.cut
selectedPatElectronsTrkIPlooseIsolation.vertexSource = selectedPatElectronsTrkIP.vertexSource
selectedPatElectronsTrkIPlooseIsolation.IpMax = selectedPatElectronsTrkIP.IpMax

patElectronSelConfiguratorLooseIsolation = objSelConfigurator(
    [ selectedPatElectronsTightId,
      selectedPatElectronsAntiCrackCut,
      selectedPatElectronsEta21,
      selectedPatElectronsPt15,
      selectedPatElectronsTrkIsoLooseIsolation,
      selectedPatElectronsEcalIsoLooseIsolation,
      selectedPatElectronsTrkLooseIsolation,
      selectedPatElectronsTrkIPlooseIsolation ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatElectronsLooseIsolation = patElectronSelConfiguratorLooseIsolation.configure(pyNameSpace = locals())

#
# select electrons for Z->electron + tau-jet analysis
#

# VBTF WP80 electron ID for pt > 20; WP70 for pt < 20
selectedPatElectronsForElecTauId.cut = cms.string('(abs(superCluster.eta) < 1.479 & pt > 20 & abs(deltaEtaSuperClusterTrackAtVtx) < 0.004 & abs(deltaPhiSuperClusterTrackAtVtx) < 0.06 & hcalOverEcal < 0.04 & sigmaIetaIeta < 0.01) | (abs(superCluster.eta) > 1.479 & pt > 20 & abs(deltaEtaSuperClusterTrackAtVtx) < 0.007 & abs(deltaPhiSuperClusterTrackAtVtx) <0.03 & hcalOverEcal < 0.025 & sigmaIetaIeta < 0.03) | (abs(superCluster.eta) < 1.479 & pt < 20 & abs(deltaEtaSuperClusterTrackAtVtx) < 0.004 & abs(deltaPhiSuperClusterTrackAtVtx) < 0.03 & hcalOverEcal < 0.025 & sigmaIetaIeta < 0.01 & (fbrem > 0.15 | (abs(superCluster.eta) < 1 & eSuperClusterOverP > 0.95) )) | (abs(superCluster.eta) > 1.479 & pt < 20 & abs(deltaEtaSuperClusterTrackAtVtx) < 0.005 & abs(deltaPhiSuperClusterTrackAtVtx) <0.02 & hcalOverEcal < 0.025 & sigmaIetaIeta < 0.03 & (fbrem > 0.15 | (abs(superCluster.eta) < 1 & eSuperClusterOverP > 0.95) ))')
selectedPatElectronsForElecTauAntiCrackCut.cut = cms.string('abs(superCluster.eta) < 1.442 | abs(superCluster.eta) > 1.566')
selectedPatElectronsForElecTauEta.cut = cms.string('abs(eta) < 2.1')
selectedPatElectronsForElecTauPt.cut = cms.string('pt > 15.')
selectedPatElectronsForElecTauIso.sumPtMaxEB = cms.double(0.01)
selectedPatElectronsForElecTauIso.sumPtMaxEE = cms.double(0.01)
selectedPatElectronsForElecTauIso.sumPtMethod = cms.string("relative")
selectedPatElectronsForElecTauConversionVeto.maxMissingInnerHits = cms.int32(0)
selectedPatElectronsForElecTauConversionVeto.invertConversionVeto = cms.bool(False)
selectedPatElectronsForElecTauTrkIP.vertexSource = selectedPatElectronsTrkIP.vertexSource
selectedPatElectronsForElecTauTrkIP.IpMax = cms.double(0.02)
#selectedPatElectronsForElecTauTrkIP.IpZMax = cms.double(0.2)

patElectronSelConfiguratorForElecTau = objSelConfigurator(
    [ selectedPatElectronsForElecTauId,
      selectedPatElectronsForElecTauAntiCrackCut,
      selectedPatElectronsForElecTauEta,
      selectedPatElectronsForElecTauPt,
      selectedPatElectronsForElecTauIso,
      selectedPatElectronsForElecTauConversionVeto,
      selectedPatElectronsForElecTauTrkIP ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatElectronsForElecTau = patElectronSelConfiguratorForElecTau.configure(pyNameSpace = locals())

#  loose isolation
selectedPatElectronsForElecTauIsoLooseIsolation.sumPtMax = cms.double(0.25)
selectedPatElectronsForElecTauTrkIPlooseIsolation.IpMax = selectedPatElectronsForElecTauTrkIP.IpMax

selectedPatElectronsForElecTauConversionVetoLooseIsolation.maxMissingInnerHits = cms.int32(2)
selectedPatElectronsForElecTauConversionVetoLooseIsolation.invertConversionVeto = cms.bool(True)

patElectronSelConfiguratorForElecTauLooseIsolation = objSelConfigurator(
    [ selectedPatElectronsForElecTauId,
      selectedPatElectronsForElecTauAntiCrackCut,
      selectedPatElectronsForElecTauEta,
      selectedPatElectronsForElecTauPt,
      selectedPatElectronsForElecTauIsoLooseIsolation,
      selectedPatElectronsForElecTauConversionVetoLooseIsolation,
      selectedPatElectronsForElecTauTrkIPlooseIsolation ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatElectronsForElecTauLooseIsolation = patElectronSelConfiguratorForElecTauLooseIsolation.configure(pyNameSpace = locals())

#
# select electrons for Z->electron + muon analysis
#

selectedPatElectronsForElecMuAntiOverlapWithMuonsVeto.dRmin = cms.double(0.3)
selectedPatElectronsForElecMuTightId.cut = cms.string('electronID("eidRobustTight") > 0')
selectedPatElectronsForElecMuAntiCrackCut.cut = selectedPatElectronsAntiCrackCut.cut
selectedPatElectronsForElecMuEta21.cut = selectedPatElectronsEta21.cut
selectedPatElectronsForElecMuPt15.cut = selectedPatElectronsPt15.cut
selectedPatElectronsForElecMuTrkIso.cut = cms.string('userIsolation("pat::TrackIso") < 2.')
selectedPatElectronsForElecMuEcalIso.cut = selectedPatElectronsEcalIso.cut
selectedPatElectronsForElecMuTrk.cut = selectedPatElectronsTrk.cut
selectedPatElectronsForElecMuTrkIP.vertexSource = selectedPatElectronsTrkIP.vertexSource
selectedPatElectronsForElecMuTrkIP.IpMax = selectedPatElectronsTrkIP.IpMax

patElectronSelConfiguratorForElecMu = objSelConfigurator(
    [ selectedPatElectronsForElecMuAntiOverlapWithMuonsVeto,
      selectedPatElectronsForElecMuTightId,
      selectedPatElectronsForElecMuAntiCrackCut,
      selectedPatElectronsForElecMuEta21,
      selectedPatElectronsForElecMuPt15,
      selectedPatElectronsForElecMuTrkIso,
      selectedPatElectronsForElecMuEcalIso,
      selectedPatElectronsForElecMuTrk,
      selectedPatElectronsForElecMuTrkIP ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatElectronsForElecMu = patElectronSelConfiguratorForElecMu.configure(pyNameSpace = locals())

selectedPatElectronsForElecMuTrkIsoLooseIsolation.cut = cms.string('userIsolation("pat::TrackIso") < 8.')
selectedPatElectronsForElecMuEcalIsoLooseIsolation.cut = cms.string('userIsolation("pat::EcalIso") < 8.')
selectedPatElectronsForElecMuTrkLooseIsolation.cut = selectedPatElectronsForElecMuTrk.cut
selectedPatElectronsForElecMuTrkIPlooseIsolation.vertexSource = selectedPatElectronsForElecMuTrkIP.vertexSource
selectedPatElectronsForElecMuTrkIPlooseIsolation.IpMax = selectedPatElectronsForElecMuTrkIP.IpMax

patElectronSelConfiguratorForElecMuLooseIsolation = objSelConfigurator(
    [ selectedPatElectronsForElecMuAntiOverlapWithMuonsVeto,
      selectedPatElectronsForElecMuTightId,
      selectedPatElectronsForElecMuAntiCrackCut,
      selectedPatElectronsForElecMuEta21,
      selectedPatElectronsForElecMuPt15,
      selectedPatElectronsForElecMuTrkIsoLooseIsolation,
      selectedPatElectronsForElecMuEcalIsoLooseIsolation,
      selectedPatElectronsForElecMuTrkLooseIsolation,
      selectedPatElectronsForElecMuTrkIPlooseIsolation ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatElectronsForElecMuLooseIsolation = patElectronSelConfiguratorForElecMuLooseIsolation.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------
# define selection criteria for pat::Muons
# (settings made here overwrite values defined in muonPatSelector_cfi)
#--------------------------------------------------------------------------------

selectedPatMuonsGlobal.cut = cms.string('isGlobalMuon()')
selectedPatMuonsEta21.cut = cms.string('abs(eta) < 2.1')
selectedPatMuonsPt15.cut = cms.string('pt > 15.')
selectedPatMuonsVbTfId.beamSpotSource = cms.InputTag("offlineBeamSpot")
selectedPatMuonsVbTfId.vertexSource = cms.InputTag("selectedPrimaryVertexHighestPtTrackSum")
selectedPatMuonsPFRelIso.chargedHadronIso.ptMin = cms.double(0.5)
selectedPatMuonsPFRelIso.chargedHadronIso.dRvetoCone = cms.double(-1.)
selectedPatMuonsPFRelIso.chargedHadronIso.dRisoCone = cms.double(0.4)
selectedPatMuonsPFRelIso.neutralHadronIso.ptMin = cms.double(1.0)
selectedPatMuonsPFRelIso.neutralHadronIso.dRvetoCone = cms.double(0.08)
selectedPatMuonsPFRelIso.neutralHadronIso.dRisoCone = cms.double(0.4)
selectedPatMuonsPFRelIso.photonIso.ptMin = cms.double(1.0)
selectedPatMuonsPFRelIso.photonIso.dRvetoCone = cms.double(0.05)
selectedPatMuonsPFRelIso.photonIso.dRisoCone = cms.double(0.4)
selectedPatMuonsPFRelIso.sumPtMax = cms.double(0.10)
selectedPatMuonsPFRelIso.sumPtMethod = cms.string("relative")
selectedPatMuonsTrk.cut = cms.string('innerTrack.isNonnull')
selectedPatMuonsTrkIP.vertexSource = cms.InputTag("selectedPrimaryVertexHighestPtTrackSum")
selectedPatMuonsTrkIP.IpMax = cms.double(0.05)

patMuonSelConfigurator = objSelConfigurator(
    [ selectedPatMuonsGlobal,
      selectedPatMuonsEta21,
      selectedPatMuonsPt15,
      selectedPatMuonsVbTfId,
      selectedPatMuonsPFRelIso,
      selectedPatMuonsTrk,
      selectedPatMuonsTrkIP ],
    src = "patMuonsMuScleFitCorrectedMomentum",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatMuons = patMuonSelConfigurator.configure(pyNameSpace = locals())

selectedPatMuonsPFRelIsoLooseIsolation.sumPtMax = cms.double(0.25)
selectedPatMuonsPFRelIsoLooseIsolation.sumPtMethod = cms.string("relative")
selectedPatMuonsTrkLooseIsolation.cut = selectedPatMuonsTrk.cut
selectedPatMuonsTrkIPlooseIsolation.vertexSource = selectedPatMuonsTrkIP.vertexSource
selectedPatMuonsTrkIPlooseIsolation.IpMax = selectedPatMuonsTrkIP.IpMax

patMuonSelConfiguratorLooseIsolation = objSelConfigurator(
    [ selectedPatMuonsGlobal,
      selectedPatMuonsEta21,
      selectedPatMuonsPt15,
      selectedPatMuonsVbTfId,
      selectedPatMuonsPFRelIsoLooseIsolation,
      selectedPatMuonsTrkLooseIsolation,
      selectedPatMuonsTrkIPlooseIsolation ],
    src = "patMuonsMuScleFitCorrectedMomentum",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatMuonsLooseIsolation = patMuonSelConfiguratorLooseIsolation.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------
# define selection criteria for pat::(PF)Taus
# (settings made here overwrite values defined in pftauPatSelector_cfi)
#--------------------------------------------------------------------------------

selectedPatTausEta23.cut = cms.string('abs(eta) < 2.3')
selectedPatTausPt20.cut = cut = cms.string('pt > 20.')
selectedPatTausLeadTrk.cut = cms.string('tauID("leadingTrackFinding") > 0.5')
selectedPatTausLeadTrkPt.cut = cms.string('tauID("leadingTrackPtCut") > 0.5')
selectedPatTausTaNCdiscr.cut = cms.string('tauID("byTaNCloose") > 0.5')
selectedPatTausProng.cut = cms.string('signalPFChargedHadrCands.size() = 1 | signalPFChargedHadrCands.size() = 3')
selectedPatTausCharge.cut = cms.string('abs(charge) > 0.5 & abs(charge) < 1.5')
selectedPatTausMuonVeto.cut = cms.string('tauID("againstMuonTight") > 0.5')
selectedPatTausElectronVeto.cut = cms.string('tauID("againstElectronLoose") > 0.5')
selectedPatTausEcalCrackVeto.cut = cms.string('abs(eta) < 1.460 | abs(eta) > 1.558')

patTauSelConfigurator = objSelConfigurator(
    [ selectedPatTausEta23,
      selectedPatTausPt20,
      selectedPatTausLeadTrk,
      selectedPatTausLeadTrkPt,
      selectedPatTausTaNCdiscr,
      selectedPatTausProng,
      selectedPatTausCharge,
      selectedPatTausMuonVeto,
      selectedPatTausElectronVeto,
      selectedPatTausEcalCrackVeto ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelCumulative = True,
    doSelIndividual = True
)

selectPatTaus = patTauSelConfigurator.configure(pyNameSpace = locals())
#
# define collections of pat::(PF)Taus used in semi-leptonic e + tau-jet channel
# (require electron and tau-jet candidates to be separated in eta-phi,
#  in order to avoid double-counting the same particle both as electron and as tau-jet;
#  apply anti-electron veto only; no need to apply anti-muon veto)
#
selectedPatTausForElecTauAntiOverlapWithElectronsVeto.dRmin = cms.double(0.3)
selectedPatTausForElecTauEta.cut = selectedPatTausEta23.cut
selectedPatTausForElecTauPt.cut = selectedPatTausPt20.cut
selectedPatTausForElecTauDecayModeFinding.cut = cms.string('tauID("decayModeFinding") > 0.5')
selectedPatTausForElecTauLeadTrkPt.cut = cms.string('tauID("decayModeFinding") > 0.5')
selectedPatTausForElecTauIso.cut = cms.string('tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5')
selectedPatTausForElecTauElectronVeto.cut = cms.string('tauID("againstElectronTight") > 0.5')
selectedPatTausForElecTauEcalCrackVeto.cut =  selectedPatTausEcalCrackVeto.cut
selectedPatTausForElecTauMuonVeto.cut = selectedPatTausMuonVeto.cut

patTauSelConfiguratorForElecTau = objSelConfigurator(
    [ selectedPatTausForElecTauAntiOverlapWithElectronsVeto,
      selectedPatTausForElecTauEta,
      selectedPatTausForElecTauPt,
      selectedPatTausForElecTauDecayModeFinding,
      selectedPatTausForElecTauLeadTrkPt,
      selectedPatTausForElecTauIso,
      selectedPatTausForElecTauElectronVeto,
      selectedPatTausForElecTauEcalCrackVeto,
      selectedPatTausForElecTauMuonVeto ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatTausForElecTau = patTauSelConfiguratorForElecTau.configure(pyNameSpace = locals())
#
# define collections of pat::(PF)Taus used in semi-leptonic mu + tau-jet channel
# (require muon and tau-jet candidates to be separated in eta-phi,
#  in order to avoid double-counting the same particle both as muon and as tau-jet;
#  apply anti-muon veto only; no need to apply anti-electron veto)
#
selectedPatTausForMuTauAntiOverlapWithMuonsVeto.dRmin = cms.double(0.3)
selectedPatTausForMuTauEta23.cut = selectedPatTausEta23.cut
selectedPatTausForMuTauPt20.cut = selectedPatTausPt20.cut
selectedPatTausForMuTauLeadTrk.cut = selectedPatTausLeadTrk.cut
selectedPatTausForMuTauLeadTrkPt.cut = selectedPatTausLeadTrkPt.cut
selectedPatTausForMuTauTaNCdiscr.cut = cms.string('tauID("byTaNCloose") > 0.5')
selectedPatTausForMuTauProng.cut = selectedPatTausProng.cut
selectedPatTausForMuTauCharge.cut = selectedPatTausCharge.cut
selectedPatTausForMuTauMuonVeto.cut = selectedPatTausMuonVeto.cut
#selectedPatTausForMuTauCaloMuonVeto.cut = cms.string('tauID("againstCaloMuon") > 0.5')
selectedPatTausForMuTauCaloMuonVeto.cut = selectedPatTausMuonVeto.cut # CV: caloMuon veto disabled for now
selectedPatTausForMuTauElectronVeto.cut = selectedPatTausElectronVeto.cut

patTauSelConfiguratorForMuTau = objSelConfigurator(
    [ selectedPatTausForMuTauAntiOverlapWithMuonsVeto,
      selectedPatTausForMuTauEta23,
      selectedPatTausForMuTauPt20,
      selectedPatTausForMuTauLeadTrk,
      selectedPatTausForMuTauLeadTrkPt,
      selectedPatTausForMuTauTaNCdiscr,
      selectedPatTausForMuTauProng,
      selectedPatTausForMuTauCharge,
      selectedPatTausForMuTauMuonVeto,
      selectedPatTausForMuTauCaloMuonVeto,
      selectedPatTausForMuTauElectronVeto ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatTausForMuTau = patTauSelConfiguratorForMuTau.configure(pyNameSpace = locals())
#
# define collections of pat::(PF)Taus used in pure hadronic tau-jet + tau-jet channel
# (no need to apply anti-electron or anti-muon vetos)
#
selectedPatTausForDiTauEta23.cut = selectedPatTausEta23.cut
selectedPatTausForDiTauPt20.cut = selectedPatTausPt20.cut
selectedPatTausForDiTauLeadTrk.cut = selectedPatTausLeadTrk.cut
selectedPatTausForDiTau1stLeadTrkPt.cut = cms.string('leadPFChargedHadrCand().isNonnull() & leadPFChargedHadrCand().pt() > 12.')
selectedPatTausForDiTau1stTaNCdiscr.cut = cms.string('tauID("byTaNCmedium") > 0.5')
selectedPatTausForDiTau1stProng.cut = selectedPatTausProng.cut
selectedPatTausForDiTau1stCharge.cut = selectedPatTausCharge.cut
selectedPatTausForDiTau1stMuonVeto.cut = selectedPatTausMuonVeto.cut
selectedPatTausForDiTau1stElectronVeto.cut = selectedPatTausElectronVeto.cut
selectedPatTausForDiTau2ndLeadTrkPt.cut = cms.string('leadPFChargedHadrCand().isNonnull() & leadPFChargedHadrCand().pt() > 8.')
selectedPatTausForDiTau2ndTaNCdiscr.cut = cms.string('tauID("byTaNCmedium") > 0.5')
selectedPatTausForDiTau2ndProng.cut = selectedPatTausProng.cut
selectedPatTausForDiTau2ndCharge.cut = selectedPatTausCharge.cut
selectedPatTausForDiTau2ndMuonVeto.cut = selectedPatTausMuonVeto.cut
selectedPatTausForDiTau2ndElectronVeto.cut = selectedPatTausElectronVeto.cut

patTauSelConfiguratorForDiTau1st = objSelConfigurator(
    [ selectedPatTausForDiTauEta23,
      selectedPatTausForDiTauPt20,
      selectedPatTausForDiTauLeadTrk,
      selectedPatTausForDiTau1stLeadTrkPt,
      selectedPatTausForDiTau1stTaNCdiscr,
      selectedPatTausForDiTau1stProng,
      selectedPatTausForDiTau1stCharge,
      selectedPatTausForDiTau1stMuonVeto,
      selectedPatTausForDiTau1stElectronVeto ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatTausForDiTau1st = patTauSelConfiguratorForDiTau1st.configure(pyNameSpace = locals())

patTauSelConfiguratorForDiTau2nd = objSelConfigurator(
    [ selectedPatTausForDiTau2ndLeadTrkPt,
      selectedPatTausForDiTau2ndTaNCdiscr,
      selectedPatTausForDiTau2ndProng,
      selectedPatTausForDiTau2ndCharge,
      selectedPatTausForDiTau2ndMuonVeto,
      selectedPatTausForDiTau2ndElectronVeto ],
    src = "selectedPatTausForDiTauLeadTrkCumulative",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatTausForDiTau2nd = patTauSelConfiguratorForDiTau2nd.configure(pyNameSpace = locals())

selectedPatTausForDiTau2ndLeadTrkPtLoose.cut = cms.string('leadPFChargedHadrCand().isNonnull() & leadPFChargedHadrCand().pt() > 5.')
selectedPatTausForDiTau2ndTaNCdiscrLoose.cut = cms.string('tauID("byTaNCloose") > -1.')
selectedPatTausForDiTau2ndProngLoose.cut = cms.string('signalPFChargedHadrCands.size() > -1')
selectedPatTausForDiTau2ndChargeLoose.cut = cms.string('abs(charge) > -1.')
selectedPatTausForDiTau2ndMuonVetoLoose.cut = selectedPatTausForDiTau2ndMuonVeto.cut
selectedPatTausForDiTau2ndElectronVetoLoose.cut = selectedPatTausForDiTau2ndElectronVeto.cut

patTauSelConfiguratorForDiTau2ndLoose = objSelConfigurator(
    [ selectedPatTausForDiTau2ndLeadTrkPtLoose,
      selectedPatTausForDiTau2ndTaNCdiscrLoose,
      selectedPatTausForDiTau2ndProngLoose,
      selectedPatTausForDiTau2ndChargeLoose,
      selectedPatTausForDiTau2ndMuonVetoLoose,
      selectedPatTausForDiTau2ndElectronVetoLoose ],
    src = "selectedPatTausForDiTauLeadTrkCumulative",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectPatTausForDiTau2ndLoose = patTauSelConfiguratorForDiTau2ndLoose.configure(pyNameSpace = locals())

selectPatTausForDiTau = cms.Sequence(selectPatTausForDiTau1st * selectPatTausForDiTau2nd * selectPatTausForDiTau2ndLoose)

# define collections of pat::(PF)Taus used in W->tau-jet + nu channel
selectedPatTausForWTauNuEta21.cut = selectedPatTausEta23.cut
selectedPatTausForWTauNuPt20.cut = cms.string("pt > 30.")
selectedPatTausForWTauNuTrkMatchVertex = cms.EDFilter("PATTauDzSelector",
                                                      vertexSource = cms.InputTag('selectedPrimaryVertexHighestPtTrackSum'),
                                                      dzMax = cms.double(0.2),
                                                      filter = cms.bool(False)
                                                      )
selectedPatTausForWTauNuLeadTrk.cut = cms.string('leadTrack().isNonnull')
selectedPatTausForWTauNuLeadTrkPt.cut = cms.string('leadTrack().isNonnull() & leadTrack().pt() > 15.')
selectedPatTausForWTauNuIso.cut = cms.string("tauID('byHPSmedium') > 0.5")
selectedPatTausForWTauNuProng.cut = selectedPatTausProng.cut
selectedPatTausForWTauNuCharge.cut = selectedPatTausCharge.cut
selectedPatTausForWTauNuMuonVeto.cut = selectedPatTausMuonVeto.cut
selectedPatTausForWTauNuElectronVeto.cut = selectedPatTausElectronVeto.cut
selectedPatTausForWTauNuEmFraction.cut = cms.string("emFraction < 0.90")
selectedPatTausForWTauNuEcalCrackVeto.cut = selectedPatTausEcalCrackVeto.cut

patTauSelConfiguratorForWTauNu = objSelConfigurator(
    [ selectedPatTausForWTauNuEta21,
      selectedPatTausForWTauNuPt20,
      selectedPatTausForWTauNuTrkMatchVertex,
      selectedPatTausForWTauNuLeadTrk,
      selectedPatTausForWTauNuLeadTrkPt,
      selectedPatTausForWTauNuMuonVeto,
      selectedPatTausForWTauNuElectronVeto,
      selectedPatTausForWTauNuEmFraction,
      selectedPatTausForWTauNuIso,
      selectedPatTausForWTauNuProng,
      selectedPatTausForWTauNuCharge,
      selectedPatTausForWTauNuEcalCrackVeto ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelIndividual = True
)
selectPatTausForWTauNu = patTauSelConfiguratorForWTauNu.configure(pyNameSpace = locals())

producePatSelLeptons = cms.Sequence (
    selectPatElectrons + selectPatElectronsLooseIsolation
   + selectPatMuons + selectPatMuonsLooseIsolation
   + selectPatTaus
   + selectPatElectronsForElecTau + selectPatElectronsForElecTauLooseIsolation
   + selectPatElectronsForElecMu + selectPatElectronsForElecMuLooseIsolation
   + selectPatTausForElecTau + selectPatTausForMuTau + selectPatTausForDiTau
   + selectPatTausForWTauNu 
)
