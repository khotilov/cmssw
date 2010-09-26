import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.Configuration.analyzeAHtoMuTau_cfi import *
from TauAnalysis.DQMTools.tools.drawJobConfigurator import *

# define template for all kinds of plots
# (specific to MSSM Higgs A/H --> mu + tau-jet analysis)
plots_AHtoMuTau = cms.PSet(
    plots = cms.PSet(  
        dqmMonitorElements = cms.vstring(''),
        processes = cms.vstring(
            'Zmumu',
            'WplusJets',
            'TTplusJets',
            'qcdSum',
            #'ZplusJets',
            'Vqq',
            'Ztautau',
            'AHsum120_tautau'
        )
    ),
    xAxis = cms.string('unlabeled'),
    #yAxis = cms.string('numEntries_linear'),
    yAxis = cms.string('numEntries_log'),
    legend = cms.string('regular'),
    labels = cms.vstring('mcNormScale'),                   
    drawOptionSet = cms.string('default'),
    stack = cms.vstring(
        'Zmumu',
        'WplusJets',
        'TTplusJets',
        'qcdSum',
        #'ZplusJets',
        'Vqq',
        'Ztautau'
    )
)

drawJobConfigurator_AHtoMuTau = drawJobConfigurator(
    template = plots_AHtoMuTau,
    dqmDirectory = ''
)

#--------------------------------------------------------------------------------
# define cut-flow control plots common to "centralJetVeto" and "centralJetBtag" analysis paths;
# show distribution of each quantity used in event selection
# (**before** quantity is cutted on)
#--------------------------------------------------------------------------------

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelPrimaryEventVertex,
    beforeCut = evtSelPrimaryEventVertexQuality,
    plot = drawJobConfigEntry(
        meName = 'VertexQuantities/VertexChi2Prob',
        title = "P(#Chi^{2}_{vtx} (after primary Event Vertex Cut)",
        xAxis = 'prob',
        name = "cutFlowControlPlots_vertexChi2Prob_afterPrimaryEventVertex"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelPrimaryEventVertexQuality,
    beforeCut = evtSelPrimaryEventVertexPosition,
    plot = drawJobConfigEntry(
        meName = 'VertexQuantities/VertexZ',
        title = "z_{vtx} (after primary Event Vertex quality Cut)",
        xAxis = 'posZ',
        name = "cutFlowControlPlots_vertexZ_afterPrimaryEventVertexQuality"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelPrimaryEventVertexPosition,
    beforeCut = evtSelGlobalMuon,
    plot = drawJobConfigEntry(
        meName = 'MuonQuantities/Muon#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Muon (after primary Event Vertex position Cut)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_muon_afterPrimaryEventVertexPosition"
    )
)    

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelGlobalMuon,
    beforeCut = evtSelMuonEta,
    plot = drawJobConfigEntry(
        meName = 'MuonQuantities/Muon#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Muon (after global Muon Cut)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_muon_afterGlobalMuon"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelMuonEta,
    beforeCut = evtSelMuonPt,
    plot = drawJobConfigEntry(
        meName = 'MuonQuantities/Muon#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Muon (after Muon #eta Cut)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_muon_afterMuonEta"
    )
)    

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelMuonPt,
    beforeCut = evtSelTauAntiOverlapWithMuonsVeto,
    plot = drawJobConfigEntry(
        meName = 'TauQuantities/Tau#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Tau (after Muon P_{T} Cut)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_tau_afterMuonPt"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelTauAntiOverlapWithMuonsVeto,
    beforeCut = evtSelTauEta,
    plots = [
        drawJobConfigEntry(
            meName = 'TauQuantities/Tau#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Tau (after Muon-Tau overlap Veto)",
            xAxis = '#PAR#',
            name = "cutFlowControlPlots_tau_afterTauAntiOverlapWithMuonsVeto"
        ),
        drawJobConfigEntry(
            meName = 'TauQuantities/TauLeadTrkPt',
            title = "Tau lead. Track (after Muon-Tau overlap Veto)",
            xAxis = 'Pt',
            name = "cutFlowControlPlots_tauLeadTrkPt_afterTauAntiOverlapWithMuonsVeto"
        )
    ]
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelTauEta,
    beforeCut = evtSelTauPt,
    plots = [
        drawJobConfigEntry(
            meName = 'TauQuantities/Tau#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Tau (after Tau #eta Cut)",
            xAxis = '#PAR#',
            name = "cutFlowControlPlots_tau_afterTauEta"
        ),
        drawJobConfigEntry(
            meName = 'TauQuantities/TauLeadTrkPt',
            title = "Tau lead. Track (after Tau #eta Cut)",
            xAxis = 'Pt',
            name = "cutFlowControlPlots_tauLeadTrkPt_afterTauEta"
        )
    ]
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelTauPt,
    beforeCut = evtSelMuonVbTfId,
    plots = [
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonTypeTracker',
            title = "Muon type 'Tracker' (after Tau P_{T} Cut)",
            xAxis = 'unlabeled',
            name = "cutFlowControlPlots_muonTypeTracker_afterTauPt"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonTypeTracker',
            title = "Muon global Track #Chi^{2} (after Tau P_{T} Cut)",
            xAxis = 'Chi2',
            name = "cutFlowControlPlots_muonGlobalTrackChi2_afterTauPt"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonGlobalTrackIPxyBeamSpot',
            title = "Muon global Track IP_{xy} wrt. BeamSpot (after Tau P_{T} Cut)",
            xAxis = 'IPxy',
            name = "cutFlowControlPlots_muonGlobalTrackIPxyBeamSpot_afterTauPt"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonNumTrackerHits',
            title = "Muon num. Tracker Hits (after Tau P_{T} Cut)",
            xAxis = 'N',
            name = "cutFlowControlPlots_muonNumTrackerHits_afterTauPt"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonNumPixelHits',
            title = "Muon num. Pixel Hits (after Tau P_{T} Cut)",
            xAxis = 'N',
            name = "cutFlowControlPlots_muonNumPixelHits_afterTauPt"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonNumberOfChambers',
            title = "Muon num. Chambers with Hits (after Tau P_{T} Cut)",
            xAxis = 'N',
            name = "cutFlowControlPlots_muonNumChamberHits_afterTauPt"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonNumberOfSegments',
            title = "Muon num. Segments (after Tau P_{T} Cut)",
            xAxis = 'N',
            name = "cutFlowControlPlots_muonNumSegments_afterTauPt"
        )
    ]
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelMuonVbTfId,
    beforeCut = evtSelMuonTrkIso,
    plots = [
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonPFChargedHadronIsoPt',
            title = "Muon Track iso. (after Muon VBTF id. Cut)",
            xAxis = 'Pt',
            name = "cutFlowControlPlots_muonAbsTrkIso_afterMuonVbTfId"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonPFChargedHadronIsoPtRel',
            title = "Muon Track rel. iso. (after Muon VBTF id. Cut)",
            xAxis = 'unlabeled',
            name = "cutFlowControlPlots_muonRelTrkIso_afterMuonVbTfId"
        )
    ]
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelMuonTrkIso,
    beforeCut = evtSelMuonEcalIso,
    plots = [
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonPFGammaIsoPt',
            title = "Muon ECAL iso. (after Muon Track iso. Cut)",
            xAxis = 'Pt',
            name = "cutFlowControlPlots_muonAbsEcalIso_afterMuonTrkIso"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonPFGammaIsoPtRel',
            title = "Muon ECAL rel. iso. (after Muon Track iso. Cut)",
            xAxis = 'unlabeled',
            name = "cutFlowControlPlots_muonRelEcalIso_afterMuonTrkIso"
        )
    ]
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelMuonEcalIso,
    beforeCut = evtSelMuonCombIso,
    plots = [
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonParticleFlowIsoPt',
            title = "Muon Combined iso. (after Muon Track iso. Cut)",
            xAxis = 'Pt',
            name = "cutFlowControlPlots_muonAbsCombIso_afterMuonEcalIso"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonParticleFlowIsoPtRel',
            title = "Muon Combined rel. iso. (after Muon Track iso. Cut)",
            xAxis = 'unlabeled',
            name = "cutFlowControlPlots_muonRelCombIso_afterMuonEcalIso"
        )
    ]
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelMuonCombIso,
    beforeCut = evtSelMuonAntiPion,
    plot = drawJobConfigEntry(
        meName = 'MuonQuantities/Muon#PAR#Compatibility',
        PAR = [ 'Calo', 'Segment' ],
        title = "Muon #PAR# compatibility (after Muon Combined iso. Cut)",
        xAxis = 'prob',
        name = "cutFlowControlPlots_muonComp_afterMuonCombIso"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelMuonAntiPion,
    beforeCut = evtSelMuonTrkIP,
    plot = drawJobConfigEntry(
        meName = 'MuonQuantities/MuonTrackIP#PAR#',
        PAR = [ 'xy', 'z' ],
        title = "Muon Track IP_{#PAR#}(after Muon #pi-Veto Cut)",
        xAxis = 'IP#PAR#',
        name = "cutFlowControlPlots_muonTrkIP_afterMuonAntiPionVeto"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelMuonTrkIP,
    beforeCut = evtSelTauLeadTrk,
    plots = [
        drawJobConfigEntry(
            meName = 'TauQuantities/Tau#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Tau (after Muon Track IP_{xy} Cut)",
            xAxis = '#PAR#',
            name = "cutFlowControlPlots_tau_afterMuonTrkIP"
        ),
        drawJobConfigEntry(
            meName = 'TauQuantities/TauLeadTrkPt',
            title = "Tau lead. Track (after Muon Track IP_{xy} Cut)",
            xAxis = 'Pt',
            name = "cutFlowControlPlots_tauLeadTrkPt_afterMuonTrkIP"
        )
    ]
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelTauLeadTrk,
    beforeCut = evtSelTauLeadTrkPt,
    plots = [
        drawJobConfigEntry(
            meName = 'TauQuantities/Tau#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Tau (after Tau lead. Track Cut)",
            xAxis = '#PAR#',
            name = "cutFlowControlPlots_tau_afterTauLeadTrk"
        ),
        drawJobConfigEntry(
            meName = 'TauQuantities/TauLeadTrkPt',
            title = "Tau lead. Track (after Tau lead. Track Cut)",
            xAxis = 'Pt',
            name = "cutFlowControlPlots_tauLeadTrkPt_afterTauLeadTrk"
       )
   ]
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelTauLeadTrkPt,
    beforeCut = evtSelTauTaNCdiscr,
    plot = drawJobConfigEntry(
        meName = 'TauQuantities/TauDiscriminatorTaNCfrQuarterPercent',
        title = "Tau TaNC discr. (after Tau lead. Track P_{T} Cut)",
        xAxis = 'Pt',
        name = "cutFlowControlPlots_tauTaNCdiscr_afterTauLeadTrkPt"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelTauTaNCdiscr,
    beforeCut = evtSelTauTrkIso,
    plot = drawJobConfigEntry(
        meName = 'TauQuantities/TauTrkIsoPt',
        title = "Tau Track iso. (after Tau TaNC discr. Cut)",
        xAxis = 'Pt',
        name = "cutFlowControlPlots_tauTrkIso_afterTauTaNCdiscr"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelTauTrkIso,
    beforeCut = evtSelTauEcalIso,
    plot = drawJobConfigEntry(
        meName = 'TauQuantities/TauEcalIsoPt',
        title = "Tau ECAL iso. (after Tau Track iso. Cut)",
        xAxis = 'Pt',
        name = "cutFlowControlPlots_tauEcalIso_afterTauTrkIso"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelTauEcalIso,
    beforeCut = evtSelTauProng,
    plot = drawJobConfigEntry(
        meName = 'TauQuantities/TauNumTracksSignalCone',
        title = "Tau Tracks in Signal Cone (after Tau ECAL iso. Cut)",
        xAxis = 'unlabeled',
        name = "cutFlowControlPlots_tauNumTracksSignalCone_afterTauEcalIso"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelTauProng,
    beforeCut = evtSelTauCharge,
    plot = drawJobConfigEntry(
        meName = 'TauQuantities/TauCharge',
        title = "Tau Charge (#Sigma Tracks in Signal Cone, after Tau 1-Prong||3-Prong Cut)",
        xAxis = 'unlabeled',
        name = "cutFlowControlPlots_tauCharge_afterTauProng"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelTauCharge,
    beforeCut = evtSelTauMuonVeto,
    plots = [
        drawJobConfigEntry(
            meName = 'TauQuantities/Tau#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Tau (after Charge(Tau) = +/-1 Cut)",
            xAxis = '#PAR#',
            name = "cutFlowControlPlots_tau_afterTauCharge"
        ),
        drawJobConfigEntry(
            meName = 'TauQuantities/TauDiscriminatorAgainstMuons',
            title = "Tau anti-Muon Discr. (after Charge(Tau) = +/-1 Cut)",
            xAxis = 'unlabeled',
            name = "cutFlowControlPlots_tauAntiMuonDiscr_afterTauCharge"
        )
    ]
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelTauMuonVeto,
    beforeCut = evtSelTauElectronVeto,
    plots = [
        drawJobConfigEntry(
            meName = 'TauQuantities/Tau#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Tau (after Tau #mu-Veto Cut)",
            xAxis = '#PAR#',
            name = "cutFlowControlPlots_tau_afterTauMuonVeto"
        ),
        drawJobConfigEntry(
            meName = 'TauQuantities/TauDiscriminatorAgainstElectrons',
            title = "Tau anti-Electron Discr. (after Tau #mu-Veto Cut)",
            xAxis = 'unlabeled',
            name = "cutFlowControlPlots_tauAntiElectronDiscr_afterTauMuonVeto"
        )
    ]
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelTauElectronVeto,
    beforeCut = evtSelDiTauCandidateForAHtoMuTauAntiOverlapVeto,
    plot = drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/DR12',
        title = "#Delta R(Muon,Tau) (after Tau e-Veto Cut)",
        xAxis = 'dR',
        name = "cutFlowControlPlots_dR12_afterTauElectronVeto"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelDiTauCandidateForAHtoMuTauAntiOverlapVeto,
    beforeCut = evtSelDiTauCandidateForAHtoMuTauZeroCharge,
    plot = drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/DiTauCandidateCharge',
        title = "Charge(Muon + Tau) (after diTau anti-Overlap Veto)",
        xAxis = 'unlabeled',
        name = "cutFlowControlPlots_diTauCharge_afterAntiOverlapVeto"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelDiTauCandidateForAHtoMuTauZeroCharge,
    beforeCut = evtSelDiTauCandidateForAHtoMuTauMt1MET,
    plot = drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/Mt1MET',
        title = "M_{T}(Muon + MET) (after opposite Charge Cut)",
        xAxis = 'Mt',
        name = "cutFlowControlPlots_mtMuonMET_afterZeroCharge"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelDiTauCandidateForAHtoMuTauMt1MET,
    beforeCut = evtSelDiTauCandidateForAHtoMuTauPzetaDiff,
    plot = drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/PzetaDiff',
        title = "P_{#zeta} - 1.5*P_{#zeta}^{vis} (after transverse Mass Cut)",
        xAxis = 'GeV',
        name = "cutFlowControlPlots_PzetaDiff_afterMt1MET"
    )
)

drawJobConfigurator_AHtoMuTau.add(
    afterCut = evtSelDiTauCandidateForAHtoMuTauPzetaDiff,
    beforeCut = evtSelDiMuPairZmumuHypothesisVeto,
    plot = drawJobConfigEntry(
	meName = 'DiMuZmumuHypothesisQuantities/VisMass',
        title = "M(Muon + Muon, Z #rightarrow #mu^{+} #mu^{-} Mass hypothesis) (after P_{#zeta} Cut)",
        xAxis = 'Mass',
        name = "cutFlowControlPlots_mZmumuHypothesis_afterPzetaDiff"
    )
)

#--------------------------------------------------------------------------------
# define cut-flow control plots specific to "non-b-tag" analysis path
#--------------------------------------------------------------------------------

drawJobConfigurator_AHtoMuTau_woBtag = copy.deepcopy(drawJobConfigurator_AHtoMuTau)
drawJobConfigurator_AHtoMuTau_woBtag.setDQMdirectory('#PROCESSDIR#/ahMuTauAnalyzer_woBtag/')

drawJobConfigurator_AHtoMuTau_woBtag.add(
    afterCut = evtSelDiMuPairZmumuHypothesisVeto,
    beforeCut = evtSelNonCentralJetEt20bTag,
    plot = drawJobConfigEntry(
	meName = 'JetQuantities/Jet#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Jet (after Z #rightarrow #mu^{+} #mu^{-} Mass hypothesis Veto)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_jet_afterZmumuMassHypothesisVeto"
    )
)

#--------------------------------------------------------------------------------
# define cut-flow control plots specific to "b-tag" analysis path
#--------------------------------------------------------------------------------

drawJobConfigurator_AHtoMuTau_wBtag = copy.deepcopy(drawJobConfigurator_AHtoMuTau)
drawJobConfigurator_AHtoMuTau_wBtag.setDQMdirectory('#PROCESSDIR#/ahMuTauAnalyzer_wBtag/')

drawJobConfigurator_AHtoMuTau_wBtag.add(
    afterCut = evtSelDiMuPairZmumuHypothesisVeto,
    beforeCut = evtSelCentralJetEt20,
    plot = drawJobConfigEntry(
	meName = 'JetQuantities/Jet#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Jet (after Z #rightarrow #mu^{+} #mu^{-} Mass hypothesis Veto)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_jet_afterZmumuMassHypothesisVeto"
    )
)

drawJobConfigurator_AHtoMuTau_wBtag.add(
    afterCut = evtSelCentralJetEt20,
    beforeCut = evtSelCentralJetEt20bTag,
    plots = [
        drawJobConfigEntry(
            meName = 'JetQuantities/Jet#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Jet (after Jet P{T} and #eta Cuts)",
            xAxis = '#PAR#',
            name = "cutFlowControlPlots_jet_afterJetPtAndEta"
        ),
        drawJobConfigEntry(
            meName = 'JetQuantities/BtagDisc_trackCountingHighEffBJetTags',
            title = "Jet b-Tag Discr. (after Jet P{T} and #eta Cuts)",
            xAxis = 'unlabeled',
            name = "cutFlowControlPlots_jetBtagDiscr_afterJetPtAndEta"
        ),
        drawJobConfigEntry(
            meName = 'JetQuantities/NumBtags_trackCountingHighEffBJetTags',
            title = "Num. Jets with b-Tag (after Jet P{T} and #eta Cuts)",
            xAxis = 'unlabeled',
            name = "cutFlowControlPlots_numBtagJets_afterJetPtAndEta"
        )
    ]
    
)

#--------------------------------------------------------------------------------
# define distributions to be plotted for events passing all event selection criteria
#--------------------------------------------------------------------------------

finalSamplePlots = [
    drawJobConfigEntry(
        meName = 'MuonQuantities/Muon#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Muon (final Event sample)",
        xAxis = '#PAR#',
        name = "finalSamplePlots_muon"
    ),
    drawJobConfigEntry(
        meName = 'MuonQuantities/MuonMatchingGenParticlePdgId',
        title = "PdgId of gen. Particle matching Muon (final Event sample)",
        xAxis = 'PdgId',
        name = "finalSamplePlots_pdgIdGenParticleMatchingMuon"
    ),
    drawJobConfigEntry(
        meName = 'TauQuantities/Tau#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Tau (final Event sample)",
        xAxis = '#PAR#',
        name = "finalSamplePlots_tau"
    ),
    drawJobConfigEntry(
        meName = 'TauQuantities/TauMatchingGenParticlePdgId',
        title = "PdgId of gen. Particle matching Tau (final Event sample)",
        xAxis = 'PdgId',
        name = "finalSamplePlots_pdgIdGenParticleMatchingTau"
    ),
    drawJobConfigEntry(
        meName = 'TauQuantities/TauLeadTrkPt',
        title = "Tau lead. Track (final Event sample)",
        xAxis = 'Pt',
        name = "finalSamplePlots_tauLeadTrkPt"
    ),
    drawJobConfigEntry(
        meName = 'TauQuantities/TauNumTracksSignalCone',
        title = "Tau Tracks in Signal Cone (final Event sample)",
        xAxis = 'unlabeled',
        name = "finalSamplePlots_tauNumTracksSignalCone"
    ),
    drawJobConfigEntry(
        meName = 'TauQuantities/TauJetRadius',
        title = "Tau Jet Radius (final Event sample)",
        xAxis = 'unlabeled',
        name = "finalSamplePlots_tauJetRadius"
    ),
    drawJobConfigEntry(
        meName = 'TauQuantities/TauDiscriminatorTaNCfrQuarterPercent',
        title = "TaNC output (fr = 0.25%) (final Event sample)",
        xAxis = 'unlabeled',
        name = "finalSamplePlots_tauDiscrTaNCfrQuarterPercent"
    ),
    drawJobConfigEntry(
        meName = 'TauQuantities/TauTaNCoutputTransform',
        title = "TaNC output (transformed) (final Event sample)",
        xAxis = 'unlabeled',
        name = "finalSamplePlots_tauTaNCtransform"
    ),
    drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/DPhi12',
        title = "#Delta#phi(Muon-Tau) (final Event sample)",
        xAxis = 'dPhi',
        name = "finalSamplePlots_dPhiMuonTau"
    ),
    drawJobConfigEntry(
        meName = 'CaloMEtQuantities/MEtPt',
        title = "CaloMET (final Event sample)",
        xAxis = 'Pt',
        name = "finalSamplePlots_caloMEt"
    ),
    drawJobConfigEntry(
        meName = 'PFMEtQuantities/MEtPt',
        title = "PFMET (final Event sample)",
        xAxis = 'Pt',
        name = "finalSamplePlots_pfMEt"
    ),
    drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/PzetaDiff',
        title = "P_{#zeta} - 1.5*P_{#zeta}^{vis} (final Event sample)",
        xAxis = 'GeV',
        name = "finalSamplePlots_PzetaDiff"
    ),
    drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/Mt1MET',
        title = "M_{T}(Muon + MET) (final Event sample)",
        xAxis = 'Mt',
        name = "finalSamplePlots_mtMuonMET"
    ),
    drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/Mt2MET',
        title = "M_{T}(Tau + MET) (final Event sample)",
        xAxis = 'Mt',
        name = "finalSamplePlots_mtTauMET"
    ),
    drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/Mt12MET',
        title = "M_{T}(Muon + Tau + MET) (final Event sample)",
        xAxis = 'Mt',
        name = "finalSamplePlots_mtMuonTauMET"
    ),
    drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/VisMass',
        title = "M_{vis}(Muon + Tau) (final Event sample)",
        xAxis = 'Mass',
        name = "finalSamplePlots_mVisible"
    ),
    drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/VisMassZllCombinedHypothesis',
        title = "M_{vis}(Muon + Tau), Z #rightarrow #ell^{+} #ell^{-} combined Hypothesis (final Event sample)",
        xAxis = 'Mass',
        name = "finalSamplePlots_mVisibleZllCombinedHypothesis"
    ),
    drawJobConfigEntry(
        meName = 'DiMuZmumuHypothesisQuantities/VisMass',
        title = "M(Muon + Muon, Z #rightarrow #mu^{+} #mu^{-} Mass hypothesis) (final Event sample)",
        xAxis = 'Mass',
        name = "finalSamplePlots_mZmumuHypothesis"
    ),
    drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/CDFmethodMass',
        title = "M(Muon + Tau + MET), CDF method (final Event sample)",
        xAxis = 'Mass',
        name = "finalSamplePlots_mCDFmethod"
    ),
    drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/CollinearApproxMass',
        title = "M(Muon + Tau), collinear Approx. (final Event sample)",
        xAxis = 'Mass',
        name = "finalSamplePlots_mCollApprox"
    ),
    drawJobConfigEntry(
        meName = 'DiTauCandidateSVfitQuantities/Mass1stSolution',
        title = "M(Muon + Tau), SV method (final Event sample)",
        xAxis = 'Mass',
        name = "finalSamplePlots_mSVmethod"
    ),
    drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/Ht12MET',
        title = "#Sigma H_{T}(Muon + Tau + MET) (final Event sample)",
        xAxis = 'Mass',
        name = "finalSamplePlots_ht"
    ),
    drawJobConfigEntry(
        meName = 'JetQuantities/Jet#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Jet (final Event sample)",
        xAxis = '#PAR#',
        name = "finalSamplePlots_jet"
    ),
    drawJobConfigEntry(
        meName = 'JetQuantities/BtagDisc_trackCountingHighEffBJetTags',
        title = "Jet b-Tag Discr. (final Event sample)",
        xAxis = 'unlabeled',
        name = "finalSamplePlots_jetBtagDiscr"
    ),
    drawJobConfigEntry(
        meName = 'JetQuantities/NumJets',
        title = "Num. Jets with b-Tag (final Event sample)",
        xAxis = 'unlabeled',
        name = "finalSamplePlots_numBtagJets"
    )
]

drawJobConfigurator_AHtoMuTau_woBtag.add(
    afterCut = evtSelNonCentralJetEt20bTag,
    plots = finalSamplePlots
)

drawJobConfigurator_AHtoMuTau_wBtag.add(
    afterCut = evtSelCentralJetEt20bTag,
    plots = finalSamplePlots
)

