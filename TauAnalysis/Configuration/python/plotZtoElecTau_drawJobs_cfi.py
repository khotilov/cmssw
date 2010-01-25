import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.Configuration.analyzeZtoElecTau_cfi import *
from TauAnalysis.DQMTools.tools.drawJobConfigurator import *

# define template for all kins of plots
# (specific to Z --> e + tau-jet analysis)
plots_ZtoElecTau = cms.PSet(
    plots = cms.PSet(  
        dqmMonitorElements = cms.vstring(''),
        processes = cms.vstring(
            #'Zee',
            'ZeePlusJets',
            'WplusJets',
            'TTplusJets',
            'qcdSum',
            'gammaPlusJetsSum',
            'Ztautau'
            #'ZtautauPlusJets'
        )
    ),
    xAxis = cms.string('unlabeled'),
    yAxis = cms.string('numEntries_linear'),
    #yAxis = cms.string('numEntries_log'),
    legend = cms.string('regular'),
    labels = cms.vstring('mcNormScale'),                   
    drawOptionSet = cms.string('default'),
    stack = cms.vstring(
         #'Zee',
         'ZeePlusJets',
         'WplusJets',
         'TTplusJets',
         'qcdSum',
         'gammaPlusJetsSum',
         'Ztautau'
         #'ZtautauPlusJets' 
    )
)

drawJobConfigurator_ZtoElecTau = drawJobConfigurator(
    template = plots_ZtoElecTau,
    dqmDirectory = '#PROCESSDIR#/zElecTauAnalyzer/'
)

#--------------------------------------------------------------------------------
# define cut-flow control plots;
# show distribution of each quantity used in event selection
# (**before** quantity is cutted on)
#--------------------------------------------------------------------------------

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelPrimaryEventVertex,
    beforeCut = evtSelPrimaryEventVertexQuality,
    plot = drawJobConfigEntry(
        meName = 'VertexQuantities/VertexChi2Prob',
        title = "P(#Chi^{2}_{vtx} (after primary Event Vertex Cut)",
        xAxis = 'prob',
        name = "cutFlowControlPlots_vertexChi2Prob_afterPrimaryEventVertex"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelPrimaryEventVertexQuality,
    beforeCut = evtSelPrimaryEventVertexPosition,
    plot = drawJobConfigEntry(
        meName = 'VertexQuantities/VertexZ',
        title = "z_{vtx} (after primary Event Vertex quality Cut)",
        xAxis = 'posZ',
        name = "cutFlowControlPlots_vertexZ_afterPrimaryEventVertexQuality"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelPrimaryEventVertexPosition,
    beforeCut = evtSelTightElectronId,
    plot = drawJobConfigEntry(
        meName = 'ElectronQuantities/Electron#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Electron (after primary Event Vertex position Cut)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_electron_afterPrimaryEventVertexPosition"
    )
)    

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelTightElectronId,
    beforeCut = evtSelElectronAntiCrack,
    plot = drawJobConfigEntry(
        meName = 'ElectronQuantities/Electron#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Electron (after Electron ID Cut)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_electron_afterTightElectronId"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelElectronAntiCrack,
    beforeCut = evtSelElectronEta,
    plot = drawJobConfigEntry(
        meName = 'ElectronQuantities/Electron#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Electron (after Electron anti-crack Cut)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_electron_afterElectronAntiCrack"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelElectronEta,
    beforeCut = evtSelElectronPt,
    plot = drawJobConfigEntry(
        meName = 'ElectronQuantities/Electron#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Electron (after Electron anti-crack Cut)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_electron_afterElectronEta"
    )
)    

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelElectronPt,
    beforeCut = evtSelTauAntiOverlapWithElectronsVeto,
    plot = drawJobConfigEntry(
        meName = 'TauQuantities/Tau#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Tau (after Electron P_{T} Cut)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_tau_afterElectronPt"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelTauAntiOverlapWithElectronsVeto,
    beforeCut = evtSelTauEta,
    plots = [
        drawJobConfigEntry(
            meName = 'TauQuantities/Tau#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Tau (after Electron-Tau overlap Veto)",
            xAxis = '#PAR#',
            name = "cutFlowControlPlots_tau_afterTauAntiOverlapWithElectronsVeto"
        ),
        drawJobConfigEntry(
            meName = 'TauQuantities/TauLeadTrkPt',
            title = "Tau lead. Track (after Electron-Tau overlap Veto)",
            xAxis = 'Pt',
            name = "cutFlowControlPlots_tauLeadTrkPt_afterTauAntiOverlapWithElectronsVeto"
            )
    ]
)

drawJobConfigurator_ZtoElecTau.add(
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

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelTauPt,
    beforeCut = evtSelElectronTrkIso,
    plot = drawJobConfigEntry(
        meName = 'ElectronQuantities/ElectronTrkIsoPt',
        title = "Electron Track iso. (after Tau P_{T} Cut)",
        xAxis = 'Pt',
        name = "cutFlowControlPlots_electronTrkIso_afterTauPt"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelElectronTrkIso,
    beforeCut = evtSelElectronEcalIso,
    plot = drawJobConfigEntry(
        meName = 'ElectronQuantities/ElectronEcalIsoPt',
        title = "Electron ECAL iso. (after Electron Track iso. Cut)",
        xAxis = 'Pt',
        name = "cutFlowControlPlots_electronEcalIso_afterElectronTrkIso"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelElectronEcalIso,
    beforeCut = evtSelElectronTrk,
    plot = drawJobConfigEntry(
        meName = 'ElectronQuantities/Electron#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Electron (after Electron ECAL iso. Cut)",
        xAxis = 'prob',
        name = "cutFlowControlPlots_electron_afterElectronEcalIso"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelElectronTrk,
    beforeCut = evtSelElectronTrkIP,
     plot = drawJobConfigEntry(
        meName = 'ElectronQuantities/ElectronTrackIP#PAR#',
        PAR = [ 'xy', 'z' ],
        title = "Electron Track IP_{#PAR#} (after Electron Track Cut)",
        xAxis = 'IP#PAR#',
        name = "cutFlowControlPlots_electronTrkIP_afterElectronTrk"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelElectronTrkIP,
    beforeCut = evtSelElectronConversionVeto,
    plot = drawJobConfigEntry(
        meName = 'ElectronQuantities/Electron#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Electron (after Electron (after Electron Track IP_{xy} Cut)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_electron_afterElectronTrkIP"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelElectronConversionVeto,
    beforeCut = evtSelTauLeadTrk,
    plots = [
        drawJobConfigEntry(
            meName = 'TauQuantities/Tau#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Tau (after Electron #gamma-Conversion Veto)",
            xAxis = '#PAR#',
            name = "cutFlowControlPlots_tau_afterElectronConversionVeto"
        ),
        drawJobConfigEntry(
            meName = 'TauQuantities/TauLeadTrkPt',
            title = "Tau lead. Track (after Electron #gamma-Conversion Veto)",
            xAxis = 'Pt',
            name = "cutFlowControlPlots_tauLeadTrkPt_afterElectronConversionVeto"
            )
    ]
)

drawJobConfigurator_ZtoElecTau.add(
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

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelTauLeadTrkPt,
    beforeCut = evtSelTauTrkIso,
    plot = drawJobConfigEntry(
        meName = 'TauQuantities/TauTrkIsoPt',
        title = "Tau Track iso. (after Tau lead. Track P_{T} Cut)",
        xAxis = 'Pt',
        name = "cutFlowControlPlots_tauTrkIso_afterTauLeadTrkPt"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelTauTrkIso,
    beforeCut = evtSelTauEcalIso,
    plot = drawJobConfigEntry(
        meName = 'TauQuantities/TauEcalIsoPt',
        title = "Tau ECAL iso. (after Tau Track iso. Cut)",
        xAxis = 'Pt',
        name = "cutFlowControlPlots_tauEcalIso_afterTauTrkIso"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelTauEcalIso,
    beforeCut = evtSelTauProng,
    plot = drawJobConfigEntry(
        meName = 'TauQuantities/TauNumTracksSignalCone',
        title = "Tau Tracks in Signal Cone (after Tau ECAL iso. Cut)",
        xAxis = 'unlabeled',
        name = "cutFlowControlPlots_tauNumTracksSignalCone_afterTauEcalIso"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelTauProng,
    beforeCut = evtSelTauCharge,
    plot = drawJobConfigEntry(
        meName = 'TauQuantities/TauCharge',
        title = "Tau Charge (#Sigma Tracks in Signal Cone, after Tau 1-Prong||3-Prong Cut)",
        xAxis = 'unlabeled',
        name = "cutFlowControlPlots_tauCharge_afterTauProng"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelTauCharge,
    beforeCut = evtSelTauElectronVeto,
    plots = [
        drawJobConfigEntry(
            meName = 'TauQuantities/Tau#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Tau (after Charge(Tau) = +/-1 Cut)",
            xAxis = '#PAR#',
            name = "cutFlowControlPlots_tau_afterTauCharge"
        ),
        drawJobConfigEntry(
            meName = 'TauQuantities/TauDiscriminatorAgainstElectrons',
            title = "Tau anti-Electron Discr. (after Charge(Tau) = +/-1 Cut)",
            xAxis = 'unlabeled',
            name = "cutFlowControlPlots_tauAntiElectronDiscr_afterTauCharge"
        )
    ]
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelTauElectronVeto,
    beforeCut = evtSelTauEcalCrackVeto,
    plot = drawJobConfigEntry(
        meName = 'TauQuantities/Tau#PAR#',
        PAR = [ 'Pt', 'Eta', 'Phi' ],
        title = "Tau (after Tau electron-Veto cut)",
        xAxis = '#PAR#',
        name = "cutFlowControlPlots_tau_afterTauElectronVeto"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelTauEcalCrackVeto,
    beforeCut = evtSelDiTauCandidateForElecTauAntiOverlapVeto,
    plot = drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/DR12',
        title = "#Delta R(Electron,Tau) (after Tau ECAL Crack Veto)",
        xAxis = 'dR',
        name = "cutFlowControlPlots_dR12_afterTauEcalCrackVeto"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelDiTauCandidateForElecTauAntiOverlapVeto,
    beforeCut = evtSelDiTauCandidateForElecTauZeroCharge,
    plot = drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/DiTauCandidateCharge',
        title = "Charge(Electron + Tau) (after diTau anti-Overlap Veto)",
        xAxis = 'unlabeled',
        name = "cutFlowControlPlots_diTauCharge_afterAntiOverlapVeto"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelDiTauCandidateForElecTauZeroCharge,
    beforeCut = evtSelDiTauCandidateForElecTauAcoplanarity12,
    plot = drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/DPhi12',
        title = "#Delta#phi(Electron-Tau) (after opposite Charge Cut)",
        xAxis = 'dPhi',
        name = "cutFlowControlPlots_dPhiElectronTau_afterZeroCharge"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelDiTauCandidateForElecTauAcoplanarity12,
    beforeCut = evtSelDiTauCandidateForElecTauMt1MET,
    plot = drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/Mt1MET',
        title = "M_{T}(Electron + MET) (after Acoplanarity(Electron-Tau) Cut)",
        xAxis = 'Mt',
        name = "cutFlowControlPlots_mtElectronMET_afterAcoplanarityElectronTau"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelDiTauCandidateForElecTauMt1MET,
    beforeCut = evtSelDiTauCandidateForElecTauPzetaDiff,
    plot = drawJobConfigEntry(
        meName = 'DiTauCandidateQuantities/PzetaDiff',
        title = "P_{#zeta} - 1.5*P_{#zeta}^{vis} (after transverse Mass Cut)",
        xAxis = 'GeV',
        name = "cutFlowControlPlots_PzetaDiff_afterMt1MET"
    )
)

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelDiTauCandidateForElecTauPzetaDiff,
    beforeCut = evtSelElecTauPairZeeHypothesisVeto,
    plot = drawJobConfigEntry(
	meName = 'DiTauCandidateZeeHypothesisQuantities/VisMassBestMach',
        title = "M_{vis}(Electron + Tau, Z #rightarrow e^{+} e^{-} Mass hypothesis) (after P_{#zeta} Cut)",
        xAxis = 'Mass',
        name = "cutFlowControlPlots_mVisibleZeeHypothesis_afterPzetaDiff"
    )
)

#--------------------------------------------------------------------------------
# define distributions to be plotted
# for events passing all event selection criteria
#--------------------------------------------------------------------------------

drawJobConfigurator_ZtoElecTau.add(
    afterCut = evtSelElecTauPairZeeHypothesisVeto,
    plots = [
        drawJobConfigEntry(
            meName = 'ElectronQuantities/Electron#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Electron (final Event sample)",
            xAxis = '#PAR#',
            name = "finalSamplePlots_electron"
        ),
        drawJobConfigEntry(
            meName = 'ElectronQuantities/ElectronMatchingGenParticlePdgId',
            title = "PdgId of gen. Particle matching Electron (final Event sample)",
            xAxis = 'PdgId',
            name = "finalSamplePlots_pdgIdGenParticleMatchingElectron"
        ),
        drawJobConfigEntry(
            meName = 'ElectronQuantities/ElectronEcalIsoPtBarrel',
            title = "Electron ECAL iso., Barrel (final Event sample)",
            xAxis = 'Pt',
            name = "finalSamplePlots_electronEcalIsoBarrel"
        ),
        drawJobConfigEntry(
            meName = 'ElectronQuantities/ElectronEcalIsoPtEndcap',
            title = "Electron ECAL iso., Endcap (final Event sample)",
            xAxis = 'Pt',
            name = "finalSamplePlots_electronEcalIsoEndcap"
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
            title = "#Delta#phi(Electron-Tau) (final Event sample)",
            xAxis = 'dPhi',
            name = "finalSamplePlots_dPhiElectronTau"
        ),
        drawJobConfigEntry(
            meName = 'CaloMEtQuantities/RAWplusJESplusMUONplusTAU_MEtPt',
            title = "MET (final Event sample)",
            xAxis = 'Pt',
            name = "finalSamplePlots_met"
        ),
        drawJobConfigEntry(
            meName = 'DiTauCandidateQuantities/PzetaDiff',
            title = "P_{#zeta} - 1.5*P_{#zeta}^{vis} (final Event sample)",
            xAxis = 'GeV',
            name = "finalSamplePlots_PzetaDiff"
        ),
        drawJobConfigEntry(
            meName = 'DiTauCandidateQuantities/Mt1MET',
            title = "M_{T}(Electron + MET) (final Event sample)",
            xAxis = 'Mt',
            name = "finalSamplePlots_mtElectronMET"
        ),
        drawJobConfigEntry(
            meName = 'DiTauCandidateQuantities/Mt2MET',
            title = "M_{T}(Tau + MET) (final Event sample)",
            xAxis = 'Mt',
            name = "finalSamplePlots_mtTauMET"
        ),
        drawJobConfigEntry(
            meName = 'DiTauCandidateQuantities/Mt12MET',
            title = "M_{T}(Electron + Tau + MET) (final Event sample)",
            xAxis = 'Mt',
            name = "finalSamplePlots_mtElectronTauMET"
        ),
        drawJobConfigEntry(
            meName = 'DiTauCandidateQuantities/VisMass',
            title = "M_{vis}(Electron + Tau) (final Event sample)",
            xAxis = 'Mass',
            name = "finalSamplePlots_mVisible"
        ),
        drawJobConfigEntry(
            meName = 'DiTauCandidateQuantities/VisMassZllCombinedHypothesis',
            title = "M_{vis}(Electron + Tau), Z #rightarrow #ell^{+} #ell^{-} combined Hypothesis (final Event sample)",
            xAxis = 'Mass',
            name = "finalSamplePlots_mVisibleZllCombinedHypothesis"
        ),
        drawJobConfigEntry(
            meName = 'DiTauCandidateZeeHypothesisQuantities/VisMassBestMach',
            title = "M_{vis}(Electron + Tau, Z #rightarrow e^{+} e^{-} Mass hypothesis) (final Event sample)",
            xAxis = 'Mass',
            name = "finalSamplePlots_mVisibleZeeHypothesis"
        ),
        drawJobConfigEntry(
            meName = 'DiTauCandidateQuantities/CDFmethodMass',
            title = "M(Electron + Tau), CDF method (final Event sample)",
            xAxis = 'Mass',
            name = "finalSamplePlots_mCDFmethod"
        ),
        drawJobConfigEntry(
            meName = 'DiTauCandidateQuantities/CollinearApproxMass',
            title = "M(Electron + Tau), collinear Approx. (final Event sample)",
            xAxis = 'Mass',
            name = "finalSamplePlots_mCollApprox"
        )
    ]
)                

