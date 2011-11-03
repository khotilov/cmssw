import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
# import configuration parameters of Z --> electron + tau-jet channel
#
# WARNING: definitions from analyzeZtoElecTau_cfi.py need to be imported before any
#          A/H --> electron + tau-jet channel specific configuration parameters are defined;
#          otherwise the import of analyzeZtoElecTau_cfi.py will **overwrite**
#          A/H --> electron + tau-jet channel specific definitions
#         (due to the fact that e.g. histogram managers of 'ZtoElecTau' channel
#          and 'AHtoElecTau' channel have the same name '...ForElecTau')
#
from TauAnalysis.Configuration.analyzeZtoElecTau_cfi import *
#--------------------------------------------------------------------------------

# import config for histogram manager filling information about phase-space simulated in Monte Carlo sample
from TauAnalysis.Core.genPhaseSpaceEventInfoHistManager_cfi import *

# import config for event weight histogram manager
from TauAnalysis.Core.eventWeightHistManager_cfi import *

# import config for electron histogram manager
from TauAnalysis.Core.electronHistManager_cfi import *
electronHistManager.src = cms.InputTag("selectedPatElectronsForElecTauTrkIPcumulative")

# import config for tau histogram manager
from TauAnalysis.Core.pftauHistManager_cfi import *
tauHistManager.src = cms.InputTag('selectedPatTausForElecTauElectronVetoCumulative')
tauHistManager.useHPSclassicAlgorithm = cms.bool(True)

# import config for di-tau histogram manager
from TauAnalysis.Core.diTauCandidateHistManager_cfi import *
diTauCandidateHistManagerForElecTau = copy.deepcopy(diTauCandidateHistManager)
diTauCandidateHistManagerForElecTau.pluginName = cms.string('diTauCandidateHistManagerForElecTau')
diTauCandidateHistManagerForElecTau.pluginType = cms.string('PATElecTauPairHistManager')
diTauCandidateHistManagerForElecTau.diTauCandidateSource = cms.InputTag('selectedElecTauPairsForAHtoElecTauZeroChargeCumulative')
diTauCandidateHistManagerForElecTau.visMassHypothesisSource = cms.InputTag('')

from TauAnalysis.Core.diTauCandidateNSVfitHistManager_cfi import *
diTauCandidateNSVfitHistManagerForElecTau = copy.deepcopy(diTauCandidateNSVfitHistManager)
diTauCandidateNSVfitHistManagerForElecTau.pluginName = cms.string('diTauCandidateNSVfitHistManagerForElecTau')
diTauCandidateNSVfitHistManagerForElecTau.pluginType = cms.string('PATElecTauPairNSVfitHistManager')
#diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = cms.InputTag('selectedElecTauPairsForAHtoElecTauZeroChargeCumulative')
diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = cms.InputTag('allElecTauPairs')
diTauCandidateNSVfitHistManagerForElecTau.nSVfitEventHypotheses = cms.PSet(
    psKine_MEt_logM_fit = cms.string('psKine_MEt_logM_fit')
    #    psKine_MEt_logM_int = cms.string('psKine_MEt_logM_int')
)

from TauAnalysis.Core.diTauCandidateZllHypothesisHistManager_cfi import *
diTauCandidateZeeHypothesisHistManagerForElecTau = copy.deepcopy(ZllHypothesisHistManager)
diTauCandidateZeeHypothesisHistManagerForElecTau.pluginName = cms.string('diTauCandidateZeeHypothesisHistManagerForElecTau')
diTauCandidateZeeHypothesisHistManagerForElecTau.pluginType = cms.string('ZllHypothesisElecTauHistManager')
diTauCandidateZeeHypothesisHistManagerForElecTau.ZllHypothesisSource = cms.InputTag('elecTauPairZeeHypothesesForAHtoElecTau')
diTauCandidateZeeHypothesisHistManagerForElecTau.dqmDirectory_store = cms.string('DiTauCandidateZeeHypothesisQuantities')

# import config for Zee veto histogram manager
elecPairHistManagerByLooseIsolation = diTauCandidateHistManager.clone(
	pluginName = cms.string('elecPairHistManagerByLooseIsolation'),
	pluginType = cms.string('PATDiElecPairHistManager'),
	diTauCandidateSource = cms.InputTag('allDiElecPairZeeHypothesesByLooseIsolation'),
	dqmDirectory_store = cms.string('DiElecZeeHypothesisByLooseIsolationQuantities')
)

# import config for central jet veto histogram manager
from TauAnalysis.Core.jetHistManager_cfi import *

# import config for missing-Et histogram managers
from TauAnalysis.Core.caloMEtHistManager_cfi import *
caloMEtHistManager.leg1Source = cms.InputTag('selectedPatElectronsForElecTauConversionVetoCumulative')
caloMEtHistManager.leg2Source = cms.InputTag('selectedPatTausForElecTauMuonVetoCumulative')
from TauAnalysis.Core.pfMEtHistManager_cfi import *
pfMEtHistManager.leg1Source = cms.InputTag('selectedPatElectronsForElecTauConversionVetoCumulative')
pfMEtHistManager.leg2Source = cms.InputTag('selectedPatTausForElecTauMuonVetoCumulative')

# import config for particle multiplicity histogram manager
from TauAnalysis.Core.particleMultiplicityHistManager_cfi import *

# import config for primary event vertex histogram manager
from TauAnalysis.Core.vertexHistManager_cfi import *
vertexHistManager.vertexSource = cms.InputTag('offlinePrimaryVerticesWithBS')

# import config for L1 & HLT histogram manager
from TauAnalysis.Core.triggerHistManager_cfi import *
triggerHistManagerForElecTau = copy.deepcopy(triggerHistManager)
triggerHistManagerForElecTau.pluginName = cms.string('triggerHistManagerForElecTau')
triggerHistManagerForElecTau.l1Bits = cms.vstring(
    'L1_SingleEG5',
    'L1_SingleEG8',
    'L1_SingleEG10',
    'L1_SingleEG12',
    'L1_SingleEG15',
    'L1_SingleIsoEG5',
    'L1_SingleIsoEG8',
    'L1_SingleIsoEG10',
    'L1_SingleIsoEG12',
    'L1_SingleIsoEG15'
)

#triggerHistManagerForElecTau.hltPaths = cms.vstring(    
#	'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v1',
#	'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2',
#	'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v4',
#    'HLT_IsoEle12_PFTau15_v3',
#    'HLT_Ele12_SW_TighterEleId_L1R_v2'
#)

# import config for binning results
# used for keeping track of number of events passing all selection criteria
from TauAnalysis.Core.dataBinner_cfi import *

# import config for binning results
# used to estimate systematic uncertainties
from TauAnalysis.Core.sysUncertaintyBinner_cfi import *
from TauAnalysis.CandidateTools.sysErrDefinitions_cfi import *
sysUncertaintyNames = [ "CENTRAL_VALUE", ]
sysUncertaintyNames.extend(
    getSysUncertaintyNames(
        [ electronSystematics,
          tauSystematics,
          elecTauPairSystematics,
          jetSystematics,
          theorySystematics ]
    )
)
sysUncertaintyBinnerForElecTauEff = sysUncertaintyBinner.clone(
    pluginName = cms.string('sysUncertaintyBinnerForElecTauEff'),
    binnerPlugins = cms.VPSet(
        dataBinner
    ),
    systematics = cms.vstring(sysUncertaintyNames)
)

sysUncertaintyHistManagerForElecTau = cms.PSet(
    pluginName = cms.string('sysUncertaintyHistManagerForElecTau'),
    pluginType = cms.string('SysUncertaintyHistManager'),
    histManagers = cms.VPSet(
        cms.PSet(
            config = diTauCandidateHistManagerForElecTau,
            systematics = cms.PSet(
                diTauCandidateSource = getSysUncertaintyParameterSets(
                    [ elecTauPairSystematics ]
                )
            )
        ),
        cms.PSet(
            config = diTauCandidateNSVfitHistManagerForElecTau,
            systematics = cms.PSet(
                diTauCandidateSource = getSysUncertaintyParameterSets(
                    [ elecTauPairSystematics ]
                )
            )
        )
    ),
    dqmDirectory_store = cms.string('sysUncertaintyHistManagerResults')
)

diTauLeg1ChargeBinGridHistManager = cms.PSet(
    pluginName = cms.string('diTauLeg1ChargeBinGridHistManager'),
    pluginType = cms.string('BinGridHistManager'),
    binning = cms.PSet(
        name = cms.string("diTauLeg1ChargeBinning"),
        config = cms.VPSet(
            cms.PSet(
                extractor = cms.PSet(
                    pluginType = cms.string("PATElecTauPairValExtractor"),
                    src = cms.InputTag('selectedElecTauPairsForAHtoElecTauZeroChargeCumulative'),
                    value = cms.string('leg1.charge')
                ),
                branchName = cms.string('diTauLeg1Charge'),
                binning = cms.PSet(
                    boundaries = cms.vdouble(0.),
                    min = cms.double(-2.),
                    max = cms.double(+2.)
                )
            )
        )
    ),
    histManagers = cms.VPSet(
        diTauCandidateHistManagerForElecTau,
        diTauCandidateNSVfitHistManagerForElecTau
    ),
    dqmDirectory_store = cms.string('diTauLeg1ChargeBinnedHistograms')
)
#--------------------------------------------------------------------------------
# define event selection criteria
#--------------------------------------------------------------------------------


# di-tau candidate selection
evtSelDiTauCandidateForAHtoElecTauAntiOverlapVeto = evtSelDiTauCandidateForElecTauAntiOverlapVeto.clone(
    pluginName = cms.string('evtSelDiTauCandidateForAHtoElecTauAntiOverlapVeto'),
    src_cumulative = cms.InputTag('diTauCandidateForAHtoElecTauAntiOverlapVeto', 'cumulative'),
    src_individual = cms.InputTag('diTauCandidateForAHtoElecTauAntiOverlapVeto', 'individual')
)
evtSelDiTauCandidateForAHtoElecTauMt1MET = evtSelDiTauCandidateForElecTauMt1MET.clone(
    pluginName = cms.string('evtSelDiTauCandidateForAHtoElecTauMt1MET'),
    src_cumulative = cms.InputTag('diTauCandidateForAHtoElecTauMt1METcut', 'cumulative'),
    src_individual = cms.InputTag('diTauCandidateForAHtoElecTauMt1METcut', 'individual')
)
evtSelDiTauCandidateForAHtoElecTauPzetaDiff = evtSelDiTauCandidateForElecTauPzetaDiff.clone(
    pluginName = cms.string('evtSelDiTauCandidateForAHtoElecTauPzetaDiff'),
    src_cumulative = cms.InputTag('diTauCandidateForAHtoElecTauPzetaDiffCut', 'cumulative'),
    src_individual = cms.InputTag('diTauCandidateForAHtoElecTauPzetaDiffCut', 'individual')
)

# opposite-sign pairs
evtSelDiTauCandidateForAHtoElecTauZeroCharge = evtSelDiTauCandidateForElecTauZeroCharge.clone(
    pluginName = cms.string('evtSelDiTauCandidateForAHtoElecTauZeroCharge'),
    src_cumulative = cms.InputTag('diTauCandidateForAHtoElecTauZeroChargeCut', 'cumulative'),
    src_individual = cms.InputTag('diTauCandidateForAHtoElecTauZeroChargeCut', 'individual')
)

# same-sign pairs
evtSelDiTauCandidateForAHtoElecTauNonZeroCharge = evtSelDiTauCandidateForElecTauNonZeroCharge.clone(
    pluginName = cms.string('evtSelDiTauCandidateForAHtoElecTauNonZeroCharge'),
    src_cumulative = cms.InputTag('diTauCandidateForAHtoElecTauNonZeroChargeCut', 'cumulative'),
    src_individual = cms.InputTag('diTauCandidateForAHtoElecTauNonZeroChargeCut', 'individual')
)

#  jet veto/b-jet candidate selection
evtSelBtagVeto = cms.PSet(
    pluginName = cms.string('evtSelBtagVeto'),
    pluginType = cms.string('BoolEventSelector'),
    src_cumulative = cms.InputTag('jetBtagVeto', 'cumulative'),
    src_individual = cms.InputTag('jetBtagVeto', 'individual'),
    systematics = cms.vstring(jetSystematics.keys())
)
evtSelJetEtCut = cms.PSet(
    pluginName = cms.string('evtSelJetEtCut'),
    pluginType = cms.string('BoolEventSelector'),
    src_cumulative = cms.InputTag('jetEtCut', 'cumulative'),
    src_individual = cms.InputTag('jetEtCut', 'individual'),
    systematics = cms.vstring(jetSystematics.keys())
)
evtSelBtagCut = cms.PSet(
    pluginName = cms.string('evtSelBtagCut'),
    pluginType = cms.string('BoolEventSelector'),
    src_cumulative = cms.InputTag('jetBtagCut', 'cumulative'),
    src_individual = cms.InputTag('jetBtagCut', 'individual'),
    systematics = cms.vstring(jetSystematics.keys())
)

#--------------------------------------------------------------------------------
# define event print-out
#--------------------------------------------------------------------------------
elecTauEventDump_woBtag = elecTauEventDump.clone(
	diTauCandidateSource = cms.InputTag('selectedElecTauPairsForAHtoElecTauZeroChargeCumulative'),
	elecTauZeeHypothesisSource = cms.InputTag('elecTauPairZeeHypothesesForAHtoElecTau'),
    #pfCandidateSource = cms.InputTag(''),
	triggerConditions = cms.vstring("evtSelDiTauCandidateForAHtoElecTauZeroCharge: passed_cumulative")
)

elecTauEventDump_wBtag = elecTauEventDump.clone(
	diTauCandidateSource = cms.InputTag('selectedElecTauPairsForAHtoElecTauZeroChargeCumulative'),
	elecTauZeeHypothesisSource = cms.InputTag('elecTauPairZeeHypothesesForAHtoElecTau'),
    #pfCandidateSource = cms.InputTag(''),
	triggerConditions = cms.vstring("evtSelDiTauCandidateForAHtoElecTauZeroCharge: passed_cumulative")
)

#--------------------------------------------------------------------------------
# define analysis sequence
# (ordered list of event selection criteria and histogram filling)
#--------------------------------------------------------------------------------

elecTauAnalysisSequenceOS_woBtag = cms.VPSet(
    # fill histograms for full event sample
    cms.PSet(
        analyzers = cms.vstring(
            'genPhaseSpaceEventInfoHistManager',
            'electronHistManager',
            'tauHistManager',
            'pfMEtHistManager',
            'vertexHistManager',
            'triggerHistManagerForElecTau'
        )
    ),

    cms.PSet(
        filter = cms.string('genPhaseSpaceCut'),
        title = cms.string('gen. Phase-Space'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'triggerHistManagerForElecTau'
        )
    ),

    # trigger selection
    cms.PSet(
        filter = cms.string('evtSelTrigger'),
        title = cms.string('Trigger'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),

    # data quality control
    cms.PSet(
        filter = cms.string('evtSelDataQuality'),
        title = cms.string('Data Quality'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager'
        )
    ),

    # electron acceptance cuts
	cms.PSet(
        filter = cms.string('evtSelElectronId'),
        title = cms.string('Electron ID'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauIdCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronAntiCrack'),
        title = cms.string('Electron crack-Veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauAntiCrackCutCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronEta'),
        title = cms.string('-2.5 < eta(Electron) < +2.5'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauEtaCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronPt'),
        title = cms.string('Pt(Electron) > 20 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauPtCumulative')
    ),

    # tau acceptance cuts
    cms.PSet(
        filter = cms.string('evtSelTauAntiOverlapWithElectronsVeto'),
        title = cms.string('Tau not overlapping with Elec.'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauPtCumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauAntiOverlapWithElectronsVetoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauEta'),
        title = cms.string('-2.3 < eta(Tau) < +2.3'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauPtCumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauEtaCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauPt'),
        title = cms.string('Pt(Tau) > 20 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'pfMEtHistManager',
            'vertexHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauPtCumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauPtCumulative')
    ),
  
    cms.PSet(
        filter = cms.string('evtSelElectronIso'),
        title = cms.string('Electron PF iso. (rel.)'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauIsoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronConversionVeto'),
        title = cms.string('Electron Track conv. veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauConversionVetoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronTrkIP'),
        title = cms.string('Electron Track IP'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'pfMEtHistManager',
            'vertexHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative')
    ),
    # 
    # selection of tau-jet candidate (id.)
    # produced in hadronic tau decay
    cms.PSet(
        filter = cms.string('evtSelTauDecayModeFinding'),
        title = cms.string('Tau decay mode find.'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
    		'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauDecayModeFindingCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauLeadTrkPt'),
        title = cms.string('Tau decay mode finding'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauLeadTrkPtCumulative'
							  )
    ),
    cms.PSet(
        filter = cms.string('evtSelTauIso'),
        title = cms.string('Tau ID (HPS loose)'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauIsoCumulative'
							  )
    ),
    cms.PSet(
        filter = cms.string('evtSelTauElectronVeto'),
        title = cms.string('Tau e-Veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauElectronVetoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauEcalCrackVeto'),
        title = cms.string('Tau ECAL crack-Veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'tauHistManager'
        ),
        replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauEcalCrackVetoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauMuonVeto'),
        title = cms.string('Tau mu-Veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'pfMEtHistManager'
        ),
        replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative')
    ),
    
    #selection of electron + tau-jet combinations
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauAntiOverlapVeto'),
        title = cms.string('dR(Electron-Tau) > 0.5'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
			'diTauCandidateHistManagerForElecTau',
            'pfMEtHistManager',
            'vertexHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
                              'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauAntiOverlapVetoCumulative')
    ),
    
    # primary event vertex selection
    cms.PSet(
        filter = cms.string('evtSelPrimaryEventVertexForElecTau'),
        title = cms.string('Common Electron+Tau Vertex'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
			'diTauCandidateHistManagerForElecTau',
            'vertexHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
                              'vertexHistManager.vertexSource = selectedPrimaryVertexForElecTau',
                              'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauAntiOverlapVetoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelPrimaryEventVertexQualityForElecTau'),
        title = cms.string('p(chi2Vertex) > 0.01'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
			'diTauCandidateHistManagerForElecTau',
            'vertexHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
                              'vertexHistManager.vertexSource = selectedPrimaryVertexQualityForElecTau',
                              'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauAntiOverlapVetoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelPrimaryEventVertexPositionForElecTau'),
        title = cms.string('-25 < zVertex < +25 cm'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
			'diTauCandidateHistManagerForElecTau',
            'vertexHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
                              'vertexHistManager.vertexSource = selectedPrimaryVertexPositionForElecTau',
                              'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauAntiOverlapVetoCumulative')
    ),
    
    # selection of electron + tau-jet combinations (continued)
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauMt1MET'),
        title = cms.string('M_{T}(Electron-MET) < 4000 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'diTauCandidateHistManagerForElecTau'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
                              'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauMt1METcumulative'
		)  
    ),
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauPzetaDiff'),
        title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} > -20 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'diTauCandidateHistManagerForElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau',
            'elecPairHistManagerByLooseIsolation'
        ),
		replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
			'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
            'diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
			'elecPairHistManagerByLooseIsolation.diTauCandidateSource = allDiElecPairZeeHypothesesByLooseIsolation'
		)
    ),
	# veto events compatible with Z --> e+ e- hypothesis
	# based on presence of a second, oppositely-charged, loosely isolated electron
	cms.PSet(
	    filter = cms.string('evtSelDiElecPairZeeHypothesisVetoByLooseIsolation'),
	    title = cms.string('charge(isoelectron + looseIsoElectron) != 0'),
	    saveRunLumiSectionEventNumbers = cms.vstring('')
	),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'diTauCandidateHistManagerForElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau',
            'diTauCandidateZeeHypothesisHistManagerForElecTau',
            'elecPairHistManagerByLooseIsolation',
            'jetHistManager'
        ),
		replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
			'jetHistManager.jetSource = selectedPatJetsForAHtoElecTauAntiOverlapWithLeptonsVetoCumulative',
			'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
            'diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
			'elecPairHistManagerByLooseIsolation.diTauCandidateSource = selectedDiElecPairZeeHypothesesByLooseIsolation',
			'diTauCandidateZeeHypothesisHistManagerForElecTau.ZllHypothesisSource = selectedElecTauPairZeeHypotheses'
		)
    ),
    # select events with no more than 1 jet with Et > 30
	cms.PSet(
		filter = cms.string('evtSelJetEtCut'),
		title = cms.string('N(jets with E_{T} > 30) < 2'),
		saveRunLumiSectionEventNumbers = cms.vstring('')
	),
	cms.PSet(
		analyzers = cms.vstring(
			'electronHistManager',
			'tauHistManager',
			'diTauCandidateHistManagerForElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau',
			'jetHistManager'
		),
		replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauElectronVetoCumulative',
			'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
			'diTauCandidateHistManagerForElecTau.visMassHypothesisSource = elecTauPairVisMassHypothesesForAHtoElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
			'jetHistManager.jetSource = selectedPatJetsForAHtoElecTauJetTagCumulative'
		)
	),
    # veto events containing b-tagged jets
    # in order to avoid overlap with "b-tag" analysis path
	cms.PSet(
		filter = cms.string('evtSelBtagVeto'),
		title = cms.string('no E_{T} > 20 GeV jet with b-Tag'),
		saveRunLumiSectionEventNumbers = cms.vstring('')
	),
	cms.PSet(
		analyzers = cms.vstring(
			'electronHistManager',
			'tauHistManager',
			'diTauCandidateHistManagerForElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau',
		),
		replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauElectronVetoCumulative',
			'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
			'diTauCandidateHistManagerForElecTau.visMassHypothesisSource = elecTauPairVisMassHypothesesForAHtoElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative'
		)
	),
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauZeroCharge'),
        title = cms.string('Charge(Electron+Tau) = 0'),
        saveRunLumiSectionEventNumbers = cms.vstring('passed_cumulative')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'genPhaseSpaceEventInfoHistManager',
            'electronHistManager',
            'tauHistManager',
            'diTauCandidateHistManagerForElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau',
            'diTauCandidateZeeHypothesisHistManagerForElecTau',
            #'diTauLeg1ChargeBinGridHistManager', breaks 4_2_X
			'jetHistManager',
            'caloMEtHistManager',
            'pfMEtHistManager',
            'particleMultiplicityHistManager',
            'vertexHistManager',
            'triggerHistManagerForElecTau',
			'sysUncertaintyHistManagerForElecTau',
			'dataBinner',
			'sysUncertaintyBinnerForElecTauEff'
        ),
		replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
			'jetHistManager.jetSource = selectedPatJetsForAHtoElecTauJetTagCumulative',
            'vertexHistManager.vertexSource = selectedPrimaryVertexHighestPtTrackSumForElecTau',
            'diTauCandidateZeeHypothesisHistManagerForElecTau.ZllHypothesisSource = elecTauPairZeeHypothesesForAHtoElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauZeroChargeCumulative',
			'diTauCandidateHistManagerForElecTau.visMassHypothesisSource = elecTauPairVisMassHypothesesForAHtoElecTau',
			'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauZeroChargeCumulative'
		)
	)
)
elecTauAnalysisSequenceOS_wBtag = cms.VPSet(
    # fill histograms for full event sample
    cms.PSet(
        analyzers = cms.vstring(
            'genPhaseSpaceEventInfoHistManager',
            'electronHistManager',
            'tauHistManager',
            'pfMEtHistManager',
            'vertexHistManager',
            'triggerHistManagerForElecTau'
        )
    ),

    cms.PSet(
        filter = cms.string('genPhaseSpaceCut'),
        title = cms.string('gen. Phase-Space'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'triggerHistManagerForElecTau'
        )
    ),

    # trigger selection
    cms.PSet(
        filter = cms.string('evtSelTrigger'),
        title = cms.string('Trigger'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'vertexHistManager'
        )
    ),
  
    # data quality control
    cms.PSet(
        filter = cms.string('evtSelDataQuality'),
        title = cms.string('Data Quality'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager'
        )
    ),

    # electron acceptance cuts
	cms.PSet(
        filter = cms.string('evtSelElectronId'),
        title = cms.string('Electron ID'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauIdCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronAntiCrack'),
        title = cms.string('Electron crack-Veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauAntiCrackCutCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronEta'),
        title = cms.string('-2.5 < eta(Electron) < +2.5'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauEtaCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronPt'),
        title = cms.string('Pt(Electron) > 20 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'triggerHistManagerForElecTau'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauPtCumulative')
    ),

    # tau acceptance cuts
    cms.PSet(
        filter = cms.string('evtSelTauAntiOverlapWithElectronsVeto'),
        title = cms.string('Tau not overlapping with Elec.'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauPtCumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauAntiOverlapWithElectronsVetoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauEta'),
        title = cms.string('-2.3 < eta(Tau) < +2.3'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauPtCumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauEtaCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauPt'),
        title = cms.string('Pt(Tau) > 20 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'pfMEtHistManager',
            'vertexHistManager',
            'triggerHistManagerForElecTau'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauPtCumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauPtCumulative')
    ),
  
    cms.PSet(
        filter = cms.string('evtSelElectronIso'),
        title = cms.string('Electron PF iso. (rel.)'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauIsoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronConversionVeto'),
        title = cms.string('Electron Track conv. veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'vertexHistManager',
            'electronHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauConversionVetoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronTrkIP'),
        title = cms.string('Electron Track IP'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'pfMEtHistManager',
            'vertexHistManager',
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative')
    ),
  
    # selection of tau-jet candidate (id.)
    # produced in hadronic tau decay
    cms.PSet(
        filter = cms.string('evtSelTauDecayModeFinding'),
        title = cms.string('Tau decay mode find.'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
    		'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauDecayModeFindingCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauLeadTrkPt'),
        title = cms.string('Tau decay mode finding'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauLeadTrkPtCumulative'
							  )
    ),
    cms.PSet(
        filter = cms.string('evtSelTauIso'),
        title = cms.string('Tau ID (HPS loose)'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauIsoCumulative'
							  )
    ),
    cms.PSet(
        filter = cms.string('evtSelTauElectronVeto'),
        title = cms.string('Tau e-Veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauElectronVetoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauEcalCrackVeto'),
        title = cms.string('Tau ECAL crack-Veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'tauHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauEcalCrackVetoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauMuonVeto'),
        title = cms.string('Tau mu-Veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager'
        ),
        replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative')
    ),
    
    #selection of electron + tau-jet combinations
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauAntiOverlapVeto'),
        title = cms.string('dR(Electron-Tau) > 0.5'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
			'diTauCandidateHistManagerForElecTau',
            'caloMEtHistManager',
            'pfMEtHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
                              'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauAntiOverlapVetoCumulative')
    ),

    # primary event vertex selection
    cms.PSet(
        filter = cms.string('evtSelPrimaryEventVertexForElecTau'),
        title = cms.string('Common Electron+Tau Vertex'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
			'diTauCandidateHistManagerForElecTau',
            'vertexHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
                              'vertexHistManager.vertexSource = selectedPrimaryVertexForElecTau',
                              'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauAntiOverlapVetoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelPrimaryEventVertexQualityForElecTau'),
        title = cms.string('p(chi2Vertex) > 0.01'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
			'diTauCandidateHistManagerForElecTau',
            'vertexHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
                              'vertexHistManager.vertexSource = selectedPrimaryVertexQualityForElecTau',
                              'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauAntiOverlapVetoCumulative')
    ),
    cms.PSet(
        filter = cms.string('evtSelPrimaryEventVertexPositionForElecTau'),
        title = cms.string('-25 < zVertex < +25 cm'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
			'diTauCandidateHistManagerForElecTau',
            'vertexHistManager'
        ),
        replace = cms.vstring('electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
                              'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
                              'vertexHistManager.vertexSource = selectedPrimaryVertexPositionForElecTau',
                              'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauAntiOverlapVetoCumulative')
    ),
    
    #selection of electron + tau-jet combinations (continued)
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauMt1MET'),
        title = cms.string('M_{T}(Electron-MET) < 4000 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'diTauCandidateHistManagerForElecTau'
        ),
        replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
            'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
            'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauMt1METcumulative'
		)  
    ),
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauPzetaDiff'),
        title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} > -20 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'diTauCandidateHistManagerForElecTau',
            'elecPairHistManagerByLooseIsolation'
        ),
		replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
			'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
            'diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
			'elecPairHistManagerByLooseIsolation.diTauCandidateSource = allDiElecPairZeeHypothesesByLooseIsolation'
		)
    ),
	# veto events compatible with Z --> e+ e- hypothesis
	# based on presence of a second, oppositely-charged, loosely isolated electron
	cms.PSet(
	    filter = cms.string('evtSelDiElecPairZeeHypothesisVetoByLooseIsolation'),
	    title = cms.string('charge(isoelectron + looseIsoElectron) != 0'),
	    saveRunLumiSectionEventNumbers = cms.vstring('')
	),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'diTauCandidateHistManagerForElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau',
            'diTauCandidateZeeHypothesisHistManagerForElecTau',
            'jetHistManager'
        ),
		replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
			'jetHistManager.jetSource = selectedPatJetsForAHtoElecTauAntiOverlapWithLeptonsVetoCumulative',
			'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
            'diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative'
		)
    ),
    # require at least one b-tagged jet with Et > 20 in the event
	cms.PSet(
		filter = cms.string('evtSelJetEtCut'),
		title = cms.string('N(jets with E_{T} > 30) < 2'),
		saveRunLumiSectionEventNumbers = cms.vstring('')
	),
	cms.PSet(
		analyzers = cms.vstring(
			'electronHistManager',
			'tauHistManager',
			'diTauCandidateHistManagerForElecTau',
			'jetHistManager'
			),
		replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauElectronVetoCumulative',
			'jetHistManager.jetSource = selectedPatJetsForAHtoElecTauJetTagCumulative',
			'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
			'diTauCandidateHistManagerForElecTau.visMassHypothesisSource = elecTauPairVisMassHypothesesForAHtoElecTau'
		)
	),
	cms.PSet(
		filter = cms.string('evtSelBtagCut'),
		title = cms.string('E_{T} > 20 GeV jet with b-Tag'),
		saveRunLumiSectionEventNumbers = cms.vstring('')
	),
	cms.PSet(
		analyzers = cms.vstring(
			'electronHistManager',
			'tauHistManager',
			'diTauCandidateHistManagerForElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau',
			'jetHistManager'
			),
		replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauElectronVetoCumulative',
			'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
			'diTauCandidateHistManagerForElecTau.visMassHypothesisSource = elecTauPairVisMassHypothesesForAHtoElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauPzetaDiffCumulative',
			'jetHistManager.jetSource = selectedPatJetsForAHtoElecTauBtagCumulative'
		)
	),
	#  finally apply opposite-sign selection
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauZeroCharge'),
        title = cms.string('Charge(Electron+Tau) = 0'),
        saveRunLumiSectionEventNumbers = cms.vstring('passed_cumulative')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'genPhaseSpaceEventInfoHistManager',
            'electronHistManager',
            'tauHistManager',
            'diTauCandidateHistManagerForElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau',
            'diTauCandidateZeeHypothesisHistManagerForElecTau',
            #'diTauLeg1ChargeBinGridHistManager',  breaks 4_2_X
			'jetHistManager',
            'caloMEtHistManager',
            'pfMEtHistManager',
            'particleMultiplicityHistManager',
            'vertexHistManager',
            'triggerHistManagerForElecTau',
			'sysUncertaintyHistManagerForElecTau',
			'dataBinner',
			'sysUncertaintyBinnerForElecTauEff'
        ),
		replace = cms.vstring(
			'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
			'tauHistManager.tauSource = selectedPatTausForElecTauMuonVetoCumulative',
			'jetHistManager.jetSource = selectedPatJetsForAHtoElecTauBtagCumulative',
            'vertexHistManager.vertexSource = selectedPrimaryVertexHighestPtTrackSumForElecTau',
            'diTauCandidateZeeHypothesisHistManagerForElecTau.ZllHypothesisSource = elecTauPairZeeHypothesesForAHtoElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauZeroChargeCumulative',
			'diTauCandidateHistManagerForElecTau.visMassHypothesisSource = elecTauPairVisMassHypothesesForAHtoElecTau',
			'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauZeroChargeCumulative'
		)
	)
)
elecTauAnalysisSequenceSS_woBtag = cms.VPSet(
    cms.PSet(
        filter = cms.string('genPhaseSpaceCut'),
        title = cms.string('gen. Phase-Space'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTrigger'),
        title = cms.string('Trigger'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelDataQuality'),
        title = cms.string('Data Quality'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronId'),
        title = cms.string('Electron ID'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronAntiCrack'),
        title = cms.string('Electron crack veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronEta'),
        title = cms.string('-2.5 < eta(Electron) < +2.5'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronPt'),
        title = cms.string('Pt(Electron) > 20 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauAntiOverlapWithElectronsVeto'),
        title = cms.string('Tau not overlapping w. Electron'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauEta'),
        title = cms.string('-2.3 < eta(Tau) < +2.3'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauPt'),
        title = cms.string('Pt(Tau) > 20 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronIso'),
        title = cms.string('Electron iso.'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronConversionVeto'),
        title = cms.string('Electron conversion veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronTrkIP'),
        title = cms.string('Electron Track IP'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauDecayModeFinding'),
        title = cms.string('Tau decay mode find.'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauLeadTrkPt'),
        title = cms.string('Tau decay mode finding'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauIso'),
        title = cms.string('Tau ID (HPS loose)'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauElectronVeto'),
        title = cms.string('Tau e-Veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
	cms.PSet(
		filter = cms.string('evtSelTauEcalCrackVeto'),
		title = cms.string('Tau ECAL crack-Veto'),
		saveRunLumiSectionEventNumbers = cms.vstring('')
	),
    cms.PSet(
        filter = cms.string('evtSelTauMuonVeto'),
        title = cms.string('Tau mu-Veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauAntiOverlapVeto'),
        title = cms.string('dR(Electron-Tau) > 0.5'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelPrimaryEventVertexForElecTau'),
        title = cms.string('Common Electron+Tau Vertex'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelPrimaryEventVertexQualityForElecTau'),
        title = cms.string('Vertex quality'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelPrimaryEventVertexPositionForElecTau'),
        title = cms.string('-24 < zVertex < +24 cm'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauMt1MET'),
        title = cms.string('M_{T}(Electron-MET) < 4000 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauPzetaDiff'),
        title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} > -20 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelDiElecPairZeeHypothesisVetoByLooseIsolation'),
        title = cms.string('charge(isoElectron+looseIsoElectron) != 0'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
	cms.PSet(
		filter = cms.string('evtSelJetEtCut'),
		title = cms.string('N(jets with E_{T} > 30) < 2'),
		saveRunLumiSectionEventNumbers = cms.vstring('')
	),
    cms.PSet(
		filter = cms.string('evtSelBtagVeto'),
		title = cms.string('no E_{T} > 20 GeV jet with b-Tag'),
		saveRunLumiSectionEventNumbers = cms.vstring('')
	),
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauNonZeroCharge'),
        title = cms.string('Charge(Muon+Tau) != 0'),
        saveRunLumiSectionEventNumbers = cms.vstring('passed_cumulative')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'diTauCandidateHistManagerForElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau',
            'diTauCandidateZeeHypothesisHistManagerForElecTau',
            'jetHistManager',
            'caloMEtHistManager',
            'pfMEtHistManager',
            'vertexHistManager',
            'triggerHistManagerForElecTau'
        ),
        replace = cms.vstring(
            'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
            'tauHistManager.tauSource = selectedPatTausForElecTauElectronVetoCumulative',
            'vertexHistManager.vertexSource = selectedPrimaryVertexHighestPtTrackSumForElecTau',
            'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauNonZeroChargeCumulative',
            'diTauCandidateHistManagerForElecTau.visMassHypothesisSource = elecTauPairVisMassHypothesesForAHtoElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauNonZeroChargeCumulative',
            'diTauCandidateZeeHypothesisHistManagerForElecTau.ZllHypothesisSource = elecTauPairZeeHypothesesForAHtoElecTau',
			'jetHistManager.jetSource = selectedPatJetsForAHtoElecTauAntiOverlapWithLeptonsVetoCumulative'
        )
    )
)

elecTauAnalysisSequenceSS_wBtag = cms.VPSet(
    cms.PSet(
        filter = cms.string('genPhaseSpaceCut'),
        title = cms.string('gen. Phase-Space'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTrigger'),
        title = cms.string('Trigger'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelDataQuality'),
        title = cms.string('Data Quality'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronId'),
        title = cms.string('Electron ID'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronAntiCrack'),
        title = cms.string('Electron crack veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronEta'),
        title = cms.string('-2.5 < eta(Electron) < +2.5'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronPt'),
        title = cms.string('Pt(Electron) > 20 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauAntiOverlapWithElectronsVeto'),
        title = cms.string('Tau not overlapping w. Electron'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauEta'),
        title = cms.string('-2.3 < eta(Tau) < +2.3'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauPt'),
        title = cms.string('Pt(Tau) > 20 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronIso'),
        title = cms.string('Electron iso.'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronConversionVeto'),
        title = cms.string('Electron conversion veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelElectronTrkIP'),
        title = cms.string('Electron Track IP'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauDecayModeFinding'),
        title = cms.string('Tau lead. Track find.'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauLeadTrkPt'),
        title = cms.string('Tau decay mode finding'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauIso'),
        title = cms.string('Tau ID (HPS loose)'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelTauElectronVeto'),
        title = cms.string('Tau e-Veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
	cms.PSet(
		filter = cms.string('evtSelTauEcalCrackVeto'),
		title = cms.string('Tau ECAL crack-Veto'),
		saveRunLumiSectionEventNumbers = cms.vstring('')
	),
    cms.PSet(
        filter = cms.string('evtSelTauMuonVeto'),
        title = cms.string('Tau mu-Veto'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauAntiOverlapVeto'),
        title = cms.string('dR(Electron-Tau) > 0.5'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelPrimaryEventVertexForElecTau'),
        title = cms.string('Common Electron+Tau Vertex'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelPrimaryEventVertexQualityForElecTau'),
        title = cms.string('Vertex quality'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelPrimaryEventVertexPositionForElecTau'),
        title = cms.string('-24 < zVertex < +24 cm'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauMt1MET'),
        title = cms.string('M_{T}(Electron-MET) < 4000 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauPzetaDiff'),
        title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} > -20 GeV'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
    cms.PSet(
        filter = cms.string('evtSelDiElecPairZeeHypothesisVetoByLooseIsolation'),
        title = cms.string('charge(isoElectron+looseIsoElectron) != 0'),
        saveRunLumiSectionEventNumbers = cms.vstring('')
    ),
	cms.PSet(
		filter = cms.string('evtSelJetEtCut'),
		title = cms.string('N(jets with E_{T} > 30) < 2'),
		saveRunLumiSectionEventNumbers = cms.vstring('')
	),
	cms.PSet(
		filter = cms.string('evtSelBtagCut'),
		title = cms.string('E_{T} > 20 GeV jet with b-Tag'),
		saveRunLumiSectionEventNumbers = cms.vstring('')
	),
    cms.PSet(
        filter = cms.string('evtSelDiTauCandidateForAHtoElecTauNonZeroCharge'),
        title = cms.string('Charge(Muon+Tau) != 0'),
        saveRunLumiSectionEventNumbers = cms.vstring('passed_cumulative')
    ),
    cms.PSet(
        analyzers = cms.vstring(
            'electronHistManager',
            'tauHistManager',
            'diTauCandidateHistManagerForElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau',
            'diTauCandidateZeeHypothesisHistManagerForElecTau',
            'jetHistManager',
            'caloMEtHistManager',
            'pfMEtHistManager',
            'vertexHistManager',
            'triggerHistManagerForElecTau'
        ),
        replace = cms.vstring(
            'electronHistManager.electronSource = selectedPatElectronsForElecTauTrkIPcumulative',
            'tauHistManager.tauSource = selectedPatTausForElecTauElectronVetoCumulative',
            'vertexHistManager.vertexSource = selectedPrimaryVertexHighestPtTrackSumForElecTau',
            'diTauCandidateHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauNonZeroChargeCumulative',
            'diTauCandidateHistManagerForElecTau.visMassHypothesisSource = elecTauPairVisMassHypothesesForAHtoElecTau',
            'diTauCandidateNSVfitHistManagerForElecTau.diTauCandidateSource = selectedElecTauPairsForAHtoElecTauNonZeroChargeCumulative',
            'diTauCandidateZeeHypothesisHistManagerForElecTau.ZllHypothesisSource = elecTauPairZeeHypothesesForAHtoElecTau',
            'jetHistManager.jetSource = selectedPatJetsForAHtoElecTauBtagCumulative'
        )
    )
)
