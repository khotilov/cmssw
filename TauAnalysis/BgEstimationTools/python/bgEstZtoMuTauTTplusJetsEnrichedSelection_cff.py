import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *

#--------------------------------------------------------------------------------
# select ttbar + jets background enriched event sample
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce collection of pat::Muons
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patMuonSelection_cfi import *

muonsBgEstTTplusJetsEnrichedVbTfId = copy.deepcopy(selectedPatMuonsVbTfId)

muonsBgEstTTplusJetsEnrichedPFRelIso = copy.deepcopy(selectedPatMuonsPFRelIso)
muonsBgEstTTplusJetsEnrichedPFRelIso.sumPtMax = cms.double(0.10)

muonSelConfiguratorBgEstTTplusJetsEnriched = objSelConfigurator(
    [ muonsBgEstTTplusJetsEnrichedVbTfId,
      muonsBgEstTTplusJetsEnrichedPFRelIso ],
    src = "selectedPatMuonsPt15Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuonsBgEstTTplusJetsEnriched = muonSelConfiguratorBgEstTTplusJetsEnriched.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of vertices compatible with background-enriched lepton collections
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.recoVertexSelectionForMuTau_cff import *

selectedPrimaryVertexForMuTauBgEstTTplusJetsEnriched = selectedPrimaryVertexForMuTau.clone(
    srcParticles = cms.VInputTag(
        'muonsBgEstTTplusJetsEnrichedPFRelIsoCumulative',
        'selectedPatTausForMuTauMuonVetoCumulative'
    )
)
selectedPrimaryVertexQualityForMuTauBgEstTTplusJetsEnriched = selectedPrimaryVertexQualityForMuTau.clone(
    src = cms.InputTag('selectedPrimaryVertexForMuTauBgEstTTplusJetsEnriched')
)
selectedPrimaryVertexPositionForMuTauBgEstTTplusJetsEnriched = selectedPrimaryVertexPositionForMuTau.clone(
    src = cms.InputTag('selectedPrimaryVertexQualityForMuTauBgEstTTplusJetsEnriched')
)
selectedPrimaryVertexHighestPtTrackSumForMuTauBgEstTTplusJetsEnriched = selectedPrimaryVertexHighestPtTrackSumForMuTau.clone(
    src = cms.InputTag('selectedPrimaryVertexPositionForMuTauBgEstTTplusJetsEnriched')
)
selectPrimaryVertexForMuTauBgEstTTplusJetsEnriched = cms.Sequence(
    selectedPrimaryVertexForMuTauBgEstTTplusJetsEnriched
    * selectedPrimaryVertexQualityForMuTauBgEstTTplusJetsEnriched
    * selectedPrimaryVertexPositionForMuTauBgEstTTplusJetsEnriched
    * selectedPrimaryVertexHighestPtTrackSumForMuTauBgEstTTplusJetsEnriched
)

#--------------------------------------------------------------------------------
# produce collection of muon + tau-jet combinations
#--------------------------------------------------------------------------------

from TauAnalysis.CandidateTools.muTauPairProduction_cff import *

muTauPairsBgEstTTplusJetsEnriched = allMuTauPairs.clone(
    srcLeg1 = cms.InputTag('muonsBgEstTTplusJetsEnrichedPFRelIsoCumulative'),
    srcLeg2 = cms.InputTag('selectedPatTausForMuTauMuonVetoCumulative'),
    verbosity = cms.untracked.int32(0)
)

produceMuTauPairsBgEstTTplusJetsEnriched = cms.Sequence(muTauPairsBgEstTTplusJetsEnriched)

muTauPairsBgEstTTplusJetsEnrichedZeroCharge = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsBgEstTTplusJetsEnriched'),
    cut = cms.string('abs(charge) < 0.5'),
    filter = cms.bool(False)
)

selectMuTauPairsBgEstTTplusJetsEnriched = cms.Sequence(muTauPairsBgEstTTplusJetsEnrichedZeroCharge)

#--------------------------------------------------------------------------------
# produce collection of pat::Jets
#--------------------------------------------------------------------------------

jetsBgEstTTplusJetsEnrichedAntiOverlapWithLeptonsVeto = cms.EDFilter("PATJetAntiOverlapSelector",
    src = cms.InputTag("patJets"),
    srcNotToBeFiltered = cms.VInputTag(
        "selectedPatElectronsTrkCumulative",
        "muonsBgEstTTplusJetsEnrichedPFRelIsoCumulative",
        "selectedPatTausForMuTauMuonVetoCumulative"
    ),
    dRmin = cms.double(0.7),
    filter = cms.bool(False)
)

jetsBgEstTTplusJetsEnrichedEta24 = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("jetsBgEstTTplusJetsEnrichedAntiOverlapWithLeptonsVeto"),
    cut = cms.string('abs(eta) < 2.4'),
    filter = cms.bool(False)
)

jetsBgEstTTplusJetsEnrichedEt40 = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("jetsBgEstTTplusJetsEnrichedEta24"),
    cut = cms.string('et > 40.'),
    filter = cms.bool(False)
)

jetsBgEstTTplusJetsEnrichedEt40bTag = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("jetsBgEstTTplusJetsEnrichedEt40"),
    cut = cms.string('bDiscriminator("trackCountingHighEffBJetTags") > 2.5'),
    filter = cms.bool(False)
)

jetsBgEstTTplusJetsEnrichedEt60 = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("jetsBgEstTTplusJetsEnrichedAntiOverlapWithLeptonsVeto"),
    cut = cms.string('et > 60.'),
    filter = cms.bool(False)
)

selectJetsBgEstTTplusJetsEnriched = cms.Sequence(
    jetsBgEstTTplusJetsEnrichedAntiOverlapWithLeptonsVeto
   + jetsBgEstTTplusJetsEnrichedEta24
   + jetsBgEstTTplusJetsEnrichedEt40 + jetsBgEstTTplusJetsEnrichedEt40bTag
   + jetsBgEstTTplusJetsEnrichedEt60
)

#--------------------------------------------------------------------------------
# produce boolean event selection flags
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.selectZtoMuTau_cff import *

cfgMuonVbTfIdCutBgEstTTplusJetsEnriched = copy.deepcopy(cfgMuonVbTfIdCut)
cfgMuonVbTfIdCutBgEstTTplusJetsEnriched.pluginName = cms.string('muonVbTfIdCutBgEstTTplusJetsEnriched')
cfgMuonVbTfIdCutBgEstTTplusJetsEnriched.src_cumulative = cms.InputTag('muonsBgEstTTplusJetsEnrichedVbTfIdCumulative')
cfgMuonVbTfIdCutBgEstTTplusJetsEnriched.systematics = cms.vstring()

cfgMuonPFRelIsoCutBgEstTTplusJetsEnriched = copy.deepcopy(cfgMuonPFRelIsoCut)
cfgMuonPFRelIsoCutBgEstTTplusJetsEnriched.pluginName = cms.string('muonPFRelIsoCutBgEstTTplusJetsEnriched')
cfgMuonPFRelIsoCutBgEstTTplusJetsEnriched.src_cumulative = cms.InputTag('muonsBgEstTTplusJetsEnrichedPFRelIsoCumulative')
cfgMuonPFRelIsoCutBgEstTTplusJetsEnriched.systematics = cms.vstring()

cfgMuTauPairBgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('muTauPairBgEstTTplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsBgEstTTplusJetsEnriched'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairZeroChargeBgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('muTauPairZeroChargeBgEstTTplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsBgEstTTplusJetsEnrichedZeroCharge'),
    minNumber = cms.uint32(1)
)

cfgJetsEt40BgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('jetsEt40BgEstTTplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('jetsBgEstTTplusJetsEnrichedEt40'),
    minNumber = cms.uint32(2)
)

cfgJetEt40bTagBgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('jetEt40bTagBgEstTTplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('jetsBgEstTTplusJetsEnrichedEt40bTag'),
    minNumber = cms.uint32(1)
)

cfgJetEt60BgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('jetEt60BgEstTTplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('jetsBgEstTTplusJetsEnrichedEt60'),
    minNumber = cms.uint32(1)
)
cfgVertexForMuTauBgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('vertexBgEstTTplusJetsEnriched'),
    pluginType = cms.string('VertexMinEventSelector'),
    src = cms.InputTag('selectedPrimaryVertexForMuTauBgEstTTplusJetsEnriched'),
    minNumber = cms.uint32(1)
)
cfgVertexQualityForMuTauBgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('vertexQualityBgEstTTplusJetsEnriched'),
    pluginType = cms.string('VertexMinEventSelector'),
    src = cms.InputTag('selectedPrimaryVertexQualityForMuTauBgEstTTplusJetsEnriched'),
    minNumber = cms.uint32(1)
)
cfgVertexPositionForMuTauBgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('vertexPositionBgEstTTplusJetsEnriched'),
    pluginType = cms.string('VertexMinEventSelector'),
    src = cms.InputTag('selectedPrimaryVertexPositionForMuTauBgEstTTplusJetsEnriched'),
    minNumber = cms.uint32(1)
)

evtSelConfiguratorBgEstTTplusJetsEnriched = eventSelFlagProdConfigurator(
    [ cfgMuonVbTfIdCutBgEstTTplusJetsEnriched,
      cfgMuonPFRelIsoCutBgEstTTplusJetsEnriched,
      cfgVertexForMuTauBgEstTTplusJetsEnriched,
      cfgVertexQualityForMuTauBgEstTTplusJetsEnriched,
      cfgVertexPositionForMuTauBgEstTTplusJetsEnriched,
      cfgMuTauPairBgEstTTplusJetsEnriched,
      cfgMuTauPairZeroChargeBgEstTTplusJetsEnriched,
      cfgJetsEt40BgEstTTplusJetsEnriched,
      cfgJetEt40bTagBgEstTTplusJetsEnriched,
      cfgJetEt60BgEstTTplusJetsEnriched ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

selectEventsBgEstTTplusJetsEnriched = evtSelConfiguratorBgEstTTplusJetsEnriched.configure()

#--------------------------------------------------------------------------------
# apply event selection criteria; fill histograms
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoMuTau_cfi import *
from TauAnalysis.BgEstimationTools.selectZtoMuTauEventVertex_cff import *

muonHistManagerBgEstTTplusJetsEnriched = copy.deepcopy(muonHistManager)
muonHistManagerBgEstTTplusJetsEnriched.pluginName = cms.string('muonHistManagerBgEstTTplusJetsEnriched')
muonHistManagerBgEstTTplusJetsEnriched.muonSource = muTauPairsBgEstTTplusJetsEnriched.srcLeg1

tauHistManagerBgEstTTplusJetsEnriched = copy.deepcopy(tauHistManager)
tauHistManagerBgEstTTplusJetsEnriched.pluginName = cms.string('tauHistManagerBgEstTTplusJetsEnriched')
tauHistManagerBgEstTTplusJetsEnriched.tauSource = muTauPairsBgEstTTplusJetsEnriched.srcLeg2

caloMEtHistManagerBgEstTTplusJetsEnriched = copy.deepcopy(caloMEtHistManager)
caloMEtHistManagerBgEstTTplusJetsEnriched.pluginName = cms.string('caloMEtHistManagerBgEstTTplusJetsEnriched')
caloMEtHistManagerBgEstTTplusJetsEnriched.leg1Source = muonHistManagerBgEstTTplusJetsEnriched.muonSource
caloMEtHistManagerBgEstTTplusJetsEnriched.leg2Source = tauHistManagerBgEstTTplusJetsEnriched.tauSource
pfMEtHistManagerBgEstTTplusJetsEnriched = copy.deepcopy(pfMEtHistManager)
pfMEtHistManagerBgEstTTplusJetsEnriched.pluginName = cms.string('pfMEtHistManagerBgEstTTplusJetsEnriched')
pfMEtHistManagerBgEstTTplusJetsEnriched.leg1Source = muonHistManagerBgEstTTplusJetsEnriched.muonSource
pfMEtHistManagerBgEstTTplusJetsEnriched.leg2Source = tauHistManagerBgEstTTplusJetsEnriched.tauSource

diTauCandidateHistManagerBgEstTTplusJetsEnriched = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerBgEstTTplusJetsEnriched.pluginName = cms.string('diTauCandidateHistManagerBgEstTTplusJetsEnriched')
diTauCandidateHistManagerBgEstTTplusJetsEnriched.diTauCandidateSource = cms.InputTag('muTauPairsBgEstTTplusJetsEnrichedZeroCharge')
diTauCandidateHistManagerBgEstTTplusJetsEnriched.visMassHypothesisSource = cms.InputTag('')

#diTauCandidateNSVfitHistManagerBgEstTTplusJetsEnriched = copy.deepcopy(diTauCandidateNSVfitHistManagerForMuTau)
#diTauCandidateNSVfitHistManagerBgEstTTplusJetsEnriched.pluginName = cms.string('diTauCandidateNSVfitHistManagerBgEstTTplusJetsEnriched')
#diTauCandidateNSVfitHistManagerBgEstTTplusJetsEnriched.diTauCandidateSource = diTauCandidateHistManagerBgEstTTplusJetsEnriched.diTauCandidateSource

jetHistManagerBgEstTTplusJetsEnriched = copy.deepcopy(jetHistManager)
jetHistManagerBgEstTTplusJetsEnriched.pluginName = cms.string('jetHistManagerBgEstTTplusJetsEnriched')
jetHistManagerBgEstTTplusJetsEnriched.jetSource = cms.InputTag('jetsBgEstTTplusJetsEnrichedEt40')

from TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauHistManager_cfi import *
tauIdEffHistManagerBgEstTTplusJetsEnriched = copy.deepcopy(tauIdEffZtoMuTauHistManager)
tauIdEffHistManagerBgEstTTplusJetsEnriched.pluginName = cms.string('tauIdEffHistManagerBgEstTTplusJetsEnriched')
tauIdEffHistManagerBgEstTTplusJetsEnriched.muonSource = muonHistManagerBgEstTTplusJetsEnriched.muonSource
tauIdEffHistManagerBgEstTTplusJetsEnriched.tauSource = tauHistManagerBgEstTTplusJetsEnriched.tauSource
tauIdEffHistManagerBgEstTTplusJetsEnriched.diTauSource = diTauCandidateHistManagerBgEstTTplusJetsEnriched.diTauCandidateSource
tauIdEffHistManagerBgEstTTplusJetsEnriched.diTauChargeSignExtractor.src = tauIdEffHistManagerBgEstTTplusJetsEnriched.diTauSource

dataBinnerBgEstTTplusJetsEnriched = copy.deepcopy(dataBinner)
dataBinnerBgEstTTplusJetsEnriched.pluginName = cms.string('dataBinnerBgEstTTplusJetsEnriched')

analyzeEventsBgEstTTplusJetsEnriched = cms.EDAnalyzer("GenericAnalyzer",

    name = cms.string('BgEstTemplateAnalyzer_TTplusJetsEnriched'),

    filters = cms.VPSet(
        evtSelGenPhaseSpace,
        evtSelTrigger,
        evtSelDataQuality,
        evtSelGlobalMuon,
        evtSelMuonEta,
        evtSelMuonPt,
        cms.PSet(
            pluginName = cms.string('muonVbTfIdCutBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muonVbTfIdCutBgEstTTplusJetsEnriched', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('muonPFRelIsoCutBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muonPFRelIsoCutBgEstTTplusJetsEnriched', 'cumulative')
        ),
        evtSelTauAntiOverlapWithMuonsVeto,
        evtSelTauEta,
        evtSelTauPt,
        evtSelTauLeadTrk,
        evtSelTauLeadTrkPt,
        evtSelTauTaNCdiscr,
        evtSelTauProng,
        evtSelTauCharge,
        evtSelTauMuonVeto,
        cms.PSet(
            pluginName = cms.string('vertexBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('vertexBgEstTTplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('vertexQualityBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('vertexQualityBgEstTTplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('vertexPositionBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('vertexPositionBgEstTTplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairBgEstTTplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairZeroChargeBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairZeroChargeBgEstTTplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('jetEt40BgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('jetsEt40BgEstTTplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('jetEt40bTagBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('jetEt40bTagBgEstTTplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('jetEt60BgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('jetEt60BgEstTTplusJetsEnriched')
        )
    ),

    analyzers = cms.VPSet(
        muonHistManagerBgEstTTplusJetsEnriched,
        tauHistManagerBgEstTTplusJetsEnriched,
        diTauCandidateHistManagerBgEstTTplusJetsEnriched,
		#diTauCandidateNSVfitHistManagerBgEstTTplusJetsEnriched,
        caloMEtHistManagerBgEstTTplusJetsEnriched,
	pfMEtHistManagerBgEstTTplusJetsEnriched,
        jetHistManagerBgEstTTplusJetsEnriched,
        #tauIdEffHistManagerBgEstTTplusJetsEnriched,
        dataBinnerBgEstTTplusJetsEnriched
    ),

    eventDumps = cms.VPSet(),

    analysisSequence = cms.VPSet(

        # generator level phase-space selection
        # (NOTE: (1) to be used in case of Monte Carlo samples
        #            overlapping in simulated phase-space only !!
        #        (2) genPhaseSpaceCut needs to be **always** the first entry in the list of cuts
        #           - otherwise the script submitToBatch.csh for submission of cmsRun jobs
        #            to the CERN batch system will not work !!)
        cms.PSet(
            filter = cms.string('genPhaseSpaceCut'),
            title = cms.string('gen. Phase-Space')
        ),
        cms.PSet(
            filter = cms.string('evtSelTrigger'),
            title = cms.string('Muon Trigger')
        ),
        cms.PSet(
            filter = cms.string('evtSelDataQuality'),
            title = cms.string('Data quality')
        ),
        cms.PSet(
            filter = cms.string('evtSelGlobalMuon'),
            title = cms.string('global Muon')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonEta'),
            title = cms.string('-2.1 < eta(Muon) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonPt'),
            title = cms.string('Pt(Muon) > 15 GeV')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauAntiOverlapWithMuonsVeto'),
            title = cms.string('Tau not overlapping w. Muon'),
        ),
        cms.PSet(
            filter = cms.string('evtSelTauEta'),
            title = cms.string('-2.1 < eta(Tau) < +2.1'),
        ),
        cms.PSet(
            filter = cms.string('evtSelTauPt'),
            title = cms.string('Pt(Tau) > 20 GeV'),
        ),
        cms.PSet(
            filter = cms.string('muonVbTfIdCutBgEstTTplusJetsEnriched'),
            title = cms.string('Muon VBTF id.'),
        ),
        cms.PSet(
            filter = cms.string('muonPFRelIsoCutBgEstTTplusJetsEnriched'),
            title = cms.string('Muon iso.')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauLeadTrk'),
            title = cms.string('Tau lead. Track find.'),
        ),
        cms.PSet(
            filter = cms.string('evtSelTauLeadTrkPt'),
            title = cms.string('Tau lead. Track Pt'),
        ),
        cms.PSet(
            filter = cms.string('evtSelTauTaNCdiscr'),
            title = cms.string('Tau TaNC discr.')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauProng'),
            title = cms.string('Tau 1||3-Prong')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauCharge'),
            title = cms.string('Charge(Tau) = +/-1')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauMuonVeto'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('vertexBgEstTTplusJetsEnriched'),
            title = cms.string('Vertex')
        ),
        cms.PSet(
            filter = cms.string('vertexQualityBgEstTTplusJetsEnriched'),
            title = cms.string('Vertex quality')
        ),
        cms.PSet(
            filter = cms.string('vertexPositionBgEstTTplusJetsEnriched'),
            title = cms.string('Vertex position')
        ),
        cms.PSet(
            filter = cms.string('muTauPairBgEstTTplusJetsEnriched'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairZeroChargeBgEstTTplusJetsEnriched'),
            title = cms.string('Charge(Muon+Tau) = 0')
        ),
        cms.PSet(
            filter = cms.string('jetEt40BgEstTTplusJetsEnriched'),
            title = cms.string('two E_{T} > 40 GeV Jets')
        ),
        cms.PSet(
            filter = cms.string('jetEt40bTagBgEstTTplusJetsEnriched'),
            title = cms.string('one E_{T} > 40 GeV Jet with b-Tag')
        ),
        cms.PSet(
            filter = cms.string('jetEt60BgEstTTplusJetsEnriched'),
            title = cms.string('one E_{T} > 60 GeV Jet')
        ),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerBgEstTTplusJetsEnriched',
                'tauHistManagerBgEstTTplusJetsEnriched',
                'diTauCandidateHistManagerBgEstTTplusJetsEnriched',
				#'diTauCandidateNSVfitHistManagerBgEstTTplusJetsEnriched',
                'caloMEtHistManagerBgEstTTplusJetsEnriched',
				'pfMEtHistManagerBgEstTTplusJetsEnriched',
                'jetHistManagerBgEstTTplusJetsEnriched',
                #'tauIdEffHistManagerBgEstTTplusJetsEnriched',
                'dataBinnerBgEstTTplusJetsEnriched'
            )
        )
    )
)

analysisSequenceBgEstTTplusJetsEnriched = cms.Sequence(analyzeEventsBgEstTTplusJetsEnriched)

#--------------------------------------------------------------------------------
# define (final) analysis sequence
#--------------------------------------------------------------------------------

bgEstTTplusJetsEnrichedAnalysisSequence = cms.Sequence(
    selectMuonsBgEstTTplusJetsEnriched
   + selectPrimaryVertexForMuTauBgEstTTplusJetsEnriched
   + produceMuTauPairsBgEstTTplusJetsEnriched + selectMuTauPairsBgEstTTplusJetsEnriched
   + selectJetsBgEstTTplusJetsEnriched
   + selectEventsBgEstTTplusJetsEnriched
   + analysisSequenceBgEstTTplusJetsEnriched
)
