#
#  SUSY-PAT configuration fragment
#
#  PAT configuration for the SUSY group - 33X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV7


import FWCore.ParameterSet.Config as cms

def addDefaultSUSYPAT(process, HLTMenu='HLT', JetMetCorrections='Summer09_7TeV_ReReco332'):
    loadPAT(process,JetMetCorrections)
    addJetMET(process)
    loadPATTriggers(process,HLTMenu)

    # Full path
    process.seqSUSYDefaultSequence = cms.Sequence( process.jpt * process.addTrackJets
                                                   *process.patDefaultSequence
                                                   * process.patTrigger*process.patTriggerEvent )

def loadPAT(process,JetMetCorrections='Summer09_7TeV_ReReco332'):
    #-- PAT standard config -------------------------------------------------------
    process.load("PhysicsTools.PatAlgos.patSequences_cff")

    #-- Changes for electron and photon ID ----------------------------------------
    # Turn off photon-electron cleaning (i.e., flag only)
    process.cleanLayer1Photons.checkOverlaps.electrons.requireNoOverlaps = False

    # Remove embedding of superClusters, will keep entire superCluster collection
    process.allLayer1Electrons.embedSuperCluster = False
    process.allLayer1Photons.embedSuperCluster   = False
    
    #-- Tuning of Monte Carlo matching --------------------------------------------
    # Also match with leptons of opposite charge
    process.electronMatch.checkCharge = False
    process.muonMatch.checkCharge     = False
    process.tauMatch.checkCharge      = False

    #-- Jet corrections -----------------------------------------------------------
    process.jetCorrFactors.corrSample = JetMetCorrections ## 'Summer09' for 10TeV, 'Summer09_7TeV' for 7TeV no ReReco
    
def loadPATTriggers(process,HLTMenu='HLT'):
    #-- Trigger matching ----------------------------------------------------------
    from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
    switchOnTrigger( process )
    process.patTriggerSequence.remove( process.patTriggerMatcher )
    process.patTriggerEvent.patTriggerMatches  = ()
    # If we have to rename the default trigger menu
    process.patTrigger.processName = HLTMenu
    process.patTriggerEvent.processName = HLTMenu

def addJetMET(process):
    #-- Extra Jet/MET collections -------------------------------------------------
    # Add a few jet collections...
    from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
    
    #-- Jet plus tracks -----------------------------------------------------------
    process.load("PhysicsTools.PatAlgos.recoLayer0.jetPlusTrack_cff")
    process.jpt = cms.Sequence( process.jptCaloJets )
	
    # CaloJets
    addJetCollection(process, cms.InputTag('iterativeCone5CaloJets'),
                     'IC5',
                     doJTA            = True,
                     doBTagging       = True,
                     jetCorrLabel     = ('IC5','Calo'),
                     doType1MET       = True,
                     doL1Cleaning     = True,
                     doL1Counters     = True,
                     doJetID          = True,
		     jetIdLabel       = "ic5",
                     genJetCollection = cms.InputTag("iterativeCone5GenJets")
                     )
    addJetCollection(process,cms.InputTag('sisCone5CaloJets'),
                     'SC5',
                     doJTA        = True,
                     doBTagging   = True,
                     jetCorrLabel = ('SC5','Calo'),
                     doType1MET   = True,
                     doL1Cleaning = True,
                     doL1Counters = True,
                     doJetID      = True,
                     jetIdLabel   = "sc5",
                     genJetCollection=cms.InputTag("sisCone5GenJets")
                     )
    # PF jets
    addJetCollection(process,cms.InputTag('ak5PFJets'),
                     'AK5PF',
                     doJTA        = True,
                     doBTagging   = True,
                     jetCorrLabel = ('AK5','PF'),
                     doType1MET   = False,
                     doL1Cleaning = True,
                     doL1Counters = True,
                     doJetID      = False,
                     genJetCollection=cms.InputTag("ak5GenJets")
                     )
    addJetCollection(process,cms.InputTag('sisCone5PFJets'),
                     'SC5PF',
                     doJTA        = True,
                     doBTagging   = True,
                     jetCorrLabel = ('SC5','PF'),
                     doType1MET   = False,
                     doL1Cleaning = True,
                     doL1Counters = True,
                     doJetID      = False,
                     genJetCollection=cms.InputTag("sisCone5GenJets")
                     )
    
    # JPT jets
    addJetCollection(process,cms.InputTag('JetPlusTrackZSPCorJetAntiKt5'),
                     'AK5JPT',
                     doJTA        = True,
                     doBTagging   = True,
                     jetCorrLabel = None,
                     doType1MET   = False,
                     doL1Cleaning = True,
                     doL1Counters = True,
                     doJetID      = False,
                     genJetCollection = cms.InputTag("ak5GenJets")
                     )
	
    #-- Track Jets ----------------------------------------------------------------
    # Select tracks for track jets
    process.load("PhysicsTools.RecoAlgos.TrackWithVertexSelector_cfi")
    process.trackWithVertexSelector.src              = cms.InputTag("generalTracks")
    process.trackWithVertexSelector.ptMax            = cms.double(500.0) 
    process.trackWithVertexSelector.normalizedChi2   = cms.double(100.0)
    process.trackWithVertexSelector.vertexTag        = cms.InputTag("offlinePrimaryVertices")
    process.trackWithVertexSelector.copyTrajectories = cms.untracked.bool(False)
    process.trackWithVertexSelector.vtxFallback      = cms.bool(False)
    process.trackWithVertexSelector.useVtx           = cms.bool(False)
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    process.tracksForJets = cms.EDProducer("ConcreteChargedCandidateProducer",
                                           src = cms.InputTag("trackWithVertexSelector"),
                                           particleType = cms.string('pi+')
                                           )
	
    #add track jet collection
    process.load('RecoJets.JetProducers.ak5TrackJets_cfi')
    process.addTrackJets = cms.Sequence(  process.trackWithVertexSelector
                                          * process.tracksForJets
                                          * process.ak5TrackJets )
    addJetCollection(process,cms.InputTag('ak5TrackJets'),
                     'AK5Track',
                     doJTA        = False,
                     doBTagging   = True,
                     jetCorrLabel = None,
                     doType1MET   = False,
                     doL1Cleaning = True,
                     doL1Counters = True,
                     doJetID      = False,
                     genJetCollection = cms.InputTag("ak5GenJets")
                     )
    
    #-- Tune contents of jet collections  -----------------------------------------
    for jetName in ( '', 'IC5', 'SC5' , 'AK5PF', 'SC5PF', 'AK5JPT', 'AK5Track'):
        module = getattr(process,'allLayer1Jets'+jetName)
        module.addTagInfos = False    # Remove tag infos
        module.embedGenJetMatch = False # Only keep reference, since we anyway keep the genJet collections
 
    # Add tcMET and PFMET
    from PhysicsTools.PatAlgos.tools.metTools import addTcMET, addPfMET
    addTcMET(process,'TC')
    addPfMET(process,'PF')

    # Rename default jet collection for uniformity
    process.cleanLayer1JetsAK5 = process.cleanLayer1Jets
    process.layer1METsAK5      = process.layer1METs

    # Modify subsequent modules
    process.cleanLayer1Hemispheres.patJets = process.cleanLayer1JetsAK5.label()
    process.countLayer1Jets.src            = process.cleanLayer1JetsAK5.label()

    # Modify counters' input
    process.allLayer1Summary.candidates.remove(cms.InputTag('layer1METs'))
    process.allLayer1Summary.candidates.append(cms.InputTag('layer1METsAK5'))
    process.cleanLayer1Summary.candidates.remove(cms.InputTag('cleanLayer1Jets'))
    process.cleanLayer1Summary.candidates.append(cms.InputTag('cleanLayer1JetsAK5'))
    # Add new jet collections to counters (MET done automatically)
    for jets in ( 'IC5', 'SC5', 'AK5PF', 'SC5PF', 'AK5JPT', 'AK5Track'):
        process.allLayer1Summary.candidates.append(cms.InputTag('allLayer1Jets'+jets))
        process.selectedLayer1Summary.candidates.append(cms.InputTag('selectedLayer1Jets'+jets))
        process.cleanLayer1Summary.candidates.append(cms.InputTag('cleanLayer1Jets'+jets))
	

def removeMCDependence( process ):
    #-- Remove MC dependence ------------------------------------------------------
    from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
    removeMCMatching(process, 'All')

def getSUSY_pattuple_outputCommands( process ):
    return [ # PAT Objects
        'keep *_cleanLayer1Photons_*_*',
        'keep *_cleanLayer1Electrons_*_*',
        'keep *_cleanLayer1Muons_*_*',
        'keep *_cleanLayer1Taus_*_*',
        'keep *_cleanLayer1Jets*_*_*',       # All Jets
        'keep *_layer1METs*_*_*',            # All METs
        'keep *_cleanLayer1Hemispheres_*_*',
        'keep *_cleanLayer1PFParticles_*_*',
        # Generator information
        'keep GenEventInfoProduct_generator_*_*',
        'keep GenRunInfoProduct_generator_*_*',
        # Generator particles/jets/MET
        'keep recoGenParticles_genParticles_*_*',
        'keep recoGenJets_iterativeCone5GenJets_*_*',
        'keep recoGenJets_sisCone5GenJets_*_*',
        'keep recoGenJets_ak5GenJets_*_*',
        'keep recoGenMETs_*_*_*',
        # Trigger information
        'keep edmTriggerResults_TriggerResults_*_HLT',
        'keep *_hltTriggerSummaryAOD_*_*',
        'keep L1GlobalTriggerObjectMapRecord_*_*_*',
        'keep L1GlobalTriggerReadoutRecord_*_*_*',
        # Others
        'keep *_muon*METValueMapProducer_*_*',   # Muon corrections to MET
        'keep *_offlinePrimaryVertices_*_*',
        'keep *_offlineBeamSpot_*_*',
        'keep *_towerMaker_*_*',                 # Keep CaloTowers for cross-cleaning
        'keep recoTracks_generalTracks_*_*',
        'keep recoSuperClusters_corrected*_*_*',
        'keep recoConversions_conversions_*_*',
        'keep recoTracks_*onversions_*_*',
        'keep HcalNoiseSummary_*_*_*', #Keep the one in RECO
        'keep recoPFCandidates_particleFlow_*_*'
        ] 

