import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.coreTools import *

from RecoTauTag.RecoTau.TauDiscriminatorTools import *
def redoPFTauDiscriminators(process,
                            oldPFTauLabel = cms.InputTag('shrinkingConePFTauProducer'),
                            newPFTauLabel = cms.InputTag('shrinkingConePFTauProducer'),
                            tauType='shrinkingConePFTau'):
    print 'Tau discriminators: ', oldPFTauLabel, '->', newPFTauLabel
    print 'Tau type: ', tauType
    tauSrc = 'PFTauProducer'

    tauDiscriminationSequence = process.patShrinkingConePFTauDiscrimination
    if tauType == 'fixedConeHighEffPFTau':
        tauDiscriminationSequence = process.patFixedConeHighEffPFTauDiscrimination
    elif tauType == 'fixedConePFTau':
        tauDiscriminationSequence = process.patFixedConePFTauDiscrimination
    elif tauType == 'shrinkingConePFTau':
        tauDiscriminationSequence = process.patShrinkingConePFTauDiscrimination
    elif tauType == 'caloTau':
        tauDiscriminationSequence = process.patCaloTauDiscrimination
        tauSrc = 'CaloTauProducer'

    process.makePatTaus.replace(process.patTaus, tauDiscriminationSequence*process.patTaus)

    massSearchReplaceParam(tauDiscriminationSequence, tauSrc, oldPFTauLabel, newPFTauLabel)

# switch to CaloTau collection
def switchToCaloTau(process,
                    pfTauLabel = cms.InputTag('shrinkingConePFTauProducer'),
                    caloTauLabel = cms.InputTag('caloRecoTauProducer')):
    process.tauMatch.src       = caloTauLabel
    process.tauGenJetMatch.src = caloTauLabel
    process.patTaus.tauSource = caloTauLabel
    process.patTaus.tauIDSources = cms.PSet(
        leadingTrackFinding = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut   = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackPtCut"),
        byIsolation         = cms.InputTag("caloRecoTauDiscriminationByIsolation"),
        againstElectron     = cms.InputTag("caloRecoTauDiscriminationAgainstElectron"),  
    )
    process.patTaus.addDecayMode = False
    ## Isolation is somewhat an issue, so we start just by turning it off
    print "NO PF Isolation will be computed for CaloTau (this could be improved later)"
    process.patTaus.isolation   = cms.PSet()
    process.patTaus.isoDeposits = cms.PSet()
    process.patTaus.userIsolation = cms.PSet()
    process.patDefaultSequence.remove(process.patPFCandidateIsoDepositSelection)
    process.patDefaultSequence.remove(process.patPFTauIsolation)
    ## adapt cleanPatTaus
    process.cleanPatTaus.preselection = 'tauID("leadingTrackFinding") > 0.5 & tauID("leadingTrackPtCut") > 0.5 & tauID("byIsolation") > 0.5 & tauID("againstElectron") > 0.5 & (signalTracks.size() = 1 | signalTracks.size() = 3)'

# internal auxiliary function to switch to **any** PFTau collection
def _switchToPFTau(process, pfTauLabelOld, pfTauLabelNew, pfTauType):

    print ' Taus: ', pfTauLabelOld, '->', pfTauLabelNew

    process.tauMatch.src       = pfTauLabelNew
    process.tauGenJetMatch.src = pfTauLabelNew
    process.tauIsoDepositPFCandidates.src = pfTauLabelNew
    process.tauIsoDepositPFCandidates.ExtractorPSet.tauSource = pfTauLabelNew
    process.tauIsoDepositPFChargedHadrons.src = pfTauLabelNew
    process.tauIsoDepositPFChargedHadrons.ExtractorPSet.tauSource = pfTauLabelNew
    process.tauIsoDepositPFNeutralHadrons.src = pfTauLabelNew
    process.tauIsoDepositPFNeutralHadrons.ExtractorPSet.tauSource = pfTauLabelNew
    process.tauIsoDepositPFGammas.src = pfTauLabelNew
    process.tauIsoDepositPFGammas.ExtractorPSet.tauSource = pfTauLabelNew
    process.patTaus.tauSource = pfTauLabelNew
    process.patTaus.tauIDSources = cms.PSet(
        leadingTrackFinding = cms.InputTag(pfTauType + "DiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut = cms.InputTag(pfTauType + "DiscriminationByLeadingTrackPtCut"),
        leadingPionPtCut = cms.InputTag(pfTauType + "DiscriminationByLeadingPionPtCut"),
        trackIsolation = cms.InputTag(pfTauType + "DiscriminationByTrackIsolation"),
        trackIsolationUsingLeadingPion = cms.InputTag(pfTauType + "DiscriminationByTrackIsolationUsingLeadingPion"),
        ecalIsolation = cms.InputTag(pfTauType + "DiscriminationByECALIsolation"),
        ecalIsolationUsingLeadingPion = cms.InputTag(pfTauType + "DiscriminationByECALIsolationUsingLeadingPion"),
        byIsolation = cms.InputTag(pfTauType + "DiscriminationByIsolation"),
        byIsolationUsingLeadingPion = cms.InputTag(pfTauType + "DiscriminationByIsolationUsingLeadingPion"),
        againstElectron = cms.InputTag(pfTauType + "DiscriminationAgainstElectron"),
        againstMuon = cms.InputTag(pfTauType + "DiscriminationAgainstMuon")
        #
        # CV: TaNC only trained for shrinkingCone PFTaus up to now,
        #     so cannot implement switch of TaNC based discriminators
        #     generically for all kinds of PFTaus yet...
        #
        #byTaNC = cms.InputTag(pfTauType + "DiscriminationByTaNC"),
        #byTaNCfrOnePercent = cms.InputTag(pfTauType + "DiscriminationByTaNCfrOnePercent"),
        #byTaNCfrHalfPercent = cms.InputTag(pfTauType + "DiscriminationByTaNCfrHalfPercent"),
        #byTaNCfrQuarterPercent = cms.InputTag(pfTauType + "DiscriminationByTaNCfrQuarterPercent"),
        #byTaNCfrTenthPercent = cms.InputTag(pfTauType + "DiscriminationByTaNCfrTenthPercent")
    )
    process.patTaus.decayModeSrc = cms.InputTag(pfTauType + "DecayModeProducer")


# switch to PFTau collection produced for fixed dR = 0.07 signal cone size
def switchToPFTauFixedCone(process,
                           pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
                           pfTauLabelNew = cms.InputTag('fixedConePFTauProducer')):
    _switchToPFTau(process, pfTauLabelOld, pfTauLabelNew, 'fixedConePFTau')
    #
    # CV: PFTauDecayMode objects produced only for shrinking cone reco::PFTaus in
    #     RecoTauTag/Configuration global_PFTau_22X_V00-02-01 and CMSSW_3_1_x tags,
    #     so need to disable embedding of PFTauDecayMode information into pat::Tau for now...
    #
    process.patTaus.addDecayMode = cms.bool(False)

# switch to PFTau collection produced for fixed dR = 0.15 signal cone size
def switchToPFTauFixedConeHighEff(process, 
                                  pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
                                  pfTauLabelNew = cms.InputTag('fixedConeHighEffPFTauProducer')):
    _switchToPFTau(process, pfTauLabelOld, pfTauLabelNew, 'fixedConeHighEffPFTau')
    #
    # CV: PFTauDecayMode objects produced only for shrinking cone reco::PFTaus in
    #     RecoTauTag/Configuration global_PFTau_22X_V00-02-01 and CMSSW_3_1_x tags,
    #     so need to disable embedding of PFTauDecayMode information into pat::Tau for now...
    #
    process.patTaus.addDecayMode = cms.bool(False)

# switch to PFTau collection produced for shrinking signal cone of size dR = 5.0/Et(PFTau)
def switchToPFTauShrinkingCone(process,
                               pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
                               pfTauLabelNew = cms.InputTag('shrinkingConePFTauProducer')):
    _switchToPFTau(process, pfTauLabelOld, pfTauLabelNew, 'shrinkingConePFTau')
    #
    # CV: TaNC only trained for shrinkingCone PFTaus up to now,
    #     so need to add TaNC based discriminators
    #     specifically for that case here...
    #
    process.patTaus.tauIDSources = cms.PSet(
        leadingTrackFinding = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackPtCut"),
        leadingPionPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingPionPtCut"),
        trackIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolation"),
        trackIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion"),
        ecalIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolation"),
        ecalIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion"),
        byIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByIsolation"),
        byIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion"),
        againstElectron = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectron"),
        againstMuon = cms.InputTag("shrinkingConePFTauDiscriminationAgainstMuon"),
        byTaNC = cms.InputTag("shrinkingConePFTauDiscriminationByTaNC"),
        byTaNCfrOnePercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrOnePercent"),
        byTaNCfrHalfPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrHalfPercent"),
        byTaNCfrQuarterPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent"),
        byTaNCfrTenthPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrTenthPercent")
    )

# Select switcher by string
def switchToPFTauByType(process, pfTauType=None, pfTauLabelNew=None,
                        pfTauLabelOld=cms.InputTag('shrinkingConePFTauProducer') ):
    mapping = { 'shrinkingConePFTau' : switchToPFTauShrinkingCone,
                'fixedConePFTau' : switchToPFTauFixedCone,
                'fixedConeHighEffPFTau' : switchToPFTauFixedConeHighEff,
                'caloTau' : switchToCaloTau }
    mapping[pfTauType](process, pfTauLabelOld=pfTauLabelOld, pfTauLabelNew=pfTauLabelNew)

# switch to PFTau collection that was default in PAT production in CMSSW_3_1_x release series
def switchTo31Xdefaults(process):
    switchToPFTauFixedCone(process)
    process.cleanPatTaus.preselection = cms.string('tauID("byIsolation") > 0')
    
# function to switch to **any** PFTau collection
# It is just to make internal function accessible externally
def switchToAnyPFTau(process,
                     pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
                     pfTauLabelNew = cms.InputTag('shrinkingConePFTauProducer'),
                     pfTauType='shrinkingConePFTau'):
    _switchToPFTau(process, pfTauLabelOld, pfTauLabelNew, pfTauType)
