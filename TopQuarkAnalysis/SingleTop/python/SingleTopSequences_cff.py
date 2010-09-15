import FWCore.ParameterSet.Config as cms

from TopQuarkAnalysis.SingleTop.SelectionCuts_top_group_control_samples_v3_cff import *

from PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi import *

from PhysicsTools.PatAlgos.patSequences_cff import *


#cFlavorHistoryProducer.matchedSrc = cms.InputTag("antikt5GenJets")
#bFlavorHistoryProducer.matchedSrc = cms.InputTag("antikt5GenJets")

#patElectrons.addElectronID = cms.bool(True)
#patElectrons.electronIDSources = electronIDSources

#patElectronIDs = cms.Sequence(simpleEleIdSequence)

#makeNewPatElectrons = cms.Sequence(patElectronIDs * patElectronIsolation * patElectrons)

patElectrons.usePV = cms.bool(False)
patMuons.usePV = cms.bool(False)


basePath = cms.Sequence(
   preselectedMETs *
   preselectedMuons *
   preselectedElectrons *
   looseElectrons
   )

flavorHistorySequence = cms.Sequence(
    cFlavorHistoryProducer *
    bFlavorHistoryProducer
    )

##Muons Sequences
baseMuonSequencePF = cms.Sequence(
    hltFilterDev *
    PVFilter *
    countLeptons *
    topElectrons *
    topMuons *
    preselectedJets *
    topJetsPF *
    countMuons 
    )


baseMuonSequence = cms.Sequence(
    hltFilterDev *
    PVFilter *
    HBHENoiseFilter *
    countLeptons *
    topMuons *
    topElectrons *
    preselectedJets *
    topJets *
    #topJetsPF *
    countMuons 
    )

baseMuonSequencePF = cms.Sequence(
    hltFilterDev *
    PVFilter *
    HBHENoiseFilter *
    countLeptons *
    topElectrons *
    topMuons *
    preselectedJets *
    #    topJets *
    topJetsPF *
    countMuonsPF 
    )

baseMuonAntiIsoSequence = cms.Sequence(
    hltFilterDev *
    PVFilter *
    vetoLeptonsIso *
    topMuonsAntiIso *
    countMuonsAntiIso *
    topElectronsAntiIso *
    vetoElectronsAntiIso 
    )
##


##Electron sequences
baseElectronSequence = cms.Sequence(
    hltFilterPhoton20 *
    PVFilter *
    HBHENoiseFilter *
    #    countLeptons *
    vetoLooseMuons *
    topElectrons *
    diElectrons *
    vetoDiElectrons *
    countElectrons *
#    electronIDIso *
    topMuons *
    preselectedJets 
    )

baseElectronSequencePF = cms.Sequence(
    hltFilterPhoton20 *
    PVFilter *
    HBHENoiseFilter *
    #    countLeptons *
    vetoLooseMuons *
    topElectrons *
    diElectrons *
    vetoDiElectrons *
    countElectrons *
#    electronIDIso *
    topMuons *
    preselectedJets 
    )



baseElectronAntiIsoSequence = cms.Sequence(
    hltFilterDev *
    PVFilter *
    vetoLeptonsIso *
    topElectronsAntiIso *
    countElectronsAntiIso *
    electronIDAntiIso *
    topMuonsAntiIso *
    vetoMuonsAntiIso 
    )


###All leptons, same path for electrons and muons

baseLeptonSequence = cms.Sequence(
    countLeptons *
    topMuons *
    topElectrons * 
    countLeptonsTight
    )

baseLeptonAntiIsoSequence = cms.Sequence(
    topElectrons *
    vetoLeptonsIso *
    countLeptonsAntiIso 
    )


baseJetMETSequence = cms.Sequence(
#Production
    topJets *
    bJets *
    forwardJets 
   )

baseJetMETSequencePF= cms.Sequence(
#Production
    topJetsPF *
    bJetsPF *
    forwardJetsPF 
   )

baseJetMETAntiIsoSequence = cms.Sequence(
#Production
    topJetsAntiIso *
    bJetsAntiIso *
    antiBJetsAntiIso *
    forwardJetsAntiIso 
   )

###Muon paths:

IsoMuons = cms.Sequence(
    baseMuonSequence +
    baseJetMETSequence
    )

IsoMuonsPF = cms.Sequence(
    baseMuonSequencePF +
    baseJetMETSequencePF
    )

AntiIsoMuons = cms.Sequence(
    baseMuonAntiIsoSequence +
    baseJetMETAntiIsoSequence
    )

###Electron paths
IsoElectrons = cms.Sequence(
    baseElectronSequence +
    baseJetMETSequence
    )

IsoElectronsPF = cms.Sequence(
    baseElectronSequence +
    baseJetMETSequencePF
    )

AntiIsoElectrons = cms.Sequence(
    baseElectronAntiIsoSequence +
    baseJetMETAntiIsoSequence
    )


PathMuonsIso = cms.Sequence(
    IsoMuons *
    countJetsNonTTBar
    )

#AntiIso Paths
PathMuonsAntiIso = cms.Sequence(
    AntiIsoMuons *
    countJetsNonTTBarAntiIso
    )



PathElectronsIso = cms.Sequence(
    IsoElectrons *
    countJetsNonTTBar
    )

#AntiIso Paths
PathElectronsAntiIso = cms.Sequence(
    AntiIsoElectrons *
    countJetsNonTTBarAntiIso
    )



###Tops Sequences
allTops= cms.Sequence(
    recoTops *
    boostedTops *
    boostedForwardJets 
    )

allAntiIsoTops= cms.Sequence(
    recoTopsAntiIso *
    boostedTopsAntiIsoTops *
    boostedForwardJetsAntiIsoTops
    )

allPseudoBJetsTops= cms.Sequence(
    recoTopsWSamplePseudoBTags *
    boostedTopsWSamplePseudoBTagsTops *
    boostedForwardJetsWSamplePseudoBTagsTops 
    )

allPseudoBJetsAntiIsoTops= cms.Sequence(
    recoTopsWSamplePseudoBTagsAntiIso *
    boostedTopsWSamplePseudoBTagsAntiIsoTops *
    boostedForwardJetsWSamplePseudoBTagsAntiIsoTops 
    )



TSampleMuonPF = cms.Sequence(
    IsoMuonsPF *
    countJetsPF *
    MTWFilterMuonsPF * 
    countBTagsPF * 
    countForwardJetsPF  
    )


TSampleMuon = cms.Sequence(
    IsoMuons *
    countJetsNonTTBar *
    MTWFilterMuons * 
    countBTags * 
    countForwardJets * 
    allTops *
    singleTopObservablesTSample #*
#    SingleTopWtransverseMassFilter
    )

WSampleMuon = cms.Sequence(
    IsoMuons *
    countJetsNonTTBar *
    pseudoBJets * 
    countWSamplePseudoBTags * 
    countForwardJets * 
    allPseudoBJetsTops *
    singleTopObservablesWSample #*
#    SingleTopWtransverseMassFilterWSample

    )

TTBarSampleMuon = cms.Sequence(
    IsoMuons *
    countBTags * 
    countForwardJets * 
    countJetsTTBar *
    allTops *
    singleTopObservablesTTBarSample #*
#    SingleTopWtransverseMassFilterTTBarSample
    )

QCDSampleMuon = cms.Sequence(
    AntiIsoMuons *
    countJetsNonTTBarAntiIso *
    bJetsAntiIso *
    forwardJetsAntiIso *
    countBTagsAntiIso *
    countForwardJetsAntiIso *
    allAntiIsoTops *
    singleTopObservablesAntiIso
    )


QCDWSampleMuon = cms.Sequence(
    AntiIsoMuons *
    countJetsNonTTBarAntiIso *
    pseudoBJetsAntiIso * 
    countWSamplePseudoBTagsAntiIso * 
    countForwardJetsAntiIso * 
    allPseudoBJetsAntiIsoTops *
    singleTopObservablesWSampleAntiIso 
    )




TSampleElectron = cms.Sequence(
    IsoElectrons *
    topMuons *
    countJetsNonTTBar *
    MTWFilterElectrons * 
    countBTags * 
    countForwardJets * 
    allTops *
    singleTopObservablesTSample #*
#    SingleTopWtransverseMassFilter
    )

TSampleElectronPF = cms.Sequence(
    IsoElectronsPF *
    countJetsPF *
    MTWFilterElectronsPF * 
    countBTagsPF * 
    countForwardJetsPF  
    )


WSampleElectron = cms.Sequence(
    IsoElectrons *
    countJetsNonTTBar *
    pseudoBJets * 
    countWSamplePseudoBTags * 
    countForwardJets * 
    allPseudoBJetsTops *
    singleTopObservablesWSample #*
#    SingleTopWtransverseMassFilterWSample

    )

TTBarSampleElectron = cms.Sequence(
    IsoElectrons *
    countBTags * 
    countForwardJets * 
    countJetsTTBar *
    allTops *
    singleTopObservablesTTBarSample #*
#    SingleTopWtransverseMassFilterTTBarSample
    )

QCDSampleElectron = cms.Sequence(
    AntiIsoElectrons *
    countJetsNonTTBarAntiIso *
    bJetsAntiIso *
    forwardJetsAntiIso *
    countBTagsAntiIso *
    countForwardJetsAntiIso *
    allAntiIsoTops *
    singleTopObservablesAntiIso
    )


QCDWSampleElectron = cms.Sequence(
    AntiIsoElectrons *
    countJetsNonTTBarAntiIso *
    pseudoBJetsAntiIso * 
    countWSamplePseudoBTagsAntiIso * 
    countForwardJetsAntiIso * 
    allPseudoBJetsAntiIsoTops *
    singleTopObservablesWSampleAntiIso 
    )






