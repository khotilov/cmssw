import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.tauTools import *

def warningIsolation():
    print "WARNING: particle based isolation must be studied"

def adaptPFMuons(process,module):

    
    print "Adapting PF Muons "
    print "***************** "
    warningIsolation()
    print 
    module.useParticleFlow = True
    module.isolation   = cms.PSet()
    module.isoDeposits = cms.PSet(
        pfChargedHadrons = cms.InputTag("isoDepMuonWithCharged"),
        pfNeutralHadrons = cms.InputTag("isoDepMuonWithNeutral"),
        pfPhotons = cms.InputTag("isoDepMuonWithPhotons")
        )
    module.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("isoValMuonWithCharged"),
        pfNeutralHadrons = cms.InputTag("isoValMuonWithNeutral"),
        pfPhotons = cms.InputTag("isoValMuonWithPhotons")
        )
    # matching the pfMuons, not the standard muons.
    process.muonMatch.src = module.pfMuonSource

    print " muon source:", module.pfMuonSource
    print " isolation  :",
    print module.isolationValues
    print " isodeposits: "
    print module.isoDeposits
    print 
    

def adaptPFElectrons(process,module):
    # module.useParticleFlow = True

    print "Adapting PF Electrons "
    print "********************* "
    warningIsolation()
    print 
    module.useParticleFlow = True
    module.isolation   = cms.PSet()
    module.isoDeposits = cms.PSet(
        pfChargedHadrons = cms.InputTag("isoDepElectronWithCharged"),
        pfNeutralHadrons = cms.InputTag("isoDepElectronWithNeutral"),
        pfPhotons = cms.InputTag("isoDepElectronWithPhotons")
        )
    module.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("isoValElectronWithCharged"),
        pfNeutralHadrons = cms.InputTag("isoValElectronWithNeutral"),
        pfPhotons = cms.InputTag("isoValElectronWithPhotons")
        )

    # COLIN: since we take the egamma momentum for pat Electrons, we must
    # match the egamma electron to the gen electrons, and not the PFElectron.  
    # -> do not uncomment the line below.
    # process.electronMatch.src = module.pfElectronSource
    # COLIN: how do we depend on this matching choice? 

    print " PF electron source:", module.pfElectronSource
    print " isolation  :"
    print module.isolationValues
    print " isodeposits: "
    print module.isoDeposits
    print 
    
    


##     print "Temporarily switching off isolation & isoDeposits for PF Electrons"
##     module.isolation   = cms.PSet()
##     module.isoDeposits = cms.PSet()
##     print "Temporarily switching off electron ID for PF Electrons"
##     module.isolation   = cms.PSet()
##     module.addElectronID = False
##     if module.embedTrack.value(): 
##         module.embedTrack = False
##         print "Temporarily switching off electron track embedding"
##     if module.embedGsfTrack.value(): 
##         module.embedGsfTrack = False
##         print "Temporarily switching off electron gsf track embedding"
##     if module.embedSuperCluster.value(): 
##         module.embedSuperCluster = False
##         print "Temporarily switching off electron supercluster embedding"

def adaptPFPhotons(process,module):
    raise RuntimeError, "Photons are not supported yet"

def adaptPFJets(process,module):
    module.embedCaloTowers   = False

def adaptPFTaus(process,tauType = 'fixedConePFTau' ): #tauType can be changed only to shrinkig cone one, otherwise request is igonred
    oldTaus = process.allLayer1Taus.tauSource
    process.allLayer1Taus.tauSource = cms.InputTag("allLayer0Taus")

    if tauType == 'shrinkingConePFTau': 
        print "PF2PAT: tauType changed from default \'fixedConePFTau\' to \'shrinkingConePFTau\'"
        process.allLayer0TausDiscrimination.PFTauProducer = cms.InputTag(tauType+"Producer")
        process.allLayer0Taus.src = cms.InputTag(tauType+"Producer")
        process.pfTauSequence.replace(process.fixedConePFTauProducer,
                                      process.shrinkingConePFTauProducer)
        
    if (tauType != 'shrinkingConePFTau' and tauType != 'fixedConePFTau'):
        print "PF2PAT: TauType \'"+tauType+"\' is not supported. Default \'fixedConePFTau\' is used instead."
        tauType = 'fixedConePFTau'
        
    redoPFTauDiscriminators(process, cms.InputTag(tauType+'Producer'),
                            process.allLayer1Taus.tauSource,
                            tauType)
    switchToAnyPFTau(process, oldTaus, process.allLayer1Taus.tauSource, tauType)


def addPFCandidates(process,src,patLabel='PFParticles',cut=""):
    from PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi import allLayer1PFParticles
    # make modules
    producer = allLayer1PFParticles.clone(pfCandidateSource = src)
    filter   = cms.EDFilter("PATPFParticleSelector", 
                    src = cms.InputTag('allLayer1' + patLabel), 
                    cut = cms.string(cut))
    counter  = cms.EDFilter("PATCandViewCountFilter",
                    minNumber = cms.uint32(0),
                    maxNumber = cms.uint32(999999),
                    src       = cms.InputTag('selectedLayer1' + patLabel))
    # add modules to process
    setattr(process, 'allLayer1'      + patLabel, producer)
    setattr(process, 'selectedLayer1' + patLabel, filter)
    setattr(process, 'countLayer1'    + patLabel, counter)
    # insert into sequence
    process.allLayer1Objects.replace(process.allLayer1Summary, producer +  process.allLayer1Summary)
    process.selectedLayer1Objects.replace(process.selectedLayer1Summary, filter +  process.selectedLayer1Summary)
    process.countLayer1Objects    += counter
    # summary tables
    process.allLayer1Summary.candidates.append(cms.InputTag('allLayer1' + patLabel))
    process.selectedLayer1Summary.candidates.append(cms.InputTag('selectedLayer1' + patLabel))

def switchToPFMET(process,input=cms.InputTag('pfMET')):
    print 'MET: using ', input
    oldMETSource = process.layer1METs.metSource
    process.layer1METs.metSource = input
    process.layer1METs.addMuonCorrections = False
    #process.patJetMETCorrections.remove(process.patMETCorrections)
    process.patDefaultSequence.remove(process.patMETCorrections)

def switchToPFJets(process,input=cms.InputTag('pfNoTau')):
    print 'Jets: using ', input
    switchJetCollection(process,
                        input,
                        doJTA=True,
                        doBTagging=True,
                        jetCorrLabel=None, 
                        doType1MET=False)  
    adaptPFJets(process, process.allLayer1Jets)

def usePF2PAT(process,runPF2PAT=True,addElectrons=False):

    # PLEASE DO NOT CLOBBER THIS FUNCTION WITH CODE SPECIFIC TO A GIVEN PHYSICS OBJECT.
    # CREATE ADDITIONAL FUNCTIONS IF NEEDED. 


    """Switch PAT to use PF2PAT instead of AOD sources. if 'runPF2PAT' is true, we'll also add PF2PAT in front of the PAT sequence"""
    """If 'addElectrons' is True electrons are added and standard PAT cleaning wrt them is performed"""

    # -------- CORE ---------------
    if runPF2PAT:
        process.load("PhysicsTools.PFCandProducer.PF2PAT_cff")

#        process.dump = cms.EDAnalyzer("EventContentAnalyzer")
        process.patDefaultSequence.replace(process.allLayer1Objects,
                                           process.PF2PAT +
                                           process.allLayer1Objects
                                           )

#    if not addElectrons:
#        removeCleaning(process)
    
    # -------- OBJECTS ------------
    # Muons
    adaptPFMuons(process,process.allLayer1Muons)

    
    # Electrons
    if addElectrons:
        adaptPFElectrons(process,process.allLayer1Electrons)
    else:
        print "PFElectrons switched off for the moment."

#COLIN: the following should be in a specific function, for the ease of understanding and maintenance.
#       please provide such a function if still needed

#    process.electronMatch.src = process.allLayer1Electrons.pfElectronSource
#    process.patDefaultSequence.remove(process.patElectronId)
#    process.patDefaultSequence.remove(process.patElectronIsolation)
#    if not addElectrons:
#        print "Temporarily switching off electrons completely"
#        removeSpecificPATObjects(process,['Electrons'])
#        process.patDefaultSequence.remove(process.patElectronId)
#        process.patDefaultSequence.remove(process.patElectronIsolation)
#        #process.countLayer1Leptons.countElectrons = False
#    else:
#        print "Standard electrons added and then PAT-cleaning wrt electrons switched on"
#        #tune PAT-cleaning
#        process.cleanLayer1Objects.remove(process.cleanLayer1Hemispheres)
#        # store only isolated electrons passing the RobustLoose-ID, without any overlap with muons (to be tuned)
#        electronPresel = 'electronID("eidRobustLoose") > 0 '
#        electronPresel += '&& trackIso < 3 '
#        electronPresel += '&& ecalIso < 5 '
#        electronPresel += '&& hcalIso < 5 '
#        #electronPresel += '&& pt > 10 '
#        print "A \"good electron\" definition for cleaner: "+electronPresel
#        process.cleanLayer1Electrons.checkOverlaps.muons.requireNoOverlaps = True
#        process.cleanLayer1Electrons.finalCut = electronPresel
#        # clean taus wrt electrons
#        process.cleanLayer1Taus.preselection = ''
#        # a trick to remove overlaps wrt other objects than electrons
#        myOverlaps = cms.PSet( electrons = cms.PSet() )
#        myOverlaps.electrons = process.cleanLayer1Taus.checkOverlaps.electrons
#        process.cleanLayer1Taus.checkOverlaps = myOverlaps
#        process.cleanLayer1Taus.checkOverlaps.electrons.preselection = electronPresel # really needed?
#        process.cleanLayer1Taus.checkOverlaps.electrons.requireNoOverlaps = True
#        # clean jets wrt electrons
#        # a trick to remove overlaps wrt other objects than electrons
#        myOverlaps.electrons = process.cleanLayer1Jets.checkOverlaps.electrons
#        process.cleanLayer1Jets.checkOverlaps = myOverlaps
#        process.cleanLayer1Jets.checkOverlaps.electrons.preselection = electronPresel # really needed?
#        process.cleanLayer1Jets.checkOverlaps.electrons.deltaR = 0.3
#        process.cleanLayer1Jets.checkOverlaps.electrons.requireNoOverlaps = True
    
    # Photons
    print "Temporarily switching off photons completely"
    removeSpecificPATObjects(process,['Photons'])
    process.patDefaultSequence.remove(process.patPhotonIsolation)
    
    # Jets
    switchToPFJets( process, cms.InputTag('pfNoTau') )
    
    # Taus
    adaptPFTaus( process ) #default (i.e. fixedConePFTau)
    #adaptPFTaus( process, tauType='shrinkingConePFTau' )
    
    # MET
    switchToPFMET(process, cms.InputTag('pfMET'))
    
    # Unmasked PFCandidates
    addPFCandidates(process,cms.InputTag('pfNoJet'),patLabel='PFParticles',cut="")
