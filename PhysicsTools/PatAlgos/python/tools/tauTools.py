import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.coreTools import *
from FWCore.GuiBrowsers.ConfigToolBase import *
from PhysicsTools.PatAlgos.tools.helpers import applyPostfix 
from RecoTauTag.RecoTau.TauDiscriminatorTools import *

def redoPFTauDiscriminators(process,
                            oldPFTauLabel = cms.InputTag('shrinkingConePFTauProducer'),
                            newPFTauLabel = cms.InputTag('shrinkingConePFTauProducer'),
                            tauType = 'shrinkingConePFTau', postfix = ""):
    print 'Tau discriminators: ', oldPFTauLabel, '->', newPFTauLabel
    print 'Tau type: ', tauType
    tauSrc = 'PFTauProducer'

    tauDiscriminationSequence = None
    if tauType == 'hpsPFTau':
        tauDiscriminationSequence =  cloneProcessingSnippet(process, process.patHPSPFTauDiscrimination, postfix)
    elif tauType == 'fixedConePFTau':
        tauDiscriminationSequence = cloneProcessingSnippet(process, process.patFixedConePFTauDiscrimination, postfix)
    elif tauType == 'shrinkingConePFTau':
        tauDiscriminationSequence = cloneProcessingSnippet(process, process.patShrinkingConePFTauDiscrimination, postfix)
    elif tauType == 'caloTau':
        tauDiscriminationSequence = cloneProcessingSnippet(process, process.patCaloTauDiscrimination, postfix)
        tauSrc = 'CaloTauProducer'
    else:
        raise StandardError, "Unkown tauType: '%s'"%tauType

    #process.makePatTaus.replace(process.patTaus, tauDiscriminationSequence*process.patTaus)
    applyPostfix(process,"makePatTaus",postfix).replace(
        applyPostfix(process,"patTaus",postfix),
        tauDiscriminationSequence*applyPostfix(process,"patTaus",postfix)
    )

    massSearchReplaceParam(tauDiscriminationSequence, tauSrc, oldPFTauLabel, newPFTauLabel)

# switch to CaloTau collection
def switchToCaloTau(process,
                    pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
                    pfTauLabelNew = cms.InputTag('caloRecoTauProducer'),
                    patTauLabel = "",
                    postfix = ""):
    print ' Taus: ', pfTauLabelOld, '->', pfTauLabelNew

    caloTauLabel = pfTauLabelNew
    applyPostfix(process, "tauMatch" + patTauLabel, postfix).src = caloTauLabel
    applyPostfix(process, "tauGenJetMatch"+ patTauLabel, postfix).src = caloTauLabel

    applyPostfix(process, "patTaus" + patTauLabel, postfix).tauSource = caloTauLabel
    applyPostfix(process, "patTaus" + patTauLabel, postfix).tauIDSources = _buildIDSourcePSet('caloRecoTau', classicTauIDSources, postfix)
#    applyPostfix(process, "patTaus" + patTauLabel, postfix).tauIDSources = cms.PSet(        
#        leadingTrackFinding = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackFinding" + postfix),
#        leadingTrackPtCut   = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackPtCut" + postfix),
#        trackIsolation      = cms.InputTag("caloRecoTauDiscriminationByTrackIsolation" + postfix),
#        ecalIsolation       = cms.InputTag("caloRecoTauDiscriminationByECALIsolation" + postfix),
#        byIsolation         = cms.InputTag("caloRecoTauDiscriminationByIsolation" + postfix),
#        againstElectron     = cms.InputTag("caloRecoTauDiscriminationAgainstElectron" + postfix),
#        againstMuon         = cms.InputTag("caloRecoTauDiscriminationAgainstMuon" + postfix)
#    )
    ## Isolation is somewhat an issue, so we start just by turning it off
    print "NO PF Isolation will be computed for CaloTau (this could be improved later)"
    applyPostfix(process, "patTaus" + patTauLabel, postfix).isolation   = cms.PSet()
    applyPostfix(process, "patTaus" + patTauLabel, postfix).isoDeposits = cms.PSet()
    applyPostfix(process, "patTaus" + patTauLabel, postfix).userIsolation = cms.PSet()

    ## adapt cleanPatTaus
    applyPostfix(process, "cleanPatTaus" + patTauLabel, postfix).preselection = \
      'tauID("leadingTrackFinding") > 0.5 & tauID("leadingTrackPtCut") > 0.5' \
     + ' & tauID("byIsolation") > 0.5 & tauID("againstElectron") > 0.5 & (signalTracks.size() = 1 | signalTracks.size() = 3)'

def _buildIDSourcePSet(pfTauType, idSources, postfix =""):
    """ Build a PSet defining the tau ID sources to embed into the pat::Tau """
    output = cms.PSet()
    for label, discriminator in idSources:
        setattr(output, label, cms.InputTag(pfTauType + discriminator + postfix))
    return output

def _switchToPFTau(process,
                   pfTauLabelOld,
                   pfTauLabelNew,
                   pfTauType,
                   idSources,
                   patTauLabel = "",
                   postfix = ""):
    """internal auxiliary function to switch to **any** PFTau collection"""  
    print ' Taus: ', pfTauLabelOld, '->', pfTauLabelNew
    
    applyPostfix(process, "tauMatch" + patTauLabel, postfix).src = pfTauLabelNew
    applyPostfix(process, "tauGenJetMatch" + patTauLabel, postfix).src = pfTauLabelNew
    
    applyPostfix(process, "tauIsoDepositPFCandidates" + patTauLabel, postfix).src = pfTauLabelNew
    applyPostfix(process, "tauIsoDepositPFCandidates" + patTauLabel, postfix).ExtractorPSet.tauSource = pfTauLabelNew
    applyPostfix(process, "tauIsoDepositPFChargedHadrons" + patTauLabel, postfix).src = pfTauLabelNew
    applyPostfix(process, "tauIsoDepositPFChargedHadrons" + patTauLabel, postfix).ExtractorPSet.tauSource = pfTauLabelNew
    applyPostfix(process, "tauIsoDepositPFNeutralHadrons" + patTauLabel, postfix).src = pfTauLabelNew
    applyPostfix(process, "tauIsoDepositPFNeutralHadrons" + patTauLabel, postfix).ExtractorPSet.tauSource = pfTauLabelNew
    applyPostfix(process, "tauIsoDepositPFGammas" + patTauLabel, postfix).src = pfTauLabelNew
    applyPostfix(process, "tauIsoDepositPFGammas" + patTauLabel, postfix).ExtractorPSet.tauSource = pfTauLabelNew
    
    applyPostfix(process, "patTaus" + patTauLabel, postfix).tauSource = pfTauLabelNew
    applyPostfix(process, "patTaus" + patTauLabel, postfix).tauIDSources = _buildIDSourcePSet(pfTauType, idSources, postfix)

    applyPostfix(process, "cleanPatTaus" + patTauLabel, postfix).preselection = \
      'tauID("leadingTrackFinding") > 0.5 & tauID("leadingPionPtCut") > 0.5 & tauID("byIsolationUsingLeadingPion") > 0.5' \
     + ' & tauID("againstMuon") > 0.5 & tauID("againstElectron") > 0.5' \
     + ' & (signalPFChargedHadrCands.size() = 1 | signalPFChargedHadrCands.size() = 3)'

# Name mapping for classic tau ID sources (present for fixed and shrinkingCones)
classicTauIDSources = [
    ("leadingTrackFinding", "DiscriminationByLeadingTrackFinding"),
    ("leadingTrackPtCut", "DiscriminationByLeadingTrackPtCut"),
    ("trackIsolation", "DiscriminationByTrackIsolation"),
    ("ecalIsolation", "DiscriminationByECALIsolation"),
    ("byIsolation", "DiscriminationByIsolation"),
    ("againstElectron", "DiscriminationAgainstElectron"),
    ("againstMuon", "DiscriminationAgainstMuon") ]

classicPFTauIDSources = [
    ("leadingPionPtCut", "DiscriminationByLeadingPionPtCut"),
    ("trackIsolationUsingLeadingPion", "DiscriminationByTrackIsolationUsingLeadingPion"),
    ("ecalIsolationUsingLeadingPion", "DiscriminationByECALIsolationUsingLeadingPion"),
    ("byIsolationUsingLeadingPion", "DiscriminationByIsolationUsingLeadingPion")]

# Tau Neural Classifier Discriminators
tancTauIDSources = [
    ("byTaNC", "DiscriminationByTaNC"),
    ("byTaNCfrOnePercent", "DiscriminationByTaNCfrOnePercent"),
    ("byTaNCfrHalfPercent", "DiscriminationByTaNCfrHalfPercent"),
    ("byTaNCfrQuarterPercent", "DiscriminationByTaNCfrQuarterPercent"),
    ("byTaNCfrTenthPercent", "DiscriminationByTaNCfrTenthPercent") ]

# Hadron-plus-strip(s) (HPS) Tau Discriminators
hpsTauIDSources = [
    ("leadingTrackFinding", "DiscriminationByDecayModeFinding"),
    ("byLooseIsolation", "DiscriminationByLooseIsolation"),
    ("byMediumIsolation", "DiscriminationByMediumIsolation"),
    ("byTightIsolation", "DiscriminationByTightIsolation"),
    ("againstElectron", "DiscriminationAgainstElectron"),
    ("againstMuon", "DiscriminationAgainstMuon")]

# Discriminators of new HPS + TaNC combined Tau id. algorithm
hpsTancTauIDSources = [
    ("leadingTrackFinding", "DiscriminationByLeadingTrackFinding"),
    ("leadingTrackPtCut", "DiscriminationByLeadingTrackPtCut"),
    ("leadingPionPtCut", "DiscriminationByLeadingPionPtCut"),
    ("byTaNCraw", "DiscriminationByTancRaw"),
    ("byTaNC", "DiscriminationByTanc"),
    ("byTaNCloose", "DiscriminationByTancLoose"),
    ("byTaNCmedium", "DiscriminationByTancMedium"),
    ("byTaNCtight", "DiscriminationByTancTight"),
    ("byDecayMode", "DiscriminationByDecayModeSelection"),
    ("byHPSloose", "DiscriminationByLooseIsolation"),
    ("byHPSmedium", "DiscriminationByMediumIsolation"),
    ("byHPStight", "DiscriminationByTightIsolation"),
    ("againstElectron", "DiscriminationAgainstElectron"),
    ("againstMuon", "DiscriminationAgainstMuon"),
    ("againstCaloMuon", "DiscriminationAgainstCaloMuon") ]

# switch to PFTau collection produced for fixed dR = 0.07 signal cone size
def switchToPFTauFixedCone(process,
                           pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
                           pfTauLabelNew = cms.InputTag('fixedConePFTauProducer'),
                           patTauLabel = "",
                           postfix = ""):
    fixedConeIDSources = copy.copy(classicTauIDSources)
    fixedConeIDSources.extend(classicPFTauIDSources)

    _switchToPFTau(process, pfTauLabelOld, pfTauLabelNew, 'fixedConePFTau', fixedConeIDSources,
                   patTauLabel = patTauLabel, postfix = postfix)

# switch to PFTau collection produced for shrinking signal cone of size dR = 5.0/Et(PFTau)
def switchToPFTauShrinkingCone(process,
                               pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
                               pfTauLabelNew = cms.InputTag('shrinkingConePFTauProducer'),
                               patTauLabel = "",
                               postfix = ""):
    shrinkingIDSources = copy.copy(classicTauIDSources)
    shrinkingIDSources.extend(classicPFTauIDSources)
    # Only shrinkingCone has associated TaNC discriminators, so add them here
    shrinkingIDSources.extend(tancTauIDSources)
    
    _switchToPFTau(process, pfTauLabelOld, pfTauLabelNew, 'shrinkingConePFTau', shrinkingIDSources,
                   patTauLabel = patTauLabel, postfix = postfix)

# switch to hadron-plus-strip(s) (HPS) PFTau collection
def switchToPFTauHPS(process, 
                     pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
                     pfTauLabelNew = cms.InputTag('hpsPFTauProducer'),
                     patTauLabel = "",
                     postfix = ""):
    _switchToPFTau(process, pfTauLabelOld, pfTauLabelNew, 'hpsPFTau', hpsTauIDSources,
                   patTauLabel = patTauLabel, postfix = postfix)
    
    ## adapt cleanPatTaus
    getattr(process, "cleanPatTaus" + patTauLabel).preselection = \
      'tauID("leadingTrackFinding") > 0.5 & tauID("byMediumIsolation") > 0.5' \
     + ' & tauID("againstMuon") > 0.5 & tauID("againstElectron") > 0.5'

# switch to hadron-plus-strip(s) (HPS) PFTau collection
def switchToPFTauHPSpTaNC(process, 
                          pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
                          pfTauLabelNew = cms.InputTag('hpsTancTaus'),
                          patTauLabel = "",
                          postfix = ""):
    _switchToPFTau(process, pfTauLabelOld, pfTauLabelNew, 'hpsTancTaus', hpsTancTauIDSources,
                   patTauLabel = patTauLabel, postfix = postfix)
    
    ## adapt cleanPatTaus
    getattr(process, "cleanPatTaus" + patTauLabel).preselection = \
      'tauID("leadingPionPtCut") > 0.5 & tauID("byHPSmedium") > 0.5' \
     + ' & tauID("againstMuon") > 0.5 & tauID("againstElectron") > 0.5'

# Select switcher by string
def switchToPFTauByType(process,
                        pfTauType = None,
                        pfTauLabelNew = None,
                        pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
                        patTauLabel = "",
                        postfix = "" ):
    mapping = { 'shrinkingConePFTau' : switchToPFTauShrinkingCone,
                'fixedConePFTau' : switchToPFTauFixedCone,
                'hpsPFTau' : switchToPFTauHPS,
                'caloRecoTau' : switchToCaloTau }
    mapping[pfTauType](process, pfTauLabelOld = pfTauLabelOld, pfTauLabelNew = pfTauLabelNew,
                       patTauLabel = patTauLabel, postfix = postfix)

# switch to PFTau collection that was default in PAT production in CMSSW_3_1_x release series
def switchTo31Xdefaults(process):
    switchToPFTauFixedCone(process)
    process.cleanPatTaus.preselection = cms.string('tauID("byIsolation") > 0')
    
class AddTauCollection(ConfigToolBase):

    """ Add a new collection of taus. Takes the configuration from the
    already configured standard tau collection as starting point;
    replaces before calling addTauCollection will also affect the
    new tau collections
    """
    _label='addTauCollection'
    _defaultParameters=dicttypes.SortedKeysDict()
    def __init__(self):
        ConfigToolBase.__init__(self)
        self.addParameter(self._defaultParameters, 'tauCollection',
                          self._defaultValue, 'Input tau collection', cms.InputTag)
        self.addParameter(self._defaultParameters, 'algoLabel',
                          self._defaultValue, "label to indicate the tau algorithm (e.g.'shrinkingCone')", str)
        self.addParameter(self._defaultParameters, 'typeLabel',
                          self._defaultValue, "label to indicate the type of constituents (either 'PFTau' or 'Tau')", str)
        self.addParameter(self._defaultParameters, 'doPFIsoDeposits',
                          True, "run sequence for computing particle-flow based IsoDeposits")
        self.addParameter(self._defaultParameters, 'standardAlgo',
                          "shrinkingCone", "standard algorithm label of the collection from which the clones " \
                         + "for the new tau collection will be taken from " \
                         + "(note that this tau collection has to be available in the event before hand)")
        self.addParameter(self._defaultParameters, 'standardType',
                          "PFTau", "standard constituent type label of the collection from which the clones " \
                         + " for the new tau collection will be taken from "\
                         + "(note that this tau collection has to be available in the event before hand)")
        
        self._parameters=copy.deepcopy(self._defaultParameters)
        self._comment = ""
        
    def getDefaultParameters(self):
        return self._defaultParameters

    def __call__(self,process,
                 tauCollection      = None,
                 algoLabel          = None,
                 typeLabel          = None,
                 doPFIsoDeposits    = None,
                 standardAlgo       = None,
                 standardType       = None):

        if tauCollection is None:
            tauCollection = self._defaultParameters['tauCollection'].value
        if algoLabel is None:
            algoLabel = self._defaultParameters['algoLabel'].value
        if typeLabel is None:
            typeLabel = self._defaultParameters['typeLabel'].value
        if doPFIsoDeposits is None:
            doPFIsoDeposits = self._defaultParameters['doPFIsoDeposits'].value
        if standardAlgo is None:
            standardAlgo = self._defaultParameters['standardAlgo'].value
        if standardType is None:
            standardType = self._defaultParameters['standardType'].value

        self.setParameter('tauCollection', tauCollection)
        self.setParameter('algoLabel', algoLabel)
        self.setParameter('typeLabel', typeLabel)
        self.setParameter('doPFIsoDeposits', doPFIsoDeposits)
        self.setParameter('standardAlgo', standardAlgo)
        self.setParameter('standardType', standardType)
   
        self.apply(process) 
        
    def toolCode(self, process):        
        tauCollection = self._parameters['tauCollection'].value
        algoLabel = self._parameters['algoLabel'].value
        typeLabel = self._parameters['typeLabel'].value
        doPFIsoDeposits = self._parameters['doPFIsoDeposits'].value
        standardAlgo = self._parameters['standardAlgo'].value
        standardType = self._parameters['standardType'].value

        ## disable computation of particle-flow based IsoDeposits
        ## in case tau is of CaloTau type
        if typeLabel == 'Tau':
#            print "NO PF Isolation will be computed for CaloTau (this could be improved later)"
            doPFIsoDeposits = False
 
        ## create old module label from standardAlgo
        ## and standardType and return
        def oldLabel(prefix = ''):
            if prefix == '':
                return "patTaus"
            else:
                return prefix + "PatTaus"

        ## capitalize first character of appended part
        ## when creating new module label
        ## (giving e.g. "patTausCaloRecoTau")
        def capitalize(label):
            return label[0].capitalize() + label[1:]    

        ## create new module label from old module
        ## label and return
        def newLabel(oldLabel):
            newLabel = oldLabel
            if ( oldLabel.find(standardAlgo) >= 0 and oldLabel.find(standardType) >= 0 ):
                oldLabel = oldLabel.replace(standardAlgo, algoLabel).replace(standardType, typeLabel)
            else:
                oldLabel = oldLabel + capitalize(algoLabel + typeLabel)
            return oldLabel

        ## clone module and add it to the patDefaultSequence
        def addClone(hook, **replaceStatements):
            ## create a clone of the hook with corresponding
            ## parameter replacements
            newModule = getattr(process, hook).clone(**replaceStatements)
            ## add the module to the sequence
            addModuleToSequence(hook, newModule)

        ## clone module for computing particle-flow IsoDeposits
        def addPFIsoDepositClone(hook, **replaceStatements):
            newModule = getattr(process, hook).clone(**replaceStatements)
            newModuleIsoDepositExtractor = getattr(newModule, "ExtractorPSet")
            setattr(newModuleIsoDepositExtractor, "tauSource", getattr(newModule, "src"))
            addModuleToSequence(hook, newModule)
            
        ## add module to the patDefaultSequence
        def addModuleToSequence(hook, newModule):
            hookModule = getattr(process, hook)
            ## add the new module with standardAlgo &
            ## standardType replaced in module label
            setattr(process, newLabel(hook), newModule)
            ## add new module to default sequence
            ## just behind the hookModule
            process.patDefaultSequence.replace( hookModule, hookModule*newModule )        

        ## add a clone of patTaus
        addClone(oldLabel(), tauSource = tauCollection)
        
        ## add a clone of selectedPatTaus    
        addClone(oldLabel('selected'), src = cms.InputTag(newLabel(oldLabel())))
        
        ## add a clone of cleanPatTaus    
        addClone(oldLabel('clean'), src=cms.InputTag(newLabel(oldLabel('selected'))))

        ## get attributes of new module
        newTaus = getattr(process, newLabel(oldLabel()))

        ## add a clone of gen tau matching
        addClone('tauMatch', src = tauCollection)
        addClone('tauGenJetMatch', src = tauCollection)

        ## add a clone of IsoDeposits computed based on particle-flow
        if doPFIsoDeposits:
            addPFIsoDepositClone('tauIsoDepositPFCandidates', src = tauCollection)
            addPFIsoDepositClone('tauIsoDepositPFChargedHadrons', src = tauCollection)
            addPFIsoDepositClone('tauIsoDepositPFNeutralHadrons', src = tauCollection)
            addPFIsoDepositClone('tauIsoDepositPFGammas', src = tauCollection)

        ## fix label for input tag
        def fixInputTag(x): x.setModuleLabel(newLabel(x.moduleLabel))

        ## provide patTau inputs with individual labels
        fixInputTag(newTaus.genParticleMatch)
        fixInputTag(newTaus.genJetMatch)
        fixInputTag(newTaus.isoDeposits.pfAllParticles)
        fixInputTag(newTaus.isoDeposits.pfNeutralHadron)
        fixInputTag(newTaus.isoDeposits.pfChargedHadron)
        fixInputTag(newTaus.isoDeposits.pfGamma)
        fixInputTag(newTaus.userIsolation.pfAllParticles.src)
        fixInputTag(newTaus.userIsolation.pfNeutralHadron.src)
        fixInputTag(newTaus.userIsolation.pfChargedHadron.src)
        fixInputTag(newTaus.userIsolation.pfGamma.src)
        
        ## set discriminators
        ## (using switchTauCollection functions)
        oldTaus = getattr(process, oldLabel())
#        if typeLabel == 'Tau':
#            switchToCaloTau(process,
#                            pfTauLabel = getattr(oldTaus, "tauSource"),
#                            caloTauLabel = getattr(newTaus, "tauSource"),
#                            patTauLabel = capitalize(algoLabel + typeLabel))
#        else:
        switchToPFTauByType(process, pfTauType = algoLabel + typeLabel,
                                pfTauLabelNew = getattr(newTaus, "tauSource"),
                                pfTauLabelOld = getattr(oldTaus, "tauSource"),
                                patTauLabel = capitalize(algoLabel + typeLabel))
       
addTauCollection=AddTauCollection()
