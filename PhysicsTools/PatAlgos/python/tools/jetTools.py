import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.helpers import *


def patchJetCorrFactors_(jetCorrFactors, newAlgo):
    """
    ------------------------------------------------------------------
    Patch to be called from:
       * switchJECSet_
       * switchJECParameters
    This function can safely be removed as soon as the L7Parton
    corrections for AK5 and AK7 are available.

    jetCorrFactors : jetCorrFactors module
    ------------------------------------------------------------------
    """
    if (newAlgo == "AK5"):
        ## voice a note to the user
        print "NOTE TO USER: L7Parton is currently taken from SC5 instead of AK5 "
        print "              This is an intermediate solution for the time being."
        ## redirect the L7Parton correction in case of AK5 or AK7
        corrLevels = getattr(jetCorrFactors, 'corrLevels').value()
        corrLevels.L7Parton = corrLevels.L7Parton.value().replace(newAlgo, 'SC5')
    if (newAlgo == "AK7"):
        ## voice a note to the user
        print "NOTE TO USER: L7Parton is currently taken from SC7 instead of AK7 "
        print "              This is an intermediate solution for the time being."
        ## redirect the L7Parton correction in case of AK5 or AK7        
        corrLevels = getattr(jetCorrFactors, 'corrLevels').value()
        corrLevels.L7Parton = corrLevels.L7Parton.value().replace(newAlgo, 'SC7')


def switchJECSet(process,
                 newName
                 ):
    """
    ------------------------------------------------------------------
    replace tags in the JetCorrFactorsProducer for end-users:

    process : process
    newName : new correction sample
    ------------------------------------------------------------------    
    """
    jetCorrFactors = getattr(process, 'jetCorrFactors')
    jetCorrFactors.corrSample = newName


def switchJECParameters(jetCorrFactors,
                        newAlgo,
                        newType="Calo",
                        oldAlgo="AK5",
                        oldType="Calo"
                        ):
    """
    ------------------------------------------------------------------    
    replace tags in the JetCorrFactorsProducer

    jetCorrFactors : jetCorrFactors module
    newAlgo        : label of new jet algo [AK5,  SC5,   KT6, ...]
    newType        : label of new jet type [Calo, Pflow, Jpt, ...]
    oldAlgo        : label of old jet alog [AK5,  SC5,   KT6, ...]
    oldType        : label of old jet type [Calo, Pflow, Jpt, ...]
    ------------------------------------------------------------------    
    """
    ## check jet correction steps; the L5Flavor step
    ## is not in the list as it is said NOT to be
    ## dependendent on the specific jet algorithm

    ## do the replacement, the first replacement is newAlgo and newType (as for 
    ## L2 and L3) the second repleacement is for newAlgo only (as for L5 and L7)
    def setCorrLevel(corrLevel):
        if (corrLevel != "none"):
            return corrLevel.value().replace(oldAlgo+oldType,newAlgo+newType).replace(oldAlgo,newAlgo)

    ## get the parameters and change it's attributes for L1 to L7
    corrLevels = getattr(jetCorrFactors, 'corrLevels').value()
    corrLevels.L1Offset   = setCorrLevel(corrLevels.L1Offset  )
    corrLevels.L2Relative = setCorrLevel(corrLevels.L2Relative)
    corrLevels.L3Absolute = setCorrLevel(corrLevels.L3Absolute)
    corrLevels.L4EMF      = setCorrLevel(corrLevels.L4EMF     )
    corrLevels.L6UE       = setCorrLevel(corrLevels.L6UE      )
    corrLevels.L7Parton   = setCorrLevel(corrLevels.L7Parton  )
    ##
    ## patch the jetCorrFactors untill the L7Parton corrections are not available yet
    ##
    patchJetCorrFactors_(jetCorrFactors, newAlgo)    


def runBTagging(process,
                jetCollection,
                label
                ) :
    """
    ------------------------------------------------------------------        
    define sequence to run b tagging on AOD input for a given jet
    collection including a JetTracksAssociatorAtVertex with name
    'jetTracksAssociatorAtVertex' + 'label'

    process       : process       
    jetCollection : input jet collection
    label         : postfix label to identify new sequence/modules

    the sequence is added to the process but not to any path; return
    value is a pair of (sequence, labels) where 'sequence' is the
    cms.Sequence, and 'labels' is a vector with the following entries:
    
     * labels['jta']      = the name of the JetTrackAssociator module
     * labels['tagInfos'] = a list of names of the TagInfo modules
     * labels['jetTags '] = a list of names of the JetTag modules
    ------------------------------------------------------------------        
    """    
    if (label == ''):
        ## label is not allowed to be empty
        raise ValueError, "label for re-running b tagging is not allowed to be empty"        

    ## import track associator & b tag configuration
    process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
    from RecoJets.JetAssociationProducers.ak5JTA_cff import ak5JetTracksAssociatorAtVertex
    process.load("RecoBTag.Configuration.RecoBTag_cff")
    import RecoBTag.Configuration.RecoBTag_cff as btag

    # add negative tag infos
    import PhysicsTools.PatAlgos.recoLayer0.bTagging_cff as nbtag
    
    ## define jetTracksAssociator; for switchJetCollection
    ## the label is 'AOD' as empty labels will lead to crashes
    ## of crab. In this case the postfix label is skiped,
    ## otherwise a postfix label is added as for the other
    ## labels
    jtaLabel = 'jetTracksAssociatorAtVertex'

    if (not label == 'AOD'):
        jtaLabel  += label
    ## define tag info labels (compare with jetProducer_cfi.py)        
    ipTILabel = 'impactParameterTagInfos'     + label
    svTILabel = 'secondaryVertexTagInfos'     + label
   #nvTILabel = 'secondaryVertexNegativeTagInfos'     + label
    seTILabel = 'softElectronTagInfos'        + label
    smTILabel = 'softMuonTagInfos'            + label
    
    ## produce tag infos
    setattr( process, ipTILabel, btag.impactParameterTagInfos.clone(jetTracks = cms.InputTag(jtaLabel)) )
    setattr( process, svTILabel, btag.secondaryVertexTagInfos.clone(trackIPTagInfos = cms.InputTag(ipTILabel)) )
   #setattr( process, nvTILabel, nbtag.secondaryVertexNegativeTagInfos.clone(trackIPTagInfos = cms.InputTag(ipTILabel)) )
    setattr( process, seTILabel, btag.softElectronTagInfos.clone(jets = jetCollection) )
    setattr( process, smTILabel, btag.softMuonTagInfos.clone(jets = jetCollection) )

    ## make VInputTag from strings
    def vit(*args) : return cms.VInputTag( *[ cms.InputTag(x) for x in args ] )
    
    ## produce btags
    setattr( process, 'jetBProbabilityBJetTags'+label, btag.jetBProbabilityBJetTags.clone(tagInfos = vit(ipTILabel)) )
    setattr( process, 'jetProbabilityBJetTags' +label, btag.jetProbabilityBJetTags.clone (tagInfos = vit(ipTILabel)) )
    setattr( process, 'trackCountingHighPurBJetTags'+label, btag.trackCountingHighPurBJetTags.clone(tagInfos = vit(ipTILabel)) )
    setattr( process, 'trackCountingHighEffBJetTags'+label, btag.trackCountingHighEffBJetTags.clone(tagInfos = vit(ipTILabel)) )
    setattr( process, 'simpleSecondaryVertexBJetTags'+label, btag.simpleSecondaryVertexBJetTags.clone(tagInfos = vit(svTILabel)) )
   #setattr( process, 'simpleSecondaryVertexNegativeBJetTags'+label, nbtag.simpleSecondaryVertexNegativeBJetTags.clone(tagInfos = vit(nvTILabel)) )
    setattr( process, 'combinedSecondaryVertexBJetTags'+label, btag.combinedSecondaryVertexBJetTags.clone(tagInfos = vit(ipTILabel, svTILabel)) )
    setattr( process, 'combinedSecondaryVertexMVABJetTags'+label, btag.combinedSecondaryVertexMVABJetTags.clone(tagInfos = vit(ipTILabel, svTILabel)) )
    setattr( process, 'softElectronByPtBJetTags'+label, btag.softElectronByPtBJetTags.clone(tagInfos = vit(seTILabel)) )
    setattr( process, 'softElectronByIP3dBJetTags'+label, btag.softElectronByIP3dBJetTags.clone(tagInfos = vit(seTILabel)) )
    setattr( process, 'softMuonBJetTags'+label, btag.softMuonBJetTags.clone(tagInfos = vit(smTILabel)) )
    setattr( process, 'softMuonByPtBJetTags'+label, btag.softMuonByPtBJetTags.clone(tagInfos = vit(smTILabel)) )
    setattr( process, 'softMuonByIP3dBJetTags'+label, btag.softMuonByIP3dBJetTags.clone(tagInfos = vit(smTILabel)) )
    
    ## define vector of (output) labels
    labels = { 'jta'      : jtaLabel, 
               'tagInfos' : (ipTILabel,svTILabel,seTILabel,smTILabel), 
               'jetTags'  : [ (x + label) for x in ('jetBProbabilityBJetTags',
                                                    'jetProbabilityBJetTags',
                                                    'trackCountingHighPurBJetTags',
                                                    'trackCountingHighEffBJetTags',
                                                    'simpleSecondaryVertexBJetTags',
                                                   #'simpleSecondaryVertexNegativeBJetTags',
                                                    'combinedSecondaryVertexBJetTags',
                                                    'combinedSecondaryVertexMVABJetTags',
                                                    'softElectronByPtBJetTags',
                                                    'softElectronByIP3dBJetTags',
                                                    'softMuonBJetTags',
                                                    'softMuonByPtBJetTags',
                                                    'softMuonByIP3dBJetTags'
                                                    )
                              ]
               }

    ## extend an existing sequence by otherLabels
    def mkseq(process, firstlabel, *otherlabels):
       seq = getattr(process, firstlabel)
       for x in otherlabels: seq += getattr(process, x)
       return cms.Sequence(seq)

    ## add tag infos to the process
    setattr( process, 'btaggingTagInfos'+label, mkseq(process, *(labels['tagInfos']) ) )
    ## add b tags to the process
    setattr( process, 'btaggingJetTags'+label,  mkseq(process, *(labels['jetTags'])  ) )
    ## add a combined sequence to the process
    seq = mkseq(process, 'btaggingTagInfos'+label, 'btaggingJetTags' + label) 
    setattr( process, 'btagging'+label, seq )
    ## return the combined sequence and the labels defined above
    return (seq, labels)


def switchJetCollection(process,
                        jetCollection,                       
                        doJTA            = True,
                        doBTagging       = True,
                        jetCorrLabel     = None,
                        doType1MET       = True,
                        genJetCollection = cms.InputTag("ak5GenJets"),
                        doJetID          = True,
                        jetIdLabel       = "ak5",
                        l1jetCollection  = cms.InputTag("allLayer1Jets")
                        ):
    """
    ------------------------------------------------------------------        
    switch the collection of jets in PAT from the default value to a
    new jet collection

    process          : process
    jetCollection    : input jet collection
    doBTagging       : run b tagging sequence for new jet collection
                       and add it to the new pat jet collection
    doJTA            : run JetTracksAssociation and JetCharge and add
                       it to the new pat jet collection (will autom.
                       be true if doBTagging is set to true)
    jetCorrLabel     : algorithm and type of JEC; use 'None' for no
                       JEC; examples are ('AK5','Calo'), ('SC7',
                       'Calo'), ('KT4','PF')
    doType1MET       : if jetCorrLabel is not 'None', set this to
                       'True' to redo the Type1 MET correction for
                       the new jet colllection; at the moment it must
                       be 'False' for non CaloJets otherwise the
                       JetMET POG module crashes.
    doJetID          : run jet id for new jet collection and add it
                       to the new pat jet collection
    genJetCollection : GenJet collection to match to

    doJetId          : add jetId variables to the added jet collection?

    jetIsLabel       : specify the label prefix of the xxxJetID object;
                       in general it is the jet collection tag like ak5,
                       kt4 sc5, aso. For more informatrino have a look
                       to SWGuidePATTools#switch_JetCollection
                       
    ------------------------------------------------------------------        
    """
    ## save label of old jet collection
#    oldLabel = process.allLayer1Jets.jetSource;

    l1jets      = getattr(process,l1jetCollection.moduleLabel)
    oldLabel    = l1jets.jetSource 
    mcjetpar    = getattr(process,l1jets.genPartonMatch.moduleLabel)
    mcjetjet    = getattr(process,l1jets.genJetMatch.moduleLabel)
    mcjetflass  = getattr(process,l1jets.JetPartonMapSource.moduleLabel)
    mcjetparass = getattr(process,mcjetflass.srcByReference.moduleLabel)
    jetcorr     = getattr(process,l1jets.jetCorrFactorsSource[0].moduleLabel)
    ## replace input jet collection for generator matches
    mcjetpar.src                = jetCollection
    mcjetjet.src                = jetCollection
    mcjetjet.matched            = genJetCollection
    mcjetparass.jets            = jetCollection

    
    ## replace input jet collection for trigger matches
    ##massSearchReplaceParam(process.patTrigMatch, 'src', oldLabel, jetCollection)

    ## replace input jet collection for pat jet production
    l1jets.jetSource = jetCollection
    
    ## make VInputTag from strings
    def vit(*args) : return cms.VInputTag( *[ cms.InputTag(x) for x in args ] )

    if (doJTA or doBTagging):
 
        ## replace jet track association
        process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
        from RecoJets.JetAssociationProducers.ak5JTA_cff import ak5JetTracksAssociatorAtVertex
        if (l1jetCollection.moduleLabel=="allLayer1Jets"):
            process.jetTracksAssociatorAtVertex = ak5JetTracksAssociatorAtVertex.clone(jets = jetCollection)
            process.makeAllLayer1Jets.replace(process.patJetCharge, process.jetTracksAssociatorAtVertex+process.patJetCharge)
            process.patJetCharge.src = 'jetTracksAssociatorAtVertex'
            l1jets.trackAssociationSource = 'jetTracksAssociatorAtVertex'
        if (l1jetCollection.moduleLabel=="pfLayer1Jets"):
            jetvtx = getattr(process,l1jets.trackAssociationSource.moduleLabel)
            jetvtx.jets = jetCollection
    else:
        ## remove the jet track association from the std
        ## sequence
        process.makeAllLayer1Jets.remove(process.patJetCharge)
        ## switch embedding of track association and jet
        ## charge estimate to 'False'
        l1jets.addAssociatedTracks = False
        l1jets.addJetCharge = False

    if (doBTagging):
        ## replace b tagging sequence; add postfix label 'AOD' as crab will
        ## crash when confronted with empy labels

        if (l1jetCollection.moduleLabel=="allLayer1Jets"):
            (btagSeq, btagLabels) = runBTagging(process, jetCollection, 'AOD')
            process.makeAllLayer1Jets.replace(process.jetTracksAssociatorAtVertex, process.jetTracksAssociatorAtVertex+btagSeq)
        if (l1jetCollection.moduleLabel=="pfLayer1Jets"):
            (btagSeq, btagLabels) = runBTagging(process, jetCollection, 'PF')
            process.makepflayer1Jets.replace(process.jetTracksAssociatorAtVertexPF,process.jetTracksAssociatorAtVertexPF+btagSeq)
        ## add b tagging sequence before running the allLayer1Jets modules
        

        ## replace corresponding tags for pat jet production
        l1jets.trackAssociationSource = btagLabels['jta']
        l1jets.tagInfoSources = cms.VInputTag( *[ cms.InputTag(x) for x in btagLabels['tagInfos'] ] )
        l1jets.discriminatorSources = cms.VInputTag( *[ cms.InputTag(x) for x in btagLabels['jetTags']  ] )
    else:
        ## remove b tagging from the std sequence
        l1jets.remove(process.secondaryVertexNegativeTagInfos)
        l1jets.remove(process.simpleSecondaryVertexNegativeBJetTags)
        ## switch embedding of b tagging for pat
        ## jet production to 'False'
        l1jets.addBTagInfo = False

    if (doJetID):
        jetIdLabelNew = jetIdLabel + 'JetID'
        l1jets.jetIDMap = cms.InputTag( jetIdLabelNew )
    else:
        l1jets.addJetID = cms.bool(False)


    if (jetCorrLabel!=None):
        ## replace jet energy corrections; catch
        ## a couple of exceptions first
        if (jetCorrLabel == False ):
            raise ValueError, "In switchJetCollection 'jetCorrLabel' must be set to 'None', not 'False'"
        if (jetCorrLabel == "None"):
            raise ValueError, "In switchJetCollection 'jetCorrLabel' must be set to 'None' (without quotes)"
        ## check for the correct format
        if (type(jetCorrLabel)!=type(('AK5','Calo'))): 
            raise ValueError, "In switchJetCollection 'jetCorrLabel' must be 'None', or of type ('Algo','Type')"

        ## switch JEC parameters to the new jet collection
        jetcorr.jetSource = jetCollection            
        switchJECParameters(jetcorr, jetCorrLabel[0], jetCorrLabel[1], oldAlgo='AK5',oldType='Calo')

        ## redo the type1MET correction for the new jet collection
        if (doType1MET):
            ## in case there is no jet correction service in the paths add it
            ## as L2L3 if possible, as combined from L2 and L3 otherwise
            if (not hasattr( process, 'L2L3JetCorrector%s%s' % jetCorrLabel )):
                setattr( process, 
                         'L2L3JetCorrector%s%s' % jetCorrLabel, 
                         cms.ESSource("JetCorrectionServiceChain",
                                      correctors = cms.vstring('L2RelativeJetCorrector%s%s' % jetCorrLabel,
                                                               'L3AbsoluteJetCorrector%s%s' % jetCorrLabel),
                                      label = cms.string('L2L3JetCorrector%s%s' % jetCorrLabel)
                                      )
                         )
            ## configure the type1MET correction the following muonMET
            ## corrections have the corMetType1Icone5 as input and are
            ## automatically correct  
            process.metJESCorAK5CaloJet.inputUncorJetsLabel = jetCollection.value()
            process.metJESCorAK5CaloJet.corrector = 'L2L3JetCorrector%s%s' % jetCorrLabel
    else:
        ## remove the jetCorrFactors from the std sequence
        process.patJetMETCorrections.remove(process.jetCorrFactors)
        ## switch embedding of jetCorrFactors off
        ## for pat jet production
        l1jets.addJetCorrFactors = False
    

def addJetCollection(process,
                     jetCollection,
                     postfixLabel,
                     doJTA        = True,
                     doBTagging   = True,
                     jetCorrLabel = None,
                     doType1MET   = True,
                     doL1Cleaning = True,                     
                     doL1Counters = False,
                     genJetCollection=cms.InputTag("ak5GenJets"),
                     doJetID      = True,
                     jetIdLabel   = "ak5"
                     ):
    """
    ------------------------------------------------------------------        
    add a new collection of jets in PAT

    process          : process
    jetCollection    : input jet collection    
    postfixLabel     : label to identify all modules that work with
                       this jet collection
    doBTagging       : run b tagging sequence for new jet collection
                       and add it to the new pat jet collection
    doJTA            : run JetTracksAssociation and JetCharge and add
                       it to the new pat jet collection (will autom.
                       be true if doBTagging is set to true)
    jetCorrLabel     : algorithm and type of JEC; use 'None' for no
                       JEC; examples are ('AK5','Calo'), ('SC7',
                       'Calo'), ('KT4','PF')
    doType1MET       : make also a new MET collection (not yet
                       implemented?)
    doL1Cleaning     : copy also the producer modules for cleanLayer1
                       will be set to 'True' automatically when
                       doL1Counters is 'True'
    doL1Counters     : copy also the filter modules that accept/reject
                       the event looking at the number of jets                       
    genJetCollection : GenJet collection to match to

    doJetId          : add jetId variables to the added jet collection?

    jetIsLabel       : specify the label prefix of the xxxJetID object;
                       in general it is the jet collection tag like ak5,
                       kt4 sc5, aso. For more informatrino have a look
                       to SWGuidePATTools#add_JetCollection

    this takes the configuration from the already-configured jets as
    starting point; replaces before calling addJetCollection will
    affect also the new jets
    ------------------------------------------------------------------                     
    """    
    ## add module as process to the makeAllLayer1Jets sequence
    def addAlso(label, value):
        existing = getattr(process, label)
        setattr( process, label+postfixLabel, value)
        process.patDefaultSequence.replace( existing, existing*value )        

    ## clone and add a module as process to the
    ## default sequence
    def addClone(label, **replaceStatements):
        new = getattr(process, label).clone(**replaceStatements)
        addAlso(label, new)
        
    ## add a clone of allLayer1Jets
    addClone('allLayer1Jets', jetSource = jetCollection)
    ## add a clone of selectedLayer1Jets    
    addClone('selectedLayer1Jets', src=cms.InputTag('allLayer1Jets'+postfixLabel))
    ## add a clone of cleanLayer1Jets    
    if (doL1Cleaning or doL1Counters):
        addClone('cleanLayer1Jets', src=cms.InputTag('selectedLayer1Jets'+postfixLabel))
    ## add a clone of countLayer1Jets    
    if (doL1Counters):
        addClone('countLayer1Jets', src=cms.InputTag('cleanLayer1Jets'+postfixLabel))

    ## attributes of allLayer1Jets
    l1Jets = getattr(process, 'allLayer1Jets'+postfixLabel)

    ## add a clone of gen jet matching
    addClone('jetPartonMatch', src = jetCollection)
    addClone('jetGenJetMatch', src = jetCollection, matched = genJetCollection)
    ## add a clone of parton and flavour associations
    addClone('jetPartonAssociation', jets = jetCollection)
    addClone('jetFlavourAssociation', srcByReference = cms.InputTag('jetPartonAssociation'+postfixLabel))

    ## fix label for input tag
    def fixInputTag(x): x.setModuleLabel(x.moduleLabel+postfixLabel)
    ## fix label for vector of input tags
    def fixVInputTag(x): x[0].setModuleLabel(x[0].moduleLabel+postfixLabel)

    ## provide allLayer1Jet inputs with individual labels
    fixInputTag(l1Jets.genJetMatch)
    fixInputTag(l1Jets.genPartonMatch)
    fixInputTag(l1Jets.JetPartonMapSource)

    ## make VInputTag from strings 
    def vit(*args) : return cms.VInputTag( *[ cms.InputTag(x) for x in args ] )

    if (doJTA or doBTagging):
        ## add clone of jet track association        
        process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
        from RecoJets.JetAssociationProducers.ak5JTA_cff import ak5JetTracksAssociatorAtVertex
        ## add jet track association module to processes
        jtaLabel = 'jetTracksAssociatorAtVertex'+postfixLabel
        setattr( process, jtaLabel, ak5JetTracksAssociatorAtVertex.clone(jets = jetCollection) )
        process.makeAllLayer1Jets.replace(process.patJetCharge, getattr(process,jtaLabel)+process.patJetCharge)
        l1Jets.trackAssociationSource = cms.InputTag(jtaLabel)
        addClone('patJetCharge', src=cms.InputTag(jtaLabel)),
        fixInputTag(l1Jets.jetChargeSource)
    else:
        ## switch embedding of track association and jet
        ## charge estimate to 'False'        
        l1Jets.addAssociatedTracks = False
        l1Jets.addJetCharge = False
    
    if (doBTagging):
        ## add b tagging sequence
        (btagSeq, btagLabels) = runBTagging(process, jetCollection, postfixLabel)
        ## add b tagging sequence before running the allLayer1Jets modules
        process.makeAllLayer1Jets.replace(getattr(process,jtaLabel), getattr(process,jtaLabel)+btagSeq)
        ## replace corresponding tags for pat jet production
        l1Jets.trackAssociationSource = cms.InputTag(btagLabels['jta'])
        l1Jets.tagInfoSources = cms.VInputTag( *[ cms.InputTag(x) for x in btagLabels['tagInfos'] ] )
        l1Jets.discriminatorSources = cms.VInputTag( *[ cms.InputTag(x) for x in btagLabels['jetTags']  ] )
    else:
        ## switch general b tagging info switch off
        l1Jets.addBTagInfo = False
        
    if (doJetID):
        l1Jets.addJetID = cms.bool(True)
        jetIdLabelNew = jetIdLabel + 'JetID'
        l1Jets.jetIDMap = cms.InputTag( jetIdLabelNew )
    else :
        l1Jets.addJetID = cms.bool(False)

    if (jetCorrLabel != None):
        ## add clone of jet energy corrections;
        ## catch a couple of exceptions first
        if (jetCorrLabel == False ):
            raise ValueError, "In addJetCollection 'jetCorrLabel' must be set to 'None', not 'False'"
        if (jetCorrLabel == "None"):
            raise ValueError, "In addJetCollection 'jetCorrLabel' must be set to 'None' (without quotes)"
        ## check for the correct format
        if type(jetCorrLabel) != type(('AK5','Calo')): 
            raise ValueError, "In switchJetCollection 'jetCorrLabel' must be 'None', or of type ('Algo','Type')"

        ## add clone of jetCorrFactors
        addClone('jetCorrFactors', jetSource = jetCollection)
        switchJECParameters( getattr(process,'jetCorrFactors'+postfixLabel), jetCorrLabel[0], jetCorrLabel[1], oldAlgo='AK5',oldType='Calo' )
        fixVInputTag(l1Jets.jetCorrFactorsSource)

        ## add a clone of the type1MET correction for the new jet collection
        if (doType1MET):
            ## in case there is no jet correction service in the paths add it
            ## as L2L3 if possible, as combined from L2 and L3 otherwise
            if not hasattr( process, 'L2L3JetCorrector%s%s' % jetCorrLabel ):
                setattr( process, 
                         'L2L3JetCorrector%s%s' % jetCorrLabel, 
                         cms.ESSource("JetCorrectionServiceChain",
                                      correctors = cms.vstring('L2RelativeJetCorrector%s%s' % jetCorrLabel,
                                                               'L3AbsoluteJetCorrector%s%s' % jetCorrLabel),
                                      label= cms.string('L2L3JetCorrector%s%s' % jetCorrLabel)
                                      )
                         )
            ## add a clone of the type1MET correction
            ## and the following muonMET correction  
            addClone('metJESCorAK5CaloJet', inputUncorJetsLabel = jetCollection.value(),
                     corrector = cms.string('L2L3JetCorrector%s%s' % jetCorrLabel)
                     )
            addClone('metJESCorAK5CaloJetMuons', uncorMETInputTag = cms.InputTag("metJESCorAK5CaloJet"+postfixLabel))
            addClone('layer1METs', metSource = cms.InputTag("metJESCorAK5CaloJetMuons"+postfixLabel))
            l1MET = getattr(process, 'layer1METs'+postfixLabel)

            ## add new met collections output to the pat summary
            process.allLayer1Summary.candidates += [ cms.InputTag('layer1METs'+postfixLabel) ]
    else:
        ## switch jetCorrFactors off
        l1Jets.addJetCorrFactors = False




def addJetID(process,
             jetSrc,
             jetIdTag
             ):
    """
    ------------------------------------------------------------------
    compute jet id for process

    process : process
    jetIdTag: Tag to append to jet id map
    ------------------------------------------------------------------    
    """
    jetIdLabel = jetIdTag + 'JetID'
    print "Making new jet ID label with label " + jetIdTag

    ## replace jet id sequence
    process.load("RecoJets.JetProducers.ak5JetID_cfi")
    setattr( process, jetIdLabel, process.ak5JetID.clone(src = jetSrc))
    process.makeAllLayer1Jets.replace( process.jetPartonMatch, getattr(process,jetIdLabel) + process.jetPartonMatch )




def setTagInfos(process,
                coll = "allLayer1Jets",
                tagInfos = cms.vstring( )
                ):
    """
    ------------------------------------------------------------------    
    replace tag infos for collection jetSrc

    process       : process
    jetCollection : jet collection to set tag infos for
    tagInfos      : tag infos to set
    ------------------------------------------------------------------    
    """
    found = False
    newTags = cms.VInputTag()
    iNewTags = 0
    for k in tagInfos :
        for j in getattr( process, coll ).tagInfoSources :
            vv = j.value();
            if ( vv.find(k) != -1 ):
                found = True
                newTags.append( j )
                
    if not found:
        raise RuntimeError,"""
        Cannot replace tag infos in jet collection""" % (coll)
    else :
        getattr(process,coll).tagInfoSources = newTags
