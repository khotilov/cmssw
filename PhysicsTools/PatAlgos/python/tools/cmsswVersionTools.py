import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.jetTools import *

def run33xOnReRecoMC( process,
                      genJets = "ak5GenJets"):
    """
    ------------------------------------------------------------------
    running GenJets for ak5 and ak7

    process : process
    genJets : which gen jets to run
    ------------------------------------------------------------------    
    """
    print "*********************************************************************"
    print "NOTE TO USER: when running on 31X samples re-recoed in 3.3.2         "
    print "              with this CMSSW version of PAT                         "
    print "              it is required to re-run the GenJet production for     "
    print "              anti-kT since that is not part of the re-reco          "
    print "*********************************************************************"
    process.load("RecoJets.Configuration.GenJetParticles_cff")
    process.load("RecoJets.JetProducers." + genJets +"_cfi")
    process.makePatJets.replace( process.patJetCharge, process.genParticlesForJets+getattr(process,genJets)+process.patJetCharge)

def run33xOn31xMC(process,
                  jetSrc = cms.InputTag("antikt5CaloJets"),
                  jetIdTag = "antikt5" ):
    """
    ------------------------------------------------------------------
    switch appropriate jet collections to run 33x on 31x MC

    process : process
    jetSrc  : jet source to use
    jetID   : jet ID to make
    ------------------------------------------------------------------    
    """
    print "*********************************************************************"
    print "NOTE TO USER: when running on 31X samples with this CMSSW version    "
    print "              of PAT the collection label for anti-kt has to be      "
    print "              switched from \'ak*\' to \'antikt*\'. This is going    "
    print "              to be done now. Also note that the *JetId collections  "
    print "              are not stored on these input files in contrary to     "
    print "              input files in 33X. Please use the _addJetId_ tool     "
    print "              as described on SWGuidePATTools, when adding new jet   "
    print "              collections! Such a line could look like this:         "
    print ""
    print "  addJetID( process, \"sisCone5CaloJets\", \"sc5\")"
    print "  from PhysicsTools.PatAlgos.tools.jetTools import *"
    print "  addJetCollection(process,cms.InputTag('sisCone5CaloJets'),"
    print "  ..."
    print "  )"
    print "*********************************************************************"
    addJetID( process, jetSrc, jetIdTag )
    # in PAT (iterativeCone5) to ak5 (anti-kt cone = 0.5)
    switchJetCollection(process, 
                        cms.InputTag('antikt5CaloJets'),   
                        doJTA            = True,            
                        doBTagging       = True,            
                        jetCorrLabel     = ('AK5','Calo'),  
                        doType1MET       = True,
                        genJetCollection = cms.InputTag("antikt5GenJets"),
                        doJetID          = True,
                        jetIdLabel       = "antikt5"
                        )
    

def restrictInputToAOD31X(process):
    """
    ------------------------------------------------------------------
    restrict input for pat tuple production to AOD 31x. Here the jet
    ID needs to be switched to false for the jets, as the information
    to produce it is not available on the 31X ADO samples.

    process : process
    ------------------------------------------------------------------    
    """
    jetProducer = getattr(process, 'patJets')
    jetProducer.addJetID = False



## ------------------------------------------------------
## CURRENTLY NOTHING IS TO BE DONE HERE
## IT REMAINS AS DUMMY FILE IN CASE OF
## LATER USE
## ------------------------------------------------------
