import FWCore.ParameterSet.Config as cms

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *

#-----------------------------------------------------------------------------------------------
# produce collections of pat::Jets passing selection criteria, specific for W-->tau nu analysis
#---------------------------------------------------------------------------------------------

# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string parser

# Select jets for W->tau nu analysis
selectedLayer1JetsAntiOverlapWithTausVetoForWTauNu = cms.EDFilter("PATJetAntiOverlapSelector",
                                                                  srcNotToBeFiltered = cms.VInputTag("selectedLayer1TausForWTauNuEcalCrackVetoCumulative"),
                                                                  dRmin = cms.double(0.7),
                                                                  filter = cms.bool(False)
                                                                  )

selectedLayer1JetsEta21ForWTauNu = cms.EDFilter("PATJetSelector",
                                                cut = cms.string('abs(eta) < 2.1'),
                                                filter = cms.bool(False)
                                                )

selectedLayer1JetsEt15ForWTauNu = cms.EDFilter("PATJetSelector",
                                               cut = cms.string('et > 15.'),
                                               filter = cms.bool(False)
                                               )

selectedLayer1JetsEt20ForWTauNu = cms.EDFilter("PATJetSelector",
                                               cut = cms.string('et > 20.'),
                                               filter = cms.bool(False)
                                               )

patJetSelConfiguratorForWTauNu = objSelConfigurator(
    [ selectedLayer1JetsAntiOverlapWithTausVetoForWTauNu,
      selectedLayer1JetsEta21ForWTauNu,
      selectedLayer1JetsEt15ForWTauNu,
      selectedLayer1JetsEt20ForWTauNu ],
    src = "cleanPatJets",
    pyModuleName = __name__,
    doSelIndividual = False
    )
selectLayer1JetsForWTauNu =patJetSelConfiguratorForWTauNu.configure(pyNameSpace = locals())

#select jets for W->taunu analysis, loose isolation
selectedLayer1JetsAntiOverlapWithTausVetoForWTauNuLooseIsolation = cms.EDFilter("PATJetAntiOverlapSelector",
                                                                                srcNotToBeFiltered = cms.VInputTag("selectedLayer1TausForWTauNuEcalCrackVetoLooseIsolationCumulative"),
                                                                                dRmin = cms.double(0.7),
                                                                                filter = cms.bool(False)
                                                                                )
selectedLayer1JetsEta21ForWTauNuLooseIsolation = cms.EDFilter("PATJetSelector",
                                                              cut = cms.string('abs(eta) < 2.1'),
                                                              filter = cms.bool(False)
                                                              )

selectedLayer1JetsEt15ForWTauNuLooseIsolation = cms.EDFilter("PATJetSelector",
                                                             cut = cms.string('et > 15.'),
                                                             filter = cms.bool(False)
                                                             )

selectedLayer1JetsEt20ForWTauNuLooseIsolation = cms.EDFilter("PATJetSelector",
                                                             cut = cms.string('et > 20.'),
                                                             filter = cms.bool(False)
                                                             )

patJetSelConfiguratorForWTauNuLooseIsolation = objSelConfigurator(
    [ selectedLayer1JetsAntiOverlapWithTausVetoForWTauNuLooseIsolation,
      selectedLayer1JetsEta21ForWTauNuLooseIsolation,
      selectedLayer1JetsEt15ForWTauNuLooseIsolation,
      selectedLayer1JetsEt20ForWTauNuLooseIsolation ],
    src = "cleanPatJets",
    pyModuleName = __name__,
    doSelIndividual = False
    )
selectLayer1JetsForWTauNuLooseIsolation =patJetSelConfiguratorForWTauNuLooseIsolation.configure(pyNameSpace = locals())

selectLayer1SelJetsForWTauNu = cms.Sequence (selectLayer1JetsForWTauNu
                                             *selectLayer1JetsForWTauNuLooseIsolation) 
