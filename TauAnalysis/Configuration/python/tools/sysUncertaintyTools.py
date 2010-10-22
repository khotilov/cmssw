import FWCore.ParameterSet.Config as cms
import re

from PhysicsTools.PatAlgos.tools.helpers import listModules

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName
from TauAnalysis.Configuration.tools.analysisSequenceTools import removeAnalyzer

from TauAnalysis.RecoTools.patLeptonSelection_cff import *
from TauAnalysis.RecoTools.patLeptonSystematics_cff import *
from TauAnalysis.CandidateTools.muTauPairProduction_cff import *
from TauAnalysis.CandidateTools.diTauPairSelectionAllKinds_cff import *
from TauAnalysis.CandidateTools.muTauPairSelectionForAHtoMuTau_cff import *
from TauAnalysis.RecoTools.patJetSelection_cff import *
from TauAnalysis.RecoTools.patJetSelectionForAHtoMuTau_cff import *
from TauAnalysis.CandidateTools.sysErrDefinitions_cfi import *
from TauAnalysis.GenSimTools.sysErrGenEventReweights_cfi import *

#--------------------------------------------------------------------------------
# generic utility functions for factorization
# usable for all channels
#--------------------------------------------------------------------------------

def removeModules(process, sequenceName, moduleNamePattern, pyNameSpace):
    # auxiliary function to remove from sequence
    # all modules matching pattern of module names passed as function argument

    moduleNamePattern_regexp = re.compile(moduleNamePattern)

    if hasattr(process, sequenceName):
        sequence = getattr(process, sequenceName)

        modules = listModules(sequence)
        
        for module in modules:
            if hasattr(module, "label"):
                moduleName = module.label()

                if moduleNamePattern_regexp.match(moduleName):
                    sequence.remove(module)

def addBoolEventSelFlagProducer(process, moduleName, expSysUncertainties, sequenceName):

    if hasattr(process, moduleName) and hasattr(process, sequenceName):

        module = getattr(process, moduleName)
        sequence = getattr(process, sequenceName)

        for expSysUncertainty in expSysUncertainties:
            module_systematic = copy.deepcopy(module)

            flags = getattr(module, "flags")
            flags_systematic = []

            for flag in flags:
                                
                flagModuleLabel_systematic = composeModuleName(flag.getModuleLabel(), expSysUncertainty)
                
                if hasattr(process, flagModuleLabel_systematic):
                    flag_systematic = copy.deepcopy(flag)
                    flag_systematic.setModuleLabel(flagModuleLabel_systematic)
                    flags_systematic.append(flag_systematic)
                else:
                    flags_systematic.append(flag)

            setattr(module_systematic, "flags", cms.VInputTag(flags_systematic))

            moduleName_systematic = composeModuleName(moduleName, expSysUncertainty)
            setattr(process, moduleName_systematic, module_systematic)

            module_systematic = getattr(process, moduleName_systematic)

            sequence._seq = sequence._seq * module_systematic
            
#--------------------------------------------------------------------------------
# functions to enable/disable estimation of systematic uncertainties
# in PAT-tuple production
#--------------------------------------------------------------------------------

def disableSysUncertainties_patTupleProduction(process):
    #print("<disableSysUncertainties_patTupleProduction>:")
    
    process.produceGenObjects.remove(process.produceSysErrGenEventReweights)

#--------------------------------------------------------------------------------
# functions to enable/disable estimation of systematic uncertainties
# specific to Z --> muon + tau-jet channel
#--------------------------------------------------------------------------------

def disableSysUncertainties_runZtoMuTau(process):
    print("<disableSysUncertainties_runZtoMuTau>:")
    print("--> **disabling** estimation of systematic uncertainties...")
    
    moduleNamePattern = "\w+Sys\w+(Up|Down)"
    pyNameSpace = None

    process.produceGenObjects.remove(process.produceSysErrGenEventReweights)
    process.producePatTupleZtoMuTauSpecific.remove(process.prodSmearedMuons)
    process.producePatTupleZtoMuTauSpecific.remove(process.prodSmearedTaus)

    if hasattr(process, "selectZtoMuTauEvents"):
        removeModules(process, "selectZtoMuTauEvents", moduleNamePattern, pyNameSpace)
    if hasattr(process, "selectZtoMuTauEventsLooseMuonIsolation"):
        removeModules(process, "selectZtoMuTauEventsLooseMuonIsolation", moduleNamePattern, pyNameSpace)

    if hasattr(process, "analyzeZtoMuTauEvents"):
        removeAnalyzer(process.analyzeZtoMuTauEvents.analysisSequence, "sysUncertaintyBinnerForMuTauAcc")
        removeAnalyzer(process.analyzeZtoMuTauEvents.analysisSequence, "sysUncertaintyBinnerForMuTauEff")
        process.analyzeZtoMuTauEvents.estimateSysUncertainties = cms.bool(False)
    if hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation"):
        removeAnalyzer(process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation.analysisSequence, "sysUncertaintyBinnerForMuTauAcc")
        removeAnalyzer(process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation.analysisSequence, "sysUncertaintyBinnerForMuTauEff")
        process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation.estimateSysUncertainties = cms.bool(False)
    if hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation"):
        removeAnalyzer(process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation.analysisSequence, "sysUncertaintyBinnerForMuTauAcc")
        removeAnalyzer(process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation.analysisSequence, "sysUncertaintyBinnerForMuTauEff")
        process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation.estimateSysUncertainties = cms.bool(False)    

def enableSysUncertainties_runZtoMuTau(process):
    print("<enableSysUncertainties_runZtoMuTau>:")
    print("--> **enabling** estimation of systematic uncertainties...")
    
    pyNameSpace = None

    setattr(patMuonSelConfigurator, "systematics", muonSystematics)
    process.selectPatMuons = patMuonSelConfigurator.configure(process = process)

    setattr(patMuonSelConfiguratorLooseIsolation, "systematics", muonSystematics)
    process.selectPatMuonsLooseIsolation = patMuonSelConfiguratorLooseIsolation.configure(process = process)

    setattr(patTauSelConfigurator, "systematics", tauSystematics)
    process.selectPatTaus = patTauSelConfigurator.configure(process = process)

    setattr(patTauSelConfiguratorForMuTau, "systematics", tauSystematics)
    process.selectPatTausForMuTau = patTauSelConfiguratorForMuTau.configure(process = process)

    setattr(patJetSelConfigurator, "systematics", jetSystematics)
    process.selectPatJets = patJetSelConfigurator.configure(process = process)

    process.smearedMET = cms.EDProducer("SmearedMETProducer",
        src = cms.InputTag('patMETs'),
        smearedParticles = cms.PSet()
    )
    smearedMETconfigurator = objProdConfigurator(
        process.smearedMET,
        pyModuleName = __name__,
        systematics = {
            "sysMuonPtUp" : {        
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatMuonsTrkIPcumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatMuonsTrkIPsysMuonPtUpCumulative')
            },
            "sysMuonPtDown" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatMuonsTrkIPcumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatMuonsTrkIPsysMuonPtDownCumulative')
            },
            "sysTauJetEnUp" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatTausForMuTauMuonVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetEnUpCumulative')
            },
            "sysTauJetEnDown" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatTausForMuTauMuonVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetEnDownCumulative')
            }, 
            "sysTauJetThetaUp" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatTausForMuTauMuonVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetThetaUpCumulative')
            }, 
            "sysTauJetThetaDown" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatTausForMuTauMuonVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetThetaDownCumulative')
            }, 
            "sysTauJetPhiUp" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatTausForMuTauMuonVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetPhiUpCumulative')
            }, 
            "sysTauJetPhiDown" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatTausForMuTauMuonVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetPhiDownCumulative')
            },
            "sysJetEnUp" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatJetsAntiOverlapWithLeptonsVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatJetsAntiOverlapWithLeptonsVetoSysJetEnUpCumulative')
            },
            "sysJetEnDown" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatJetsAntiOverlapWithLeptonsVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatJetsAntiOverlapWithLeptonsVetoSysJetEnDownCumulative')
            }
        }
    )
    process.prodSmearedMET = smearedMETconfigurator.configure(process = process)

    setattr(muTauPairProdConfigurator, "systematics", {
        "sysMuonPtUp" : {
            "srcLeg1" : cms.InputTag('selectedPatMuonsTrkIPsysMuonPtUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysMuonPtUp')
        },
        "sysMuonPtDown" : {
            "srcLeg1" : cms.InputTag('selectedPatMuonsTrkIPsysMuonPtDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysMuonPtDown')           
        },
        "sysTauJetEnUp" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetEnUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetEnUp')
        },
        "sysTauJetEnDown" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetEnDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetEnDown')
        }, 
        "sysTauJetThetaUp" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetThetaUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetThetaUp')
        }, 
        "sysTauJetThetaDown" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetThetaDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetThetaDown')
        }, 
        "sysTauJetPhiUp" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetPhiUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetPhiUp')
        }, 
        "sysTauJetPhiDown" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetPhiDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetPhiDown')                        
        },
        "sysJetEnUp" : {
            "srcMET"  : cms.InputTag('smearedMETsysJetEnUp')
        },
        "sysJetEnDown" : {
            "srcMET"  : cms.InputTag('smearedMETsysJetEnDown')
        }
    })
    process.produceMuTauPairs = muTauPairProdConfigurator.configure(process = process)
    process.produceMuTauPairs._seq = process.prodSmearedMET * process.produceMuTauPairs._seq

    setattr(muTauPairProdConfiguratorLooseMuonIsolation, "systematics", {
        "sysMuonPtUp" : {
            "srcLeg1" : cms.InputTag('selectedPatMuonsTrkIPlooseIsolationSysMuonPtUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysMuonPtUp')
        },
        "sysMuonPtDown" : {
            "srcLeg1" : cms.InputTag('selectedPatMuonsTrkIPlooseIsolationSysMuonPtDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysMuonPtDown')           
        },
        "sysTauJetEnUp" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetEnUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetEnUp')
        },
        "sysTauJetEnDown" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetEnDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetEnDown')
        }, 
        "sysTauJetThetaUp" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetThetaUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetThetaUp')
        }, 
        "sysTauJetThetaDown" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetThetaDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetThetaDown')
        }, 
        "sysTauJetPhiUp" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetPhiUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetPhiUp')
        }, 
        "sysTauJetPhiDown" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauMuonVetoSysTauJetPhiDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetPhiDown')            
            
        },
        "sysJetEnUp" : {
            "srcMET"  : cms.InputTag('smearedMETsysJetEnUp')
        },
        "sysJetEnDown" : {
            "srcMET"  : cms.InputTag('smearedMETsysJetEnDown')
        }
    })
    process.produceMuTauPairsLooseMuonIsolation = muTauPairProdConfiguratorLooseMuonIsolation.configure(process = process)

    setattr(patMuTauPairSelConfigurator, "systematics", muTauPairSystematics)
    process.selectMuTauPairs = patMuTauPairSelConfigurator.configure(process = process)
    
    setattr(patMuTauPairSelConfiguratorLooseMuonIsolation, "systematics", muTauPairSystematics)
    process.selectMuTauPairsLooseMuonIsolation = patMuTauPairSelConfiguratorLooseMuonIsolation.configure(process = process)

    if hasattr(process, "isRecZtoMuTau"):        
        expSysUncertainties = getSysUncertaintyNames(
            [ muonSystematics,
              tauSystematics,
              jetSystematics ]
        )
        addBoolEventSelFlagProducer(process, "isRecZtoMuTau", expSysUncertainties, "selectZtoMuTauEvents")

    if hasattr(process, "analyzeZtoMuTauEvents"):
        process.analyzeZtoMuTauEvents.estimateSysUncertainties = cms.bool(True)
    if hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation"):
        process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation.estimateSysUncertainties = cms.bool(True)
    if hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation"):
        process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation.estimateSysUncertainties = cms.bool(True)    

#--------------------------------------------------------------------------------
# functions to enable/disable estimation of systematic uncertainties
# specific to Z --> electron + tau-jet channel
#--------------------------------------------------------------------------------

def disableSysUncertainties_runZtoElecTau(process):
    #print("<disableSysUncertainties_runZtoElecTau>:")
    
    moduleNamePattern = "\w+Sys\w+(Up|Down)"
    pyNameSpace = None

    process.produceGenObjects.remove(process.produceSysErrGenEventReweights)

    removeModules(process, "selectZtoElecTauEvents", moduleNamePattern, pyNameSpace)
    if hasattr(process, "selectZtoElecTauEventsLooseElectronIsolation"):
        removeModules(process, "selectZtoElecTauEventsLooseElectronIsolation", moduleNamePattern, pyNameSpace)

#--------------------------------------------------------------------------------
# functions to enable/disable estimation of systematic uncertainties
# specific to Z --> electron + muon channel
#--------------------------------------------------------------------------------

def disableSysUncertainties_runZtoElecMu(process):
    #print("<disableSysUncertainties_runZtoElecMu>:")
    
    moduleNamePattern = "\w+Sys\w+(Up|Down)"
    pyNameSpace = None

    process.produceGenObjects.remove(process.produceSysErrGenEventReweights)

    removeModules(process, "selectZtoElecMuEvents", moduleNamePattern, pyNameSpace)
    if hasattr(process, "selectZtoElecMuEventsLooseElectronIsolation"):
        removeModules(process, "selectZtoElecMuEventsLooseElectronIsolation", moduleNamePattern, pyNameSpace)

#--------------------------------------------------------------------------------
# functions to enable/disable estimation of systematic uncertainties
# specific to Z --> tau-jet + tau-jet channel
#--------------------------------------------------------------------------------

def disableSysUncertainties_runZtoDiTau(process):
    #print("<disableSysUncertainties_runZtoDiTau>:")
    
    moduleNamePattern = "\w+Sys\w+(Up|Down)"
    pyNameSpace = None

    process.produceGenObjects.remove(process.produceSysErrGenEventReweights)

    removeModules(process, "selectZtoDiTauEvents", moduleNamePattern, pyNameSpace)

#--------------------------------------------------------------------------------
# functions to enable/disable estimation of systematic uncertainties
# specific to W --> tau-jet + nu channel
#--------------------------------------------------------------------------------

def disableSysUncertainties_runWtoTauNu(process):
    #print("<disableSysUncertainties_runWtoTauNu>:")
    
    moduleNamePattern = "\w+Sys\w+(Up|Down)"
    pyNameSpace = None

    process.produceGenObjects.remove(process.produceSysErrGenEventReweights)

    removeModules(process, "selectWtoTauNuEvents", moduleNamePattern, pyNameSpace)
    if hasattr(process, "selectWtoTauNuEventsLooseIsolation"):
        removeModules(process, "selectWtoTauNuEventsLooseIsolation", moduleNamePattern, pyNameSpace)

#--------------------------------------------------------------------------------
# functions to enable/disable estimation of systematic uncertainties
# specific to MSSM Higgs --> muon + tau-jet channel
#--------------------------------------------------------------------------------

def disableSysUncertainties_runAHtoMuTau(process):
    print("<disableSysUncertainties_runAHtoMuTau>:")
    print("--> **disabling** estimation of systematic uncertainties...")

    disableSysUncertainties_runZtoMuTau(process)
    
    moduleNamePattern = "\w+Sys\w+(Up|Down)"
    pyNameSpace = None

    if hasattr(process, "selectAHtoMuTauEvents"):
        removeModules(process, "selectAHtoMuTauEvents", moduleNamePattern, pyNameSpace)
    if hasattr(process, "selectAHtoMuTauEventsLooseMuonIsolation"):
        removeModules(process, "selectAHtoMuTauEventsLooseMuonIsolation", moduleNamePattern, pyNameSpace)

    if hasattr(process, "analyzeAHtoMuTauEvents_woBtag"):
        removeAnalyzer(process.analyzeAHtoMuTauEvents_woBtag.analysisSequence,
                       "sysUncertaintyBinnerForMuTauAccCentralJetVeto")
        removeAnalyzer(process.analyzeAHtoMuTauEvents_woBtag.analysisSequence,
                       "sysUncertaintyBinnerForMuTauEff")
        removeAnalyzer(process.analyzeAHtoMuTauEvents_woBtag.analysisSequence,
                       "sysUncertaintyHistManagerForMuTau")
        process.analyzeAHtoMuTauEvents_woBtag.estimateSysUncertainties = cms.bool(False)
    if hasattr(process, "analyzeAHtoMuTauEvents_wBtag"):
        removeAnalyzer(process.analyzeAHtoMuTauEvents_wBtag.analysisSequence,
                       "sysUncertaintyBinnerForMuTauAccCentralJetBtag")
        removeAnalyzer(process.analyzeAHtoMuTauEvents_wBtag.analysisSequence,
                       "sysUncertaintyBinnerForMuTauEff")
        removeAnalyzer(process.analyzeAHtoMuTauEvents_wBtag.analysisSequence,
                       "sysUncertaintyHistManagerForMuTau")
        process.analyzeAHtoMuTauEvents_wBtag.estimateSysUncertainties = cms.bool(False)
    if hasattr(process, "analyzeAHtoMuTauEvents_woBtag_factorizedWithMuonIsolation"):
        removeAnalyzer(process.analyzeAHtoMuTauEvents_woBtag_factorizedWithMuonIsolation.analysisSequence,
                       "sysUncertaintyBinnerForMuTauAccCentralJetVeto")
        removeAnalyzer(process.analyzeAHtoMuTauEvents_woBtag_factorizedWithMuonIsolation.analysisSequence,
                       "sysUncertaintyBinnerForMuTauEff")
        removeAnalyzer(process.analyzeAHtoMuTauEvents_woBtag_factorizedWithMuonIsolation.analysisSequence,
                       "sysUncertaintyHistManagerForMuTau")                
        process.analyzeAHtoMuTauEvents_woBtag_factorizedWithMuonIsolation.estimateSysUncertainties = cms.bool(False)
    if hasattr(process, "analyzeAHtoMuTauEvents_wBtag_factorizedWithMuonIsolation"):
        removeAnalyzer(process.analyzeAHtoMuTauEvents_wBtag_factorizedWithMuonIsolation.analysisSequence,
                       "sysUncertaintyBinnerForMuTauAccCentralJetBtag")
        removeAnalyzer(process.analyzeAHtoMuTauEvents_wBtag_factorizedWithMuonIsolation.analysisSequence,
                       "sysUncertaintyBinnerForMuTauEff")
        removeAnalyzer(process.analyzeAHtoMuTauEvents_wBtag_factorizedWithMuonIsolation.analysisSequence,
                       "sysUncertaintyHistManagerForMuTau")                              
        process.analyzeAHtoMuTauEvents_wBtag_factorizedWithMuonIsolation.estimateSysUncertainties = cms.bool(False)    
    if hasattr(process, "analyzeAHtoMuTauEvents_woBtag_factorizedWithoutMuonIsolation"):
        removeAnalyzer(process.analyzeAHtoMuTauEvents_woBtag_factorizedWithoutMuonIsolation.analysisSequence,
                       "sysUncertaintyBinnerForMuTauAccCentralJetVeto")
        removeAnalyzer(process.analyzeAHtoMuTauEvents_woBtag_factorizedWithoutMuonIsolation.analysisSequence,
                       "sysUncertaintyBinnerForMuTauEff")
        removeAnalyzer(process.analyzeAHtoMuTauEvents_woBtag_factorizedWithoutMuonIsolation.analysisSequence,
                       "sysUncertaintyHistManagerForMuTau")                
        process.analyzeAHtoMuTauEvents_woBtag_factorizedWithoutMuonIsolation.estimateSysUncertainties = cms.bool(False)
    if hasattr(process, "analyzeAHtoMuTauEvents_wBtag_factorizedWithoutMuonIsolation"):
        removeAnalyzer(process.analyzeAHtoMuTauEvents_wBtag_factorizedWithoutMuonIsolation.analysisSequence,
                       "sysUncertaintyBinnerForMuTauAccCentralJetBtag")
        removeAnalyzer(process.analyzeAHtoMuTauEvents_wBtag_factorizedWithoutMuonIsolation.analysisSequence,
                       "sysUncertaintyBinnerForMuTauEff")
        removeAnalyzer(process.analyzeAHtoMuTauEvents_wBtag_factorizedWithoutMuonIsolation.analysisSequence,
                       "sysUncertaintyHistManagerForMuTau")               
        process.analyzeAHtoMuTauEvents_wBtag_factorizedWithoutMuonIsolation.estimateSysUncertainties = cms.bool(False)    

def enableSysUncertainties_runAHtoMuTau(process):
    print("<enableSysUncertainties_runAHtoMuTau>:")
    print("--> **enabling** estimation of systematic uncertainties...")

    enableSysUncertainties_runZtoMuTau(process)

    setattr(patMuTauPairSelConfiguratorForAHtoMuTau, "systematics", muTauPairSystematics)
    process.selectMuTauPairsForAHtoMuTau = patMuTauPairSelConfiguratorForAHtoMuTau.configure(process = process)
    
    setattr(patMuTauPairSelConfiguratorForAHtoMuTauLooseMuonIsolation, "systematics", muTauPairSystematics)
    process.selectMuTauPairsForAHtoMuTauLooseMuonIsolation = \
      patMuTauPairSelConfiguratorForAHtoMuTauLooseMuonIsolation.configure(process = process)

    setattr(patJetSelConfiguratorForAHtoMuTau, "systematics", jetSystematics)
    process.selectPatJetsForAHtoMuTau = patJetSelConfiguratorForAHtoMuTau.configure(process = process)

    setattr(patJetSelConfiguratorForAHtoMuTauLooseMuonIsolation, "systematics", jetSystematics)
    process.selectPatJetsForAHtoMuTauLooseMuonIsolation = \
      patJetSelConfiguratorForAHtoMuTauLooseMuonIsolation.configure(process = process)

    expSysUncertainties = getSysUncertaintyNames(
        [ muonSystematics,
          tauSystematics,
          jetSystematics ]
    )
    if hasattr(process, "isRecAHtoMuTauCentralJetVeto"):        
        addBoolEventSelFlagProducer(process, "isRecAHtoMuTauCentralJetVeto", expSysUncertainties, "selectAHtoMuTauEvents")
    if hasattr(process, "isRecAHtoMuTauCentralJetBtag"):        
        addBoolEventSelFlagProducer(process, "isRecAHtoMuTauCentralJetBtag", expSysUncertainties, "selectAHtoMuTauEvents")    

    if hasattr(process, "analyzeAHtoMuTauEvents_woBtag"):
        process.analyzeAHtoMuTauEvents_woBtag.estimateSysUncertainties = cms.bool(True)
    if hasattr(process, "analyzeAHtoMuTauEvents_wBtag"):
        process.analyzeAHtoMuTauEvents_wBtag.estimateSysUncertainties = cms.bool(True)    
    if hasattr(process, "analyzeAHtoMuTauEvents_woBtag_factorizedWithMuonIsolation"):
        process.analyzeAHtoMuTauEvents_woBtag_factorizedWithMuonIsolation.estimateSysUncertainties = cms.bool(True)
    if hasattr(process, "analyzeAHtoMuTauEvents_wBtag_factorizedWithMuonIsolation"):
        process.analyzeAHtoMuTauEvents_wBtag_factorizedWithMuonIsolation.estimateSysUncertainties = cms.bool(True)    
    if hasattr(process, "analyzeAHtoMuTauEvents_woBtag_factorizedWithoutMuonIsolation"):
        process.analyzeAHtoMuTauEvents_woBtag_factorizedWithoutMuonIsolation.estimateSysUncertainties = cms.bool(True)
    if hasattr(process, "analyzeAHtoMuTauEvents_wBtag_factorizedWithoutMuonIsolation"):
        process.analyzeAHtoMuTauEvents_wBtag_factorizedWithoutMuonIsolation.estimateSysUncertainties = cms.bool(True)    









