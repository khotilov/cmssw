import FWCore.ParameterSet.Config as cms

process = cms.Process('makeZtoElecMuPlots')

process.load("TauAnalysis.Configuration.dumpZtoElecMu_cff")
process.load("TauAnalysis.Configuration.plotZtoElecMu_cff")

process.DQMStore = cms.Service("DQMStore")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

# import utility function to enable factorization
#from TauAnalysis.Configuration.factorizationTools import enableFactorization_makeZtoElecMuPlots
#enableFactorization_makeZtoElecMuPlots(process)

process.makeZtoElecMuPlots = cms.Sequence( process.loadZtoElecMu
                                          +process.addZtoElecMu
                                          +process.saveZtoElecMu
                                          +process.dumpZtoElecMu 
                                          +process.plotZtoElecMu )

process.p = cms.Path(process.makeZtoElecMuPlots)

# print-out all python configuration parameter information
#print process.dumpPython()
