## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

## load photon sequencesup to selectedPatPhotons
process.load("PhysicsTools.PatAlgos.producersLayer1.photonProducer_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.photonSelector_cfi")

## make sure to keep the created objects
process.out.outputCommands = ['keep *_selectedPat*_*_*',]

## to run in scheduled mode uncomment the following lines
#process.p = cms.Path(
#    process.makePatPhotons *
#    process.selectedPatPhotons
#)

## to run in un-scheduled mode uncomment the following lines
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")
process.p = cms.Path(
    process.selectedPatPhotons
    )

## ------------------------------------------------------
#  In addition you usually want to change the following
#  parameters:
## ------------------------------------------------------
#
#   process.GlobalTag.globaltag =  ...    ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
#                                         ##
#   process.source.fileNames =  ...       ##  (e.g. 'file:AOD.root')
#                                         ##
process.maxEvents.input = 100
#                                         ##
#   process.out.outputCommands = [ ... ]  ##  (e.g. taken from PhysicsTools/PatAlgos/python/patEventContent_cff.py)
#                                         ##
process.out.fileName = 'patTuple_onlyPhotons.root'
#                                         ##
#   process.options.wantSummary = False   ##  (to suppress the long output at the end of the job)
