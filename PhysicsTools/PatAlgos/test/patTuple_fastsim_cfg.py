## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

## let it run
process.p = cms.Path(
        process.patDefaultSequence
    )



## ------------------------------------------------------
#  In addition you usually want to change the following
#  parameters:
## ------------------------------------------------------
process.source.fileNames    = [ '/store/relval/CMSSW_3_6_0_pre6/RelValTTbar/GEN-SIM-RECO/MC_36Y_V4-v1/0011/060537EE-4E45-DF11-AE40-003048679012.root']                                       
process.GlobalTag.globaltag = 'MC_3XY_V14::All'

#                                         ##
#   process.maxEvents.input = ...         ##  (e.g. -1 to run on all events)
#                                         ##
#   process.out.outputCommands = [ ... ]  ##  (e.g. taken from PhysicsTools/PatAlgos/python/patEventContent_cff.py)
#                                         ##
#   process.out.fileName = ...            ##  (e.g. 'myTuple.root')
#                                         ##
#   process.options.wantSummary = True    ##  (to suppress the long output at the end of the job) 
