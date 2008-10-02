#G.Benelli
#This fragment is used to add the SimpleMemoryCheck
#and Timing services output to the log of the simulation
#performance candles.
#It is meant to be used with the cmsDriver.py option
#--customise in the following fashion:
#E.g.
#./cmsDriver.py HZZLLLL -e 190 -n 50 --step=GEN --customise=Simulation.py >& HZZLLLL_190_GEN.log&
#or
#./cmsDriver.py MINBIAS -n 50 --step=GEN --customise=Simulation.py >& MINBIAS_GEN.log&

import FWCore.ParameterSet.Config as cms
def customise(process):
    #Renaming the process
    process.__dict__['_Process__name']='DIGIPILEUP'
    #Adding SimpleMemoryCheck service:
    process.SimpleMemoryCheck=cms.Service("SimpleMemoryCheck",
                                          ignoreTotal=cms.untracked.int32(1),
                                          oncePerEventMode=cms.untracked.bool(True))
    #Adding Timing service:
    process.Timing=cms.Service("Timing")

    #Add these 3 lines to put back the summary for timing information at the end of the logfile
    #(needed for TimeReport report)
    process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
        )
            
    #Overwriting the fileNames to be used by the MixingModule
    #when invoking cmsDriver.py with the --PU option
    process.mix.input.fileNames = cms.untracked.vstring('file:../MinBias_TimeSize/MINBIAS__GEN,SIM.root')

    return(process)
