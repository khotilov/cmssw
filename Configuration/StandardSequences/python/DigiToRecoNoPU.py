import FWCore.ParameterSet.Config as cms

def customise(process):
    REDIGIInputEventSkimming= cms.PSet(
        inputCommands=cms.untracked.vstring('drop *')
        )

    REDIGIInputEventSkimming.inputCommands.extend(process.SimG4CoreRAW.outputCommands) 
    REDIGIInputEventSkimming.inputCommands.extend(process.GeneratorInterfaceRAW.outputCommands) 
    REDIGIInputEventSkimming.inputCommands.extend(process.IOMCRAW.outputCommands) 

    process.source.inputCommands = REDIGIInputEventSkimming.inputCommands
    
    if hasattr(process,"RandomNumberGeneratorService"):
        del process.RandomNumberGeneratorService.theSource
    else:    
        process.load("IOMC/RandomEngine/IOMC_cff")
        del process.RandomNumberGeneratorService.theSource

    process.RandomNumberGeneratorService.restoreStateLabel = cms.untracked.string('randomEngineStateProducer')

    # Output definition for RAW
    process.outputRaw = cms.OutputModule("PoolOutputModule",
       outputCommands = process.RAWSIMEventContent.outputCommands,
       fileName = cms.untracked.string('New_RAWSIM.root'),
       dataset = cms.untracked.PSet(
           dataTier = cms.untracked.string(''),
           filterName = cms.untracked.string('')
       )
    )

    process.out_step_raw = cms.EndPath(process.outputRaw)
    process.schedule.append(process.out_step_raw)
    
    return(process)
