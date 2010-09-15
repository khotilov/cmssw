import FWCore.ParameterSet.Config as cms

def customise(process):
    REDIGIInputEventSkimming= cms.PSet(
        inputCommands=cms.untracked.vstring('drop *')
        )

    #HLTCleaning= cms.PSet(
    #    inputCommands=cms.untracked.vstring('drop FEDRawDataCollection_*_*_*')
    #    )

    #REDIGIInputEventSkimming.inputCommands.extend(process.RecoGenJetsFEVT.outputCommands)
    #REDIGIInputEventSkimming.inputCommands.extend(process.RecoGenMETFEVT.outputCommands)
    REDIGIInputEventSkimming.inputCommands.extend(process.SimG4CoreRAW.outputCommands) 
    REDIGIInputEventSkimming.inputCommands.extend(process.GeneratorInterfaceRAW.outputCommands) 
    REDIGIInputEventSkimming.inputCommands.extend(process.IOMCRAW.outputCommands) 

    #REDIGIInputEventSkimming.inputCommands.extend(process.HLTriggerRAW.outputCommands) 
    #REDIGIInputEventSkimming.inputCommands.extend(HLTCleaning.inputCommands)
    
    process.source.inputCommands = REDIGIInputEventSkimming.inputCommands
    process.source.dropDescendantsOfDroppedBranches=cms.untracked.bool(False)
    
    process.RandomNumberGeneratorService.restoreStateLabel = cms.untracked.string('randomEngineStateProducer')

    # Remove the old RNGState product on output
    RNGStateCleaning= cms.PSet(
        outputCommands=cms.untracked.vstring('drop RandomEngineStates_*_*_*',
                                             'keep RandomEngineStates_*_*_'+process.name_())
        )
    process.outputModules_().iteritems().next().outputCommands.extend(RNGStateCleaning.outputCommands)
    #process.output.outputCommands.extend(RNGStateCleaning.outputCommands)

    # REDO the GenJets etc. in case labels have been changed
    process.load('Configuration/StandardSequences/Generator_cff')
    process.fixGenInfo = cms.Path(process.genJetMET)
    process.schedule.append(process.fixGenInfo)
    
    # Output definition for RAW
    #process.outputRaw = cms.OutputModule("PoolOutputModule",
    #                                     outputCommands = process.RAWSIMEventContent.outputCommands,
    #                                     fileName = cms.untracked.string('New_RAWSIM.root'),
    #                                     dataset = cms.untracked.PSet(dataTier = cms.untracked.string(''),
    #                                                                  filterName = cms.untracked.string('')
    #                                                                  )
    #                                     )
    
    #process.out_step_raw = cms.EndPath(process.outputRaw)
    #process.schedule.append(process.out_step_raw)
                
    return(process)
