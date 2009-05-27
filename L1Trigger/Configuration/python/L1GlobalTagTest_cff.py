# cff file defining sequences to print the L1 content of a global tag 
# for L1 sub-systems


import FWCore.ParameterSet.Config as cms

# Regional Calorimeter Trigger
#

l1RCTParametersTest = cms.EDAnalyzer('L1RCTParametersTester')
l1RCTChannelMaskTest = cms.EDAnalyzer('L1RCTChannelMaskTester')

printGlobalTagL1Rct = cms.Sequence(l1RCTParametersTest*l1RCTChannelMaskTest)

# Global Calorimeter Trigger
#

import L1TriggerConfig.GctConfigProducers.l1GctConfigDump_cfi
printGlobalTagL1Gct = L1TriggerConfig.GctConfigProducers.l1GctConfigDump_cfi.l1GctConfigDump.clone()

# DT TPG
#

# MISSING
#printGlobalTagL1DtTPG = cms.Sequence()

# DT TF
#

DTExtLutTester = cms.EDAnalyzer('DTExtLutTester')
DTPhiLutTester = cms.EDAnalyzer('DTPhiLutTester')
DTPtaLutTester = cms.EDAnalyzer('DTPtaLutTester')
DTEtaPatternLutTester = cms.EDAnalyzer('DTEtaPatternLutTester')
DTQualPatternLutTester = cms.EDAnalyzer('DTQualPatternLutTester')
DTTFParametersTester = cms.EDAnalyzer('DTTFParametersTester')
DTTFMasksTester = cms.EDAnalyzer('DTTFMasksTester')

printGlobalTagL1DtTF = cms.Sequence(DTExtLutTester
                                    *DTPhiLutTester
                                    *DTPtaLutTester
                                    *DTEtaPatternLutTester
                                    *DTQualPatternLutTester
                                    *DTTFParametersTester
                                    *DTTFMasksTester
                                    )

# CSC TF
#

# MISSING
#printGlobalTagL1CscTF = cms.Sequence()

# RPC TRigger
#
dumpL1RPCConfig = cms.EDAnalyzer('DumpL1RPCConfig',
          fileName = cms.string('PrintGlobalTag_L1RPCConfig.log'))
dumpConeDefinition = cms.EDAnalyzer('DumpConeDefinition',
          fileName = cms.string('PrintGlobalTag_L1RPCConeDefinition.log'))


printGlobalTagL1Rpc = cms.Sequence(dumpL1RPCConfig*dumpConeDefinition)

# Global Muon Trigger
#

printL1GmtParameters = cms.EDProducer('L1MuGlobalMuonTrigger',
    Debug = cms.untracked.int32(9),
    BX_min = cms.int32(-1),
    BX_max = cms.int32(1),
    BX_min_readout = cms.int32(-1),
    BX_max_readout = cms.int32(1),
    DTCandidates = cms.InputTag('none'),
    RPCbCandidates = cms.InputTag('none'),
    CSCCandidates = cms.InputTag('none'),
    RPCfCandidates = cms.InputTag('none'),
    MipIsoData = cms.InputTag('none'),
    WriteLUTsAndRegs = cms.untracked.bool(False)
)

printL1GmtMuScales = cms.EDAnalyzer('L1MuScalesTester')

printGlobalTagL1Gmt = cms.Sequence(printL1GmtParameters*printL1GmtMuScales)


# Global Trigger
#

l1GtStableParametersTest = cms.EDAnalyzer('L1GtStableParametersTester')
l1GtParametersTest = cms.EDAnalyzer('L1GtParametersTester')
l1GtBoardMapsTest = cms.EDAnalyzer('L1GtBoardMapsTester')
l1GtPsbSetupTest = cms.EDAnalyzer('L1GtPsbSetupTester')
l1GtPrescaleFactorsAndMasksTest = cms.EDAnalyzer('L1GtPrescaleFactorsAndMasksTester')
l1GtTriggerMenuTest = cms.EDAnalyzer('L1GtTriggerMenuTester')

printGlobalTagL1Gt = cms.Sequence(l1GtStableParametersTest
                                  *l1GtParametersTest
                                  *l1GtBoardMapsTest
                                  *l1GtPsbSetupTest
                                  *l1GtPrescaleFactorsAndMasksTest
                                  *l1GtTriggerMenuTest
                                  )


# all L1 records
printGlobalTagL1 = cms.Sequence(printGlobalTagL1Rct
                                *printGlobalTagL1Gct
#                                *printGlobalTagL1DtTPG
                                *printGlobalTagL1DtTF
#                                *printGlobalTagL1CscTF
                                *printGlobalTagL1Rpc
                                *printGlobalTagL1Gmt
                                *printGlobalTagL1Gt
                                )


