import FWCore.ParameterSet.Config as cms

# This is the generic cfi file for CSC unpacking

muonCSCDigis = cms.EDFilter("CSCDCCUnpacker",
    # Define input to the unpacker
    InputObjects = cms.InputTag("source"),
    # Use CSC examiner to check for corrupt or semi-corrupt data & avoid unpacker crashes
    UseExaminer = cms.bool(True),
    # This mask is needed by the examiner 
    ExaminerMask = cms.uint32(0x1FEBF3F6),
    # Use Examiner to unpack good chambers and skip only bad ones
    UseSelectiveUnpacking = cms.bool(True),
    # This mask simply reduces error reporting
    ErrorMask = cms.uint32(0x0),
    # Unpack general status digis?
    UnpackStatusDigis = cms.bool(False),
    # Unpack FormatStatus digi?
    UseFormatStatus = cms.bool(True),                        
    # Turn on lots of output                            
    Debug = cms.untracked.bool(False),
    PrintEventNumber = cms.untracked.bool(False)
)


