import FWCore.ParameterSet.Config as cms

##################################################################

# useful options
gtDigisExist=0  # =1 use existing gtDigis on the input file, =0 extract gtDigis from the RAW data collection
isData=1 # =1 running on real data, =1 running on MC

OUTPUT_HIST='hltbits.root'
NEVTS=1000

##################################################################

process = cms.Process("ANALYSIS")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/express/BeamCommissioning09/OfflineMonitor/FEVTHLTALL/v2/000/122/314/F8F49C1B-60D8-DE11-AB34-001D09F28F0C.root' # file from run 122314    
    #'/store/relval/CMSSW_3_3_2/RelValZTT/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V8-v2/0000/CC203608-59C8-DE11-B1E6-0018F3D096EA.root' # MC TTBAR
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( NEVTS ),
    skipBadFiles = cms.bool(True)
    )

process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# Which AlCa condition for what. Available from pre11
# * DESIGN_31X_V1 - no smearing, alignment and calibration constants = 1.  No bad channels.
# * MC_31X_V1 (was IDEAL_31X) - conditions intended for 31X physics MC production: no smearing,
#   alignment and calibration constants = 1.  Bad channels are masked.
# * STARTUP_31X_V1 (was STARTUP_31X) - conditions needed for HLT 8E29 menu studies: As MC_31X_V1 (including bad channels),
#   but with alignment and calibration constants smeared according to knowledge from CRAFT.
# * CRAFT08_31X_V1 (was CRAFT_31X) - conditions for CRAFT08 reprocessing.
# * CRAFT_31X_V1P, CRAFT_31X_V1H - initial conditions for 2009 cosmic data taking - as CRAFT08_31X_V1 but with different
#   tag names to allow append IOV, and DT cabling map corresponding to 2009 configuration (10 FEDs).
# Meanwhile...:
#process.GlobalTag.globaltag = 'MC_31X_V9::All'
process.GlobalTag.globaltag = 'STARTUP31X_V8::All'

process.load('Configuration/StandardSequences/SimL1Emulator_cff')

# OpenHLT specificss
# Define the HLT reco paths
process.load("HLTrigger.Configuration.HLT_FULL_cff")
# Remove the PrescaleService which, in 31X, it is expected once HLT_XXX_cff is imported

process.DQM = cms.Service( "DQM",)
process.DQMStore = cms.Service( "DQMStore",)

# AlCa OpenHLT specific settings

# Define the analyzer modules
process.load("HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi")
process.hltbitanalysis.hltresults = cms.InputTag( 'TriggerResults','','HLT' )
process.hltbitanalysis.RunParameters.HistogramFile=OUTPUT_HIST

if (gtDigisExist):
    process.analyzeThis = cms.Path( process.hltbitanalysis )
else:
    process.analyzeThis = cms.Path(process.HLTBeginSequence + process.hltbitanalysis )
    process.hltbitanalysis.l1GtReadoutRecord = cms.InputTag( 'hltGtDigis','',process.name_() )


# pdt
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# Schedule the whole thing
process.schedule = cms.Schedule( 
    process.analyzeThis )

#########################################################################################
#
#nc=0
if (isData):  # replace all instances of "rawDataCollector" with "source" in InputTags
    from FWCore.ParameterSet import Mixins
    for module in process.__dict__.itervalues():
        if isinstance(module, Mixins._Parameterizable):
            for parameter in module.__dict__.itervalues():
                if isinstance(parameter, cms.InputTag):
                    if parameter.moduleLabel == 'rawDataCollector':
                        parameter.moduleLabel = 'source'
                        #print "Replacing in module: ", module
                        #nc=nc+1
#print "Number of replacements: ", nc
