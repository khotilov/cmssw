import FWCore.ParameterSet.Config as cms

##################################################################

# useful options
isData=1 # =1 running on real data, =0 running on MC


OUTPUT_HIST='openhlt.root'
NEVTS=100
MENU="LUMI8e29" # LUMI8e29 or LUMI1e31
isRelval=1 # =1 for running on MC RelVals, =0 for standard production MC, no effect for data 

WhichHLTProcess="HLT"

####   MC cross section weights in pb, use 1 for real data  ##########

XS_7TeV_MinBias=7.126E10  # from Summer09 production
XS_10TeV_MinBias=7.528E10
XS_900GeV_MinBias=5.241E10
XSECTION=1.
FILTEREFF=1.              # gen filter efficiency

if (isData):
    XSECTION=1.         # cross section weight in pb
    FILTEREFF=1.
    MENU="LUMI8e29"

#####  Global Tag ###############################################
    
# Which AlCa condition for what. 
# See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions#Global_t
# for the complete list of global tags for a given release
# Meanwhile...:

if (isData):
    # GLOBAL_TAG='GR09_H_V6OFF::All' # collisions 2009
    # GLOBAL_TAG='GR10_H_V6A::All' # collisions2010, CMSSW_3_6_X
    GLOBAL_TAG='GR_R_38X_V9A::All' # collisions2010, preliminary tag for CMSSW_3_8_X
else:
    GLOBAL_TAG='MC_31X_V2::All'
    if (MENU == "LUMI8e29"): GLOBAL_TAG= 'STARTUP3X_V15::All'
    
    
##################################################################

process = cms.Process("ANALYSIS")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/data/Run2010A/JetMETTau/RAW/v1/000/138/564/B2D1DA16-0180-DF11-AE0D-001D09F295FB.root' # file from run 138564
    #'/store/mc/Summer09/MinBias900GeV/GEN-SIM-RAW/MC_31X_V3-v1/0029/EA9124F5-A28C-DE11-B015-0016367B47AF.root' # 900 GeV MinBias MC
    #'/store/mc/Summer09/MinBias/GEN-SIM-RAW/MC_31X_V3_7TeV-v1/0077/9E4D6E88-87A1-DE11-9D46-00E081329956.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( NEVTS ),
    skipBadFiles = cms.bool(True)
    )

process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = GLOBAL_TAG
process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')


process.load('Configuration/StandardSequences/SimL1Emulator_cff')

# OpenHLT specificss
# Define the HLT reco paths
process.load("HLTrigger.HLTanalyzers.HLTopen_cff")

# Remove the PrescaleService which, in 31X, it is expected once HLT_XXX_cff is imported
# del process.PrescaleService ## ccla no longer needed in for releases in 33x+?

process.DQM = cms.Service( "DQM",)
process.DQMStore = cms.Service( "DQMStore",)

# AlCa OpenHLT specific settings

# Define the analyzer modules
process.load("HLTrigger.HLTanalyzers.HLTAnalyser_cfi")
process.analyzeThis = cms.Path( process.hltanalysis )

process.hltanalysis.RunParameters.HistogramFile=OUTPUT_HIST
process.hltanalysis.xSection=XSECTION
process.hltanalysis.filterEff=FILTEREFF
process.hltanalysis.l1GtReadoutRecord = cms.InputTag( 'hltGtDigis','',process.name_() ) # get gtDigis extract from the RAW
process.hltanalysis.hltresults = cms.InputTag( 'TriggerResults','',WhichHLTProcess)
process.hltanalysis.HLTProcessName = WhichHLTProcess
process.hltanalysis.ht = "hltJet30Ht"
process.hltanalysis.genmet = "genMetTrue"

if (MENU == "LUMI8e29"):
    # get the objects associated with the 8e29 menu
    process.hltanalysis.recjets = "hltMCJetCorJetIcone5HF07"    
    process.hltanalysis.ht = "hltJet15UHt"    
    process.hltanalysis.IsoPixelTracksL3 = "hltHITIPTCorrector8E29"
    process.hltanalysis.IsoPixelTracksL2 = "hltIsolPixelTrackProd8E29"
    if (isData == 0):
        if(isRelval == 0):
            process.hltTrigReport.HLTriggerResults = "TriggerResults::HLT8E29"
            process.hltanalysis.l1GtObjectMapRecord = "hltL1GtObjectMap::HLT8E29"
            process.hltanalysis.hltresults = "TriggerResults::HLT8E29"
        else:
            process.hltTrigReport.HLTriggerResults = "TriggerResults::HLT"
            process.hltanalysis.l1GtObjectMapRecord = "hltL1GtObjectMap::HLT"
            process.hltanalysis.hltresults = "TriggerResults::HLT"

# pdt
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# Schedule the whole thing
if (MENU == "LUMI8e29"):
    process.schedule = cms.Schedule(
        process.DoHLTJetsU,
        process.DoHltMuon,
        process.DoHLTPhoton,
        ##    process.DoHLTElectron,
        process.DoHLTElectronStartUpWindows,
        process.DoHLTElectronLargeWindows,
        process.DoHLTElectronSiStrip,
        process.DoHLTTau,
        process.DoHLTBTag,
        process.DoHLTAlCaECALPhiSym,
        process.DoHLTAlCaPi0Eta8E29,
        # process.DoHLTIsoTrack8E29, 
        process.DoHLTMinBiasPixelTracks,
        process.analyzeThis)
else:
    process.schedule = cms.Schedule( 
        process.DoHLTJets, 
        process.DoHltMuon, 
        process.DoHLTPhoton, 
        ##    process.DoHLTElectron, 
        process.DoHLTElectronStartUpWindows, 
        process.DoHLTElectronLargeWindows,
        process.DoHLTElectronSiStrip,
        process.DoHLTTau, 
        process.DoHLTBTag,
        process.DoHLTAlCaECALPhiSym,
        process.DoHLTAlCaPi0Eta1E31,
        # process.DoHLTIsoTrack,
        process.DoHLTMinBiasPixelTracks,
        process.analyzeThis)

#########################################################################################
#
if (isData):  # replace all instances of "rawDataCollector" with "source" in InputTags
    from FWCore.ParameterSet import Mixins
    for module in process.__dict__.itervalues():
        if isinstance(module, Mixins._Parameterizable):
            for parameter in module.__dict__.itervalues():
                if isinstance(parameter, cms.InputTag):
                    if parameter.moduleLabel == 'rawDataCollector':
                        parameter.moduleLabel = 'source'
else:
    if (MENU == "LUMI8e29"):
        from FWCore.ParameterSet import Mixins
        for module in process.__dict__.itervalues():
            if isinstance(module, Mixins._Parameterizable):
                for parameter in module.__dict__.itervalues():
                    if isinstance(parameter, cms.InputTag):
                        if parameter.moduleLabel == 'rawDataCollector':
                            if(isRelval == 0):
                                parameter.moduleLabel = 'rawDataCollector::HLT8E29'
