import FWCore.ParameterSet.Config as cms

from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi import recoSampleDefinitionsZtoMuTau_7TeV
from TauAnalysis.Configuration.makePlots2_grid import makePlots
from TauAnalysis.Configuration.userRegistry import getHarvestingFilePath, getJobId

process = cms.Process('prodBgEstTemplatesZtoMuTau')

channel = 'ZtoMuTau_bgEstTemplate'
inputFilePath = getHarvestingFilePath(channel)
jobId = getJobId(channel)

makePlots(process, channel = channel, samples = recoSampleDefinitionsZtoMuTau_7TeV,
          inputFilePath = inputFilePath, jobId = jobId,
          analyzer_drawJobConfigurator_indOutputFileName_sets = [
              [ 'BgEstTemplateAnalyzer_ZmumuJetMisIdEnriched' ],
              [ 'BgEstTemplateAnalyzer_ZmumuMuonMisIdEnriched' ],
              [ 'BgEstTemplateAnalyzer_WplusJetsEnriched' ],
              [ 'BgEstTemplateAnalyzer_TTplusJetsEnriched' ],
              [ 'BgEstTemplateAnalyzer_QCDenriched' ]
          ],
          dqmDirectoryFilterStatistics = {
              'factorizationDisabled' : 'BgEstTemplateAnalyzer_QCDenriched/FilterStatistics/',
              'factorizationEnabled' : 'BgEstTemplateAnalyzer_QCDenriched/FilterStatistics/'
          },
          dumpDQMStore = False)

#--------------------------------------------------------------------------------
# Print-out cut-flow information for background enriched samples
# from which template histograms are obtained
#
def configureDQMDumpFilterStatisticsTable(dqmDirectoriesProcess, dqmDirectoryFilterStatisticsTable):
    module = cms.EDAnalyzer("DQMDumpFilterStatisticsTables",
        dqmDirectories = cms.PSet(),
        columnsSummaryTable = cms.vstring("Passed"),
        printSummaryTableOnly = cms.bool(True)                                                    
    )

    for processName, dqmDirectoryProcess in dqmDirectoriesProcess.items():
        dqmDirectory_full = '/harvested' + '/' + dqmDirectoryProcess + '/' + dqmDirectoryFilterStatisticsTable + '/'
        setattr(module.dqmDirectories, processName, cms.string(dqmDirectory_full))

    return module

dqmDirectoriesProcess = {
    'Ztautau'    : 'ZtautauSum',
    'Zmumu'      : 'ZmumuSum',
    'WplusJets'  : 'WplusJets_madgraph_T3',
    'QCD'        : 'qcdSum',
    'TTplusJets' : 'TTplusJets_madgraph_T3',
    'diBoson'    : 'VVsum',
    'smSum'      : 'smSum',
    'Data'       : 'data'
}    

process.dumpZtoMuTau_bgEstZmumuJetMisIdEnriched = configureDQMDumpFilterStatisticsTable(dqmDirectoriesProcess,
    'BgEstTemplateAnalyzer_ZmumuJetMisIdEnriched/FilterStatistics')
process.dumpZtoMuTau_bgEstZmumuMuonMisIdEnriched = configureDQMDumpFilterStatisticsTable(dqmDirectoriesProcess,
    'BgEstTemplateAnalyzer_ZmumuMuonMisIdEnriched/FilterStatistics')
process.dumpZtoMuTau_bgEstWplusJetsEnriched = configureDQMDumpFilterStatisticsTable(dqmDirectoriesProcess,
    'BgEstTemplateAnalyzer_WplusJetsEnriched/FilterStatistics')
process.dumpZtoMuTau_bgEstTTplusJetsEnriched = configureDQMDumpFilterStatisticsTable(dqmDirectoriesProcess,
    'BgEstTemplateAnalyzer_TTplusJetsEnriched/FilterStatistics')
process.dumpZtoMuTau_bgEstQCDenriched = configureDQMDumpFilterStatisticsTable(dqmDirectoriesProcess,
    'BgEstTemplateAnalyzer_QCDenriched/FilterStatistics')

process.dumpDQMStore = cms.EDAnalyzer("DQMStoreDump")

process.dumpZtoMuTauBgEstTemplates = cms.Sequence(
    #process.dumpDQMStore
    process.dumpZtoMuTau_bgEstZmumuJetMisIdEnriched + process.dumpZtoMuTau_bgEstZmumuMuonMisIdEnriched
   + process.dumpZtoMuTau_bgEstWplusJetsEnriched
   + process.dumpZtoMuTau_bgEstTTplusJetsEnriched
   + process.dumpZtoMuTau_bgEstQCDenriched
)    
#--------------------------------------------------------------------------------

process.q = cms.Path(process.dumpZtoMuTauBgEstTemplates)

process.schedule = cms.Schedule(process.p, process.q)

# print-out all python configuration parameter information
#print process.dumpPython()
