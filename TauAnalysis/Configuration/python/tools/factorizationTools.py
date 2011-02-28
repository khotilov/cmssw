import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.DQMTools.tools.composeSubDirectoryName import composeSubDirectoryName

import TauAnalysis.Configuration.tools.factorization as factorizer

#--------------------------------------------------------------------------------
# generic utility functions for factorization
# usable for all channels
#--------------------------------------------------------------------------------

def replaceEventSelections(analyzer, evtSel_replacements):
    # auxiliary function to replace in configuration of GenericAnalyzer
    # "tight" by "loose" cuts for factorization purposes

    for evtSel_replacement in evtSel_replacements:

        # check that all entries in evtSel_replacements list contain exactly two entries
        # (one for the "tight" cut to be replaced and one for the "loose" cut used as replacement)
        if len(evtSel_replacement) != 2:
            raise ValueError("Invalid 'evtSel_replacements' Parameter !!")

        evtSel_tight = evtSel_replacement[0]
        evtSelName_tight = getattr(evtSel_tight, "pluginName").value()
        evtSel_loose = evtSel_replacement[1]
        evtSelName_loose = getattr(evtSel_loose, "pluginName").value()

        for evtSel_i in analyzer.filters:
            if getattr(evtSel_i, "pluginName").value() == evtSelName_tight:

                analyzer.filters.remove(evtSel_i)
                analyzer.filters.append(evtSel_loose)

                if evtSelName_tight != evtSelName_loose:
                    analysisSequence = getattr(analyzer, "analysisSequence")
                    for analysisSequenceEntry in analysisSequence:
                        if hasattr(analysisSequenceEntry, "filter") and getattr(analysisSequenceEntry, "filter").value() == evtSelName_tight:
                            setattr(analysisSequenceEntry, "filter", cms.string(evtSelName_loose))

                print("Replaced in " + getattr(analyzer, "name").value() + ": "
                      + getattr(evtSel_tight, "pluginName").value() + " by " + getattr(evtSel_loose, "pluginName").value()
                      + " (version with factorization enabled)")

def _replaceAnyAnalyzerModules(analyzer, analyzersAttributeName, analyzerModule_replacements):
    # auxiliary function to replace analyzer/systematic analyzer modules
    # in configuration of GenericAnalyzer

    for analyzerModule_replacement in analyzerModule_replacements:

        # check that all entries in analyzerModule_replacements list contain exactly two entries
        # (one for the "old" module to be replaced and one for the "new" module used as replacement)
        if len(analyzerModule_replacement) != 2:
            raise ValueError("Invalid 'analyzerModule_replacements' Parameter !!")

        analyzerModule_old = analyzerModule_replacement[0]
        analyzerModuleName_old = getattr(analyzerModule_old, "pluginName").value()
        analyzerModule_new = analyzerModule_replacement[1]
        analyzerModuleName_new = getattr(analyzerModule_new, "pluginName").value()

        if hasattr(analyzer, analyzersAttributeName):
            analyzers = getattr(analyzer, analyzersAttributeName)

            for analyzerModule_i in analyzers:
                analyzerModuleName_i = getattr(analyzerModule_i, "pluginName").value()

                if analyzerModuleName_i == analyzerModuleName_old:

                    analyzers.remove(analyzerModule_i)
                    analyzers.append(analyzerModule_new)

                    if analyzerModuleName_old != analyzerModuleName_new:
                        analysisSequence = getattr(analyzer, "analysisSequence")
                        for analysisSequenceEntry in analysisSequence:
                            if hasattr(analysisSequenceEntry, "analyzers"):
                                numAnalyzerModules = len(getattr(analysisSequenceEntry, "analyzers"))
                                for j in range(numAnalyzerModules):
                                    analyzerModuleName_j = getattr(analysisSequenceEntry, "analyzers")[j]
                                    if analyzerModuleName_j == analyzerModuleName_old:
                                        getattr(analysisSequenceEntry, "analyzers")[j] = analyzerModuleName_new
 
                    print("Replaced in " + getattr(analyzer, "name").value() + ": "
                          + analyzerModuleName_old + " by " + analyzerModuleName_new
                          + " (version with factorization enabled)")
        else:
            raise ValueError("GenericAnalyzer module %s has not attribute %s !!" % (getattr(analyzer, "name"), analyzersAttributeName))

def replaceAnalyzerModules(analyzer, analyzerModule_replacements):
    # auxiliary function to replace analyzer modules
    # in configuration of GenericAnalyzer

    _replaceAnyAnalyzerModules(analyzer, "analyzers", analyzerModule_replacements)

def replaceSysAnalyzerModules(analyzer, sysAnalyzerModule_replacements):
    # auxiliary function to replace systematics analyzer modules
    # in configuration of GenericAnalyzer

    _replaceAnyAnalyzerModules(analyzer, "analyzers_systematic", sysAnalyzerModule_replacements)

#
#--------------------------------------------------------------------------------
#

def composeDirectoryName(dqmDirectory, factorizationLabel):
    if dqmDirectory.rfind("_") == -1:
        return dqmDirectory + '_' + factorizationLabel + '/'
    else:
        return dqmDirectory[:dqmDirectory.rindex("_")] + '_' + factorizationLabel + dqmDirectory[dqmDirectory.rindex("_"):] + '/'

def composeSubDirectoryNames_plots(evtSelList):
    # auxiliary function to compose names of dqmSubDirectories
    # in which histograms are stored

    dqmSubDirectoryNames = []
    for iEvtSel in range(len(evtSelList) - 1):
        afterCut = evtSelList[iEvtSel]
        beforeCut = evtSelList[iEvtSel + 1]

        dqmSubDirectoryNames.append(composeSubDirectoryName(afterCut = afterCut, beforeCut = beforeCut))

    return dqmSubDirectoryNames

def composeSubDirectoryNames_filterStatistics(evtSelList):
    # auxiliary function to compose names of dqmSubDirectories
    # in which FilterStatistics objects are stored

    dqmSubDirectoryNames = []
    for evtSel in evtSelList:
        dqmSubDirectoryNames.append(getattr(evtSel, "pluginName").value())

    return dqmSubDirectoryNames

def composeFactorizationSequence(process,
                                 processName,
                                 dqmDirectoryIn_factorizedTightEvtSel, evtSel_factorizedTight,
                                 dqmDirectoryIn_factorizedLooseEvtSel, evtSel_factorizedLoose,
                                 meName_numerator, meName_denominator,
                                 dqmDirectoryOut,
                                 dropInputDirectories = True):
    # compose sequence applying factorization
    # to histograms and FilterStatistics objects

    # configure EDAnalyzer that copies histograms filled **before**
    # cuts used for factorization are applied
    dqmHistScaler_plotsTightEvtSel = cms.EDAnalyzer("DQMHistScaler",
        dqmDirectory_input = cms.string(dqmDirectoryIn_factorizedTightEvtSel),
        dqmSubDirectories_input = cms.vstring(
            composeSubDirectoryNames_plots([ None ] + evtSel_factorizedTight + [ evtSel_factorizedLoose[0] ])
        ),
        scaleFactor = cms.double(1.),
        dqmDirectory_output = cms.string(dqmDirectoryOut)
    )

    # configure EDAnalyzer that copies FilterStatistics objects filled **before**
    # cuts used for factorization are applied
    dqmHistScaler_filterStatTightEvtSel = cms.EDAnalyzer("DQMHistScaler",
        dqmDirectory_input = cms.string(dqmDirectoryIn_factorizedTightEvtSel + "FilterStatistics" + "/"),
        dqmSubDirectories_input = cms.vstring(
            composeSubDirectoryNames_filterStatistics(evtSel_factorizedTight)
        ),
        scaleFactor = cms.double(1.),
        dqmDirectory_output = cms.string(dqmDirectoryOut + "FilterStatistics" + "/")
    )

    # configure EDAnalyzer that copies histograms filled **after**
    # cuts used for factorization are applied
    dqmHistScaler_plotsLooseEvtSel = cms.EDAnalyzer("DQMHistScaler",
        dqmDirectory_input = cms.string(dqmDirectoryIn_factorizedLooseEvtSel),
        dqmSubDirectories_input = cms.vstring(
            composeSubDirectoryNames_plots(evtSel_factorizedLoose + [ None ])
        ),
        dqmDirectory_output = cms.string(dqmDirectoryOut)
    )

    # configure EDAnalyzer that copies FilterStatistics objects filled **after**
    # cuts used for factorization are applied
    dqmHistScaler_filterStatLooseEvtSel = cms.EDAnalyzer("DQMHistScaler",
        dqmDirectory_input = cms.string(dqmDirectoryIn_factorizedLooseEvtSel + "FilterStatistics" + "/"),
        dqmSubDirectories_input = cms.vstring(
            composeSubDirectoryNames_filterStatistics(evtSel_factorizedLoose)
        ),
        dqmDirectory_output = cms.string(dqmDirectoryOut + "FilterStatistics" + "/")
    )

    # automatically add to meNames
    # suffix indicating how MonitorElements of type float
    # need to get scaled and added
    #
    # NOTE: definitions of meOptionsSeparator and meOptionsNumWeighted
    #       need to match those in TauAnalysis/Core/src/FilterStatisticsService.cc
    #
    meOptionsSeparator = "#"
    meOptionsNumWeighted = "".join([meOptionsSeparator, "a1", meOptionsSeparator, "s1"])

    if meName_numerator is not None and meName_denominator is not None:
        dqmDirectory_factorizedLooseSel = cms.string(dqmDirectoryIn_factorizedLooseEvtSel + "FilterStatistics" + "/")
        dqmDirectory_factorizedTightSel = cms.string(dqmDirectoryIn_factorizedTightEvtSel + "FilterStatistics" + "/")
        meType = cms.string("real")

        setattr(dqmHistScaler_plotsLooseEvtSel, "meName_numerator", cms.string("".join([meName_numerator, meOptionsNumWeighted])))
        setattr(dqmHistScaler_plotsLooseEvtSel, "meName_denominator", cms.string("".join([meName_denominator, meOptionsNumWeighted])))
        setattr(dqmHistScaler_plotsLooseEvtSel, "dqmDirectory_factorizedLooseSel", dqmDirectory_factorizedLooseSel)
        setattr(dqmHistScaler_plotsLooseEvtSel, "dqmDirectory_factorizedTightSel", dqmDirectory_factorizedTightSel)
        setattr(dqmHistScaler_plotsLooseEvtSel, "meType", meType)

        setattr(dqmHistScaler_filterStatLooseEvtSel, "meName_numerator", cms.string("".join([meName_numerator, meOptionsNumWeighted])))
        setattr(dqmHistScaler_filterStatLooseEvtSel, "meName_denominator", cms.string("".join([meName_denominator, meOptionsNumWeighted])))
        setattr(dqmHistScaler_filterStatLooseEvtSel, "dqmDirectory_factorizedLooseSel", dqmDirectory_factorizedLooseSel)
        setattr(dqmHistScaler_filterStatLooseEvtSel, "dqmDirectory_factorizedTightSel", dqmDirectory_factorizedTightSel)
        setattr(dqmHistScaler_filterStatLooseEvtSel, "meType", meType)
    else:
        setattr(dqmHistScaler_plotsLooseEvtSel, "scaleFactor", cms.double(1.))

        setattr(dqmHistScaler_filterStatLooseEvtSel, "scaleFactor", cms.double(1.))

    # delete original histograms and FilterStatistics objects
    # after applying factorization
    # (add drop command to last DQMHistScaler module put in sequence)
    if dropInputDirectories:
        setattr(dqmHistScaler_filterStatLooseEvtSel, "drop", cms.vstring(
            [ dqmDirectoryIn_factorizedLooseEvtSel,
              dqmDirectoryIn_factorizedTightEvtSel ] ))

    # add EDAnalyzers copying histograms and FilterStatistics objects
    # to process object
    setattr(process, "dqmHistScaler_plotsFactorizedTightEvtSel" + "_" + processName, dqmHistScaler_plotsTightEvtSel)
    setattr(process, "dqmHistScaler_filterStatFactorizedTightEvtSel" + "_" + processName, dqmHistScaler_filterStatTightEvtSel)
    setattr(process, "dqmHistScaler_plotsFactorizedLooseEvtSel" + "_" + processName, dqmHistScaler_plotsLooseEvtSel)
    setattr(process, "dqmHistScaler_filterStatFactorizedLooseEvtSel" + "_" + processName, dqmHistScaler_filterStatLooseEvtSel)

    # return sequence of all EDAnalyzers
    factorizationSequence = cms.Sequence(
        dqmHistScaler_plotsTightEvtSel + dqmHistScaler_filterStatTightEvtSel
       + dqmHistScaler_plotsLooseEvtSel + dqmHistScaler_filterStatLooseEvtSel
    )

    return factorizationSequence

#--------------------------------------------------------------------------------
# utility functions specific to factorization
# of muon isolation efficiencies in Z --> mu + tau-jet channel
#--------------------------------------------------------------------------------

def enableFactorization_runZtoMuTau(process):
    process.load("TauAnalysis.Configuration.selectZtoMuTau_factorized_cff")
    process.selectZtoMuTauEvents_factorized = cms.Sequence(
        process.selectZtoMuTauEvents
       * process.selectZtoMuTauEventsLooseMuonIsolation
    )
    process.p.replace(process.selectZtoMuTauEvents, process.selectZtoMuTauEvents_factorized)
    process.load("TauAnalysis.Configuration.analyzeZtoMuTau_factorized_cff")
    process.analyzeZtoMuTauSequence_factorized = cms.Sequence(
        process.analyzeZtoMuTauSequence_factorizedWithoutMuonIsolation
       * process.analyzeZtoMuTauSequence_factorizedWithMuonIsolation
    )
    process.p.replace(process.analyzeZtoMuTauSequence, process.analyzeZtoMuTauSequence_factorized)

##############################################
# EK: Temporary factorization update
##############################################

def enableFactorization_makeZtoMuTauPlots_grid2(
    process,
    factorizationSequenceName = "loadAndFactorizeZtoMuTauSamples",
    samplesToFactorize = [ 'InclusivePPmuX', 'PPmuXptGt20Mu10', 'PPmuXptGt20Mu15' ],
    relevantMergedSamples = [ 'qcdSum', 'smBgSum', 'smSum' ],
    mergedToRecoSampleDict = {},
    mergedSampleAdderModule = lambda sample: 'addZtoMuTau_%s' % (sample)):

    process.load("TauAnalysis.Configuration.analyzeZtoMuTau_cfi")

    # define list of event selection criteria on "tight" muon isolation branch
    # of the analysis, **before** applying factorization of muon track + ECAL
    # isolation efficiencies
    evtSelZtoMuTau_factorizedTight = [
        'evtSelGenPhaseSpace',
        'evtSelTrigger',
        'evtSelDataQuality',
        'evtSelPrimaryEventVertex',
        'evtSelPrimaryEventVertexQuality',
        'evtSelPrimaryEventVertexPosition',
        'evtSelGlobalMuon',
        'evtSelMuonEta',
        'evtSelMuonPt',
        'evtSelTauAntiOverlapWithMuonsVeto',
        'evtSelTauEta',
        'evtSelTauPt',
        'evtSelMuonVbTfId',
        'evtSelMuonPFRelIso',
    ]

    # define list of event selection criteria on "loose" muon isolation branch
    # of the analysis, **after** applying factorization of muon track + ECAL
    # isolation efficiencies
    evtSelZtoMuTau_factorizedLoose = [
        'evtSelMuonTrkIP',
        'evtSelTauLeadTrk',
        'evtSelTauLeadTrkPt',
        'evtSelTauTaNCdiscr',
        'evtSelTauProng',
        'evtSelTauCharge',
        'evtSelTauMuonVeto',
        'evtSelTauElectronVeto',
        'evtSelDiTauCandidateForMuTauAntiOverlapVeto',
        'evtSelDiTauCandidateForMuTauMt1MET',
        'evtSelDiTauCandidateForMuTauPzetaDiff',
        'evtSelDiMuPairZmumuHypothesisVetoByLooseIsolation'
    ]

    evtSelZtoMuTau_factorizedLooseOS = evtSelZtoMuTau_factorizedLoose \
      + [ 'evtSelDiTauCandidateForMuTauZeroCharge' ]
    evtSelZtoMuTau_factorizedLooseSS = evtSelZtoMuTau_factorizedLoose \
      + [ 'evtSelDiTauCandidateForMuTauNonZeroCharge' ]  

    analyzers = {
        'zMuTauAnalyzerOS' : {
            'tight_cuts' : evtSelZtoMuTau_factorizedTight,
            'loose_cuts' : evtSelZtoMuTau_factorizedLooseOS,
            'loose_analyzer' :
            'zMuTauAnalyzerOS_factorizedWithoutMuonIsolation',
            'tight_analyzer' :
            'zMuTauAnalyzerOS_factorizedWithMuonIsolation',
        },
	'zMuTauAnalyzerSS' : {
            'tight_cuts' : evtSelZtoMuTau_factorizedTight,
            'loose_cuts' : evtSelZtoMuTau_factorizedLooseSS,
            'loose_analyzer' :
            'zMuTauAnalyzerSS_factorizedWithoutMuonIsolation',
            'tight_analyzer' :
            'zMuTauAnalyzerSS_factorizedWithMuonIsolation',
        }
    }

    # Loop over the samples and create sequences
    # for each of the factorization jobs and add them to the factorization
    # sequence
    factorizationSequence = getattr(process, factorizationSequenceName)
    for sample in samplesToFactorize:
        factorizer.factorize(
            process,
            input_dir= '/harvested/%s' % sample,
            output_dir = '/harvested/%s_factorized' % sample,
	    sequence = factorizationSequence,
            analyzers = analyzers
        )

    dqmDirectoryOut = lambda sample:'/harvested/%s_factorized'% sample
    dqmDirectoryOutUnfactorized = lambda sample:'/harvested/%s'% sample

    # Now update any of the relevant mergers
    for mergedSample in relevantMergedSamples:
        # Get the module that is doing the merging, if it exists
        if not hasattr(process.mergeSamplesZtoMuTau, "merge_%s" % (mergedSample)): continue
        merger = getattr(process.mergeSamplesZtoMuTau, "merge_%s" % (mergedSample))

        # Get the subsamples associated with this merged sample
        subsamples = mergedToRecoSampleDict[mergedSample]['samples']
        # Set the adder to use our new factorized inputs
        def merge_directories(_list):
            for sample in _list:
                if sample in samplesToFactorize:
                    yield dqmDirectoryOut(sample)
                else:
                    yield dqmDirectoryOutUnfactorized(sample)

        merger.dqmDirectories_input = cms.vstring(list(merge_directories(subsamples)))

    # Update the plot sources in the plot jobs.  Note that we don't need to do
    # this for the merged samples, since we have replaced the HistAdder sources
    for scale in [ 'linear', 'log' ]:
        for analyzer in analyzers.keys():
            plotterModuleName = 'plot' + analyzer + '_' + scale
            print "Trying to update plotter:", plotterModuleName
            if hasattr(process, plotterModuleName):
                print " - found it in the process, modifying"
                plotterModuleProcesses = getattr(process, plotterModuleName).processes
                for sample in samplesToFactorize:
                    print "Plotter is plotting %s, changing dqm directory in!" %  sample
                    if hasattr(plotterModuleProcesses, sample):
                        getattr(plotterModuleProcesses, sample).dqmDirectory = \
                                cms.string("/harvested/%s_factorized" % sample)

#--------------------------------------------------------------------------------
# utility functions specific to factorization
# of muon isolation efficiencies in Z --> e + mu channel
#--------------------------------------------------------------------------------

def enableFactorization_runZtoElecMu(process):
    process.load("TauAnalysis.Configuration.selectZtoElecMu_factorized_cff")
    process.selectZtoElecMuEvents_factorized = cms.Sequence(
        process.selectZtoElecMuEvents
       * process.selectZtoElecMuEventsLooseElectronIsolation
    )
    process.p.replace(process.selectZtoElecMuEvents, process.selectZtoElecMuEvents_factorized)
    process.load("TauAnalysis.Configuration.analyzeZtoElecMu_factorized_cff")
    process.analyzeZtoElecMuEvents_factorized = cms.Sequence(
        process.analyzeZtoElecMuEvents_factorizedWithoutElectronIsolation
       * process.analyzeZtoElecMuEvents_factorizedWithElectronIsolation
    )
    process.p.replace(process.analyzeZtoElecMuEvents, process.analyzeZtoElecMuEvents_factorized)

def enableFactorization_makeZtoElecMuPlots(process,
        dqmDirectoryIn_InclusivePPmuX = '/harvested/InclusivePPmuX/zElecMuAnalyzer',
        dqmDirectoryOut_InclusivePPmuX = '/harvested/InclusivePPmuX_factorized/zElecMuAnalyzer',
        dqmDirectoryIn_PPmuXptGt20 = '/harvested/PPmuXptGt20/zElecMuAnalyzer',
        dqmDirectoryOut_PPmuXptGt20 = '/harvested/PPmuXptGt20_factorized/zElecMuAnalyzer',
        modName_addZtoElecMu_qcdSum = "addZtoElecMu_qcdSum",
        modName_addZtoElecMu_smSum = "addZtoElecMu_smSum",
        seqName_addZtoElecMu = "addZtoElecMu",
        pyObjectLabel = ""):
    process.load("TauAnalysis.Configuration.analyzeZtoElecMu_cfi")

    # define list of event selection criteria on "tight" muon isolation branch of the analysis,
    # **before** applying factorization of muon track + ECAL isolation efficiencies
    evtSelZtoElecMu_factorizedTight = [
        process.genPhaseSpaceCut,
        process.evtSelTrigger,
        process.evtSelPrimaryEventVertex,
        process.evtSelPrimaryEventVertexQuality,
        process.evtSelPrimaryEventVertexPosition,
        process.evtSelTightElectronId,
        process.evtSelElectronAntiCrack,
        process.evtSelElectronEta,
        process.evtSelElectronPt,
        process.evtSelGlobalMuon,
        process.evtSelMuonEta,
        process.evtSelMuonPt,
        process.evtSelElectronTrkIso,
        process.evtSelElectronEcalIso
    ]

    # define list of event selection criteria on "loose" muon isolation branch of the analysis,
    # **after** applying factorization of muon track + ECAL isolation efficiencies
    evtSelZtoElecMu_factorizedLoose = [
        process.evtSelElectronTrk,
        process.evtSelElectronTrkIP,
        process.evtSelMuonPFRelIso,
        process.evtSelMuonAntiPion,
        process.evtSelMuonTrkIP,
        process.evtSelDiTauCandidateForElecMuAntiOverlapVeto,
        process.evtSelDiTauCandidateForElecMuZeroCharge,
        process.evtSelDiTauCandidateForElecMuAcoplanarity12,
        process.evtSelDiTauCandidateForElecMuMt1MET,
        process.evtSelDiTauCandidateForElecMuMt2MET,
        process.evtSelDiTauCandidateForElecMuPzetaDiff
    ]

    # defines names of MonitorElements used as numerator and denominator
    # to compute factorization scale-factor
    meNameZtoElecMu_numerator = "evtSelElectronEcalIso/passed_cumulative_numWeighted"
    meNameZtoElecMu_denominator = "evtSelElectronTrkIso/processed_cumulative_numWeighted"

    # configure sequence for applying factorization to "InclusivePPmuX" process
    # (QCD background sample for Pt(hat) < 20 GeV region in phase-space)
    scaleZtoElecMu_InclusivePPmuX = composeFactorizationSequence(
        process = process,
        processName = "InclusivePPmuX" + "_" + pyObjectLabel,
        dqmDirectoryIn_factorizedTightEvtSel = composeDirectoryName(dqmDirectoryIn_InclusivePPmuX, "factorizedWithElectronIsolation"),
        evtSel_factorizedTight = evtSelZtoElecMu_factorizedTight,
        dqmDirectoryIn_factorizedLooseEvtSel = composeDirectoryName(dqmDirectoryIn_InclusivePPmuX, "factorizedWithoutElectronIsolation"),
        evtSel_factorizedLoose = evtSelZtoElecMu_factorizedLoose,
        meName_numerator = meNameZtoElecMu_numerator,
        meName_denominator = meNameZtoElecMu_denominator,
        dqmDirectoryOut = dqmDirectoryOut_InclusivePPmuX + '/'
    )

    scaleZtoElecMuName_InclusivePPmuX = "scaleZtoElecMu_InclusivePPmuX" + "_" + pyObjectLabel
    setattr(process, scaleZtoElecMuName_InclusivePPmuX, scaleZtoElecMu_InclusivePPmuX)

    # configure sequence for applying factorization to "PPmuXPPmuXptGt20" process
    # (QCD background sample for Pt(hat) > 20 GeV region in phase-space)
    scaleZtoElecMu_PPmuXptGt20 = composeFactorizationSequence(
        process = process,
        processName = "PPmuXptGt20" + "_" + pyObjectLabel,
        dqmDirectoryIn_factorizedTightEvtSel = composeDirectoryName(dqmDirectoryIn_PPmuXptGt20, "factorizedWithElectronIsolation"),
        evtSel_factorizedTight = evtSelZtoElecMu_factorizedTight,
        dqmDirectoryIn_factorizedLooseEvtSel = composeDirectoryName(dqmDirectoryIn_PPmuXptGt20, "factorizedWithoutElectronIsolation"),
        evtSel_factorizedLoose = evtSelZtoElecMu_factorizedLoose,
        meName_numerator = meNameZtoElecMu_numerator,
        meName_denominator = meNameZtoElecMu_denominator,
        dqmDirectoryOut = dqmDirectoryOut_PPmuXptGt20 + '/'
    )

    scaleZtoElecMuName_PPmuXptGt20 = "scaleZtoElecMu_PPmuXptGt20" + "_" + pyObjectLabel
    setattr(process, scaleZtoElecMuName_PPmuXptGt20, scaleZtoElecMu_PPmuXptGt20)

    # compute QCD background sum using factorized histograms and FilterStatistics objects
    addZtoElecMu_qcdSum = getattr(process, modName_addZtoElecMu_qcdSum)
    addZtoElecMu_qcdSum.qcdSum.dqmDirectories_input = cms.vstring(
        dqmDirectoryOut_InclusivePPmuX + '/',
        dqmDirectoryOut_PPmuXptGt20 + '/'
    )
    addZtoElecMu = cms.Sequence(
        getattr(process, scaleZtoElecMuName_InclusivePPmuX)
       + getattr(process, scaleZtoElecMuName_PPmuXptGt20)
    )
    addZtoElecMu._seq = addZtoElecMu._seq * getattr(process, modName_addZtoElecMu_qcdSum)
    if hasattr(process, modName_addZtoElecMu_smSum):
        addZtoElecMu._seq = addZtoElecMu._seq * getattr(process, modName_addZtoElecMu_smSum)
    setattr(process, seqName_addZtoElecMu, addZtoElecMu)

    if hasattr(process, "plotZtoElecMu"):
        process.plotZtoElecMu.processes.InclusivePPmuX.dqmDirectory = cms.string('/harvested/InclusivePPmuX_factorized')
        process.plotZtoElecMu.processes.PPmuXptGt20.dqmDirectory = cms.string('/harvested/PPmuXptGt20_factorized')

#--------------------------------------------------------------------------------
# utility functions specific to factorization
# of muon isolation efficiencies in Z --> e + tau channel
#--------------------------------------------------------------------------------

def enableFactorization_runZtoElecTau(process):
    process.load("TauAnalysis.Configuration.selectZtoElecTau_factorized_cff")
    process.selectZtoElecTauEvents_factorized = cms.Sequence( process.selectZtoElecTauEvents
                                                            *process.selectZtoElecTauEventsLooseElectronIsolation )
    process.p.replace(process.selectZtoElecTauEvents, process.selectZtoElecTauEvents_factorized)
    process.load("TauAnalysis.Configuration.analyzeZtoElecTau_factorized_cff")
    process.analyzeZtoElecTauSequence_factorized = cms.Sequence( process.analyzeZtoElecTauSequence_factorizedWithoutElectronIsolation
                                                             *process.analyzeZtoElecTauSequence_factorizedWithElectronIsolation )
    process.p.replace(process.analyzeZtoElecTauSequence, process.analyzeZtoElecTauSequence_factorized)



def enableFactorization_makeZtoElecTauPlots_grid2(
    process,
    factorizationSequenceName = "loadAndFactorizeZtoElecTauSamples",
    samplesToFactorize = [ ],
    relevantMergedSamples = [ 'qcdSum','photonPlusJetsSum' ],
    mergedToRecoSampleDict = {},
    mergedSampleAdderModule = lambda sample: 'addZtoElecTau_%s' % (sample)):

    process.load("TauAnalysis.Configuration.analyzeZtoElecTau_cfi")

    # define list of event selection criteria on "tight" electron isolation branch
    # of the analysis, **before** applying factorization of electron isolation efficiency
    evtSelZtoElecTau_factorizedTight = [
        'evtSelGenPhaseSpace',
        'evtSelTrigger',
        'evtSelPrimaryEventVertex',
        'evtSelPrimaryEventVertexQuality',
        'evtSelPrimaryEventVertexPosition',
		'evtSelElectronId',
        'evtSelElectronAntiCrack',
        'evtSelElectronEta',
        'evtSelElectronPt',
        'evtSelTauAntiOverlapWithElectronsVeto',
        'evtSelTauEta',
        'evtSelTauPt',
        'evtSelElectronIso',
        'evtSelElectronConversionVeto'
    ]

    # define list of event selection criteria on "loose" electron isolation branch
    # of the analysis, **after** applying factorization of electron isolation efficiency
    evtSelZtoElecTau_factorizedLoose = [
        'evtSelElectronTrkIP',
        'evtSelTauLeadTrk',
        'evtSelTauLeadTrkPt',
        'evtSelTauTaNCdiscr',
        'evtSelTauProng',
        'evtSelTauCharge',
        'evtSelTauElectronVeto',
		'evtSelTauEcalCrackVeto',
        'evtSelTauMuonVeto',
        'evtSelDiTauCandidateForElecTauAntiOverlapVeto',
        'evtSelDiTauCandidateForElecTauMt1MET',
        'evtSelDiTauCandidateForElecTauPzetaDiff',
        'evtSelDiElecPairZeeHypothesisVetoByLooseIsolation'
    ]

    evtSelZtoElecTau_factorizedLooseOS = evtSelZtoElecTau_factorizedLoose \
      + [ 'evtSelDiTauCandidateForElecTauZeroCharge' ]
    evtSelZtoElecTau_factorizedLooseSS = evtSelZtoElecTau_factorizedLoose \
      + [ 'evtSelDiTauCandidateForElecTauNonZeroCharge' ]  

    analyzers = {
        'zElecTauAnalyzerOS' : {
            'tight_cuts' : evtSelZtoElecTau_factorizedTight,
            'loose_cuts' : evtSelZtoElecTau_factorizedLooseOS,
            'loose_analyzer' :
            'zElecTauAnalyzerOS_factorizedWithoutElectronIsolation',
            'tight_analyzer' :
            'zElecTauAnalyzerOS_factorizedWithElectronIsolation',
        },
		'zElecTauAnalyzerSS' : {
            'tight_cuts' : evtSelZtoElecTau_factorizedTight,
            'loose_cuts' : evtSelZtoElecTau_factorizedLooseSS,
            'loose_analyzer' :
            'zElecTauAnalyzerSS_factorizedWithoutElectronIsolation',
            'tight_analyzer' :
            'zElecTauAnalyzerSS_factorizedWithElectronIsolation',
        }
    }

    # Loop over the samples and create sequences
    # for each of the factorization jobs and add them to the factorization
    # sequence
    factorizationSequence = getattr(process, factorizationSequenceName)
    for sample in samplesToFactorize:
        factorizer.factorize(
            process,
            input_dir= '/harvested/%s' % sample,
            output_dir = '/harvested/%s_factorized' % sample,
	        sequence = factorizationSequence,
            analyzers = analyzers
        )

    dqmDirectoryOut = lambda sample:'/harvested/%s_factorized'% sample
    dqmDirectoryOutUnfactorized = lambda sample:'/harvested/%s'% sample

    # Now update any of the relevant mergers
    for mergedSample in relevantMergedSamples:
        # Get the module that is doing the merging, if it exists
        if not hasattr(process.mergeSamplesZtoElecTau, "merge_%s" % (mergedSample)): continue
        merger = getattr(process.mergeSamplesZtoElecTau, "merge_%s" % (mergedSample))

        # Get the subsamples associated with this merged sample
        subsamples = mergedToRecoSampleDict[mergedSample]['samples']
        # Set the adder to use our new factorized inputs
        def merge_directories(_list):
            for sample in _list:
                if sample in samplesToFactorize:
                    yield dqmDirectoryOut(sample)
                else:
                    yield dqmDirectoryOutUnfactorized(sample)

        merger.dqmDirectories_input = cms.vstring(list(merge_directories(subsamples)))

    # Update the plot sources in the plot jobs.  Note that we don't need to do
    # this for the merged samples, since we have replaced the HistAdder sources
    for scale in [ 'linear', 'log' ]:
        for analyzer in analyzers.keys():
            plotterModuleName = 'plot' + analyzer + '_' + scale
            print "Trying to update plotter:", plotterModuleName
            if hasattr(process, plotterModuleName):
                print " - found it in the process, modifying"
                plotterModuleProcesses = getattr(process, plotterModuleName).processes
                for sample in samplesToFactorize:
                    print "Plotter is plotting %s, changing dqm directory in!" %  sample
                    if hasattr(plotterModuleProcesses, sample):
                        getattr(plotterModuleProcesses, sample).dqmDirectory = \
                                cms.string("/harvested/%s_factorized" % sample)

#--------------------------------------------------------------------------------
# utility functions specific to factorization
# of muon isolation efficiencies in A/H --> e + tau channel
#--------------------------------------------------------------------------------

def enableFactorization_runAHtoElecTau(process):
    process.load("TauAnalysis.Configuration.selectAHtoElecTau_factorized_cff")
    process.selectAHtoElecTauEvents_factorized = cms.Sequence( process.selectAHtoElecTauEvents
                                                            *process.selectAHtoElecTauEventsLooseElectronIsolation )
    process.p.replace(process.selectAHtoElecTauEvents, process.selectAHtoElecTauEvents_factorized)
    process.load("TauAnalysis.Configuration.analyzeAHtoElecTau_factorized_cff")
    process.analyzeAHtoElecTauSequence_factorized = cms.Sequence( process.analyzeAHtoElecTauSequence_factorizedWithoutElectronIsolation
                                                             *process.analyzeAHtoElecTauSequence_factorizedWithElectronIsolation )
    process.p.replace(process.analyzeAHtoElecTauSequence, process.analyzeAHtoElecTauSequence_factorized)


def enableFactorization_makeAHtoElecTauPlots_grid2(
    process,
    factorizationSequenceName = "loadAndFactorizeAHtoElecTauSamples",
    samplesToFactorize = [  ],
    relevantMergedSamples = [ 'qcdSum', 'photonPlusJetsSum' ],
    mergedToRecoSampleDict = {},
    mergedSampleAdderModule = lambda sample, btag: 'addAHtoElecTau_%s_%s' % (btag, sample)):

    process.load("TauAnalysis.Configuration.analyzeAHtoElecTau_cfi")

    # define list of event selection criteria on "tight" electron isolation branch
    # of the analysis, **before** applying factorization of electron isolation efficiency
    tight_cuts = [
        'evtSelGenPhaseSpace',
        'evtSelTrigger',
        'evtSelPrimaryEventVertex',
        'evtSelPrimaryEventVertexQuality',
        'evtSelPrimaryEventVertexPosition',
        'evtSelElectronId',
		'evtSelElectronAntiCrack',
        'evtSelElectronEta',
        'evtSelElectronPt',
        'evtSelTauAntiOverlapWithElectronsVeto',
        'evtSelTauEta',
        'evtSelTauPt',
        'evtSelElectronIso',
        'evtSelElectronConversionVeto',
    ]

    # define list of event selection criteria on "loose" electron isolation branch
    # of the analysis, **after** applying factorization of electron isolation efficiency
    loose_cuts_base = [
        'evtSelElectronTrkIP',
        'evtSelTauLeadTrk',
        'evtSelTauLeadTrkPt',
        'evtSelTauTaNCdiscr',
        'evtSelTauProng',
        'evtSelTauCharge',
        'evtSelTauElectronVeto',
        'evtSelTauEcalCrackVeto',
        'evtSelTauMuonVeto',
        'evtSelDiTauCandidateForAHtoElecTauAntiOverlapVeto',
        'evtSelDiTauCandidateForAHtoElecTauMt1MET',
        'evtSelDiTauCandidateForAHtoElecTauPzetaDiff',
        'evtSelDiElecPairZeeHypothesisVetoByLooseIsolation',
    ]

    loose_cuts_woBtag = loose_cuts_base \
      + [ 'evtSelNonCentralJetEt20bTag', 'evtSelDiTauCandidateForAHtoElecTauZeroCharge' ]
    loose_cuts_wBtag = loose_cuts_base \
      + [ 'evtSelCentralJetEt20', 'evtSelCentralJetEt20bTag', 'evtSelDiTauCandidateForAHtoElecTauZeroCharge' ]

    # Same sign verions of final cuts
    loose_cuts_woBtag_ss = loose_cuts_base \
      + [ 'evtSelNonCentralJetEt20bTag', 'evtSelDiTauCandidateForAHtoElecTauNonZeroCharge' ]
    loose_cuts_wBtag_ss = loose_cuts_base \
      + [ 'evtSelCentralJetEt20', 'evtSelCentralJetEt20bTag', 'evtSelDiTauCandidateForAHtoElecTauNonZeroCharge' ]

    analyzers = {
        'ahElecTauAnalyzerOS_wBtag' : {
            'tight_cuts' : tight_cuts,
            'loose_cuts' : loose_cuts_wBtag,
            'loose_analyzer' :
            'ahElecTauAnalyzerOS_wBtag_factorizedWithoutElectronIsolation',
            'tight_analyzer' :
            'ahElecTauAnalyzerOS_wBtag_factorizedWithElectronIsolation'
        },
        'ahElecTauAnalyzerOS_woBtag' : {
            'tight_cuts' : tight_cuts,
            'loose_cuts' : loose_cuts_woBtag,
            'loose_analyzer' :
            'ahElecTauAnalyzerOS_woBtag_factorizedWithoutElectronIsolation',
            'tight_analyzer' :
            'ahElecTauAnalyzerOS_woBtag_factorizedWithElectronIsolation'
        },
        'ahElecTauAnalyzerSS_wBtag' : {
            'tight_cuts' : tight_cuts,
            'loose_cuts' : loose_cuts_wBtag_ss,
            'loose_analyzer' :
            'ahElecTauAnalyzerSS_wBtag_factorizedWithoutElectronIsolation',
            'tight_analyzer' :
            'ahElecTauAnalyzerSS_wBtag_factorizedWithElectronIsolation'
        },
        'ahElecTauAnalyzerSS_woBtag' : {
            'tight_cuts' : tight_cuts,
            'loose_cuts' : loose_cuts_woBtag_ss,
            'loose_analyzer' :
            'ahElecTauAnalyzerSS_woBtag_factorizedWithoutElectronIsolation',
            'tight_analyzer' :
            'ahElecTauAnalyzerSS_woBtag_factorizedWithElectronIsolation'
        },
    }

    # Loop over the samples and create sequences
    # for each of the factorization jobs and add them to the factorization
    # sequence
    factorizationSequence = getattr(process, factorizationSequenceName)
    for sample in samplesToFactorize:
        factorizer.factorize(
            process,
            input_dir= '/harvested/%s' % sample,
            output_dir = '/harvested/%s_factorized' % sample,
            sequence = factorizationSequence,
            analyzers = analyzers
        )

    dqmDirectoryOut = lambda sample:'/harvested/%s_factorized' % sample
    dqmDirectoryOutUnfactorized = lambda sample:'/harvested/%s'% sample

    # Now update any of the relevant mergers
    for mergedSample in relevantMergedSamples:
        # Get the module that is doing the merging, if it exists
        merger_name = "merge_%s" % mergedSample
        if not hasattr(process.mergeSamplesAHtoElecTau, merger_name):
            print "factorizationTools: Expected to update ",\
                    merger_name, "but it's not in the process! skipping.."
            continue
        merger = getattr(process.mergeSamplesAHtoElecTau, merger_name)

        # Get the subsamples associated with this merged sample
        subsamples = mergedToRecoSampleDict[mergedSample]['samples']
        # Set the adder to use our new factorized inputs
        def merge_directories(_list):
            for sample in _list:
                if sample in samplesToFactorize:
                    yield dqmDirectoryOut(sample)
                else:
                    yield dqmDirectoryOutUnfactorized(sample)

        merger.dqmDirectories_input = cms.vstring(list(merge_directories(subsamples)))

    # Update the plot sources in the plot jobs.  Note that we don't need to do
    # this for the merged samples, since we have replaced the HistAdder sources
    for scale in ['linear', 'log']:
        for analyzer in analyzers.keys():
            plotterModuleName = 'plot' + analyzer + '_' + scale
            print "Trying to update plotter:", plotterModuleName
            if hasattr(process, plotterModuleName):
                print " - found it in the process, modifying"
                plotterModuleProcesses = getattr(process, plotterModuleName).processes
                for sample in samplesToFactorize:
                    if hasattr(plotterModuleProcesses, sample):
                        print "Plotter is plotting %s, changing dqm directory in!" %  sample
                        getattr(plotterModuleProcesses, sample).dqmDirectory = \
                          cms.string("/harvested/%s_factorized" % sample)

#--------------------------------------------------------------------------------
# utility functions specific to factorization
# of tau id. efficiencies in Z --> tau-jet + tau-jet channel
#--------------------------------------------------------------------------------

def enableFactorization_runZtoDiTau(process):
    process.load("TauAnalysis.Configuration.selectZtoDiTau_factorized_cff")
    process.selectZtoDiTauEvents_factorized = cms.Sequence(
        process.selectZtoDiTauEvents
       * process.selectZtoDiTauEventsLoose2ndTau
    )
    process.p.replace(process.selectZtoDiTauEvents, process.selectZtoDiTauEvents_factorized)
    process.load("TauAnalysis.Configuration.analyzeZtoDiTau_factorized_cff")
    process.analyzeZtoDiTauEvents_factorized = cms.Sequence(
        process.analyzeZtoDiTauEvents_factorizedLoose2ndTau
       * process.analyzeZtoDiTauEvents_factorizedTight2ndTau
    )
    process.p.replace(process.analyzeZtoDiTauEvents, process.analyzeZtoDiTauEvents_factorized)

def enableFactorization_makeZtoDiTauPlots_grid(
    process,
    factorizationSequenceName = "loadAndFactorizeZtoDiTauSamples",
    samplesToFactorize = [ 'qcdDiJet' ],
    relevantMergedSamples = [ 'qcdSum' ],
    mergedToRecoSampleDict = {},
    mergedSampleAdderModule = lambda sample: 'addZtoDiTau_%s' % (sample),
    dqmDirectoryOut =
    lambda sample:'/harvested/%s_factorized/zDiTauAnalyzer/'% (sample),
    dqmDirectoryOutUnfactorized =
    lambda sample:'/harvested/%s/zDiTauAnalyzer/'% (sample),
    dqmDirectoryTight =
    lambda sample:'/harvested/%s/zDiTauAnalyzer_factorizedTight2ndTau/' % (sample),
    dqmDirectoryLoose =
    lambda sample:'/harvested/%s/zDiTauAnalyzer_factorizedLoose2ndTau/' % (sample),
    pyObjectLabel = ""):

    process.load("TauAnalysis.Configuration.analyzeZtoDiTau_cfi")

    # define list of event selection criteria on "tight" tau id. branch
    # of the analysis, **before** applying factorization of
    # lead. track Pt, track isolation and ECAL isolation efficiencies
    evtSelZtoDiTau_factorizedTight = [
        process.evtSelGenPhaseSpace,
        ##process.evtSelTrigger,
        process.evtSelDataQuality,
        process.evtSelPrimaryEventVertex,
        process.evtSelPrimaryEventVertexQuality,
        process.evtSelPrimaryEventVertexPosition,
        process.evtSelFirstTauEta,
        process.evtSelFirstTauPt,
        process.evtSelSecondTauEta,
        process.evtSelSecondTauPt,
        process.evtSelFirstTauLeadTrk,
        process.evtSelFirstTauLeadTrkPt,
        process.evtSelFirstTauTaNCdiscr,
        process.evtSelFirstTauTrkIso,
        process.evtSelFirstTauEcalIso,
        process.evtSelFirstTauProng,
        process.evtSelFirstTauCharge,
        process.evtSelFirstTauMuonVeto,
        process.evtSelFirstTauElectronVeto,
        process.evtSelSecondTauLeadTrk,
        process.evtSelSecondTauLeadTrkPt,
        process.evtSelSecondTauTaNCdiscr,
        process.evtSelSecondTauTrkIso,
        process.evtSelSecondTauEcalIso,
        process.evtSelSecondTauProng,
        process.evtSelSecondTauCharge
    ]

    # define list of event selection criteria on "loose" tau id. branch
    # of the analysis, **after** applying factorization of
    # lead. track Pt, track isolation and ECAL isolation efficiencies
    evtSelZtoDiTau_factorizedLoose = [
        process.evtSelSecondTauMuonVeto,
        process.evtSelSecondTauElectronVeto,
        process.evtSelDiTauCandidateForDiTauAntiOverlapVeto,
        process.evtSelDiTauCandidateForDiTauZeroCharge,
        process.evtSelDiTauCandidateForDiTauAcoplanarity,
        process.evtSelDiTauCandidateForDiTauPzetaDiff,
        #process.evtSelCentralJetVeto,
        process.evtSelTrigger
    ]

    # defines names of MonitorElements used as numerator and denominator
    # to compute factorization scale-factor
    meNameZtoDiTau_numerator = "evtSelDiTauCandidateForDiTauPzetaDiff/passed_cumulative_numWeighted"
    meNameZtoDiTau_denominator = "evtSelSecondTauLeadTrkPt/processed_cumulative_numWeighted"

    # Loop over the samples and create sequences
    # for each of the factorization jobs and add them to the factorization
    # sequence
    factorizationSequence = getattr(process, factorizationSequenceName)
    for sample in samplesToFactorize:
        new_factorization_sequence = composeFactorizationSequence(
            process = process,
            processName = sample + "_" + pyObjectLabel,
            dqmDirectoryIn_factorizedTightEvtSel = dqmDirectoryTight(sample),
            evtSel_factorizedTight = evtSelZtoDiTau_factorizedTight,
            dqmDirectoryIn_factorizedLooseEvtSel = dqmDirectoryLoose(sample),
            evtSel_factorizedLoose = evtSelZtoDiTau_factorizedLoose,
            meName_numerator = meNameZtoDiTau_numerator,
            meName_denominator = meNameZtoDiTau_denominator,
            dqmDirectoryOut = dqmDirectoryOut(sample),
            dropInputDirectories = False
        )
        new_factorization_seq_name = "scaleZtoDiTau_%s_%s" % (sample, pyObjectLabel)
        setattr(process, new_factorization_seq_name, new_factorization_sequence)
        factorizationSequence += new_factorization_sequence

    # Now update any of the relevant mergers
    for mergedSample in relevantMergedSamples:
        # Get the module that is doing the merging, if it exists
        if not hasattr(process.mergeSamplesZtoDiTau, "merge_%s_zDiTauAnalyzer" % (mergedSample)): continue
        merger = getattr(process.mergeSamplesZtoDiTau, "merge_%s_zDiTauAnalyzer" % (mergedSample))

        # Get the subsamples associated with this merged sample
        subsamples = mergedToRecoSampleDict[mergedSample]['samples']
        # Set the adder to use our new factorized inputs
        def merge_directories(_list):
            for sample in _list:
                if sample in samplesToFactorize:
                    yield dqmDirectoryOut(sample)
                else:
                    yield dqmDirectoryOutUnfactorized(sample)

        merger.dqmDirectories_input = cms.vstring(list(merge_directories(subsamples)))

    # Update the plot sources in the plot jobs.  Note that we don't need to do
    # this for the merged samples, since we have replaced the HistAdder sources
    for plotterModuleName in [ 'plotZtoDiTau_log', 'plotZtoDiTau_linear' ]:
        plotterModuleProcesses = getattr(process, plotterModuleName).processes
        for sample in samplesToFactorize:
            if hasattr(plotterModuleProcesses, sample):
                getattr(plotterModuleProcesses, sample).dqmDirectory = \
                        cms.string("/harvested/%s_factorized" % sample)

#--------------------------------------------------------------------------------
# utility functions specific to factorization
# of tau isolation efficiencies in W --> tau-jet + nu channel
#--------------------------------------------------------------------------------

def enableFactorization_runWtoTauNu(process):
    process.load("TauAnalysis.Configuration.selectWtoTauNu_factorized_cff")
    process.selectWtoTauNuEvents_factorized = cms.Sequence( process.selectWtoTauNuEvents
                                                           *process.selectWtoTauNuEventsLooseTauIsolation )
    process.p.replace(process.selectWtoTauNuEvents, process.selectWtoTauNuEvents_factorized)
    process.load("TauAnalysis.Configuration.analyzeWtoTauNu_factorized_cff")
    process.analyzeWtoTauNuEvents_factorized = cms.Sequence( process.analyzeWtoTauNuEvents_factorizedWithoutTauIsolation
                                                            *process.analyzeWtoTauNuEvents_factorizedWithTauIsolation )
    process.p.replace(process.analyzeWtoTauNuEvents, process.analyzeWtoTauNuEvents_factorized)

def enableFactorization_makeWtoTauNuPlots(process,
      dqmDirectoryIn_qcd_W = '/harvested/qcd_W/wTauNuAnalyzer',
      dqmDirectoryOut_qcd_W = '/harvested/qcd_W_factorized/wTauNuAnalyzer',
      modName_addWtoTauNu_qcd = "addWtoTauNu_qcd",
     # modName_addWtoTauNu_smSum = "addWtoTauNu_smSum",
      seqName_addWtoTauNu = "addWtoTauNu",
      pyObjectLabel = ""):
    process.load("TauAnalysis.Configuration.analyzeWtoTauNu_cfi")
    # define list of event selection criteria on "tight" tau isolation branch of the analysis,
    # **before** applying factorization of tauTaNC+prong+charge efficiencies
    evtSelWtoTauNu_factorizedTight = [
        process.evtSelPrimaryEventVertex,
        process.evtSelPrimaryEventVertexQuality,
        process.evtSelPrimaryEventVertexPosition,
        process.evtSelTauEta,
        process.evtSelTauPt,
        process.evtSelPFMetPt,
        process.evtSelMetPt,
        process.evtSelTauLeadTrk,
        process.evtSelTauLeadTrkPt,
        process.evtSelTauIso,
        process.evtSelTauTaNC,
        process.evtSelTauProng,
        process.evtSelTauCharge,
        process.evtSelTauMuonVeto,
        process.evtSelTauElectronVeto,
        process.evtSelTauEcalCrackVeto
        ]

    # define list of event selection criteria on "loose" muon isolation branch of the analysis,
    # **after** applying factorization of muon track + ECAL isolation efficiencies
    evtSelWtoTauNu_factorizedLoose = [
        process.evtSelCentralJetVeto,
        process.evtSelRecoilEnergyFromCaloTowers,
        process.evtSelMetTopology
    ]

  # defines names of MonitorElements used as numerator and denominator
    # to compute factorization scale-factor
    meNameWtoTauNu_numerator = "evtSelTauEcalCrackVeto/passed_cumulative_numWeighted"
    meNameWtoTauNu_denominator = "evtSelTauLeadTrkPt/processed_cumulative_numWeighted"


   # configure sequence for applying factorization to "qcd_W" process (QCD background sample for Pt(hat) > 15 GeV)
    scaleWtoTauNu_qcd_W = composeFactorizationSequence(
        process = process,
        processName = "qcd_W" + "_" + pyObjectLabel,
        dqmDirectoryIn_factorizedTightEvtSel = "/harvested/qcd_W/wTauNuAnalyzer_factorizedWithTauIsolation/",
        evtSel_factorizedTight = evtSelWtoTauNu_factorizedTight,
        dqmDirectoryIn_factorizedLooseEvtSel = "/harvested/qcd_W/wTauNuAnalyzer_factorizedWithoutTauIsolation/",
        evtSel_factorizedLoose = evtSelWtoTauNu_factorizedLoose,
        meName_numerator = meNameWtoTauNu_numerator,
        meName_denominator = meNameWtoTauNu_denominator,
        dqmDirectoryOut = dqmDirectoryOut_qcd_W + '/'
    )

    scaleWtoTauNuName_qcd_W = "scaleWtoTauNu_qcd_W" + "_" + pyObjectLabel
    setattr(process,scaleWtoTauNuName_qcd_W, scaleWtoTauNu_qcd_W)

    addWtoTauNu_qcd = getattr(process, modName_addWtoTauNu_qcd)
    addWtoTauNu_qcd.qcd.dqmDirectories_input = cms.vstring(
        dqmDirectoryOut_qcd_W + '/'
        )
    addWtoTauNu = cms.Sequence(
        getattr(process, scaleWtoTauNuName_qcd_W)
        )
    addWtoTauNu._seq = addWtoTauNu._seq * getattr(process, modName_addWtoTauNu_qcd)
#    if hasattr(process, modName_addWtoTauNu_smSum):
#        addWtoTauNu._seq = addWtoTauNu._seq * getattr(process,modName_addWtoTauNu_smSum)
    setattr(process,seqName_addWtoTauNu, addWtoTauNu)

    process.plotWtoTauNu.processes.qcd_W.dqmDirectory = cms.string('/harvested/qcd_W_factorized')

#--------------------------------------------------------------------------------
# utility functions specific to factorization
# of muon isolation efficiencies in MSSM Higgs A/H --> mu + tau-jet channel
#--------------------------------------------------------------------------------

def enableFactorization_runAHtoMuTau(process):
    process.load("TauAnalysis.Configuration.selectAHtoMuTau_factorized_cff")
    process.selectAHtoMuTauEvents_factorized = cms.Sequence(
        process.selectAHtoMuTauEvents
       * process.selectAHtoMuTauEventsLooseMuonIsolation
    )
    process.p.replace(process.selectAHtoMuTauEvents, process.selectAHtoMuTauEvents_factorized)
    process.load("TauAnalysis.Configuration.analyzeAHtoMuTau_factorized_cff")
    process.analyzeAHtoMuTauSequence_factorized = cms.Sequence(
        process.analyzeAHtoMuTauSequence_factorizedWithoutMuonIsolation
       * process.analyzeAHtoMuTauSequence_factorizedWithMuonIsolation
    )
    process.p.replace(process.analyzeAHtoMuTauSequence, process.analyzeAHtoMuTauSequence_factorized)

def enableFactorization_makeAHtoMuTauPlots_grid2(
    process,
    factorizationSequenceName = "loadAndFactorizeAHtoMuTauSamples",
    samplesToFactorize = [ 'InclusivePPmuX', 'PPmuXptGt20Mu10', 'PPmuXptGt20Mu15' ],
    relevantMergedSamples = [ 'qcdSum', ],
    mergedToRecoSampleDict = {},
    mergedSampleAdderModule = lambda sample, btag: 'addAHtoMuTau_%s_%s' % (btag, sample)):

    process.load("TauAnalysis.Configuration.analyzeAHtoMuTau_cfi")

    # define list of event selection criteria on "tight" muon isolation branch
    # of the analysis, **before** applying factorization of muon track + ECAL
    # isolation efficiencies
    tight_cuts = [
        'evtSelGenPhaseSpace',
        'evtSelTrigger',
        'evtSelDataQuality',
        'evtSelPrimaryEventVertex',
        'evtSelPrimaryEventVertexQuality',
        'evtSelPrimaryEventVertexPosition',
        'evtSelGlobalMuon',
        'evtSelMuonEta',
        'evtSelMuonPt',
        'evtSelTauAntiOverlapWithMuonsVeto',
        'evtSelTauEta',
        'evtSelTauPt',
        'evtSelMuonVbTfId',
        'evtSelMuonPFRelIso',
    ]

    # define list of event selection criteria on "loose" muon isolation branch
    # of the analysis, **after** applying factorization of muon track + ECAL
    # isolation efficiencies
    loose_cuts_base = [
        'evtSelMuonTrkIP',
        'evtSelTauLeadTrk',
        'evtSelTauLeadTrkPt',
        'evtSelTauTaNCdiscr',
        'evtSelTauProng',
        'evtSelTauCharge',
        'evtSelTauMuonVeto',
        'evtSelTauElectronVeto',
        'evtSelDiTauCandidateForAHtoMuTauAntiOverlapVeto',
        'evtSelDiTauCandidateForAHtoMuTauMt1MET',
        'evtSelDiTauCandidateForAHtoMuTauPzetaDiff',
        'evtSelDiMuPairZmumuHypothesisVetoByLooseIsolation',
    ]

    loose_cuts_woBtag = loose_cuts_base \
      + [ 'evtSelNonCentralJetEt20bTag', 'evtSelDiTauCandidateForAHtoMuTauZeroCharge' ]
    loose_cuts_wBtag = loose_cuts_base \
      + [ 'evtSelCentralJetEt20', 'evtSelCentralJetEt20bTag', 'evtSelDiTauCandidateForAHtoMuTauZeroCharge' ]

    # Same sign verions of final cuts
    loose_cuts_woBtag_ss = loose_cuts_base \
      + [ 'evtSelNonCentralJetEt20bTag', 'evtSelDiTauCandidateForAHtoMuTauNonZeroCharge' ]
    loose_cuts_wBtag_ss = loose_cuts_base \
      + [ 'evtSelCentralJetEt20', 'evtSelCentralJetEt20bTag', 'evtSelDiTauCandidateForAHtoMuTauNonZeroCharge' ]

    analyzers = {
        'ahMuTauAnalyzerOS_wBtag' : {
            'tight_cuts' : tight_cuts,
            'loose_cuts' : loose_cuts_wBtag,
            'loose_analyzer' :
            'ahMuTauAnalyzerOS_wBtag_factorizedWithoutMuonIsolation',
            'tight_analyzer' :
            'ahMuTauAnalyzerOS_wBtag_factorizedWithMuonIsolation',
        },
        'ahMuTauAnalyzerOS_woBtag' : {
            'tight_cuts' : tight_cuts,
            'loose_cuts' : loose_cuts_woBtag,
            'loose_analyzer' :
            'ahMuTauAnalyzerOS_woBtag_factorizedWithoutMuonIsolation',
            'tight_analyzer' :
            'ahMuTauAnalyzerOS_woBtag_factorizedWithMuonIsolation',
        },
        'ahMuTauAnalyzerSS_wBtag' : {
            'tight_cuts' : tight_cuts,
            'loose_cuts' : loose_cuts_wBtag_ss,
            'loose_analyzer' :
            'ahMuTauAnalyzerSS_wBtag_factorizedWithoutMuonIsolation',
            'tight_analyzer' :
            'ahMuTauAnalyzerSS_wBtag_factorizedWithMuonIsolation',
        },
        'ahMuTauAnalyzerSS_woBtag' : {
            'tight_cuts' : tight_cuts,
            'loose_cuts' : loose_cuts_woBtag_ss,
            'loose_analyzer' :
            'ahMuTauAnalyzerSS_woBtag_factorizedWithoutMuonIsolation',
            'tight_analyzer' :
            'ahMuTauAnalyzerSS_woBtag_factorizedWithMuonIsolation',
        },
    }

    # Loop over the samples and create sequences
    # for each of the factorization jobs and add them to the factorization
    # sequence
    factorizationSequence = getattr(process, factorizationSequenceName)
    for sample in samplesToFactorize:
        factorizer.factorize(
            process,
            input_dir= '/harvested/%s' % sample,
            output_dir = '/harvested/%s_factorized' % sample,
            sequence = factorizationSequence,
            analyzers = analyzers
        )

    dqmDirectoryOut = lambda sample:'/harvested/%s_factorized' % sample
    dqmDirectoryOutUnfactorized = lambda sample:'/harvested/%s'% sample

    # Now update any of the relevant mergers
    for mergedSample in relevantMergedSamples:
        # Get the module that is doing the merging, if it exists
        merger_name = "merge_%s" % mergedSample
        if not hasattr(process.mergeSamplesAHtoMuTau, merger_name):
            print "factorizationTools: Expected to update ",\
                    merger_name, "but it's not in the process! skipping.."
            continue
        merger = getattr(process.mergeSamplesAHtoMuTau, merger_name)

        # Get the subsamples associated with this merged sample
        subsamples = mergedToRecoSampleDict[mergedSample]['samples']
        # Set the adder to use our new factorized inputs
        def merge_directories(_list):
            for sample in _list:
                if sample in samplesToFactorize:
                    yield dqmDirectoryOut(sample)
                else:
                    yield dqmDirectoryOutUnfactorized(sample)

        merger.dqmDirectories_input = cms.vstring(list(merge_directories(subsamples)))

    # Update the plot sources in the plot jobs.  Note that we don't need to do
    # this for the merged samples, since we have replaced the HistAdder sources
    for scale in ['linear', 'log']:
        for analyzer in analyzers.keys():
            plotterModuleName = 'plot' + analyzer + '_' + scale
            print "Trying to update plotter:", plotterModuleName
            if hasattr(process, plotterModuleName):
                print " - found it in the process, modifying"
                plotterModuleProcesses = getattr(process, plotterModuleName).processes
                for sample in samplesToFactorize:
                    if hasattr(plotterModuleProcesses, sample):
                        print "Plotter is plotting %s, changing dqm directory in!" %  sample
                        getattr(plotterModuleProcesses, sample).dqmDirectory = \
                          cms.string("/harvested/%s_factorized" % sample)


def enableFactorization_makeAHtoMuTauPlots_grid(
    process,
    factorizationSequenceName = "loadAndFactorizeAHtoMuTauSamples",
    samplesToFactorize = [ 'InclusivePPmuX', 'PPmuXptGt20Mu10', 'PPmuXptGt20Mu15' ],
    relevantMergedSamples = [ 'qcdSum', ],
    mergedToRecoSampleDict = {},
    mergedSampleAdderModule = lambda sample, btag: 'addAHtoMuTau_%s_%s' % (btag, sample),
    dqmDirectoryOut =
    lambda sample, btag:'/harvested/%s_factorized/ahMuTauAnalyzer_%s/'% (sample, btag),
    dqmDirectoryOutUnfactorized =
    lambda sample, btag:'/harvested/%s/ahMuTauAnalyzer_%s/'% (sample, btag),
    dqmDirectoryTight =
    lambda sample, btag:'/harvested/%s/ahMuTauAnalyzer_%s_factorizedWithMuonIsolation/' % (sample, btag),
    dqmDirectoryLoose =
    lambda sample, btag:'/harvested/%s/ahMuTauAnalyzer_%s_factorizedWithoutMuonIsolation/' % (sample, btag),
    pyObjectLabel = ""):

    process.load("TauAnalysis.Configuration.analyzeAHtoMuTau_cfi")

    # define list of event selection criteria on "tight" muon isolation branch
    # of the analysis, **before** applying factorization of muon track + ECAL
    # isolation efficiencies
    evtSelAHtoMuTau_factorizedTight = [
        process.evtSelGenPhaseSpace,
        process.evtSelTrigger,
        process.evtSelDataQuality,
        process.evtSelPrimaryEventVertex,
        process.evtSelPrimaryEventVertexQuality,
        process.evtSelPrimaryEventVertexPosition,
        process.evtSelGlobalMuon,
        process.evtSelMuonEta,
        process.evtSelMuonPt,
        process.evtSelTauAntiOverlapWithMuonsVeto,
        process.evtSelTauEta,
        process.evtSelTauPt,
        process.evtSelMuonVbTfId,
        process.evtSelMuonPFRelIso
    ]

    # define list of event selection criteria on "loose" muon isolation branch
    # of the analysis, **after** applying factorization of muon track + ECAL
    # isolation efficiencies
    evtSelAHtoMuTau_factorizedLoose = [
        process.evtSelMuonTrkIP,
        process.evtSelTauLeadTrk,
        process.evtSelTauLeadTrkPt,
        process.evtSelTauTaNCdiscr,
        process.evtSelTauTrkIso,
        process.evtSelTauEcalIso,
        process.evtSelTauProng,
        process.evtSelTauCharge,
        process.evtSelTauMuonVeto,
        process.evtSelTauElectronVeto,
        process.evtSelDiTauCandidateForAHtoMuTauAntiOverlapVeto,
        process.evtSelDiTauCandidateForAHtoMuTauMt1MET,
        process.evtSelDiTauCandidateForAHtoMuTauPzetaDiff,
        process.evtSelDiMuPairZmumuHypothesisVetoByLooseIsolation
    ]

    # Make specialized cases for the w/ w/o btag cases
    evtSelAHtoMuTau_factorizedLoose_specialized = {
        'woBtag' : copy.deepcopy(evtSelAHtoMuTau_factorizedLoose),
        'wBtag' : copy.deepcopy(evtSelAHtoMuTau_factorizedLoose)
    }
    evtSelAHtoMuTau_factorizedLoose_specialized['woBtag'].append(
        process.evtSelNonCentralJetEt20bTag)
    evtSelAHtoMuTau_factorizedLoose_specialized['woBtag'].append(
        process.evtSelDiTauCandidateForAHtoMuTauZeroCharge)
    evtSelAHtoMuTau_factorizedLoose_specialized['wBtag'].append(
        process.evtSelCentralJetEt20)
    evtSelAHtoMuTau_factorizedLoose_specialized['wBtag'].append(
        process.evtSelCentralJetEt20bTag)
    evtSelAHtoMuTau_factorizedLoose_specialized['wBtag'].append(
        process.evtSelDiTauCandidateForAHtoMuTauZeroCharge)

    # defines names of MonitorElements used as numerator and denominator
    # to compute factorization scale-factor
    meNameAHtoMuTau_numerator = {
        'woBtag' : "evtSelDiTauCandidateForAHtoMuTauZeroCharge/passed_cumulative_numWeighted",
        'wBtag' : "evtSelDiTauCandidateForAHtoMuTauZeroCharge/passed_cumulative_numWeighted"
    }
    meNameAHtoMuTau_denominator = "evtSelMuonPFRelIso/processed_cumulative_numWeighted"

    # Loop over the samples and btag options and create sequences
    # for each of the factorization jobs and add them to the factorization
    # sequence
    factorizationSequence = getattr(process, factorizationSequenceName)
    for sample, bTagOption in [(sample, bTagOption)
                               for bTagOption in ['woBtag', 'wBtag']
                               for sample in samplesToFactorize]:
        ##if bTagOption == 'wBtag':
        ##    print "BTag analysis chain disabled - not modifying factorization"\
        ##            " sequence for", sample
        print "Adding sample:", sample, " to factorization:", \
                factorizationSequenceName
        new_factorization_sequence = composeFactorizationSequence(
            process = process,
            processName = sample + "_" + bTagOption + "_" + pyObjectLabel,
            dqmDirectoryIn_factorizedTightEvtSel = dqmDirectoryTight(
                sample, bTagOption),
            evtSel_factorizedTight = evtSelAHtoMuTau_factorizedTight,
            dqmDirectoryIn_factorizedLooseEvtSel = dqmDirectoryLoose(
                sample, bTagOption),
            evtSel_factorizedLoose = evtSelAHtoMuTau_factorizedLoose_specialized[bTagOption],
            meName_numerator = meNameAHtoMuTau_numerator[bTagOption],
            meName_denominator = meNameAHtoMuTau_denominator,
            dqmDirectoryOut = dqmDirectoryOut(sample, bTagOption),
            dropInputDirectories = False
        )
        new_factorization_seq_name = "scaleAHtoMuTau_%s_%s_%s" % (
            bTagOption, sample, pyObjectLabel)
        setattr(process, new_factorization_seq_name, new_factorization_sequence)
        factorizationSequence += new_factorization_sequence

    # Now update any of the relevant mergers
    for btag in ['woBtag', 'wBtag']:
        for mergedSample in relevantMergedSamples:
            # Get the module that is doing the merging, if it exists
            merger_name = "merge_%s_ahMuTauAnalyzer_%s" % (mergedSample, btag)
            if not hasattr(process.mergeSamplesAHtoMuTau, merger_name):
                print "factorizationTools: Expected to update ",\
                        merger_name, "but it's not in the process! skipping.."
                continue
            merger = getattr(process.mergeSamplesAHtoMuTau, merger_name)

            # Get the subsamples associated with this merged sample
            subsamples = mergedToRecoSampleDict[mergedSample]['samples']
            # Set the adder to use our new factorized inputs
            def merge_directories(_list):
                for sample in _list:
                    if sample in samplesToFactorize:
                        yield dqmDirectoryOut(sample, btag)
                    else:
                        yield dqmDirectoryOutUnfactorized(sample, btag)

            merger.dqmDirectories_input = cms.vstring(list(merge_directories(subsamples)))

    # Update the plot sources in the plot jobs.  Note that we don't need to do
    # this for the merged samples, since we have replaced the HistAdder sources
    for plotterModuleName in [ 'plotAHtoMuTau_woBtag_log', 'plotAHtoMuTau_woBtag_linear',
                               'plotAHtoMuTau_wBtag_log',  'plotAHtoMuTau_wBtag_linear' ]:
        if hasattr(process, plotterModuleName):
            plotterModuleProcesses = getattr(process, plotterModuleName).processes
            for sample in samplesToFactorize:
                if hasattr(plotterModuleProcesses, sample):
                    getattr(plotterModuleProcesses, sample).dqmDirectory = \
                      cms.string("/harvested/%s_factorized" % sample)

