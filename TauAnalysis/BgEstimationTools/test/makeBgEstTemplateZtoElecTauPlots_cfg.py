import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.DQMTools.plotterStyleDefinitions_cfi import *
from TauAnalysis.BgEstimationTools.templateHistDefinitions_cfi import \
  drawJobTemplateHist, drawJobTemplateHistIntegrated, drawJobAnalysisHistData, drawJobAnalysisHistMC, drawJobAnalysisHistZmumuEmbedding, \
  plotBgEstData, plotBgEstMC_pure, plotBgEstMC_smSum, plotAnalysisMC_pure, plotAnalysis_ZmumuEmbedding, \
  plotHistZtoElecTau, configurePlotZtoElecTauIntegrated, plotZtoElecTauStacked_template, plotHistZtoElecTauStacked
from TauAnalysis.BgEstimationTools.tools.drawTemplateHistConfigurator import drawTemplateHistConfigurator
from TauAnalysis.DQMTools.tools.drawJobConfigurator import *
from TauAnalysis.Configuration.userRegistry import getHarvestingFilePath, getJobId

process = cms.Process('makeBgEstTemplateZtoElecTauPlots')

process.DQMStore = cms.Service("DQMStore")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

dqmDirectories_processes = {
    'Ztautau' : {
	'template' : 'ZtautauSum',
        'analysis' : 'ZtautauSum'    },
    'Ztautau_from_ZeeEmbedding' : {
	'template' : 'ZtautauSum',
        'analysis' : 'ZtautauSum'
    },    
    'Zee' : {
	'template' : 'ZeeSum',
        'analysis' : 'ZeeSum'
    },
    'WplusJets' : {
	'template' : 'WplusJets_madgraph_Pat',
        'analysis' : 'WplusJets_madgraph_Pat'
    },
    'QCD' : {
	'template' : 'qcdSum',
        'analysis' : 'qcdSum'
    },
    'PhotonPlusJets' : {
	'template' : 'photonPlusJetsSum',
	'analysis' : 'photonPlusJetsSum'
    },
    
    'smSum' : {
	'template' : 'smSum',
         'analysis' : 'smSum'
    },
    'Data' : {
	'template' : 'data',
        'analysis' : 'data'
    }
}

dqmDirectory_analysis = 'BgEstTemplateAnalyzer_ZtautauEnriched/afterDiTauCandidateZeeHypothesisVetoBgEstZtautauEnriched/'

dqmDirectories_control = {
	'ZtautauEnriched'           : dqmDirectory_analysis,
	#'ZeeJetMisIdEnriched'     : 'BgEstTemplateAnalyzer_ZeeJetMisIdEnriched/afterDiTauCandidateZeeHypothesisVetoBgEstZeeJetMisIdEnriched/',
	'ZeeElectronMisIdEnriched' : 'BgEstTemplateAnalyzer_ZeeElectronMisIdEnriched/afterDiTauCandidateZeeHypothesisVetoBgEstZeeElectronMisIdEnriched/',
	'WplusJetsEnriched'        : 'BgEstTemplateAnalyzer_WplusJetsEnriched/afterDiTauCandidateZeeHypothesisVetoBgEstWplusJetsEnriched/',
	#'TTplusJetsEnriched'      : 'BgEstTemplateAnalyzer_TTplusJetsEnriched/afterDiTauCandidateZeeHypothesisVetoBgEstTTplusJetsEnriched/',
	'PhotonPlusJetsEnriched'   : 'BgEstTemplateAnalyzer_PhotonPlusJetsEnriched/afterDiTauCandidateZeeHypothesisVetoBgEstPhotonPlusJetsEnriched/',
	'QCDenriched'              : 'BgEstTemplateAnalyzer_QCDenriched/afterDiTauCandidateZeeHypothesisVetoBgEstQCDenriched/'
}

plotTitles = {
    'visMass'   : "M_{vis}(Electron + Tau)",
    'SVfitMass' : "M(Electron + Tau), SVfit method",
    'pfMEt'     : "Particle-Flow MEt",
    'visPt'     : "P_{T}(Electron + Tau)",
    'electronPt'    : "Electron P_{T}", 
    'electronEta'   : "Electron #eta", 
    'electronPhi'   : "Electron #phi", 
    'tauPt'     : "Tau--jet P_{T}", 
    'tauEta'    : "Tau--jet #eta", 
    'tauPhi'    : "Tau--jet #phi" 
}

xAxes = {
    'visMass'   : 'Mass',
    'SVfitMass' : 'Mass',
    'pfMEt'     : 'Pt',
    'visPt'     : 'Pt',
    'electronPt'    : 'Pt',
    'electronEta'   : 'Eta', 
    'electronPhi'   : 'Phi', 
    'tauPt'     : 'Pt',
    'tauEta'    : 'Eta', 
    'tauPhi'    : 'Phi'
}

rebinning_control = {
    'ZtautauEnriched' : {
	'visMass'   : 1,
#        'SVfitMass' : 1
    }, 
    'ZeeJetMisIdEnriched' : {
        'visMass'   : 1,
        #'SVfitMass' : 5,
        #'pfMEt'     : 3,
	#'visPt'     : 2,
        #'electronPt'    : 3, 
        #'electronEta'   : 3, 
        #'electronPhi'   : 3, 
        #'tauPt'     : 3, 
        #'tauEta'    : 3, 
        #'tauPhi'    : 3
    },
    'ZeeElectronMisIdEnriched' : {
        'visMass'   : 1,
#        'SVfitMass' : 1,
##         'pfMEt'     : 3,
## 	'visPt'     : 2,
##         'electronPt'    : 3, 
##         'electronEta'   : 3, 
##         'electronPhi'   : 3, 
##         'tauPt'     : 3, 
##         'tauEta'    : 3, 
##         'tauPhi'    : 3
    },
    'WplusJetsEnriched' : {
        'visMass'   : 1,
#        'SVfitMass' : 2,
##         'pfMEt'     : 3,
## 	'visPt'     : 2,
##         'electronPt'    : 3, 
##         'electronEta'   : 3, 
##         'electronPhi'   : 3, 
##         'tauPt'     : 3, 
##         'tauEta'    : 3, 
##         'tauPhi'    : 3
    },
##     'TTplusJetsEnriched' : {
##         'visMass'   : 5,
##         'SVfitMass' : 5,
##         'pfMEt'     : 5,
##     	'visPt'     : 5,
##         'electronPt'    : 3, 
##         'electronEta'   : 3, 
##         'electronPhi'   : 3, 
##         'tauPt'     : 3, 
##         'tauEta'    : 3, 
##         'tauPhi'    : 3	
##     },
    'QCDenriched'   : {
        'visMass'   : 1,
#        'SVfitMass' : 2,
##         'pfMEt'     : 3,
## 	'visPt'     : 2,
##         'electronPt'    : 3, 
##         'electronEta'   : 3, 
##         'electronPhi'   : 3, 
##         'tauPt'     : 3, 
##         'tauEta'    : 3, 
##         'tauPhi'    : 3
    },
    'PhotonPlusJetsEnriched'   : {
        'visMass'   : 1,
#        'SVfitMass' : 2,
##         'pfMEt'     : 3,
## 	'visPt'     : 2,
##         'electronPt'    : 3, 
##         'electronEta'   : 3, 
##         'electronPhi'   : 3, 
##         'tauPt'     : 3, 
##         'tauEta'    : 3, 
##         'tauPhi'    : 3
    }    
}    

pureProcessBgEnrichedSelections = {
	#'ZeeJetMisIdEnriched'  : 'Zee',
    'ZeeElectronMisIdEnriched' : 'Zee',
    'WplusJetsEnriched'        : 'WplusJets',
    'PhotonPlusJetsEnriched'   : 'PhotonPlusJets',
    'QCDenriched'              : 'QCD'
}

meNames = {
    'visMass'   : 'DiTauCandidateQuantities/VisMass',
#    'SVfitMass' : 'DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass',
##     'pfMEt'     : 'PFMEtQuantities/MEtPt',
##     'visPt'     : 'DiTauCandidateQuantities/VisPt',
##     'electronPt'    : 'ElectronQuantities/ElectronPt',
##     'electronEta'   : 'ElectronQuantities/ElectronEta', 
##     'electronPhi'   : 'ElectronQuantities/ElectronPhi', 
##     'tauPt'     : 'TauQuantities/TauPt',
##     'tauEta'    : 'TauQuantities/TauEta',
##     'tauPhi'    : 'TauQuantities/TauPhi'
}

#--------------------------------------------------------------------------------
# build dictionary, mapping
#  [process][selection]['template'/'analysis'][distName][distType] --> meName
#--------------------------------------------------------------------------------

meNameMapping = {}

for processName in dqmDirectories_processes.keys():
	#print "processName: ",processName
	meNameMapping[processName] = {}
	for selectionName in dqmDirectories_control.keys():
		#print "\tselectionName: ",selectionName
		if processName == 'Ztautau_from_ZeeEmbedding' and selectionName != 'ZtautauEnriched':
			continue
		meNameMapping[processName][selectionName] = {}
		for selectionType in [ "analysis", "template" ]:
			#print "\t\tselectionType: ",selectionType
			if selectionName == 'ZtautauEnriched' and selectionType != 'analysis':
				continue
			if processName == 'diBoson' and selectionType == 'analysis':
				continue
			meNameMapping[processName][selectionName][selectionType] = {}
			for distName in meNames.keys():
				#print "\t\t\tdistName: ",distName
				# check if histogram exists;
				# CV: a histogram is assumed to exist it has an entry in the rebinning options dictionary
				if rebinning_control[selectionName].get(distName) is not None:
					meNameMapping[processName][selectionName][selectionType][distName] = {}
					meName_distribution = None
					if selectionType == 'analysis':
						meName_distribution = \
								'/analysis/harvested/' + dqmDirectories_processes[processName]['analysis'] \
								+ '/' + dqmDirectory_analysis + meNames[distName]
					else:
						meName_distribution = \
								'/template/harvested/' + dqmDirectories_processes[processName]['template'] \
								+ '/' + dqmDirectories_control[selectionName] + meNames[distName]
					meNameMapping[processName][selectionName][selectionType][distName]['distribution'] = meName_distribution

					#print "\t\t\t\tmeName_distribution: ", meName_distribution
           
# overwrite meNames for mass distributions in ZeeJetMisIdEnriched and ZeeElectronMisIdEnriched selections
# (require generator level matching to quark/gluon jets and electrons on analysis level)
#meNameMapping['Zee']['ZeeJetMisIdEnriched']['analysis']['visMass']['distribution'] = \
		#  '/analysis/DQMData/harvested/' + dqmDirectories_processes['Zee']['analysis'] \
		# + '/' + dqmDirectory_analysis + 'DiTauCandidateQuantities/VisMassGenLeg2Jet'
## meNameMapping['Zee']['ZeeElectronMisIdEnriched']['analysis']['visMass']['distribution'] = \
##   '/analysis/DQMData/harvested/' + dqmDirectories_processes['Zee']['analysis'] \
##  + '/' + dqmDirectory_analysis + 'DiTauCandidateQuantities/VisMassGenLeg2Electron'
##meNameMapping['Zmumu']['ZmumuJetMisIdEnriched']['analysis']['SVfitMass']['distribution'] = \
##  '/analysis/DQMData/harvested/' + dqmDirectories_processes['Zmumu']['analysis'] \
## + '/' + dqmDirectory_analysis + 'DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/MassGenLeg2Jet'
##meNameMapping['Zmumu']['ZmumuElectronMisIdEnriched']['analysis']['SVfitMass']['distribution'] = \
##  '/analysis/DQMData/harvested/' + dqmDirectories_processes['Zmumu']['analysis'] \
## + '/' + dqmDirectory_analysis + 'DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/MassGenLeg2Electron'

for processName, meNameMapping_process in meNameMapping.items():
	for selectionName, meNameMapping_selectionName in meNameMapping_process.items():
		for selectionType, meNameMapping_selectionType in meNameMapping_selectionName.items():
			for distName, meNameMapping_distName in meNameMapping_selectionType.items():
				meName_distribution = meNameMapping_distName['distribution']
				meNameMapping_distName['rebinned']   = meName_distribution + "Rebinned_%s" % selectionName
				meNameMapping_distName['template']   = meName_distribution + "Shape_%s" % selectionName
				meNameMapping_distName['integrated'] = meName_distribution + "IntegratedShape_%s" % selectionName

#print("meNameMapping:")
#print(meNameMapping)

#--------------------------------------------------------------------------------
# load template histogram for Z --> tau+ tau- signal events,
# produced by MCEmbeddingTools from Z --> mu+ mu- events selected in (pseudo)data
#--------------------------------------------------------------------------------

process.loadTemplateHistZtoElecTau_Ztautau = cms.EDAnalyzer("DQMFileLoader",
    Ztautau = cms.PSet(
        inputFileNames = cms.vstring(
            'file:/afs/crc.nd.edu/user/l/lantonel/CMSSW_4_1_3/src/TauAnalysis/BgEstimationTools/test/Run01/plotsZtoElecTau_bgEstTemplate_all.root'
        ),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string(
#            '/DQMData/harvested/ZtautauSum' + '/' + dqmDirectories_control['ZtautauEnriched']
            '/analysis/harvested/ZtautauSum' + '/' + dqmDirectories_control['ZtautauEnriched']	    
        )
    )
)

#--------------------------------------------------------------------------------
# load template histogram of for different types of background events,
# obtained from background enriched control samples
#--------------------------------------------------------------------------------

process.loadTemplateHistZtoElecTau = cms.EDAnalyzer("DQMFileLoader",
    Ztautau = cms.PSet(
        inputFileNames = cms.vstring(
	     'file:/afs/crc.nd.edu/user/l/lantonel/CMSSW_4_1_3/src/TauAnalysis/BgEstimationTools/test/Run01/plotsZtoElecTau_bgEstTemplate_all.root'
        ),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('/template')
    )
)

#--------------------------------------------------------------------------------
# load distributions expected for different signal/background processes
# and observed in (pseudo)data in events passing final analysis selection criteria
#--------------------------------------------------------------------------------

process.loadAnalysisHistZtoElecTau = cms.EDAnalyzer("DQMFileLoader",
    data = cms.PSet(
        inputFileNames = cms.vstring(
	     'file:/afs/crc.nd.edu/user/l/lantonel/CMSSW_4_1_3/src/TauAnalysis/BgEstimationTools/test/Run01/plotsZtoElecTau_bgEstTemplate_all.root'
        ),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('/analysis')
    )
)

#--------------------------------------------------------------------------------
# combine two bins into one,
# in order to mitigate effect of small event statistics
#--------------------------------------------------------------------------------

process.rebinHistZtoElecTau = cms.Sequence()

jobsHistRebinning = {} # keys = [selectionName][distName]

for processName, meNameMapping_process in meNameMapping.items():
    for selectionName, meNameMapping_selectionName in meNameMapping_process.items():
        if jobsHistRebinning.get(selectionName) is None:
	    jobsHistRebinning[selectionName] = {}
        for selectionType, meNameMapping_selectionType in meNameMapping_selectionName.items():
            for distName, meNameMapping_distName in meNameMapping_selectionType.items():
                if jobsHistRebinning[selectionName].get(distName) is None:
	            jobsHistRebinning[selectionName][distName] = {}
                    jobsHistRebinning[selectionName][distName]['config'] = []
                    jobsHistRebinning[selectionName][distName]['binning'] = cms.PSet(
                        x = cms.PSet(
                           combineBins = cms.uint32(rebinning_control[selectionName][distName])
                        )
                    )
                jobsHistRebinning[selectionName][distName]['config'].append(cms.PSet(
                    meName_original = cms.string(meNameMapping_distName['distribution']),
                    meName_rebinned = cms.string(meNameMapping_distName['rebinned'])
                ))

for selectionName, jobsHistRebinning_selectionName in jobsHistRebinning.items():
    for distName, jobsHistRebinning_distName in jobsHistRebinning_selectionName.items():
        rebinningModule = cms.EDAnalyzer("DQMHistRebinner",
	    config = cms.VPSet(jobsHistRebinning[selectionName][distName]['config']),
            binning = jobsHistRebinning[selectionName][distName]['binning']
        )
	rebinningModuleName = "rebinTemplateHistZtoElecTau%s%s" % (selectionName, distName)
        setattr(process, rebinningModuleName, rebinningModule)

        process.rebinHistZtoElecTau += rebinningModule

#--------------------------------------------------------------------------------
# correct for "bias" of W + jets template histograms
# caused by cuts on Mt(electron + tau-jet) transverse mass and CDF (Pzeta - 1.5*PzetaVis) variable
# (cf. CMS AN-2010/088)
#--------------------------------------------------------------------------------

#meNameMapping['Data']['WplusJetsEnriched']['template']['visMass']['corrDistribution'] = \
#  meNameMapping['Data']['WplusJetsEnriched']['template']['visMass']['distribution'] + "Corrected"
#meNameMapping['Data']['WplusJetsEnriched']['template']['SVfitMass']['corrDistribution'] = \
#  meNameMapping['Data']['WplusJetsEnriched']['template']['SVfitMass']['distribution'] + "Corrected"

## process.correctTemplateHistZtoElecTau = cms.EDAnalyzer("DQMHistBiasCorrection",
##     config = cms.VPSet(
##         cms.PSet(
##             meName_uncorrected = cms.string(meNameMapping['Data']['WplusJetsEnriched']['template']['visMass']['rebinned']),
##             meName_corrected = cms.string(meNameMapping['Data']['WplusJetsEnriched']['template']['visMass']['corrDistribution']),
##             meName_corrNumerator = cms.string(meNameMapping['WplusJets']['WplusJetsEnriched']['analysis']['visMass']['rebinned']),
##             meName_corrDenominator = cms.string(meNameMapping['WplusJets']['WplusJetsEnriched']['template']['visMass']['rebinned'])
##         ),
##         cms.PSet(
##             meName_uncorrected = cms.string(meNameMapping['Data']['WplusJetsEnriched']['template']['SVfitMass']['rebinned']),
##             meName_corrected = cms.string(meNameMapping['Data']['WplusJetsEnriched']['template']['SVfitMass']['corrDistribution']),
##             meName_corrNumerator = cms.string(meNameMapping['WplusJets']['WplusJetsEnriched']['analysis']['SVfitMass']['rebinned']),
##             meName_corrDenominator = cms.string(meNameMapping['WplusJets']['WplusJetsEnriched']['template']['SVfitMass']['rebinned'])
##         )
##     )
## ) 


#--------------------------------------------------------------------------------
# normalize to unit area distribution of visible electron + tau-jet mass for Z --> tau+ tau- signal
# and different types of background events
#--------------------------------------------------------------------------------

## meNameMapping['Data']['WplusJetsEnriched']['template']['visMass']['corrTemplate'] = \
##   meNameMapping['Data']['WplusJetsEnriched']['template']['visMass']['distribution'] + "CorrectedShape"
## meNameMapping['Data']['WplusJetsEnriched']['template']['SVfitMass']['corrTemplate'] = \
##   meNameMapping['Data']['WplusJetsEnriched']['template']['SVfitMass']['distribution'] + "CorrectedShape"

jobsHistNormalization = []

for processName, meNameMapping_process in meNameMapping.items():
    for selectionName, meNameMapping_selectionName in meNameMapping_process.items():
        for selectionType, meNameMapping_selectionType in meNameMapping_selectionName.items():
            for distName, meNameMapping_distName in meNameMapping_selectionType.items():
                jobsHistNormalization.append(cms.PSet(
                    meName_input = cms.string(meNameMapping_distName['rebinned']),
                    meName_output = cms.string(meNameMapping_distName['template'])
                ))

## for distName in [ "visMass", "SVfitMass"]:
##     jobsHistNormalization.append(cms.PSet(
##         meName_input = cms.string(meNameMapping['Data']['WplusJetsEnriched']['template'][distName]['corrDistribution']),
##         meName_output = cms.string(meNameMapping['Data']['WplusJetsEnriched']['template'][distName]['corrTemplate'])
##     ))
    
process.normalizeHistZtoElecTau = cms.EDAnalyzer("DQMHistNormalizer",
    config = cms.VPSet(jobsHistNormalization),
    norm = cms.double(1.)
)

#--------------------------------------------------------------------------------
# produce cumulative distributions
#--------------------------------------------------------------------------------

jobsHistIntegration = []

for processName, meNameMapping_process in meNameMapping.items():
    for selectionName, meNameMapping_selectionName in meNameMapping_process.items():
        for selectionType, meNameMapping_selectionType in meNameMapping_selectionName.items():
            for distName, meNameMapping_distName in meNameMapping_selectionType.items():
	        jobsHistIntegration.append(cms.PSet(
                    meName_input = cms.string(meNameMapping_distName['template']),
                    meName_output = cms.string(meNameMapping_distName['integrated']),
                    integrateFrom = cms.string("right")
                ))
    
process.integrateHistZtoElecTau = cms.EDAnalyzer("DQMHistIntegrator",
    config = cms.VPSet(jobsHistIntegration)
)

#--------------------------------------------------------------------------------
# compute Kolmogorov-Smirnov probabilities for distributions observed in Data
# to agree with Monte Carlo predictions in the background enriched regions
#--------------------------------------------------------------------------------

## jobsKolmogorovTest = []

## meNameMapping_mc = meNameMapping['smSum']
## meNameMapping_data = meNameMapping['Data']
## for selectionName, meNameMapping_selectionName_data in meNameMapping_data.items():
##     selectionType = 'template'
##     if meNameMapping_selectionName_data.get(selectionType) is not None:
##         meNameMapping_selectionType_data = meNameMapping_selectionName_data[selectionType]
##         for distName in [ "visMass", "SVfitMass", "pfMEt" ]:
##             if meNameMapping_selectionType_data.get(distName) is not None:
##                 meNameMapping_distName_data = meNameMapping_selectionType_data[distName]
                   
##                 meName_distribution_data = meNameMapping_distName_data['distribution']
## 	        meName_compatibility = meName_distribution + "Compatibility"
##                 meNameMapping_distName_data['compatibility'] = meName_compatibility

##                 jobsKolmogorovTest.append(cms.PSet(
##                     meName_test = cms.string(meNameMapping_distName_data['template']),
##                     meName_reference = cms.string(meNameMapping_mc[selectionName][selectionType][distName]['template']),
##                     meName_compatibility = cms.string(meName_compatibility)
##     	        ))
    
## process.compKolmogorovProbZtoElecTau = cms.EDAnalyzer("DQMHistKolmogorovTest",
##     config = cms.VPSet(jobsKolmogorovTest)
## )

#--------------------------------------------------------------------------------
# plot template histograms obtained from Monte Carlo
# compared to the shapes determined by background enriched regions in (pseudo)Data
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorZtoElecTau = drawTemplateHistConfigurator(
    template = drawJobTemplateHist
)

drawAnalysisHistConfiguratorZtoElecTauData = drawTemplateHistConfigurator(
    template = drawJobAnalysisHistData
)

drawAnalysisHistConfiguratorZtoElecTauMC = drawTemplateHistConfigurator(
    template = drawJobAnalysisHistMC
)

drawAnalysisHistConfiguratorZtoElecTauZmumuEmbedding = drawTemplateHistConfigurator(
    template = drawJobAnalysisHistZmumuEmbedding
)

process.plotTemplateHistZtoElecTauIntegrated = cms.Sequence()

meNameMapping_mcSMsum = meNameMapping['smSum']
meNameMapping_data = meNameMapping['Data']
for selectionName, meNameMapping_selectionName_data in meNameMapping_data.items():
    print "selectionName: ",selectionName
    if pureProcessBgEnrichedSelections.get(selectionName) is not None:
        meNameMapping_mcPureProcess = meNameMapping[pureProcessBgEnrichedSelections[selectionName]]
        selectionType = 'template'
        if meNameMapping_selectionName_data.get(selectionType) is not None:
            meNameMapping_selectionType_data = meNameMapping_selectionName_data[selectionType]
	    for distName, meNameMapping_distName_data in meNameMapping_selectionType_data.items():
                # configure histograms for MC(SM|bgEst) vs. MC(pure|bgEst) vs. Data|bgEst comparisson
                # (showing Data to Monte Carlo comparisson and effect of impurities in background enriched control regions)
		name_check = ("%s_%s" % (selectionName, distName))
		print "name_check: ", name_check
		#print meNameMapping_mcPureProcess[selectionName][selectionType][distName]['template']
		#print meNameMapping_distName_data['template']
                drawTemplateHistConfiguratorZtoElecTau.add(
                    meNames = [
                        meNameMapping_mcSMsum[selectionName][selectionType][distName]['template'],
                        meNameMapping_mcPureProcess[selectionName][selectionType][distName]['template'],
                        meNameMapping_distName_data['template']
                    ],
                    name = ("%s_%s" % (selectionName, distName)),
                    title = ("%s: %s" % (selectionName, plotTitles[distName]))
                )

                if meNameMapping_distName_data.get('compatibility') is not None:
                    plotHistZtoElecTauIntegrated = configurePlotZtoElecTauIntegrated(
    	                meNameMapping_distName_data['template'],
    	                meNameMapping_mcSMsum[selectionName][selectionType][distName]['template'],
                        meNameMapping_distName_data['compatibility'],
                        "%s_%s" % (selectionName, distName),
	                "%s: %s" % (selectionName, plotTitles[distName]),
	                "plotBgEstTemplateZtoElecTauIntegrated_#PLOT#.pdf"
                    )
                    plotHistZtoElecTauIntegratedName = "plotHistZtoElecTauIntegrated%s%s" % (selectionName, distName)
                    setattr(process, plotHistZtoElecTauIntegratedName, plotHistZtoElecTauIntegrated)
	            process.plotTemplateHistZtoElecTauIntegrated += plotHistZtoElecTauIntegrated

                # configure histograms for Data|bgEst vs. MC(pure|analysis) comparisson
                # (showing possible bias on mass shape caused by event selection criteria applied in analysis)
                drawAnalysisHistConfiguratorZtoElecTauData.add(
                    meNames = [
                        meNameMapping_distName_data['template'],
                        meNameMapping_mcPureProcess[selectionName]['analysis'][distName]['template']
                    ],
                    name = ("%s_%s" % (selectionName, distName)),
                    title = ("%s: %s" % (selectionName, plotTitles[distName]))
                )

                # configure histograms for MC(SM|bgEst) vs. MC(pure|bgEst)
                # (showing effect of impurities in background enriched control regions)
                drawAnalysisHistConfiguratorZtoElecTauMC.add(
                    meNames = [
                        meNameMapping_mcSMsum[selectionName][selectionType][distName]['template'],
                        meNameMapping_mcPureProcess[selectionName]['analysis'][distName]['template']
                    ],
                    name = ("%s_%s" % (selectionName, distName)),
                    title = ("%s: %s" % (selectionName, plotTitles[distName]))
                )

## for distName in [ "visMass", "SVfitMass"]:
##     drawAnalysisHistConfiguratorZtoElecTauData.add(
##         meNames = [
##             meNameMapping_data['WplusJetsEnriched']['template'][distName]['corrTemplate'],
##             meNameMapping['WplusJets']['WplusJetsEnriched']['analysis'][distName]['template']
##         ],
##         name = ("WplusJetsEnriched_%s_corrected" % distName),
##         title = ("WplusJetsEnriched: %s, corrected" % plotTitles[distName])
##     )

	#drawAnalysisHistConfiguratorZtoElecTauZmumuEmbedding.add(
	#    meNames = [
	#        meNameMapping['Ztautau']['ZmumuEmbedding']['analysis'][distName]['template'],
	#        meNameMapping['Ztautau_from_ZmumuEmbedding']['ZmumuEmbedding']['analysis'][distName]['template']
	#    ],
	#    name = ("ZmumuEmbedding_%s" % distName),
	#    title = ("ZmumuEmbedding: %s" % plotTitles[distName])
	#)

process.plotTemplateHistZtoElecTau = plotHistZtoElecTau.clone(
    drawJobs = drawTemplateHistConfiguratorZtoElecTau.configure(),
    indOutputFileName = cms.string('plotBgEstTemplateZtoElecTau_#PLOT#.pdf')
)

process.plotAnalysisHistZtoElecTauData = plotHistZtoElecTau.clone(
    drawJobs = drawAnalysisHistConfiguratorZtoElecTauData.configure(),
    indOutputFileName = cms.string('plotBgEstTemplateData_vs_AnalysisZtoElecTau_#PLOT#.pdf')
)

process.plotAnalysisHistZtoElecTauMC = plotHistZtoElecTau.clone(
    drawJobs = drawAnalysisHistConfiguratorZtoElecTauMC.configure(),
    indOutputFileName = cms.string('plotBgEstTemplateMC_vs_AnalysisZtoElecTau_#PLOT#.pdf')
)

#process.plotAnalysisHistZtoElecTauZmumuEmbedding = plotHistZtoElecTau.clone(
#    drawJobs = drawAnalysisHistConfiguratorZtoElecTauZmumuEmbedding.configure(),
#    indOutputFileName = cms.string('plotZmumuEmbedding_vs_AnalysisZtoElecTau_#PLOT#.pdf')
#)

process.plotAnalysisHistZtoElecTauStacked = cms.Sequence()

meNameMapping_data = meNameMapping['Data']
for selectionName, meNameMapping_selectionName_data in meNameMapping_data.items():
    selectionType = 'template'
    if meNameMapping_selectionName_data.get(selectionType) is not None:
        meNameMapping_selectionType_data = meNameMapping_selectionName_data[selectionType]

        drawJobConfiguratorZtoElecTauStacked = drawJobConfigurator(
            template = plotZtoElecTauStacked_template,
            dqmDirectory = '#PROCESSDIR#' + '/' + dqmDirectories_control[selectionName]
        )

        plots = []

        for distName in meNameMapping_selectionType_data.keys():
	    plots.append(drawJobConfigEntry(
                meName = meNames[distName] + "Rebinned_%s" % selectionName,
                title = plotTitles[distName],
                xAxis = xAxes[distName],
                name = ("%s" % distName)
            ))

        drawJobConfiguratorZtoElecTauStacked.add(plots = plots)

        dqmHistPlotterModule = plotHistZtoElecTauStacked.clone(
            processes = cms.PSet(
                Ztautau = plotHistZtoElecTauStacked.processes.Ztautau.clone(
                    dqmDirectory = cms.string('/template/harvested/' + dqmDirectories_processes['Ztautau']['template'])
                ),
                Zee = plotHistZtoElecTauStacked.processes.Zee.clone(
                    dqmDirectory = cms.string('/template/harvested/' + dqmDirectories_processes['Zee']['template'])
                ),
                WplusJets = plotHistZtoElecTauStacked.processes.WplusJets.clone(
                    dqmDirectory = cms.string('/template/harvested/' + dqmDirectories_processes['WplusJets']['template'])
                ),
				#TTplusJets = plotHistZtoElecTauStacked.processes.TTplusJets.clone(
				#    dqmDirectory = cms.string('/template/harvested/' + dqmDirectories_processes['TTplusJets']['template'])
				#),
				#diBoson = plotHistZtoElecTauStacked.processes.diBoson.clone(
				#    dqmDirectory = cms.string('/template/harvested/' + dqmDirectories_processes['diBoson']['template'])
				#),
                QCD = plotHistZtoElecTauStacked.processes.QCD.clone(
                    dqmDirectory = cms.string('/template/harvested/' + dqmDirectories_processes['QCD']['template'])
                ),
		PhotonPlusJets = plotHistZtoElecTauStacked.processes.PhotonPlusJets.clone(
		    dqmDirectory = cms.string('/template/harvested/' + dqmDirectories_processes['PhotonPlusJets']['template'])
		),
                Data = plotHistZtoElecTauStacked.processes.Data.clone(
                    dqmDirectory = cms.string('/template/harvested/' + dqmDirectories_processes['Data']['template'])
                )
            ),
            drawJobs = drawJobConfiguratorZtoElecTauStacked.configure(),
            indOutputFileName = cms.string('bgEstControlZtoElecTau_%s_#PLOT#.pdf' % selectionName)
        )
        dqmHistPlotterModuleName = "plotZtoElecTauStacked%s" % selectionName
        setattr(process, dqmHistPlotterModuleName, dqmHistPlotterModule)

        process.plotAnalysisHistZtoElecTauStacked += dqmHistPlotterModule
    
process.saveBgEstTemplateHistZtoElecTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string(
        getHarvestingFilePath('ZtoElecTau_bgEstTemplate') + '/' + 'bgEstTemplateHistZtoElecTau_skimmed.root'
        #'/data1/veelken/CMSSW_3_8_x/plots/ZtoElecTau_bgEstTemplate/2011Feb01_HPSloose/bgEstTemplateHistZtoElecTau_skimmed.root'
	#'skim.root'
    ),
    outputCommands = cms.vstring(
#        'drop *',
        'keep *',
	
##         'keep /analysis/harvested/*/zElecTauAnalyzerOS/afterEvtSelDiTauCandidateForElecTauZeroCharge/*',
##         'keep /analysis/harvested/ZtoElecTau_from_ZmumuEmbedding/*',                                      
##         'keep /template/harvested/ZtautauSum/*',
##         'keep /template/harvested/Ztautau_from_ZmumuEmbedding/*',
##         'keep /template/harvested/Zmumu/*',
##         'keep /template/harvested/Zmumu_powheg/*',
##         'keep /template/harvested/qcdSum/*',
##         'keep /template/harvested/WplusJets_madgraph/*',
##         'keep /template/harvested/WplusJetsSum/*',
##         'keep /template/harvested/TTplusJets/*',
##         'keep /template/harvested/TTplusJets_madgraph/*',
##         'keep /template/harvested/VVsum/*',
##         'keep /template/harvested/data/*'
    )
)

process.dumpDQMStore = cms.EDAnalyzer("DQMStoreDump")

process.p = cms.Path(
	#process.loadTemplateHistZtoElecTau_Ztautau
	process.loadTemplateHistZtoElecTau
	+ process.loadAnalysisHistZtoElecTau
   #+ process.dumpDQMStore
   + process.rebinHistZtoElecTau
   #+ process.correctTemplateHistZtoElecTau
   + process.normalizeHistZtoElecTau
   + process.integrateHistZtoElecTau 
   #+ process.compKolmogorovProbZtoElecTau
   ##+ process.dumpDQMStore 
    + process.plotTemplateHistZtoElecTau
   #+ process.plotTemplateHistZtoElecTauIntegrated
   + process.plotAnalysisHistZtoElecTauData
   + process.plotAnalysisHistZtoElecTauMC
   #+ process.plotAnalysisHistZtoElecTauZmumuEmbedding    
   + process.plotAnalysisHistZtoElecTauStacked
   + process.saveBgEstTemplateHistZtoElecTau 
)

# print-out all python configuration parameter information
#print process.dumpPython()


  
