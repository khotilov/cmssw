import FWCore.ParameterSet.Config as cms
import copy

def modifyPlots(process):
	TZtZeWGQ = cms.vstring(
		'Data',
		'TTplusJets', 
		'Ztautau',
		'Zee', 
		'WplusJets', 
		'gammaPlusJetsSum',
		'qcdEMenrichedSum',
		'qcdBCtoESum'
	)
	ZtTZeWGQ = cms.vstring(
		'Data',
		'Ztautau',
		'TTplusJets', 
		'Zee', 
		'WplusJets', 
		'gammaPlusJetsSum',
		'qcdEMenrichedSum',
		'qcdBCtoESum'
	)
	TZtWGQZe = cms.vstring(
		'Data',
		'TTplusJets', 
		'Ztautau',
		'WplusJets', 
		'gammaPlusJetsSum',
		'qcdEMenrichedSum',
		'qcdBCtoESum',
		'Zee'
	)
	ZtTWGQZe = cms.vstring(
		'Data',
		'Ztautau',
		'TTplusJets', 
		'WplusJets', 
		'gammaPlusJetsSum',
		'qcdEMenrichedSum',
		'qcdBCtoESum',
		'Zee'
	)
	TZtGWZeQ = cms.vstring(
		'Data',
		'TTplusJets', 
		'Ztautau',
		'gammaPlusJetsSum',
		'WplusJets', 
		'Zee', 
		'qcdEMenrichedSum',
		'qcdBCtoESum'
	)
	TGQZtWZe = cms.vstring(
		'Data',
		'TTplusJets', 
		'gammaPlusJetsSum',
		'qcdEMenrichedSum',
		'qcdBCtoESum',
		'Ztautau',
		'WplusJets', 
		'Zee'
	)
	ZtTQGZeW = cms.vstring(
		'Data',
		'Ztautau',
		'TTplusJets', 
		'qcdEMenrichedSum',
		'qcdBCtoESum',
		'gammaPlusJetsSum',
		'Zee',
		'WplusJets'
	)
	ZtTGWZeQ = cms.vstring(
		'Data',
		'Ztautau',
		'TTplusJets', 
		'gammaPlusJetsSum',
		'WplusJets',
		'Zee',
		'qcdEMenrichedSum',
		'qcdBCtoESum'
	)
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electron_afterElectronEta.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electron_afterElectronEta.plots.processes = ZtTZeWGQ

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electron_afterElectronId.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electron_afterElectronId.plots.processes = ZtTZeWGQ

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electronTrkIP_afterElectronConversionVeto.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electronTrkIP_afterElectronConversionVeto.plots.processes = ZtTWGQZe

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electron_afterElectronEcalIso.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electron_afterElectronEcalIso.plots.processes = ZtTWGQZe

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electronTrkIso_afterTauPt.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electronTrkIso_afterTauPt.plots.processes = ZtTZeWGQ

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electronEcalIso_afterElectronTrkIso.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electronEcalIso_afterElectronTrkIso.plots.processes = ZtTZeWGQ

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electronEcalIsoBarrel_afterElectronTrkIso.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electronEcalIsoBarrel_afterElectronTrkIso.plots.processes = ZtTWGQZe

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electronEcalIsoEndcap_afterElectronTrkIso.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electronEcalIsoEndcap_afterElectronTrkIso.plots.processes = ZtTWGQZe

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electron_afterElectronAntiCrack.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electron_afterElectronAntiCrack.plots.processes = ZtTZeWGQ

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electron_afterPrimaryEventVertexPosition.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electron_afterPrimaryEventVertexPosition.plots.processes = ZtTZeWGQ

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tauLeadTrkPt_afterElectronTrkIP.plots.processes = ZtTQGZeW

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tauLeadTrkPt_afterTauAntiOverlapWithElectronsVeto.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tauLeadTrkPt_afterTauAntiOverlapWithElectronsVeto.plots.processes = ZtTGWZeQ

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tauLeadTrkPt_afterTauEta.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tauLeadTrkPt_afterTauEta.plots.processes = ZtTGWZeQ

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tauLeadTrkPt_afterTauLeadTrk.plots.processes = ZtTQGZeW

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tau_afterElectronTrkIP.plots.processes = ZtTQGZeW
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_electron_afterElectronTrkIP.plots.processes = ZtTQGZeW

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tau_afterElectronPt.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tau_afterElectronPt.plots.processes = ZtTZeWGQ

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tau_afterTauAntiOverlapWithElectronsVeto.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tau_afterTauAntiOverlapWithElectronsVeto.plots.processes = ZtTZeWGQ

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tauAntiElectronDiscr_afterTauCharge.plots.processes = TGQZtWZe

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tau_afterTauEta.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_tau_afterTauEta.plots.processes = ZtTZeWGQ

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_vertexChi2Prob_afterPrimaryEventVertex.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_vertexChi2Prob_afterPrimaryEventVertex.plots.processes = ZtTZeWGQ

	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_vertexZ_afterPrimaryEventVertexQuality.yAxis = cms.string('numEntries_log')
	process.plotZtoElecTau.drawJobs.cutFlowControlPlots_vertexZ_afterPrimaryEventVertexQuality.plots.processes = ZtTZeWGQ

