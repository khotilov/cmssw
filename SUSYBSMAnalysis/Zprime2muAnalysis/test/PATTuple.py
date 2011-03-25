#!/usr/bin/env python

# Parameters you usually want to change are here at the top of the
# file.

# Change the following to False to run on data.
input_is_MC = True

# Change the following parameters appropriately for the samples over
# which you wish to run, most especially the globalTag (since the list
# of files will be configured by CRAB itself so it doesn't matter what
# gets put here for it).
globalTag = 'START311_V2::All'
files = ['file:/uscms/home/tucker/nobackup/Spring11/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/AODSIM/PU_S1_START311_V1G1-v1/0002/181B6610-494E-E011-A311-001A4BE1C5D4.root']
maxEvents = 100

########################################################################

# Load and configure the process appropriately.

import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTuple_cfg import process
process.GlobalTag.globaltag = globalTag
process.source.fileNames = files
process.maxEvents.input = maxEvents

# Do the right thing depending on whether MC truth is supposed to be
# available.
if not input_is_MC:
    from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import removeMCUse
    removeMCUse(process)

# Our path. The name 'p' is expected by the configuration of the
# OutputModule in PATTuple_cfg.
process.p = cms.Path(process.patDefaultSequence)

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import AODOnly
AODOnly(process)

########################################################################

# Changing things/uncommenting below here can be done by slightly more
# advanced users :-) You can use a bunch of core tools of PAT and some
# Exotica tools to tailor your PAT configuration; for a few, examples
# are commented out in the next several lines.

#from PhysicsTools.PatAlgos.tools.coreTools import restrictInputToAOD, removeAllPATObjectsBut, removeSpecificPATObjects
#restrictInputToAOD(process)
#removeAllPATObjectsBut(process, ['Muons', 'Electrons'])
#removeSpecificPATObjects(process, ['Taus'])

# Some MC samples have the HLT process name different from "HLT".
from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import switchHLTProcessName
switchHLTProcessName(process, "REDIGI311X")

# You can change the default muon cut from "isGlobalMuon ||
# isTrackerMuon" to your taste. Some examples (only uncomment one!):
#process.selectedPatMuons.cut = 'isTrackerMuon && pt > 1 && p > 2.5 && innerTrack.hitPattern.numberOfValidTrackerHits > 12 && innerTrack.normalizedChi2 < 5 && abs(dB) < 0.5 && abs(dZ) < 5 && muonID("TMLastStationAngTight")'
#process.selectedPatMuons.cut = 'muonID("GlobalMuonPromptTight") && muonID("TMOneStationLoose") && (globalTrack.hitPattern.numberOfValidMuonCSCHits + globalTrack.hitPattern.numberOfValidMuonDTHits) >= 1 && innerTrack.hitPattern.trackerLayersWithMeasurement >= 6'

# You may only want to keep electrons that pass the HEEP cuts, but
# does this screw up the PAT jet cleaning? Needs to be studied.
#process.selectedPatElectrons.cut = 'userInt("HEEPId") == 0'

# If you're studying e.g. e/mu dileptons, you may want to change the
# requirement from one/two muons to two leptons (= electrons + muons).
#process.countPatMuons.minNumber = 0
#process.countPatLeptons.minNumber = 2
