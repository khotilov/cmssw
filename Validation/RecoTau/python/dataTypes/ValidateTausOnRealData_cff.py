import FWCore.ParameterSet.Config as cms
from Validation.RecoTau.RecoTauValidation_cfi import *
import copy

from RecoJets.Configuration.RecoPFJets_cff import *
import PhysicsTools.PatAlgos.tools.helpers as helpers

kinematicSelectedTauValDenominatorRealData = cms.EDFilter( ##FIXME: this should be a filter
   "TauValPFJetSelector", #"GenJetSelector"
   src = cms.InputTag("ak5PFJets"),
   cut = kinematicSelectedTauValDenominatorCut,#cms.string('pt > 5. && abs(eta) < 2.5'), #Defined: Validation.RecoTau.RecoTauValidation_cfi 
   filter = cms.bool(False)
)

procAttributes = dir(proc) #Takes a snapshot of what there in the process
helpers.cloneProcessingSnippet( proc, proc.TauValNumeratorAndDenominator, 'RealData') #clones the sequence inside the process with RealData postfix
helpers.cloneProcessingSnippet( proc, proc.TauEfficiencies, 'RealData') #clones the sequence inside the process with RealData postfix
helpers.massSearchReplaceAnyInputTag(proc.TauValNumeratorAndDenominatorRealData, 'kinematicSelectedTauValDenominator', 'kinematicSelectedTauValDenominatorRealData') #sets the correct input tag

#adds to TauValNumeratorAndDenominator modules in the sequence RealData to the extention name
zttLabeler = lambda module : SetValidationExtention(module, 'RealData')
zttModifier = ApplyFunctionToSequence(zttLabeler)
proc.TauValNumeratorAndDenominatorRealData.visit(zttModifier)

#-----------------------------------------Sets binning
binning = cms.PSet(
    pt = cms.PSet( nbins = cms.int32(8), min = cms.double(0.), max = cms.double(300.) ), #hinfo(75, 0., 150.)
    eta = cms.PSet( nbins = cms.int32(4), min = cms.double(-3.), max = cms.double(3.) ), #hinfo(60, -3.0, 3.0);
    phi = cms.PSet( nbins = cms.int32(4), min = cms.double(-180.), max = cms.double(180.) ), #hinfo(36, -180., 180.);
    pileup = cms.PSet( nbins = cms.int32(16), min = cms.double(0.), max = cms.double(80.) ),#hinfo(25, 0., 25.0);
    )
zttModifier = ApplyFunctionToSequence(lambda m: setBinning(m,binning))
proc.TauValNumeratorAndDenominatorRealData.visit(zttModifier)
#-----------------------------------------

#Sets the correct naming to efficiency histograms
proc.efficienciesRealData.plots = Utils.SetPlotSequence(proc.TauValNumeratorAndDenominatorRealData)

#checks what's new in the process (the cloned sequences and modules in them)
newProcAttributes = filter( lambda x: (x not in procAttributes) and (x.find('RealData') != -1), dir(proc) )

#spawns a local variable with the same name as the proc attribute, needed for future process.load
for newAttr in newProcAttributes:
    locals()[newAttr] = getattr(proc,newAttr)


produceDenominatorRealData = cms.Sequence(
      kinematicSelectedTauValDenominatorRealData
      )

produceDenominator = produceDenominatorRealData

runTauValidationBatchMode = cms.Sequence(
      produceDenominator
      +TauValNumeratorAndDenominatorRealData
      )

runTauValidation = cms.Sequence(
      runTauValidationBatchMode
      +TauEfficienciesRealData
      )
