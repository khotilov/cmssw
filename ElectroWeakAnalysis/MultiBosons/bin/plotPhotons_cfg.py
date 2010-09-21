import FWCore.ParameterSet.Config as cms

process = cms.Process("FWLitePlots")


process.inputs = cms.PSet (
    fileNames = cms.vstring(
        # Your data goes here:
#         '/Volumes/MyBook/Data/TTbar/shyft_35x_v3/ljmet_1.root'
#         'file:zmm_pat.root'
#         'file:/tmp/veverka/tmp.YPPUW29322/CMSSW_3_6_2/src/PhysicsTools/SelectorUtils/bin/wm_pat.root'
  "file:VGammaSkim_numEvent10.root"
        )
)

process.outputs = cms.PSet (
    outputName = cms.string('photonPlots.root')
)

# process.load("ElectroWeakAnalysis.MultiBosons.Histogramming.leafKineHistos_cff")

## Convert cms.EDAnalyzer to cms.PSet
from ElectroWeakAnalysis.MultiBosons.Histogramming.photonHistos_cfi import photonHistos
process.photonHistos = cms.PSet(
  src = cms.InputTag(photonHistos.src.value() ),
  histograms = cms.VPSet() + #photonHistos.histograms.copy() +
    [
      cms.PSet(
        itemsToPlot = cms.untracked.int32(-1),
        min = cms.untracked.double(-5),
        max = cms.untracked.double(5),
        nbins = cms.untracked.int32(100),
        name = cms.untracked.string("seedTime"),
        description = cms.untracked.string("Photons:seed time (ns):events / bin"),
        plotquantity = cms.untracked.string('userFloat("photonUserData:seedTime")')
      ),
      cms.PSet(
        itemsToPlot = cms.untracked.int32(-1),
        min = cms.untracked.double(0),
        max = cms.untracked.double(1),
        nbins = cms.untracked.int32(100),
        name = cms.untracked.string("seedOutOfTimeChi2"),
        description = cms.untracked.string("Photons:seed out of time #chi^{2}:events / bin"),
        plotquantity = cms.untracked.string('userFloat("photonUserData:seedOutOfTimeChi2")')
      ),
      cms.PSet(
        itemsToPlot = cms.untracked.int32(-1),
        min = cms.untracked.double(0),
        max = cms.untracked.double(1),
        nbins = cms.untracked.int32(100),
        name = cms.untracked.string("seedChi2"),
        description = cms.untracked.string("Photons:seed #chi^{2}:events / bin"),
        plotquantity = cms.untracked.string('userFloat("photonUserData:seedChi2")')
      ),
      cms.PSet(
        itemsToPlot = cms.untracked.int32(-1),
        min = cms.untracked.double(0),
        max = cms.untracked.double(1),
        nbins = cms.untracked.int32(100),
        name = cms.untracked.string("seedSwissCross"),
        description = cms.untracked.string("Photons:seed swiss cross 1-S4/S1, events / bin"),
        plotquantity = cms.untracked.string('userFloat("photonUserData:seedSwissCross")')
      ),
      cms.PSet(
        itemsToPlot = cms.untracked.int32(-1),
        min = cms.untracked.double(-0.5),
        max = cms.untracked.double(20.5),
        nbins = cms.untracked.int32(21),
        name = cms.untracked.string("seedRecoFlag"),
        description = cms.untracked.string("Photons:seed reco flag:events per flag"),
        plotquantity = cms.untracked.string('userInt("photonUserData:seedRecoFlag")')
      ),
      cms.PSet(
        itemsToPlot = cms.untracked.int32(-1),
        min = cms.untracked.double(-0.5),
        max = cms.untracked.double(20.5),
        nbins = cms.untracked.int32(21),
        name = cms.untracked.string("seedSeverityLevel"),
        description = cms.untracked.string("Photons:Seed severity level:events"),
        plotquantity = cms.untracked.string('userInt("photonUserData:seedSeverityLevel")')
      ),
      cms.PSet(
        itemsToPlot = cms.untracked.int32(-1),
        min = cms.untracked.double(0),
        max = cms.untracked.double(1),
        nbins = cms.untracked.int32(100),
        name = cms.untracked.string("seedE1OverE9"),
        description = cms.untracked.string("Phtons:Seed E1/E9 from severity level:events"),
        plotquantity = cms.untracked.string('userFloat("photonUserData:seedE1OverE9")')
      ),
      cms.PSet(
        itemsToPlot = cms.untracked.int32(-1),
        min = cms.untracked.double(0),
        max = cms.untracked.double(1),
        nbins = cms.untracked.int32(100),
        name = cms.untracked.string("r19"),
        description = cms.untracked.string("eMax / e3x3"),
        plotquantity = cms.untracked.string('userFloat("photonUserData:eMax")/' +
          'userFloat("photonUserData:e3x3")'
          )
      ),
      cms.PSet(
        itemsToPlot = cms.untracked.int32(-1),
        min = cms.untracked.double(0),
        max = cms.untracked.double(100),
        nbins = cms.untracked.int32(100),
        name = cms.untracked.string("e3x3User"),
        description = cms.untracked.string("e3x3 user/cluster shape tools"),
        plotquantity = cms.untracked.string('userFloat("photonUserData:e3x3")')
      ),
      cms.PSet(
        itemsToPlot = cms.untracked.int32(-1),
        min = cms.untracked.double(0),
        max = cms.untracked.double(100),
        nbins = cms.untracked.int32(100),
        name = cms.untracked.string("e3x3Photon"),
        description = cms.untracked.string("e3x3 photon"),
        plotquantity = cms.untracked.string('e3x3')
      ),
      cms.PSet(
        itemsToPlot = cms.untracked.int32(-1),
        min = cms.untracked.double(0),
        max = cms.untracked.double(1),
        nbins = cms.untracked.int32(100),
        name = cms.untracked.string("covIEtaIEta"),
        description = cms.untracked.string("cov(i#eta,i#eta)"),
        plotquantity = cms.untracked.string('userFloat("photonUserData:covIEtaIEta")')
      ),
      cms.PSet(
        itemsToPlot = cms.untracked.int32(-1),
        min = cms.untracked.double(0),
        max = cms.untracked.double(2),
        nbins = cms.untracked.int32(100),
        name = cms.untracked.string("pi0Discriminator"),
        description = cms.untracked.string("pi_{0} discriminator"),
        plotquantity = cms.untracked.string(
          'userFloat("pi0Discriminator")'
          )
      ),
    ],
  eventWeight = cms.double(1.0),
  outputDirectory = cms.string('test')
)


# process.wplusjetsAnalysis.muonSrc = "cleanPatMuons"
# process.wplusjetsAnalysis.electronSrc = "cleanPatElectrons"
# process.wplusjetsAnalysis.jetSrc = "cleanPatJets"
#
# process.wplusjetsAnalysis.minJets = 2

