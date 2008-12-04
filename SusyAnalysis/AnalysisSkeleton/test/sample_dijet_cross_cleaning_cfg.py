import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")
process.load("FWCore.MessageLogger.MessageLogger_cfi")







# Include PAT Layer 0 & 1 if not running on pattified data
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V4::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

# CaloTowerConstituentsMap needed for Electron/Photon-Jet cleaning
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)

# Cross-cleaner setup
process.load("SusyAnalysis.PatCrossCleaner.patCrossCleaner_cfi")
# Switch on/off some components
process.patcrosscleaner.doMuonJetCC        = True
process.patcrosscleaner.doElectronJetCC    = True
process.patcrosscleaner.doPhotonJetCC      = False
process.patcrosscleaner.doElectronPhotonCC = False
# Change the jet energy corrections
process.patcrosscleaner.L1JetCorrector      = 'none'
process.patcrosscleaner.L2JetCorrector      = 'L2RelativeJetCorrectorIC5Calo'
process.patcrosscleaner.L3JetCorrector      = 'L3AbsoluteJetCorrectorIC5Calo'
process.patcrosscleaner.L4JetCorrector      = 'none'
process.patcrosscleaner.L5udsJetCorrector   = 'none'
process.patcrosscleaner.L5gluonJetCorrector = 'none'
process.patcrosscleaner.L5cJetCorrector     = 'none'
process.patcrosscleaner.L5bJetCorrector     = 'none'
process.patcrosscleaner.L6JetCorrector      = 'none'
process.patcrosscleaner.L7udsJetCorrector   = 'L7PartonJetCorrectorIC5qJet'
process.patcrosscleaner.L7gluonJetCorrector = 'L7PartonJetCorrectorIC5gJet'
process.patcrosscleaner.L7cJetCorrector     = 'L7PartonJetCorrectorIC5cJet'
process.patcrosscleaner.L7bJetCorrector     = 'L7PartonJetCorrectorIC5bJet'

# Parameters for electron-jet cross-cleaning
process.patcrosscleaner.ElectronJetCrossCleaning.SusyAnalyzerCleaning = True
process.patcrosscleaner.ElectronJetCrossCleaning.deltaR_min        = 0.5
process.patcrosscleaner.ElectronJetCrossCleaning.SharedEtoJetE     = 0.7
process.patcrosscleaner.ElectronJetCrossCleaning.SharedEForNIsoEle = -1.
process.patcrosscleaner.ElectronJetCrossCleaning.IsolationKey  = 'TrackerIso'
process.patcrosscleaner.ElectronJetCrossCleaning.IsoValueCut   = 1.
# Parameters for photon-jet cross-cleaning
process.patcrosscleaner.PhotonJetCrossCleaning.deltaR_min   = 0.5
process.patcrosscleaner.PhotonJetCrossCleaning.IsoValueCut  = 0.3
process.patcrosscleaner.PhotonJetCrossCleaning.IsolationKey = 'CaloIso'
# Parameters for muon-jet cross-cleaning
process.patcrosscleaner.MuonJetCrossCleaning.deltaR_min   = 0.2
process.patcrosscleaner.MuonJetCrossCleaning.caloIso_max  = 10.0
process.patcrosscleaner.MuonJetCrossCleaning.trackIso_max = 10.0




process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('QCD250to500-madgraph.root')
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    

   '/store/relval/CMSSW_2_1_9/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_Tauola_v1/0002/008F0E5C-5C8E-DD11-A113-001617C3B6DC.root'

    )

)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.dijet = cms.EDFilter("SusyDiJetAnalysis",
    genTag = cms.InputTag("genParticles"),
    tauTag = cms.InputTag("selectedLayer1Taus"),
   elecTag = cms.InputTag("patcrosscleaner:ccElectrons"),                  
    photTag = cms.InputTag("selectedLayer1Photons"),
   # jetTag = cms.InputTag("selectedLayer1Jets") ,
    jetTag = cms.InputTag("patcrosscleaner:ccJets"),             
   #muonTag = cms.InputTag("selectedLayer1Muons"),
     muonTag = cms.InputTag("patcrosscleaner:ccMuons"),                  
    #metTag = cms.InputTag("selectedLayer1METs"),                          
   metTag = cms.InputTag("patcrosscleaner:ccMETs"),

 selections = cms.PSet(
        selectionSequence = cms.vstring('Preselection', 
            'FinalJet', 
            'FinalDirectLeptonVeto', 
            'FinalMaxNumJetsSelector', 
            'DPhi', 
            'Alpha', 
            'Hemisphere'),
        selectors = cms.PSet(
            FinalMaxNumJetsSelector = cms.PSet(
                maxEt = cms.double(30.0),
              #   jetTag = cms.InputTag("selectedLayer1Jets"),
                jetTag = cms.InputTag("patcrosscleaner:ccJets"),
                maxNumJets = cms.uint32(100),
                selector = cms.string('MaxNumJetsEventSelector')
            ),
            DPhi = cms.PSet(
                maxDPhi = cms.double(3.15),
               #  jetTag = cms.InputTag("selectedLayer1Jets"),
               jetTag = cms.InputTag("patcrosscleaner:ccJets"),
                selector = cms.string('DPhiEventSelector')
            ),
            Hemisphere = cms.PSet(
                dPhiHemispheresMin = cms.double(0.0),
                dPhiHemispheresMax = cms.double(3.2),
                hemisphereTag = cms.InputTag("selectedLayer2Hemispheres"),
                selector = cms.string('HemisphereSelector')
            ),
            FinalDirectLeptonVeto = cms.PSet(
                #electronTag = cms.InputTag("selectedLayer1Electrons"),
                electronTag = cms.InputTag("patcrosscleaner:ccElectrons"),
                tauTag = cms.InputTag("selectedLayer1Taus"),
                minMuonEt = cms.double(30000.0),
                tauIsolation = cms.double(0.5),
                selector = cms.string('DirectLeptonVetoSelector'),
               muonTag = cms.InputTag("patcrosscleaner:ccMuons"),
               # muonTag = cms.InputTag("selectedLayer1Muons"),
                minTauEt = cms.double(30000.0),
                minElectronEt = cms.double(30000.0),
                muonIsolation = cms.double(0.5),
                electronIsolation = cms.double(0.5)
            ),
            FinalJet = cms.PSet(
                maxEMFraction = cms.vdouble(999.0, 999.0),
                maxEta = cms.vdouble(10.0, 10.0),
                correction = cms.string('HAD'),
                flavour = cms.string('GLU'),
                selector = cms.string('JetEventSelector'),
                # jetTag = cms.InputTag("selectedLayer1Jets"),
                jetTag = cms.InputTag("patcrosscleaner:ccJets"),
                minEt = cms.vdouble(50.0, 50.0)
            ),
            PreJet = cms.PSet(
                maxEMFraction = cms.vdouble(999.0, 999.0),
                maxEta = cms.vdouble(10.0, 10.0),
                correction = cms.string('HAD'),
                 flavour = cms.string('GLU'),
                selector = cms.string('JetEventSelector'),
                # jetTag = cms.InputTag("selectedLayer1Jets"),
                jetTag = cms.InputTag("patcrosscleaner:ccJets"),
                minEt = cms.vdouble(50.0, 50.0)
            ),
            Preselection = cms.PSet(
                components = cms.vstring('PreJet'),
                selector = cms.string('EventSelectorAND')
            ),
            Alpha = cms.PSet(
                minAlpha = cms.double(0.0),
                jetTag = cms.InputTag("patcrosscleaner:ccJets"),
               # jetTag = cms.InputTag("selectedLayer1Jets"),
                selector = cms.string('AlphaSelector')
            )
        )
    ),
                             
                  

    pathNames = cms.vstring('HLT_Jet180','HLT_DiJetAve130','HLT_MET50','HLT_Mu9')          ,    
                   
    eventWeight = cms.double(0.298),
  
    triggerResults = cms.InputTag("TriggerResults","","HLT"),
    plotSelection = cms.vstring('Preselection'),
                             
         
    
  
)

#process.genParticles.abortOnUnknownPDGCode = False

process.selectedLayer2Hemispheres = cms.EDProducer("PATHemisphereProducer",
    patJets = cms.InputTag("patcrosscleaner","ccJets"),
    maxTauEta = cms.double(-1.0),
    maxPhotonEta = cms.double(5.0),
    minMuonEt = cms.double(10.0),
    patMuons = cms.InputTag("patcrosscleaner","ccMuons"),
    maxElectronEta = cms.double(5.0),
    patElectrons = cms.InputTag("patcrosscleaner","ccElectrons"),
    patMets = cms.InputTag("patcrosscleaner","ccMETs"),
    maxMuonEta = cms.double(5.0),
    minTauEt = cms.double(1000000.0),
    minPhotonEt = cms.double(200000.0),
    minElectronEt = cms.double(10.0),
    minJetEt = cms.double(50.0),
    combinationMethod = cms.int32(3),
    maxJetEta = cms.double(5.0),
    seedMethod = cms.int32(3),
    patPhotons = cms.InputTag("selectedLayer1Photons"),
    patTaus = cms.InputTag("selectedLayer1Taus")
)


process.MessageLogger.categories.extend(['SelectorSequence', 'EventSelectorAND', 'EventSelectorOR', 'JetEventSelector', 'BJetEventSelector', 
    'MetEventSelector', 'AlphaSelector', 'DPhiEventSelector', 'SusyDiJetEvent'])
process.MessageLogger.cerr.threshold = 'Verbatim'
process.MessageLogger.cerr.default.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 1000



#process.p = cms.Path(process.patLayer0*process.patLayer1*process.patcrosscleaner*process.selectedLayer2Hemispheres*process.dijet)
#process.p = cms.Path(process.patLayer0*process.patLayer1*process.patcrosscleaner*process.selectedLayer2Hemispheres*process.dijet)
process.p = cms.Path(process.patLayer0*process.patLayer1*process.patcrosscleaner*process.selectedLayer2Hemispheres*process.dijet)


## Necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)
