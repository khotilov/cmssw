import copy
from ElectroWeakAnalysis.MultiBosons.Skimming.options import options as defaultOptions

def applyJobOptions(options):
  """
  Set multiple options defined by the jobType argument.
  To be called after options.parseArguments()
  """
  options.parseArguments()
  jobOptions = copy.deepcopy(defaultOptions)

  ## Set the default trigger match paths
  jobOptions.electronTriggerMatchPaths = """
    HLT_Ele10_LW_L1R
    HLT_Ele12_SW_EleIdIsol_L1R
    HLT_Ele15_LW_L1R
    HLT_Ele15_SW_EleId_L1R
    HLT_Ele15_SW_L1R
    HLT_Ele15_SW_LooseTrackIso_L1R
    HLT_Ele17_SW_CaloEleId_L1R
    HLT_Ele17_SW_EleIdIsol_L1R
    HLT_Ele17_SW_LEleId_L1R
    HLT_Ele20_SW_L1R
    HLT_Photon10_L1R
    HLT_Photon15_L1R
    """.split()
  jobOptions.muonTriggerMatchPaths = """
    HLT_Mu9
    HLT_Mu11
    """.split()
  jobOptions.electronTriggerMatchPaths = """
    HLT_Ele10_LW_L1R
    HLT_Ele12_SW_EleIdIsol_L1R
    HLT_Ele15_LW_L1R
    HLT_Ele15_SW_EleId_L1R
    HLT_Ele15_SW_L1R
    HLT_Ele15_SW_LooseTrackIso_L1R
    HLT_Ele17_SW_CaloEleId_L1R
    HLT_Ele17_SW_EleIdIsol_L1R
    HLT_Ele17_SW_LEleId_L1R
    HLT_Ele20_SW_L1R
    HLT_Photon10_L1R
    HLT_Photon15_L1R
    """.split()


  if options.jobType == "testMC":
    jobOptions.maxEvents = 100
    jobOptions.inputFiles = [
      "/store/relval/CMSSW_3_5_7/" +
      "RelValZMM/GEN-SIM-RECO/START3X_V26-v1/0012/" +
      file for file in """
        10B71379-4549-DF11-9D80-003048D15D22.root
        34FD3B1D-6949-DF11-9529-0018F3D09612.root
        4C8D7358-4449-DF11-86BB-003048678A6C.root
      """.split()
    ]
    jobOptions.globalTag = "START38_V10::All"
    jobOptions.reportEvery = 1
    jobOptions.isRealData = False
    jobOptions.use35XInput = True
    jobOptions.isMaxEventsOutput = True
    jobOptions.wantSummary = False
    jobOptions.hltPaths = ["HLT_Mu9"]
    jobOptions.electronTriggerMatchPaths = """
      HLT_Ele10_LW_L1R
      HLT_Ele12_SW_EleIdIsol_L1R
      HLT_Ele15_LW_L1R
      HLT_Ele15_SW_EleId_L1R
      HLT_Ele15_SW_L1R
      HLT_Ele15_SW_LooseTrackIso_L1R
      HLT_Ele17_SW_CaloEleId_L1R
      HLT_Ele17_SW_EleIdIsol_L1R
      HLT_Ele17_SW_LEleId_L1R
      HLT_Ele20_SW_L1R
      HLT_Photon10_L1R
      HLT_Photon15_L1R
      """.split()
    jobOptions.muonTriggerMatchPaths = """
      HLT_Mu9
      HLT_Mu11
      """.split()

  # end of testMC options <-------------------------------------

  elif options.jobType == "testRealData":
    jobOptions.maxEvents = 1
    jobOptions.inputFiles = [
      "/store/data/Run2010A/Mu/RECO/v4/000/" +
      file for file in """
        142/933/C043478A-79A7-DF11-BB4C-001D09F2426D.root
        142/933/C277D905-95A7-DF11-A70D-003048F024F6.root
        142/933/CC18C413-9CA7-DF11-91B8-003048F118AA.root
        142/933/D433F278-9FA7-DF11-8403-003048F11114.root
        142/933/DE260EAC-A3A7-DF11-B4A8-0030487CD7B4.root
        143/657/A8367E04-73AE-DF11-9989-001D09F297EF.root
        143/657/B81CBEDC-C3AE-DF11-8909-003048F118D2.root
        143/657/BE540563-9AAE-DF11-BF6D-003048D3750A.root
        143/657/BEF1F8EB-DFAE-DF11-AF55-0030487A1FEC.root
        143/657/C0905609-92AE-DF11-809B-0030487CD77E.root
        143/657/D2E2D0AB-8EAE-DF11-A325-0030487D0D3A.root
        143/657/D689519B-D2AE-DF11-8A48-003048F024C2.root
        143/657/D8EB995B-98AE-DF11-B1CB-0030487A3232.root
        143/657/DC26FBFB-A2AE-DF11-8459-0030487C2B86.root
        143/657/DE4FFA93-95AE-DF11-AFAA-001D09F24259.root
        143/657/E07818F5-A9AE-DF11-B252-0030487CD6F2.root
      """.split()
    ]
    jobOptions.globalTag = "GR10_P_V9::All"
    jobOptions.reportEvery = 1
    jobOptions.isRealData = True
    jobOptions.use35XInput = False
    jobOptions.isMaxEventsOutput = True
    jobOptions.wantSummary = False
#     jobOptions.hltPaths = ["HLT_Mu%d" % i for i in [0, 3, 5, 7, 9, 11]]
    jobOptions.hltPaths = ["HLT_Mu9"]
    jobOptions.muonTriggerMatchPaths = """
      HLT_L1Mu14_L1ETM30
      HLT_L1Mu14_L1SingleJet6U
      HLT_L1Mu14_L1SingleEG10
      HLT_L1Mu20
      HLT_DoubleMu3
      HLT_Mu3
      HLT_Mu5
      HLT_Mu9
      HLT_L2Mu9
      HLT_L2Mu11
      HLT_L1Mu30
      HLT_Mu7
      HLT_L2Mu15
      """.split()  ## Definition of MU PD for run 142933
  # end of testRealData options <-----------------------------------

  elif options.jobType == "testSummer10":
    jobOptions.inputFiles = [
        "rfio:/castor/cern.ch/user/z/zkliu/TestRECO/36x_ZJet_Madgraph_tauola_TestRECO.root"
        ]
    jobOptions.skimType = "ElectronPhoton"
    jobOptions.globalTag = "START38_V10::All"
    jobOptions.reportEvery = 1
    jobOptions.isRealData = False
    jobOptions.use35XInput = False
    jobOptions.isMaxEventsOutput = True
    jobOptions.wantSummary = True
    jobOptions.hltProcessName = "REDIGI36X"
    jobOptions.hltPaths = ["HLT_Ele15_LW_L1R", "HLT_Ele15_SW_L1R"]

  elif options.jobType == "PromptReco36X":
    jobOptions.globalTag = "GR10_P_V7::All"
    jobOptions.isRealData = True
    jobOptions.use35XInput = False
    jobOptions.maxEvents = -1

  elif options.jobType == "PromptReco38X":
    jobOptions.globalTag = "GR10_P_V9::All"
    jobOptions.isRealData = True
    jobOptions.use35XInput = False
    jobOptions.maxEvents = -1

  elif options.jobType == "ReReco36X":
    jobOptions.globalTag = "GR_R_36X_V12B::All"
    jobOptions.isRealData = True
    jobOptions.use35XInput = False
    jobOptions.maxEvents = -1

  elif options.jobType == "ReReco38X":
    jobOptions.globalTag = "GR_R_38X_V9A::All"
    jobOptions.isRealData = True
    jobOptions.use35XInput = False
    jobOptions.maxEvents = -1

  elif options.jobType == "MC36X":
    jobOptions.globalTag = "START3X_V26::All"
    jobOptions.isRealData = False
    jobOptions.use35XInput = True
    jobOptions.maxEvents = -1

  elif options.jobType != "":
    raise RuntimeError, "Unknown jobType option `%s'" % options.jobType

  ## Set all the non-default option values equal to jobOptions to preserve options set
  ##+ in the cfg file
  for name, value in options._singletons.items() + options._lists.items():
    if value != getattr(defaultOptions, name):
      setattr(jobOptions, name, value)

  jobOptions.parseArguments()

  ## Set all the options equal to jobOptions
  for name, value in jobOptions._singletons.items() + jobOptions._lists.items():
    setattr(options, name, value)

#   return jobOptions

# applyMultiOptionTag(options) <----------------------------------------------
