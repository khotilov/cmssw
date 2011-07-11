
# import the definition of the steps and input files:
from  Configuration.PyReleaseValidation.relval_steps import *

# here only define the workflows as a combination of the steps defined above:
workflows = {}

# each workflow defines a name and a list of steps to be done. 
# if no explicit name/label given for the workflow (first arg),
# the name of step1 will be used

## production tests
workflows[1] = ['', ['ProdMinBias','DIGIPROD1','RECOPROD1']]
workflows[2] = ['', ['ProdTTbar','DIGIPROD1','RECOPROD1']]
workflows[3] = ['', ['ProdQCD_Pt_3000_3500','DIGIPROD1','RECOPROD1']]

### data ###
workflows[4.5]  = ['', ['RunCosmicsA','RECOCOSD','ALCACOSD']]
#workflows[4.45] = ['', ['RunCosmicsA','RECOD']]
workflows[4.6]  = ['', ['MinimumBias2010A','RECOSKIM','HARVESTD']]
workflows[4.7]  = ['', ['MinimumBias2010B','RECOSKIMALCA']]
workflows[4.8]  = ['', ['WZMuSkim2010A','RECOSKIM']]
workflows[4.9]  = ['', ['WZEGSkim2010A','RECOSKIM']]
workflows[4.10] = ['', ['WZMuSkim2010B','RECOSKIM']]
workflows[4.11] = ['', ['WZEGSkim2010B','RECOSKIM']]

workflows[4.12] = ['', ['RunMinBias2010B','RECOD']]
workflows[4.13] = ['', ['RunMu2010B','RECOD']]
workflows[4.14] = ['', ['RunElectron2010B','RECOD']]
workflows[4.15] = ['', ['RunPhoton2010B','RECOD']]
workflows[4.16] = ['', ['RunJet2010B','RECOD']]


workflows[4.17] = ['', ['RunMinBias2011A','RECOD','HARVESTD','SKIMD']]
workflows[4.18] = ['', ['RunMu2011A','RECOD']]
workflows[4.19] = ['', ['RunElectron2011A','RECOD']]
workflows[4.20] = ['', ['RunPhoton2011A','RECOD']]
workflows[4.21] = ['', ['RunJet2011A','RECOD']]

workflows[4.22] = ['', ['RunCosmics2011A','RECOCOSD','ALCACOSD','SKIMCOSD']]

#workflows[4.23] = ['',[]
#workflows[4.24] = ['',['WMuSkim2011A','RECOSKIM']]
#workflows[4.25] = ['',['WElSkim2011A','RECOSKIM']]
#workflows[4.26] = ['',['ZMuSkim2011A','RECOSKIM']]
#workflows[4.27] = ['',['ZElSkim2011A','RECOSKIM']]
#workflows[4.28] = ['',['HighMet2011A','RECOSKIM']]

workflows[4.51] = ['',['RunHI2010','RECOHID']]

### fastsim ###
#workflows[5.1] = ['TTbar', ['TTbarFS1']]
workflows[6.3] = ['TTbar', ['TTbarFS']]
workflows[5.2] = ['SingleMuPt10', ['SingleMuPt10FS']]
workflows[5.3] = ['SingleMuPt100', ['SingleMuPt100FS']]
#workflows[6.1] = ['ZEE', ['ZEEFS1']]
workflows[6.2] = ['ZEE', ['ZEEFS']]
workflows[39]  = ['QCDFlatPt153000', ['QCDFlatPt153000FS']]
workflows[6.4] = ['H130GGgluonfusion', ['H130GGgluonfusionFS']]

### standard set ###
#workflows[10] = ['', ['MinBias','DIGI','RECO']]
#workflows[12] = ['', ['QCD_Pt_3000_3500','DIGI','RECO']]
#workflows[14] = ['', ['QCD_Pt_80_120','DIGI','RECO']]
workflows[16] = ['', ['SingleElectronPt10','DIGI','RECO']]
workflows[17] = ['', ['SingleElectronPt35','DIGI','RECO']]
workflows[18] = ['', ['SingleGammaPt10','DIGI','RECO']]
workflows[19] = ['', ['SingleGammaPt35','DIGI','RECO']]
workflows[20] = ['', ['SingleMuPt10','DIGI','RECO']]
workflows[21] = ['', ['SingleMuPt100','DIGI','RECO']]
workflows[22] = ['', ['SingleMuPt1000','DIGI','RECO']]
#workflows[24] = ['', ['TTbar','DIGI','RECO','HARVEST']]
#workflows[28] = ['', ['ZEE','DIGI','RECO']]
workflows[35] = ['', ['Wjet_Pt_80_120','DIGI','RECO']]
workflows[36] = ['', ['Wjet_Pt_3000_3500','DIGI','RECO']]
workflows[37] = ['', ['LM1_sfts','DIGI','RECO']]
workflows[38] = ['', ['QCD_FlatPt_15_3000','DIGI','RECO']]

workflows[9]  = ['', ['Higgs200ChargedTaus','DIGI','RECO']]
workflows[13] = ['', ['QCD_Pt_3000_3500','DIGI','RECO']]
workflows[23] = ['', ['JpsiMM','DIGI','RECO']]
workflows[25] = ['', ['TTbar','DIGI','RECO','ALCATT2']]
workflows[26] = ['', ['WE','DIGI','RECO','HARVEST2']]
workflows[29] = ['', ['ZEE','DIGI','RECO','ALCAELE']]
workflows[31] = ['', ['ZTT','DIGI','RECO']]
workflows[32] = ['', ['H130GGgluonfusion','DIGI','RECO']]
workflows[33] = ['', ['PhotonJets_Pt_10','DIGI','RECO']]
workflows[34] = ['', ['QQH1352T_Tauola','DIGI','RECO']]

workflows[7]  = ['', ['Cosmics','DIGICOS','RECOCOS','ALCACOS']]
workflows[8]  = ['', ['BeamHalo','DIGICOS','RECOCOS','ALCABH']]
workflows[11] = ['', ['MinBias','DIGI','RECOMIN','ALCAMIN']]
workflows[15] = ['', ['QCD_Pt_80_120','DIGI','RECOQCD']]
workflows[27] = ['', ['WM','DIGI','RECOMU']]
workflows[30] = ['', ['ZMM','DIGI','RECOMU']]


### HI test ###
workflows[40] = ['',['HydjetQ_MinBias_2760GeV','DIGIHI','RECOHI']]
workflows[41] = ['',['HydjetQ_B0_2760GeV','DIGIHI','RECOHI']]


