
# import the definition of the steps and input files:
from  Configuration.PyReleaseValidation.relval_steps import *

# here only define the workflows as a combination of the steps defined above:
workflows = {}

# each workflow defines a name and a list of steps to be done. 
# if no explicit name/label given for the workflow (first arg),
# the name of step1 will be used

workflows[501]=['',['MinBias_TuneZ2_7TeV_pythia6','HARVGEN']]
workflows[502]=['',['QCD_Pt-30_TuneZ2_7TeV_pythia6','HARVGEN']]
workflows[503]=['',['TT_TuneZ2_7TeV_pythia6-evtgen','HARVGEN']]
workflows[504]=['',['DYToLL_M-50_TuneZ2_7TeV_pythia6-tauola','HARVGEN']]
workflows[505]=['',['WToLNu_TuneZ2_7TeV_pythia6-tauola','HARVGEN']]
workflows[506]=['',['MinBias_7TeV_pythia8','HARVGEN']]
workflows[507]=['',['QCD_Pt-30_7TeV_pythia8','HARVGEN']]
workflows[508]=['',['QCD_Pt-30_7TeV_herwig6','HARVGEN']]
workflows[509]=['',['QCD_Pt-30_7TeV_herwigpp','HARVGEN']]
workflows[510]=['',['GluGluTo2Jets_M-100_7TeV_exhume','HARVGEN']]
workflows[511]=['',['DYToMuMu_M-20_7TeV_mcatnlo','HARVGEN']]
workflows[512]=['',['TT_7TeV_mcatnlo','HARVGEN']]
workflows[513]=['',['WminusToENu_7TeV_mcatnlo','HARVGEN']]
workflows[514]=['',['WminusToMuNu_7TeV_mcatnlo','HARVGEN']]
workflows[515]=['',['WplusToENu_7TeV_mcatnlo','HARVGEN']]
workflows[516]=['',['WplusToMuNu_7TeV_mcatnlo','HARVGEN']]
workflows[517]=['',['QCD_Ht-100To250_TuneZ2_7TeV_madgraph-tauola','HARVGEN']]
workflows[518]=['',['QCD_Ht-250To500_TuneZ2_7TeV_madgraph-tauola','HARVGEN']]
workflows[519]=['',['QCD_Ht-500To1000_TuneZ2_7TeV_madgraph-tauola','HARVGEN']]
workflows[520]=['',['TTJets_TuneZ2_7TeV_madgraph-tauola','HARVGEN']]
workflows[521]=['',['WJetsLNu_TuneZ2_7TeV_madgraph-tauola','HARVGEN']]
workflows[522]=['',['ZJetsLNu_TuneZ2_7TeV_madgraph-tauola','HARVGEN']]
workflows[523]=['',['QCD2Jets_Pt-40To120_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[524]=['',['QCD3Jets_Pt-40To120_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[525]=['',['QCD4Jets_Pt-40To120_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[526]=['',['QCD5Jets_Pt-40To120_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[527]=['',['TT0Jets_Et-40_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[528]=['',['TT1Jets_Et-40_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[529]=['',['TT2Jets_Et-40_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[530]=['',['TT3Jets_Et-40_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[531]=['',['W0Jets_Pt-0To100_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[532]=['',['W1Jets_Pt-0To100_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[533]=['',['W2Jets_Pt-0To100_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[534]=['',['W3Jets_Pt-0To100_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[535]=['',['Z0Jets_Pt-0To100_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[536]=['',['Z1Jets_Pt-0To100_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[537]=['',['Z2Jets_Pt-0To100_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[538]=['',['Z3Jets-Pt_0To100_TuneZ2_7TeV_alpgen_tauola','HARVGEN']]
workflows[539]=['',['ZJetsLNu_Tune4C_7TeV_madgraph-pythia8','HARVGEN']]

