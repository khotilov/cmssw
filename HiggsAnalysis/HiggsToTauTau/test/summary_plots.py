from optparse import OptionParser, OptionGroup

## set up the option parser
parser = OptionParser(usage="usage: %prog [options]",
                      description="Script to combine 2011+2012 and high and low pt event categories. This script requires that the script run_macros.py has been executed beforehand and produced resultrs w/o issues.")
## direct options
parser.add_option("-a", "--analysis", dest="analysis", default="sm", type="choice", help="Type of analysis (sm or mssm). Lower case is required. [Default: sm]", choices=["sm", "mssm"])
parser.add_option("-t", "--type", dest="type", default="rescaled", type="string", help="Type of plots, unscaled or rescaled. [Default: \"rescaled\"]")
## check number of arguments; in case print usage
(options, args) = parser.parse_args()


## a set of pre-defined lists
channels   = [
    "emu",
    "eleTau",
    "muTau",
    "mumu",
    ]

categories_sm = [
    "0jet_low",
    "0jet_high",
    "0jet",
    "boost_low",
    "boost_high",
   #"boost",  ## these are not combined at the moment due to different binning in low and high pt
    "vbf",
    ]

categories_mssm = [
    #"0jet_low",
    #"0jet_high",
    #"0jet",
    #"boost_low",
    #"boost_high",
    #"boost",
    "nobtag",
    "btag",
    ] 

extra = {
    "emu"   : "#tau_{e}#tau_{#mu}",
    "eleTau"  : "#tau_{e}#tau_{h}",
    "muTau" : "#tau_{#mu}#tau_{h}",
    "mumu"  : "#tau_{#mu}#tau_{#mu}",
    }

log = {
    ("emu"    , "0jet_low"  ) : ["false",],
    ("emu"    , "0jet_high" ) : ["false",],
    ("emu"    , "0jet"      ) : ["false",],
    ("emu"    , "boost_low" ) : ["false",],
    ("emu"    , "boost_high") : ["false",],
    ("emu"    , "boost"     ) : ["false",],
    ("emu"    , "btag_low"  ) : ["false",],
    ("emu"    , "btag_high" ) : ["false",],
    ("emu"    , "nobtag"    ) : ["true", "false"],
    ("emu"    , "btag"      ) : ["true", "false"],
    ("emu"    , "vbf"       ) : ["false",],
    ("muTau"  , "0jet_low"  ) : ["false",],
    ("muTau"  , "0jet_high" ) : ["false",],
    ("muTau"  , "0jet"      ) : ["false",],
    ("muTau"  , "boost_low" ) : ["false",],
    ("muTau"  , "boost_high") : ["false",],
    ("muTau"  , "boost"     ) : ["false",],
    ("muTau"  , "btag_low"  ) : ["false",],
    ("muTau"  , "btag_high" ) : ["false",],
    ("muTau"  , "nobtag"    ) : ["true", "false"],
    ("muTau"  , "btag"      ) : ["true", "false"],
    ("muTau"  , "vbf"       ) : ["false",],
    ("eleTau" , "0jet_low"  ) : ["false",],
    ("eleTau" , "0jet_high" ) : ["false",],
    ("eleTau" , "0jet"      ) : ["false",],
    ("eleTau" , "boost_low" ) : ["false",],
    ("eleTau" , "boost_high") : ["false",],
    ("eleTau" , "boost"     ) : ["false",],
    ("eleTau" , "btag_low"  ) : ["false",],
    ("eleTau" , "btag_high" ) : ["false",],
    ("eleTau" , "nobtag"    ) : ["true", "false"],
    ("eleTau" , "btag"      ) : ["true", "false"],
    ("eleTau" , "vbf"       ) : ["false",],
    ("mumu"   , "0jet_low"  ) : ["false", ],
    ("mumu"   , "0jet_high" ) : ["false", ],
    ("mumu"   , "0jet"      ) : ["false", ],
    ("mumu"   , "boost_low" ) : ["false", ],
    ("mumu"   , "boost_high") : ["false", ],
    ("mumu"   , "boost"     ) : ["false", ],
    ("mumu"   , "btag_low"  ) : ["false", ],
    ("mumu"   , "btag_high" ) : ["false", ],
    ("mumu"   , "nobtag"    ) : ["true", "false"],
    ("mumu"   , "btag"      ) : ["true", "false"],
    ("mumu"   , "vbf"       ) : ["false",],
    }

max = {
    ("emu"    , "0jet_low"  ) : ["-1",],
    ("emu"    , "0jet_high" ) : ["-1",],
    ("emu"    , "0jet"      ) : ["-1",],
    ("emu"    , "boost_low" ) : ["-1",],
    ("emu"    , "boost_high") : ["-1",],
    ("emu"    , "boost"     ) : ["-1",],
    ("emu"    , "btag_low"  ) : ["-1",],
    ("emu"    , "btag_high" ) : ["-1",],
    ("emu"    , "nobtag"    ) : ["3000","2500"],
    ("emu"    , "btag"      ) : ["40","45"],
    ("emu"    , "vbf"       ) : ["-1",],
    ("muTau"  , "0jet_low"  ) : ["-1",],
    ("muTau"  , "0jet_high" ) : ["-1",],
    ("muTau"  , "0jet"      ) : ["-1",],
    ("muTau"  , "boost_low" ) : ["-1",],
    ("muTau"  , "boost_high") : ["-1",],
    ("muTau"  , "boost"     ) : ["-1",],
    ("muTau"  , "btag_low"  ) : ["-1",],
    ("muTau"  , "btag_high" ) : ["-1",],
    ("muTau"  , "nobtag"    ) : ["8000","5000"],
    ("muTau"  , "btag"      ) : ["100","50"],
    ("muTau"  , "vbf"       ) : ["-1",],
    ("eleTau" , "0jet_low"  ) : ["-1",],
    ("eleTau" , "0jet_high" ) : ["-1",],
    ("eleTau" , "0jet"      ) : ["-1",],
    ("eleTau" , "boost_low" ) : ["-1",],
    ("eleTau" , "boost_high") : ["-1",],
    ("eleTau" , "boost"     ) : ["-1",],
    ("eleTau" , "btag_low"  ) : ["-1",],
    ("eleTau" , "btag_high" ) : ["-1",],
    ("eleTau" , "nobtag"    ) : ["2000","1500"],
    ("eleTau" , "btag"      ) : ["20","20"],
    ("eleTau" , "vbf"       ) : ["-1",],
    ("mumu"   , "0jet_low"  ) : ["-1",],
    ("mumu"   , "0jet_high" ) : ["-1",],
    ("mumu"   , "0jet"      ) : ["-1",],
    ("mumu"   , "boost_low" ) : ["-1",],
    ("mumu"   , "boost_high") : ["-1",],
    ("mumu"   , "boost"     ) : ["-1",],
    ("mumu"   , "btag_low"  ) : ["-1",],
    ("mumu"   , "btag_high" ) : ["-1",],
    ("mumu"   , "nobtag"    ) : ["200000","120000"],
    ("mumu"   , "btag"      ) : ["500","300"],
    ("mumu"   , "vbf"       ) : ["-1",],
    }

min = {
    ("emu"    , "0jet_low"  ) : ["0",  ],
    ("emu"    , "0jet_high" ) : ["0",  ],
    ("emu"    , "0jet"      ) : ["0",  ],
    ("emu"    , "boost_low" ) : ["0",  ],
    ("emu"    , "boost_high") : ["0",  ],
    ("emu"    , "boost"     ) : ["0",  ],
    ("emu"    , "btag_low"  ) : ["0",  ],
    ("emu"    , "btag_high" ) : ["0",  ],
    ("emu"    , "nobtag"    ) : ["1e-1", "0"],
    ("emu"    , "btag"      ) : ["1e-1", "0"],
    ("emu"    , "vbf"       ) : ["0",  ],
    ("muTau"  , "0jet_low"  ) : ["0",  ],
    ("muTau"  , "0jet_high" ) : ["0",  ],
    ("muTau"  , "0jet"      ) : ["0",  ],
    ("muTau"  , "boost_low" ) : ["0",  ],
    ("muTau"  , "boost_high") : ["0",  ],
    ("muTau"  , "boost"     ) : ["0",  ],
    ("muTau"  , "btag_low"  ) : ["0",  ],
    ("muTau"  , "btag_high" ) : ["0",  ],
    ("muTau"  , "nobtag"    ) : ["1e-3", "0"],
    ("muTau"  , "btag"      ) : ["1e-3", "0"],
    ("muTau"  , "vbf"       ) : ["0",  ],
    ("eleTau" , "0jet_low"  ) : ["0",  ],
    ("eleTau" , "0jet_high" ) : ["0",  ],
    ("eleTau" , "0jet"      ) : ["0",  ],
    ("eleTau" , "boost_low" ) : ["0",  ],
    ("eleTau" , "boost_high") : ["0",  ],
    ("eleTau" , "boost"     ) : ["0",  ],
    ("eleTau" , "btag_low"  ) : ["0",  ],
    ("eleTau" , "btag_high" ) : ["0",  ],
    ("eleTau" , "nobtag"    ) : ["1e-3", "0"],
    ("eleTau" , "btag"      ) : ["1e-3", "0"],
    ("eleTau" , "vbf"       ) : ["0",  ],
    ("mumu"   , "0jet_low"  ) : ["0",  ],
    ("mumu"   , "0jet_high" ) : ["0",  ],
    ("mumu"   , "0jet"      ) : ["0",  ],
    ("mumu"   , "boost_low" ) : ["0",  ],
    ("mumu"   , "boost_high") : ["0",  ],
    ("mumu"   , "boost"     ) : ["0",  ],
    ("mumu"   , "btag_low"  ) : ["0",  ],
    ("mumu"   , "btag_high" ) : ["0",  ],
    ("mumu"   , "nobtag"    ) : ["1e-2", "0"],
    ("mumu"   , "btag"      ) : ["1e-3", "0"],
    ("mumu"   , "vbf"       ) : ["0",  ],  
    }

import os

type = options.type
categories = categories_sm if options.analysis == "sm" else categories_mssm

print "C R E A T I N G   S U M M A R Y   P L O T S :"
print "CATEGORIES:", categories
print "CHANNELS  :", channels

## combine 2011+2012
for chn in channels :
    for cat in categories :
        ## combine high and low pt categories, make sure in your 
        ## list that {CAT}_low and {CAT}_high are run beforehand
        if cat == "0jet" : ## or cat == "boost" : ## in boost low and high pt have a different binning
            #print "hadd {CHN}_{CAT}_{TYPE}_7+8TeV.root {CHN}_{CAT}_low_{TYPE}_7+8TeV.root {CHN}_{CAT}_high_{TYPE}_7+8TeV.root".format(CHN=chn, CAT=cat, TYPE=type)
            os.system("hadd -f {CHN}_{CAT}_{TYPE}_7+8TeV.root {CHN}_{CAT}_low_{TYPE}_7+8TeV.root {CHN}_{CAT}_high_{TYPE}_7+8TeV.root".format(CHN=chn, CAT=cat, TYPE=type))
        else :
            #print "hadd {CHN}_{CAT}_{TYPE}_7+8TeV.root {CHN}_{CAT}_{TYPE}_7TeV_{LOG}.root {CHN}_{CAT}_{TYPE}_8TeV_{LOG}.root".format(
            #    CHN=chn, CAT=cat, TYPE=type, LOG="LOG" if log[(chn, cat)] else "")
            os.system("hadd -f {CHN}_{CAT}_{TYPE}_7+8TeV.root {CHN}_{CAT}_{TYPE}_7TeV_*.root {CHN}_{CAT}_{TYPE}_8TeV_*.root".format(
                CHN=chn, CAT=cat, TYPE=type))#, LOG="LOG" if log[(chn, cat)]==True else ""))
## make plots
for chn in channels :
    for cat in categories :
        print chn, cat
        for i in range(len(log[chn,cat])) :
            os.system("root -l -q -b {CMSSW_BASE}/src/HiggsAnalysis/HiggsToTauTau/macros/postfit.C+\\(\\\"{CHN}_{CAT}_{TYPE}_7+8TeV.root\\\",\\\"{ANA}\\\",\\\"{LABEL}\\\",\\\"{EXTRA}\\\",{MIN},{MAX},{LOG}\)".format(CMSSW_BASE=os.environ['CMSSW_BASE'], CHN=chn, CAT=cat, TYPE=type, ANA=options.analysis.upper(), LABEL="2011+2012", EXTRA=extra[chn], MIN=min[chn,cat][i], MAX=max[chn,cat][i], LOG=log[chn,cat][i]))
        
