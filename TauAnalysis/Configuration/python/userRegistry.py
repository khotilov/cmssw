import os

userSettings = {
    'squires' : {
        # The job ID to use per default.
        'current' : 'Run29',
        # List of runs
        'jobs' : {
        },
        'global' : {
            'drawOptions' : {
                'canvasSizeX' : 640,
                'canvasSizeY' : 800
            }
        }
    },
    'friis' : {
        # The job ID to use per default.
        'current' : 'RunSVTestApr01',
        # List of runs
        'jobs' : {
            'RunSVTestApr01' : {
                'AHtoMuTau' : {
                    'analysisFilePath' : "/castor/cern.ch/user/f/friis/RunSVTestApr01/",
                    'harvestingFilePath' : "/data2/friis/RunSVTestApr01",
                    'batchHarvest' : "/castor/cern.ch/user/f/friis/RunSVTestApr01_harvest/"
                }
            },
            'RunSVTestApr04' : {
                'AHtoMuTau' : {
                    'analysisFilePath' : "/castor/cern.ch/user/f/friis/RunSVTestApr04/",
                    'harvestingFilePath' : "/data2/friis/",
                    'batchHarvest' : "/castor/cern.ch/user/f/friis/RunSVTestApr04_harvest/",
                    'skimSource' : "RunSVTestApr01"
                }
            }
        },
        'global' : {
            'drawOptions' : {
                'canvasSizeX' : 640,
                'canvasSizeY' : 800
            },
            'harvestScripts' : '/tmp/friis/harvest_scripts',
        }
    },
    'veelken': {
        'current' : {
            'ZtoMuTau'          : '2011Jun12',
            'ZtoMuTau_tauIdEff' : '2011Jul23'
            #'ZtoMuTau_tauIdEff' : '2011Jul06_mauro'
        },
        'jobs' : {
            '2011Jul23' : {
                'ZtoMuTau_tauIdEff' : {
                    'analysisFilePath' : "/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/TauIdEffMeas/",
                    'harvestingFilePath' : "/data2/veelken/CMSSW_4_2_x/skims/ZtoMuTau_tauIdEff/",
                    'tmpFilePath' : "/data1/veelken/tmp/ZtoMuTau_tauIdEff/",
                    'batchHarvest' : "/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/harvesting/ZtoMuTau_tauIdEff/2011Jul23/"
                }
            },
            '2011Jul06_mauro' : {
                'ZtoMuTau_tauIdEff' : {
                    'analysisFilePath' : "/castor/cern.ch/user/m/mverzett/tagprobe/skims/cms42x/TauIdEffMeas_2011Jul06/",
                    'harvestingFilePath' : "/data2/veelken/CMSSW_4_2_x/skims/ZtoMuTau_tauIdEff/",
                    'tmpFilePath' : "/data1/veelken/tmp/ZtoMuTau_tauIdEff/",
                    'batchHarvest' : "/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/harvesting/ZtoMuTau_tauIdEff/2011Jul06_mauro/"
                }
            },
            '2011Jul06v2' : {
                'ZtoMuTau_tauIdEff' : {
                    'analysisFilePath' : "/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/TauIdEffMeas/",
                    'harvestingFilePath' : "/data2/veelken/CMSSW_4_2_x/skims/ZtoMuTau_tauIdEff/",
                    'tmpFilePath' : "/data1/veelken/tmp/ZtoMuTau_tauIdEff/",
                    'batchHarvest' : "/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/harvesting/ZtoMuTau_tauIdEff/2011Jul06v2/"
                }
            },
            '2011Jul01_mauro' : {
                'ZtoMuTau_tauIdEff' : {
                    'analysisFilePath' : "/castor/cern.ch/user/m/mverzett/tagprobe/skims/cms42x/TauIdEffMeas_2011Jul01_v2",
                    'harvestingFilePath' : "/data2/veelken/CMSSW_4_2_x/skims/ZtoMuTau_tauIdEff/",
                    'tmpFilePath' : "/data1/veelken/tmp/ZtoMuTau_tauIdEff/",
                    'batchHarvest' : "/castor/cern.ch/user/m/mverzett/tagprobe/skims/cms42x/TauIdEffMeas_Harvest_2011Jul01_v2",
                }
            },
            '2011Jun25' : {
                'ZtoMuTau_tauIdEff' : {
                    'analysisFilePath' : "/castor/cern.ch/user/v/veelken/CMSSW_4_1_x/skims/TauIdEffMeas/",
                    'harvestingFilePath' : "/data2/veelken/CMSSW_4_1_x/skims/ZtoMuTau_tauIdEff/",
                    'tmpFilePath' : "/data1/veelken/tmp/ZtoMuTau_tauIdEff/",
                    'batchHarvest' : "/castor/cern.ch/user/m/mverzett/tagprobe/skims/TauIdEffMeas_Harvested_Jun06/"
                }
            }, 
            '2011Jun12' : {
                'ZtoMuTau' : {
                    'analysisFilePath' : "/castor/cern.ch/user/v/veelken/CMSSW_4_1_x/plots/ZtoMuTau/",
                    'harvestingFilePath' : "/data1/veelken/CMSSW_4_1_x/plots/ZtoMuTau/",
                    'tmpFilePath' : "/data2/veelken/tmp/ZtoMuTau/",
                    'jobId' : "2011Jun12",
                    'batchHarvest' : "/castor/cern.ch/user/v/veelken/CMSSW_4_1_x/harvesting/ZtoMuTau/2011Jun12/"
                }
            }
        },
        'global' : {
            'drawOptions' : {
                'canvasSizeX' : 800,
                'canvasSizeY' : 640
            },
            'harvestScripts' : '/tmp/veelken/harvest_scripts',
        }
    },
    'jkolb': {
        'current' : 'Run31',
        'jobs' : {
            # Summer 42X MC + 42X ReReco and PR data: HPS tau ID eff: byLooseCombinedIsolationDeltaBetaCorr --> running at T3
            'Run32' : {
                'ZtoMuTau_bgEstTemplate' : {
                    'analysisFilePath' : "/afs/crc.nd.edu/user/j/jkolb/CMSSW_4_2_4_patch1/src/TauAnalysis/BgEstimationTools/test/",
                    'harvestingFilePath' : "/pscratch/ndcms/bestman/jkolb/TauResults",
                    'tmpFilePath' : "/pscratch/ndcms/bestman/jkolb/TauResults/tmp"
                },
               	'ZtoMuTau' : {
                    'analysisFilePath' : "/afs/crc.nd.edu/user/j/jkolb/CMSSW_4_2_4_patch1/src/TauAnalysis/Configuration/test/",
                    'harvestingFilePath' : "/pscratch/ndcms/bestman/jkolb/TauResults",
                    'tmpFilePath' : "/pscratch/ndcms/bestman/jkolb/TauResults/tmp"
                }
			},
            # Spring11 PR v1/2 data + MC: HPS tau ID eff: byLooseCombinedIsolationDeltaBetaCorr --> running at T3
            'Run31' : {
                'ZtoMuTau_bgEstTemplate' : {
                    'analysisFilePath' : "/afs/crc.nd.edu/user/j/jkolb/CMSSW_4_1_3/src/TauAnalysis/BgEstimationTools/test/",
                    'harvestingFilePath' : "/pscratch/ndcms/bestman/jkolb/TauResults",
                    'tmpFilePath' : "/pscratch/ndcms/bestman/jkolb/TauResults/tmp"
                },
               	'ZtoMuTau' : {
                    'analysisFilePath' : "/afs/crc.nd.edu/user/j/jkolb/CMSSW_4_1_3/src/TauAnalysis/Configuration/test/",
                    'harvestingFilePath' : "/pscratch/ndcms/bestman/jkolb/TauResults",
                    'tmpFilePath' : "/pscratch/ndcms/bestman/jkolb/TauResults/tmp"
                }
			},
            # Spring11 PR v1/2 data + MC: HPS tau ID eff: byLooseIsolation --> running at T3
            'Run30' : {
                'ZtoMuTau_bgEstTemplate' : {
                    'analysisFilePath' : "/afs/crc.nd.edu/user/j/jkolb/CMSSW_4_1_3/src/TauAnalysis/BgEstimationTools/test/",
                    'harvestingFilePath' : "/pscratch/ndcms/bestman/jkolb/TauResults",
                    'tmpFilePath' : "/pscratch/ndcms/bestman/jkolb/TauResults/tmp"
                },
               	'ZtoMuTau' : {
                    'analysisFilePath' : "/afs/crc.nd.edu/user/j/jkolb/CMSSW_4_1_3/src/TauAnalysis/Configuration/test/",
                    'harvestingFilePath' : "/pscratch/ndcms/bestman/jkolb/TauResults",
                    'tmpFilePath' : "/pscratch/ndcms/bestman/jkolb/TauResults/tmp"
                }
			},
            # Spring11 PR v1/2 data + MC: HPS tau ID eff: byLooseIsolation
            'Run29' : {
                'ZtoMuTau_bgEstTemplate' : {
                    'analysisFilePath' : "/user/j/jkolb/TauResults/Run29/",
                    'harvestingFilePath' : "/data/ndpc0/c/jkolb/TauResults",
                    'tmpFilePath' : "/data/ndpc0/c/jkolb/TauResults/tmp/"
                },
                'ZtoMuTau' : {
                    'analysisFilePath' : "/user/j/jkolb/TauResults/Run29/",
                    'harvestingFilePath' : "/data/ndpc0/c/jkolb/TauResults",
                    'tmpFilePath' : "/data/ndpc0/c/jkolb/TauResults/tmp/"
                }
            },
            # (retry) Spring11 PR v1/2 data + MC: HPS tau ID eff: byLooseCombinedIsolationDeltaBetaCorr
            'Run28' : {
                'ZtoMuTau_bgEstTemplate' : {
                    'analysisFilePath' : "/user/j/jkolb/TauResults/Run28/",
                    'harvestingFilePath' : "/data/ndpc0/c/jkolb/TauResults",
                    'tmpFilePath' : "/data/ndpc0/c/jkolb/TauResults/tmp/"
                },
                'ZtoMuTau' : {
                    'analysisFilePath' : "/user/j/jkolb/TauResults/Run28/",
                    'harvestingFilePath' : "/data/ndpc0/c/jkolb/TauResults",
                    'tmpFilePath' : "/data/ndpc0/c/jkolb/TauResults/tmp/"
                }
            },
            # Spring11 PR v1/2 data + MC: HPS tau ID eff: byLooseCombinedIsolationDeltaBetaCorr
            'Run27' : {
                'ZtoMuTau_bgEstTemplate' : {
                    'analysisFilePath' : "/user/j/jkolb/TauResults/Run27/",
                    'harvestingFilePath' : "/data/ndpc0/c/jkolb/TauResults",
                    'tmpFilePath' : "/data/ndpc0/c/jkolb/TauResults/tmp/"
                },
                'ZtoMuTau' : {
                    'analysisFilePath' : "/user/j/jkolb/TauResults/Run27/",
                    'harvestingFilePath' : "/data/ndpc0/c/jkolb/TauResults",
                    'tmpFilePath' : "/data/ndpc0/c/jkolb/TauResults/tmp/"
                }
            },
            # Spring11 data + MC: may10 rereco + PR -- move to 4_2_X
            'Run25' : {
                'AHtoElecTau' : {
                    'analysisFilePath' : "/user/j/jkolb/TauResults/Run25/",
                    'harvestingFilePath' : "/data/ndpc0/c/jkolb/TauResults",
                    'tmpFilePath' : "/data/ndpc0/c/jkolb/TauResults/tmp/"
                },
                'ZtoMuTau' : {
                    'analysisFilePath' : "/user/j/jkolb/TauResults/Run25/",
                    'harvestingFilePath' : "/data/ndpc0/c/jkolb/TauResults",
                    'tmpFilePath' : "/data/ndpc0/c/jkolb/TauResults/tmp/"
                },
                'ZtoMuTau_bgEstTemplate' : {
                    'analysisFilePath' : "/user/j/jkolb/TauResults/Run25/",
                    'harvestingFilePath' : "/data/ndpc0/c/jkolb/TauResults",
                    'tmpFilePath' : "/data/ndpc0/c/jkolb/TauResults/tmp/"
                },

            },
            # Spring11 data + MC: 165/pb  -- electron pt cut to 18 GeV
            'Run24' : {
                'AHtoElecTau' : {
                    'analysisFilePath' : "/user/j/jkolb/TauResults/Run24/",
                    'harvestingFilePath' : "/data/ndpc0/c/jkolb/TauResults",
                    'tmpFilePath' : "/data/ndpc0/c/jkolb/TauResults/tmp/"
                }
            },
            # Spring11 data + MC: xxx/pb  -- add deltaBeta PU correction, removed Mt(e+MET) cut, added PzetaDiff, added vertexMult reweight with custom LUT
            'Run23' : {
                'AHtoElecTau' : {
                    'analysisFilePath' : "/user/j/jkolb/TauResults/Run23/AHtoElecTau/",
                    'harvestingFilePath' : "/data/ndpc0/c/jkolb/TauResults",
                    'tmpFilePath' : "/data/ndpc0/c/jkolb/TauResults/tmp/"
                }
            }
        },
        'global' : {
            'drawOptions' : {
                'canvasSizeX' : 800,
                'canvasSizeY' : 640
            },
            'harvestScripts' : '/tmp/jkolb/harvest_scripts',
        }
    },
    'liis': {
        'current' : 'RunPileup',
        'jobs' : {
            'RunPileup' : {
                'WtoTauNu' : {
                    'analysisFilePath' : "/castor/cern.ch/user/l/liis/CMSSW_38X/Histograms/Pileup/run2",
                    'harvestingFilePath' : "/afs/cern.ch/user/l/liis/scratch0",
                    'tmpFilePath' : "/tmp/liis",
                    'jobId' : "RunTest",
                    }
                },
            'Run01' : {
                'WtoTauNu' : {
                    'analysisFilePath' : "/castor/cern.ch/user/l/liis/CMSSW_38X/Histograms/HPSTau",
                    'harvestingFilePath' : "/tmp/liis",
                    'tmpFilePath' : "/tmp/liis",
                    'jobID' : "Run01",
                    }
                }
            },
        'global' : {
            'drawOptions' : {
                'canvasSizeX' : 800,
                'canvasSizeY' : 640
                },
#            'harvestScripts' : 'tmp/liis/harvest_scripts'
            }
        },
        'lusito': {
        'current' : {
            'ZtoMuTau'               : '2011Feb15',
            'ZtoMuTau_bgEstTemplate' : '2011Feb15'
        },
        'jobs' : {
            '2011Feb15' : {
                'ZtoMuTau' : {
                    'analysisFilePath' : "/castor/cern.ch/user/l/lusito/CMSSW_3_8_x/plots15feb/ZtoMuTau/",
                    'harvestingFilePath' : "/tmp/lusito/CMSSW_3_8_x/plots15feb/ZtoMuTau_tmp/",
                    'tmpFilePath' : "/tmp/lusito/CMSSW_3_8_x/plots15feb/ZtoMuTau/",
                     # Directory containing selected events
                    #'pickevents' : '/castor/cern.ch/user/f/friis/fixme_delete/'
                },
                'ZtoMuTau_bgEstTemplate' : {
                    'analysisFilePath' : "/castor/cern.ch/user/l/lusito/CMSSW_3_8_x/plots15feb/ZtoMuTau_bgEstTemplate/",
                    'harvestingFilePath' : "/tmp/lusito/CMSSW_3_8_x/plots15feb/ZtoMuTau_bgEstTemplate_tmp/",
                    'tmpFilePath' : "/tmp/lusito/CMSSW_3_8_x/plots15feb/ZtoMuTau_bgEstTemplate/",
                     # Directory containing selected events
                    #'pickevents' : '/castor/cern.ch/user/f/friis/fixme_delete/'
                }
            }
          }
        },
        'lantonel': {
        'current' : 'Run01',
        
        'jobs' : {
           'Run01' : {
                'ZtoElecTau_bgEstTemplate' : {
                    'analysisFilePath' : "/afs/crc.nd.edu/user/l/lantonel/CMSSW_3_8_7_patch1/src/TauAnalysis/BgEstimationTools/test/",
                    'harvestingFilePath' : "/afs/crc.nd.edu/user/l/lantonel/CMSSW_3_8_7_patch1/src/TauAnalysis/BgEstimationTools/test/",
                    'tmpFilePath' : "/afs/crc.nd.edu/user/l/lantonel/CMSSW_3_8_7_patch1/src/TauAnalysis/BgEstimationTools/test/"
                 },
            }
         },
        
        'global' : {
            'drawOptions' : {
                'canvasSizeX' : 800,
                'canvasSizeY' : 640
            },
            'harvestScripts' : '/tmp/lantonel/harvest_scripts',
        }
        
        },    

    

    }

def mine():
    return userSettings[os.environ['LOGNAME']]

def getJobId(channel):
    # Get the current job id for this channel
    jobIdInfo = mine()['current']
    # Check if we need to remap the job ids for different channels.
    if isinstance(jobIdInfo, dict):
        return jobIdInfo[channel]
    else:
        return jobIdInfo

def check_slash(dir):
    " Make sure a directory has a trailing slash "
    if not dir.endswith('/'):
        dir += '/'
    return dir

def getInfo(channel, jobid = None):
    # Get the current job id for this channel
    if jobid is None:
        jobid = getJobId(channel)
    # Get the channel information for the jobId
    return mine()['jobs'][jobid][channel]

def overrideJobId(channel, jobId):
    # Set the job ID to something else temporarily
    mine()['current'] = jobId

def getPickEventsPath(channel):
    return getInfo(channel)['pickevents']

def getHarvestScriptLocation():
    " Get location of harvest scripts "
    return mine()['global']['harvestScripts']

def getSkimEvents(channel, jobid = None):
    " Get the directory on castor that containing the skimmed files "
    skimJobId = getInfo(channel, jobid)['skimSource']
    # Get the location of the harvested & merged skims
    #output_path = os.path.normpath(
    #    '/castor/cern.ch/' + getInfo(channel, skimJobId)['batchHarvest'])
    output_path = os.path.normpath(
         getInfo(channel, skimJobId)['batchHarvest'])
    return check_slash(output_path)

def getAnalysisFilePath(channel):
    return getInfo(channel)['analysisFilePath']

def getHarvestingFilePath(channel, jobid=None):
    user_settings = getInfo(channel, jobid)
    if jobid is None:
        jobid = getJobId(channel)
    harvestingFilePath = os.path.join(
        user_settings['harvestingFilePath'], jobid)
    if not os.path.exists(harvestingFilePath):
        try:
            os.makedirs(harvestingFilePath)
        except:
            print "Failed to make harvesting file path: %s" % harvestingFilePath
    if not os.path.isdir(harvestingFilePath):
        print "WARNING: Harvesting file path %s is not a directory!" % harvestingFilePath
    return check_slash(harvestingFilePath)

def makeSkimStatFileMapper(channel, jobid=None):
    skimJobId = getInfo(channel, jobid)['skimSource']
    print skimJobId
    skim_path = getHarvestingFilePath(channel, skimJobId)
    print skim_path
    def mapper(sample):
        return os.path.join(skim_path, 'harvested_%s_%s_%s.root' %
                            (channel, sample, skimJobId))
    return mapper

def getLocalHarvestingFilePath(channel):
    return os.path.join(getHarvestingFilePath(channel), 'local')

def getTmpFilePath(channel):
    return getInfo(channel)['tmpFilePath']

def getBatchHarvestLocation(channel):
    " Where to store the output histograms when harvesting on LXBatch "
    output = getInfo(channel)['batchHarvest']
    # Add a trailing slash if we don't have one
    return check_slash(output)

def getConfigFileName(channel):
    return getInfo(channel)['configFileName']
