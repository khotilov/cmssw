#!/usr/bin/env python

import urllib
import string
import os
import sys
import LaunchOnCondor
import glob

def ComputeLimits(InputPattern, syst):
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino300_2C"' , '"Gluino300"', 0.0    / 0.3029 , 0.0    / 0.4955 , 1.0    / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino400_2C"' , '"Gluino400"', 0.0    / 0.3029 , 0.0    / 0.4955 , 1.0    / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino500_2C"' , '"Gluino500"', 0.0    / 0.3029 , 0.0    / 0.4955 , 1.0    / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino600_2C"' , '"Gluino600"', 0.0    / 0.3029 , 0.0    / 0.4955 , 1.0    / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino700_2C"' , '"Gluino700"', 0.0    / 0.3029 , 0.0    / 0.4955 , 1.0    / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino800_2C"' , '"Gluino800"', 0.0    / 0.3029 , 0.0    / 0.4955 , 1.0    / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino900_2C"' , '"Gluino900"', 0.0    / 0.3029 , 0.0    / 0.4955 , 1.0    / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1000_2C"', '"Gluino1000"',0.0    / 0.3029 , 0.0    / 0.4955 , 1.0    / 0.2015   , syst])

#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino300_f0"' , '"Gluino300"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino400_f0"' , '"Gluino400"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino500_f0"' , '"Gluino500"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino600_f0"' , '"Gluino600"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino700_f0"' , '"Gluino700"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino800_f0"' , '"Gluino800"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino900_f0"' , '"Gluino900"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1000_f0"', '"Gluino1000"',0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])

        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino300_f1"' , '"Gluino300"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino400_f1"' , '"Gluino400"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino500_f1"' , '"Gluino500"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino600_f1"' , '"Gluino600"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino700_f1"' , '"Gluino700"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino800_f1"' , '"Gluino800"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino900_f1"' , '"Gluino900"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1000_f1"', '"Gluino1000"',0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1100_f1"', '"Gluino1100"',0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])

        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino300_f5"' , '"Gluino300"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino400_f5"' , '"Gluino400"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino500_f5"' , '"Gluino500"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino600_f5"' , '"Gluino600"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino700_f5"' , '"Gluino700"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino800_f5"' , '"Gluino800"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino900_f5"' , '"Gluino900"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1000_f5"', '"Gluino1000"',0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1100_f5"', '"Gluino1100"',0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])

#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino300N_f0"', '"Gluino300N"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino400N_f0"', '"Gluino400N"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino500N_f0"', '"Gluino500N"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino600N_f0"', '"Gluino600N"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino700N_f0"', '"Gluino700N"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino800N_f0"', '"Gluino800N"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino900N_f0"', '"Gluino900N"', 0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])
#        launchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1000N_f0"','"Gluino1000N"',0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015   , syst])

        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino300N_f1"', '"Gluino300N"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino400N_f1"', '"Gluino400N"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino500N_f1"', '"Gluino500N"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino600N_f1"', '"Gluino600N"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino700N_f1"', '"Gluino700N"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino800N_f1"', '"Gluino800N"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino900N_f1"', '"Gluino900N"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1000N_f1"','"Gluino1000N"',0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1100N_f1"','"Gluino1100N"',0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])

#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino300N_f5"', '"Gluino300N"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino400N_f5"', '"Gluino400N"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino500N_f5"', '"Gluino500N"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
#        launchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino600N_f5"', '"Gluino600N"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino700N_f5"', '"Gluino700N"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino800N_f5"', '"Gluino800N"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino900N_f5"', '"Gluino900N"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1000N_f5"','"Gluino1000N"',0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015   , syst])

#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino600Z_f1"' , '"Gluino600Z"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino700Z_f1"' , '"Gluino700Z"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino800Z_f1"' , '"Gluino800Z"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015   , syst])

#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop130_2C"'   , '"Stop130"', 0.0    / 0.1705 , 0.0    / 0.4868 , 1.0    / 0.3427   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop200_2C"'   , '"Stop200"', 0.0    / 0.1705 , 0.0    / 0.4868 , 1.0    / 0.3427   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop300_2C"'   , '"Stop300"', 0.0    / 0.1705 , 0.0    / 0.4868 , 1.0    / 0.3427   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop400_2C"'   , '"Stop400"', 0.0    / 0.1705 , 0.0    / 0.4868 , 1.0    / 0.3427   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop500_2C"'   , '"Stop500"', 0.0    / 0.1705 , 0.0    / 0.4868 , 1.0    / 0.3427   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop600_2C"'   , '"Stop600"', 0.0    / 0.1705 , 0.0    / 0.4868 , 1.0    / 0.3427   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop700_2C"'   , '"Stop700"', 0.0    / 0.1705 , 0.0    / 0.4868 , 1.0    / 0.3427   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop800_2C"'   , '"Stop800"', 0.0    / 0.1705 , 0.0    / 0.4868 , 1.0    / 0.3427   , syst])

        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop130"'      , '"Stop130"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop200"'      , '"Stop200"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop300"'      , '"Stop300"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop400"'      , '"Stop400"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop500"'      , '"Stop500"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop600"'      , '"Stop600"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop700"'      , '"Stop700"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop800"'      , '"Stop800"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])

        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop130N"'     , '"Stop130N"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop200N"'     , '"Stop200N"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop300N"'     , '"Stop300N"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop400N"'     , '"Stop400N"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop500N"'     , '"Stop500N"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop600N"'     , '"Stop600N"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop700N"'     , '"Stop700N"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop800N"'     , '"Stop800N"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])

#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop300Z"'      , '"Stop300Z"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop400Z"'      , '"Stop400Z"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Stop500Z"'      , '"Stop500Z"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   , syst])

        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"GMStau100"', '"GMStau100"' ,-1, -1, -1  , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"GMStau126"', '"GMStau126"' ,-1, -1, -1  , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"GMStau156"', '"GMStau156"' ,-1, -1, -1  , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"GMStau200"', '"GMStau200"' ,-1, -1, -1  , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"GMStau247"', '"GMStau247"' ,-1, -1, -1  , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"GMStau308"', '"GMStau308"' ,-1, -1, -1  , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"GMStau370"', '"GMStau370"' ,-1, -1, -1  , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"GMStau432"', '"GMStau432"' ,-1, -1, -1  , syst])
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"GMStau494"', '"GMStau494"' ,-1, -1, -1  , syst])

#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"PPStau100"', '"PPStau100"' ,-1, -1, -1  , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"PPStau126"', '"PPStau126"' ,-1, -1, -1  , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"PPStau156"', '"PPStau156"' ,-1, -1, -1  , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"PPStau200"', '"PPStau200"' ,-1, -1, -1  , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"PPStau247"', '"PPStau247"' ,-1, -1, -1  , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"PPStau308"', '"PPStau308"' ,-1, -1, -1  , syst])

#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"DCStau121"', '"DCStau121"' ,-1, -1, -1  , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"DCStau182"', '"DCStau182"' ,-1, -1, -1  , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"DCStau242"', '"DCStau242"' ,-1, -1, -1  , syst])
#        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"DCStau302"', '"DCStau302"' ,-1, -1, -1  , syst])



if len(sys.argv)==1:
	print "Please pass in argument a number between 0 and 2"
        print "  0 - Submit the Core of the (TkOnly+TkTOF) Analysis     --> submitting 2x 1 jobs"
        print "  1 - Run the control plot macro                         --> submitting    0 jobs"
        print "  2 - Run the Optimization macro based on best Exp Limit --> submitting 2x75 jobs"
        print "  3 - Run the exclusion plot macro                       --> submitting    0 jobs"
	sys.exit()

elif sys.argv[1]=='0':	
        print 'ANALYSIS'
        FarmDirectory = "FARM"
        JobName = "HscpAnalysis"
	LaunchOnCondor.Jobs_RunHere = 1
	LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)	
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_DATA"'  , 0, '"dedxASmi"'  ,'"dedxHarm2"'  , '"combined"', 0.0, 0.0, 0.0, 35, 1.5]) #TkOnly
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', 0, '"dedxASmi"'  ,'"dedxHarm2"'  , '"combined"', 0.0, 0.0, 0.0, 35, 1.5]) #TkOnly
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_MC"'    , 0, '"dedxASmi"'  ,'"dedxHarm2"'  , '"combined"', 0.0, 0.0, 0.0, 35, 1.5]) #TkOnly
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_DATA"'  , 2, '"dedxASmi"'  ,'"dedxHarm2"'  , '"combined"', 0.0, 0.0, 0.0, 35, 1.5]) #TkTOF 
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', 2, '"dedxASmi"'  ,'"dedxHarm2"'  , '"combined"', 0.0, 0.0, 0.0, 35, 1.5]) #TkTOF
        LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_MC"'    , 2, '"dedxASmi"'  ,'"dedxHarm2"'  , '"combined"', 0.0, 0.0, 0.0, 35, 1.5]) #TkTOF
	LaunchOnCondor.SendCluster_Submit()

elif sys.argv[1]=='1':
        print 'PLOTTING'
	os.system('root Analysis_Step5.C++ -l -b -q')


elif sys.argv[1]=='2':
        print 'OPTIMIZATION'
        FarmDirectory = "FARM"
        JobName = "HscpLimits"
        LaunchOnCondor.Jobs_RunHere = 1
        LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
#        ComputeLimits('"Results/dedxASmi/combined/Eta25/PtMin35/Type0/"', '""')
#        ComputeLimits('"Results/dedxASmi/combined/Eta25/PtMin35/Type2/"', '""')
#        ComputeLimits('"Results/dedxASmi/combined/Eta20/PtMin35/Type0/"', '""')
#        ComputeLimits('"Results/dedxASmi/combined/Eta20/PtMin35/Type2/"', '""')
        ComputeLimits('"Results/dedxASmi/combined/Eta15/PtMin35/Type0/"', '""')
        ComputeLimits('"Results/dedxASmi/combined/Eta15/PtMin35/Type2/"', '""')

#        ComputeLimits('"Results/dedxASmi/combined/Eta25/PtMin35/Type0/"', '"_SystP"')
#        ComputeLimits('"Results/dedxASmi/combined/Eta25/PtMin35/Type2/"', '"_SystP"')
#        ComputeLimits('"Results/dedxASmi/combined/Eta25/PtMin35/Type0/"', '"_SystI"')
#        ComputeLimits('"Results/dedxASmi/combined/Eta25/PtMin35/Type2/"', '"_SystI"')
#        ComputeLimits('"Results/dedxASmi/combined/Eta25/PtMin35/Type0/"', '"_SystM"')
#        ComputeLimits('"Results/dedxASmi/combined/Eta25/PtMin35/Type2/"', '"_SystM"')
#        ComputeLimits('"Results/dedxASmi/combined/Eta25/PtMin35/Type0/"', '"_SystT"')
#        ComputeLimits('"Results/dedxASmi/combined/Eta25/PtMin35/Type2/"', '"_SystT"')
        LaunchOnCondor.SendCluster_Submit()


elif sys.argv[1]=='3':
        print 'EXCLUSION'
#        os.system('root Analysis_Step6.C++\'(\"tmp\")\' -l -b -q')
        os.system('sh Analysis_Step6.sh')
else:
	print 'Unknwon case: use an other argument or no argument to get help'



