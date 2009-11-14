#! /usr/bin/env python

import os, sys, optparse, math

copyargs = sys.argv[:]
for i in range(len(copyargs)):
    if copyargs[i] == "":
        copyargs[i] = "\"\""
commandline = " ".join(copyargs)

usage = """./%prog DIRNAME ITERATIONS INITIALGEOM INPUTFILES [options]

Creates (overwrites) a directory for each of the iterations and creates (overwrites)
submitJobs.sh with the submission sequence and dependencies.

DIRNAME        directories will be named DIRNAME01, DIRNAME02, etc.
ITERATIONS     number of iterations
INITIALGEOM    SQLite file containing muon geometry with tag names
               DTAlignmentRcd, DTAlignmentErrorRcd, CSCAlignmentRcd, CSCAlignmentErrorRcd
INPUTFILES     Python file defining 'fileNames', a list of input files as
               strings (create with findQualityFiles.py)"""

parser = optparse.OptionParser(usage)
parser.add_option("-j", "--jobs",
                   help="approximate number of \"gather\" subjobs",
                   type="int",
                   default=50,
                   dest="subjobs")
parser.add_option("-s", "--submitJobs",
                   help="alternate name of submitJobs.sh script (please include .sh extension); a file with this name will be OVERWRITTEN",
                   type="string",
                   default="submitJobs.sh",
                   dest="submitJobs")
parser.add_option("-b", "--big",
                  help="if invoked, subjobs will also be run on cmscaf1nd",
                  action="store_true",
                  dest="big")
parser.add_option("--mapplots",
                  help="if invoked, draw \"map plots\"",
                  action="store_true",
                  dest="mapplots")
parser.add_option("--segdiffplots",
                  help="if invoked, draw \"segment-difference plots\"",
                  action="store_true",
                  dest="segdiffplots")
parser.add_option("--globalTag",
                  help="GlobalTag for alignment/calibration conditions (typically all conditions except muon and tracker alignment)",
                  type="string",
                  default="CRAFT0831X_V1::All",
                  dest="globaltag")
parser.add_option("--trackerconnect",
                  help="connect string for tracker alignment (frontier://... or sqlite_file:...)",
                  type="string",
                  default="",
                  dest="trackerconnect")
parser.add_option("--trackeralignment",
                  help="name of TrackerAlignmentRcd tag",
                  type="string",
                  default="Alignments",
                  dest="trackeralignment")
parser.add_option("--trackerAPEconnect",
                  help="connect string for tracker APEs (frontier://... or sqlite_file:...)",
                  type="string",
                  default="",
                  dest="trackerAPEconnect")
parser.add_option("--trackerAPE",
                  help="name of TrackerAlignmentErrorRcd tag (tracker APEs)",
                  type="string",
                  default="AlignmentErrors",
                  dest="trackerAPE")
parser.add_option("--iscosmics",
                  help="if invoked, use cosmic track refitter instead of the standard one",
                  action="store_true",
                  dest="iscosmics")
parser.add_option("--station123params",
                  help="alignable parameters for DT stations 1, 2, 3 (see SWGuideAlignmentAlgorithms#Selection_of_what_to_align)",
                  type="string",
                  default="111111",
                  dest="station123params")
parser.add_option("--station4params",
                  help="alignable parameters for DT station 4",
                  type="string",
                  default="100011",
                  dest="station4params")
parser.add_option("--cscparams",
                  help="alignable parameters for CSC chambers",
                  type="string",
                  default="100011",
                  dest="cscparams")
parser.add_option("--minTrackPt",
                  help="minimum allowed track transverse momentum (in GeV)",
                  type="string",
                  default="100",
                  dest="minTrackPt")
parser.add_option("--maxTrackPt",
                  help="maximum allowed track transverse momentum (in GeV)",
                  type="string",
                  default="200",
                  dest="maxTrackPt")
parser.add_option("--minTrackerHits",
                  help="minimum number of tracker hits",
                  type="int",
                  default=15,
                  dest="minTrackerHits")
parser.add_option("--maxTrackerRedChi2",
                  help="maximum tracker chi^2 per degrees of freedom",
                  type="string",
                  default="10",
                  dest="maxTrackerRedChi2")
parser.add_option("--allowTIDTEC",
                  help="if invoked, allow tracks that pass through the tracker's !TID/!TEC region (recommended)",
                  action="store_true",
                  dest="allowTIDTEC")
parser.add_option("--twoBin",
                  help="if invoked, apply the \"two-bin method\" to control charge-antisymmetric errors",
                  action="store_true",
                  dest="twoBin")
parser.add_option("--weightAlignment",
                  help="if invoked, segments will be weighted by ndf/chi^2 in the alignment",
                  action="store_true",
                  dest="weightAlignment")
parser.add_option("--minAlignmentHits",
                  help="minimum number of hits required to align a chamber",
                  type="int",
                  default=30,
                  dest="minAlignmentHits")

if len(sys.argv) < 5:
    raise SystemError, "Too few arguments.\n\n"+parser.format_help()

DIRNAME = sys.argv[1]
ITERATIONS = int(sys.argv[2])
INITIALGEOM = sys.argv[3]
INPUTFILES = sys.argv[4]

options, args = parser.parse_args(sys.argv[5:])
mapplots = options.mapplots
segdiffplots = options.segdiffplots
globaltag = options.globaltag
trackerconnect = options.trackerconnect
trackeralignment = options.trackeralignment
trackerAPEconnect = options.trackerAPEconnect
trackerAPE = options.trackerAPE
iscosmics = str(options.iscosmics)
station123params = options.station123params
station4params = options.station4params
cscparams = options.cscparams
minTrackPt = options.minTrackPt
maxTrackPt = options.maxTrackPt
minTrackerHits = str(options.minTrackerHits)
maxTrackerRedChi2 = options.maxTrackerRedChi2
allowTIDTEC = str(options.allowTIDTEC)
twoBin = str(options.twoBin)
weightAlignment = str(options.weightAlignment)
minAlignmentHits = str(options.minAlignmentHits)

execfile(INPUTFILES)
stepsize = int(math.ceil(1.*len(fileNames)/options.subjobs))
pwd = str(os.getcwdu())

bsubfile = ["#!/bin/sh", ""]
bsubnames = []
last_align = None

directory = ""
for iteration in range(1, ITERATIONS+1):
    if iteration == 1:
        inputdb = INITIALGEOM
        inputdbdir = directory[:]
    else:
        inputdb = director + ".db"
        inputdbdir = directory[:]

    directory = "%s%02d/" % (DIRNAME, iteration)
    director = directory[:-1]
    os.system("rm -rf %s; mkdir %s" % (directory, directory))
    os.system("cp gather_cfg.py %s" % directory)
    os.system("cp align_cfg.py %s" % directory)

    bsubfile.append("cd %s" % directory)

    for jobnumber in range(options.subjobs):
        gather_fileName = "%sgather%03d.sh" % (directory, jobnumber)
        inputfiles = " ".join(fileNames[jobnumber*stepsize:(jobnumber+1)*stepsize])

        if mapplots or segdiffplots: copyplots = "plotting*.root"
        else: copyplots = ""

        copytrackerdb = ""
        if trackerconnect[0:12] == "sqlite_file:": copytrackerdb += "%s " % trackerconnect[12:]
        if trackerAPEconnect[0:12] == "sqlite_file:": copytrackerdb += "%s " % trackerAPEconnect[12:]

        if len(inputfiles) > 0:
            file(gather_fileName, "w").write("""#/bin/sh
# %(commandline)s

export ALIGNMENT_CAFDIR=`pwd`

cd %(pwd)s
eval `scramv1 run -sh`
export ALIGNMENT_AFSDIR=`pwd`

export ALIGNMENT_INPUTFILES='%(inputfiles)s'
export ALIGNMENT_ITERATION=%(iteration)d
export ALIGNMENT_JOBNUMBER=%(jobnumber)d
export ALIGNMENT_MAPPLOTS=%(mapplots)s
export ALIGNMENT_SEGDIFFPLOTS=%(segdiffplots)s
export ALIGNMENT_GLOBALTAG=%(globaltag)s
export ALIGNMENT_INPUTDB=%(inputdb)s
export ALIGNMENT_TRACKERCONNECT=%(trackerconnect)s
export ALIGNMENT_TRACKERALIGNMENT=%(trackeralignment)s
export ALIGNMENT_TRACKERAPECONNECT=%(trackerAPEconnect)s
export ALIGNMENT_TRACKERAPE=%(trackerAPE)s
export ALIGNMENT_ISCOSMICS=%(iscosmics)s
export ALIGNMENT_STATION123PARAMS=%(station123params)s
export ALIGNMENT_STATION4PARAMS=%(station4params)s
export ALIGNMENT_CSCPARAMS=%(cscparams)s
export ALIGNMENT_MINTRACKPT=%(minTrackPt)s
export ALIGNMENT_MAXTRACKPT=%(maxTrackPt)s
export ALIGNMENT_MINTRACKERHITS=%(minTrackerHits)s
export ALIGNMENT_MAXTRACKERREDCHI2=%(maxTrackerRedChi2)s
export ALIGNMENT_ALLOWTIDTEC=%(allowTIDTEC)s
export ALIGNMENT_TWOBIN=%(twoBin)s
export ALIGNMENT_WEIGHTALIGNMENT=%(weightAlignment)s
export ALIGNMENT_MINALIGNMENTHITS=%(minAlignmentHits)s

cp -f %(directory)sgather_cfg.py %(inputdbdir)s%(inputdb)s %(copytrackerdb)s $ALIGNMENT_CAFDIR/
cd $ALIGNMENT_CAFDIR/
ls -l
cmsRun gather_cfg.py
ls -l
cp -f *.tmp %(copyplots)s $ALIGNMENT_AFSDIR/%(directory)s
""" % vars())
            os.system("chmod +x %s" % gather_fileName)
            bsubfile.append("echo %sgather%03d.sh" % (directory, jobnumber))

            if last_align is None: waiter = ""
            else: waiter = "-w \"ended(%s)\"" % last_align            
            if options.big: queue = "cmscaf1nd"
            else: queue = "cmscaf1nh"

            bsubfile.append("bsub -q %s -J \"%s_gather%03d\" %s gather%03d.sh" % (queue, director, jobnumber, waiter, jobnumber))

            bsubnames.append("ended(%s_gather%03d)" % (director, jobnumber))

    file("%sconvert-db-to-xml_cfg.py" % directory, "w").write("""from Alignment.MuonAlignment.convertSQLitetoXML_cfg import *
process.PoolDBESSource.connect = \"sqlite_file:%(directory)s%(director)s.db\"
process.MuonGeometryDBConverter.outputXML.fileName = \"%(directory)s%(director)s.xml\"
process.MuonGeometryDBConverter.outputXML.relativeto = \"ideal\"
process.MuonGeometryDBConverter.outputXML.suppressDTChambers = False
process.MuonGeometryDBConverter.outputXML.suppressDTSuperLayers = True
process.MuonGeometryDBConverter.outputXML.suppressDTLayers = True
process.MuonGeometryDBConverter.outputXML.suppressCSCChambers = False
process.MuonGeometryDBConverter.outputXML.suppressCSCLayers = True

process.MuonGeometryDBConverter.getAPEs = True
process.PoolDBESSource.toGet = cms.VPSet(
    cms.PSet(record = cms.string(\"DTAlignmentRcd\"), tag = cms.string(\"DTAlignmentRcd\")),
    cms.PSet(record = cms.string(\"DTAlignmentErrorRcd\"), tag = cms.string(\"DTAlignmentErrorRcd\")),
    cms.PSet(record = cms.string(\"CSCAlignmentRcd\"), tag = cms.string(\"CSCAlignmentRcd\")),
    cms.PSet(record = cms.string(\"CSCAlignmentErrorRcd\"), tag = cms.string(\"CSCAlignmentErrorRcd\")),
      )
""" % vars())

    copytrackerdb = ""
    if trackerconnect[0:12] == "sqlite_file:": copytrackerdb += "%s " % trackerconnect[12:]
    if trackerAPEconnect[0:12] == "sqlite_file:": copytrackerdb += "%s " % trackerAPEconnect[12:]

    file("%salign.sh" % directory, "w").write("""#!/bin/sh
# %(commandline)s

export ALIGNMENT_CAFDIR=`pwd`

cd %(pwd)s
eval `scramv1 run -sh`
export ALIGNMENT_AFSDIR=`pwd`
export ALIGNMENT_INPUTDB=%(inputdb)s
export ALIGNMENT_ITERATION=%(iteration)d
export ALIGNMENT_GLOBALTAG=%(globaltag)s
export ALIGNMENT_INPUTDB=%(inputdb)s
export ALIGNMENT_TRACKERCONNECT=%(trackerconnect)s
export ALIGNMENT_TRACKERALIGNMENT=%(trackeralignment)s
export ALIGNMENT_TRACKERAPECONNECT=%(trackerAPEconnect)s
export ALIGNMENT_TRACKERAPE=%(trackerAPE)s
export ALIGNMENT_ISCOSMICS=%(iscosmics)s
export ALIGNMENT_STATION123PARAMS=%(station123params)s
export ALIGNMENT_STATION4PARAMS=%(station4params)s
export ALIGNMENT_CSCPARAMS=%(cscparams)s
export ALIGNMENT_MINTRACKPT=%(minTrackPt)s
export ALIGNMENT_MAXTRACKPT=%(maxTrackPt)s
export ALIGNMENT_MINTRACKERHITS=%(minTrackerHits)s
export ALIGNMENT_MAXTRACKERREDCHI2=%(maxTrackerRedChi2)s
export ALIGNMENT_ALLOWTIDTEC=%(allowTIDTEC)s
export ALIGNMENT_TWOBIN=%(twoBin)s
export ALIGNMENT_WEIGHTALIGNMENT=%(weightAlignment)s
export ALIGNMENT_MINALIGNMENTHITS=%(minAlignmentHits)s

cp -f %(directory)salign_cfg.py %(directory)sconvert-db-to-xml_cfg.py %(inputdbdir)s%(inputdb)s %(directory)s*.tmp  %(copytrackerdb)s $ALIGNMENT_CAFDIR/
cd $ALIGNMENT_CAFDIR/
export ALIGNMENT_ALIGNMENTTMP=`ls alignment*.tmp`

ls -l
cmsRun align_cfg.py
cp -f MuonAlignmentFromReference_report.py $ALIGNMENT_AFSDIR/%(directory)s%(director)s_report.py
cp -f MuonAlignmentFromReference_outputdb.db $ALIGNMENT_AFSDIR/%(directory)s%(director)s.db
cp -f MuonAlignmentFromReference_plotting.root $ALIGNMENT_AFSDIR/%(directory)s%(director)s.root

cd $ALIGNMENT_AFSDIR
cmsRun %(directory)sconvert-db-to-xml_cfg.py
""" % vars())
    os.system("chmod +x %salign.sh" % directory)

    bsubfile.append("echo %salign.sh" % directory)
    bsubfile.append("bsub -q cmscaf1nd -J \"%s_align\" -w \"%s\" align.sh" % (director, " && ".join(bsubnames)))
    bsubfile.append("cd ..")
    bsubfile.append("")
    bsubnames = []
    last_align = "%s_align" % director

file(options.submitJobs, "w").write("\n".join(bsubfile))
os.system("chmod +x %s" % options.submitJobs)
