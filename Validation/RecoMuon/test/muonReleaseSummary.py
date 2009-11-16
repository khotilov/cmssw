#! /usr/bin/env python

import os
import sys
import fileinput
import string

NewVersion='3_3_3'
RefVersion='3_3_2'
NewRelease='CMSSW_'+NewVersion
RefRelease='CMSSW_'+RefVersion
#NewRelease='Summer09'
#RefRelease='Summer09_pre1'

NewCondition='MC'
RefCondition='MC'
#NewCondition='STARTUP'
#RefCondition='STARTUP'

NewFastSim=False
RefFastSim=False

if (NewCondition=='MC'):
    samples= ['RelValSingleMuPt10','RelValSingleMuPt100','RelValSingleMuPt1000','RelValTTbar']
    if (NewFastSim|RefFastSim):
        samples= ['RelValSingleMuPt10','RelValSingleMuPt100']
elif (NewCondition=='STARTUP'):
    samples= ['RelValTTbar','RelValZMM','RelValJpsiMM']
    if (NewFastSim|RefFastSim):
        samples= ['']
# These are some of the (pre)production samples, to be included by hand:
#samples= ['ppMuXLoose', 'InclusiveMu5_Pt50', 'InclusiveMu5_Pt250', 'ZmumuJet_Pt0to15', 'ZmumuJet_Pt300toInf', 'ZmumuJet_Pt80to120']
#samples= ['InclusiveMu5_Pt50', 'ZmumuJet_Pt0to15', 'ZmumuJet_Pt300toInf', 'ZmumuJet_Pt80to120']

Submit=True
Publish=False

GetFilesFromCastor=True
GetRefsFromCastor=True
CastorRepository = '/castor/cern.ch/cms/store/temp/dqm/offline/harvesting_output/mc/relval'
#CastorRepository = '/castor/cern.ch/user/n/nuno/relval/harvest'
#CastorRepository = '/castor/cern.ch/user/n/nuno/preproduction/harvest'
#CastorRepository = '/castor/cern.ch/user/j/jhegeman/preproduction_summer09/3_1_2'

ValidateHLT=True


if (NewFastSim):
    NewTag = NewCondition+'_noPU_ootb_FSIM'
    NewLabel=NewCondition+'_31X_V9_FastSim'
    NewFormat='GEN-SIM-DIGI-RECO'
else:
    NewTag = NewCondition+'_noPU_ootb'
    NewLabel=NewCondition+'_31X_V9'
#    NewLabel=NewCondition+'31X_V3_preproduction_312'
    if (NewCondition=='STARTUP'):
        NewLabel=NewCondition+'31X_V8'
    NewFormat='GEN-SIM-RECO'

if (RefFastSim):
    RefTag = RefCondition+'_noPU_ootb_FSIM'
    RefLabel=RefCondition+'_31X_FastSim'
    RefFormat='GEN-SIM-DIGI-RECO'
else:
    RefTag = RefCondition+'_noPU_ootb'
    RefLabel=RefCondition+'_31X_V8'
#    RefLabel=RefCondition+'_31X_V2_preproduction_311'
    if (RefCondition=='STARTUP'):
        RefLabel=RefCondition+'31X_V8'
    RefFormat='GEN-SIM-RECO'

NewLabel=NewLabel+'-v1'
RefLabel=RefLabel+'-v1'

if (NewFastSim):
    NewCondition=NewCondition+'_FSIM'
if (RefFastSim):
    RefCondition=RefCondition+'_FSIM'


NewRepository = '/afs/cern.ch/cms/Physics/muon/CMSSW/Performance/RecoMuon/Validation/val'
RefRepository = '/afs/cern.ch/cms/Physics/muon/CMSSW/Performance/RecoMuon/Validation/val'
CastorRefRepository = '/castor/cern.ch/user/a/aperrott/ValidationRecoMuon'


macro='macro/TrackValHistoPublisher.C'

def replace(map, filein, fileout):
    replace_items = map.items()
    while 1:
        line = filein.readline()
        if not line: break
        for old, new in replace_items:
            line = string.replace(line, old, new)
        fileout.write(line)
    fileout.close()
    filein.close()

###############################

for sample in samples :

    newdir=NewRepository+'/'+NewRelease+'/'+NewTag+'/'+sample 


    if(os.path.exists(NewRelease+'/'+NewTag+'/'+sample)==False):
        os.makedirs(NewRelease+'/'+NewTag+'/'+sample)

    if(os.path.exists(RefRelease+'/'+RefTag+'/'+sample)==False):
        os.makedirs(RefRelease+'/'+RefTag+'/'+sample)

    checkFile = NewRelease+'/'+NewTag+'/'+sample+'/globalMuons_tpToGlbAssociation.pdf'
    if (NewFastSim):
        checkFile = NewRelease+'/'+NewTag+'/'+sample+'/globalMuons_tpToGlbAssociationFS.pdf'
    if (os.path.isfile(checkFile)==True):
        print "Files of type "+checkFile+' exist alredy: delete them first, if you really want to overwrite them'
    else:
        newSample=NewRepository+'/'+NewRelease+'/'+NewTag+'/'+sample+'/'+'val.'+sample+'.root'
        refSample=RefRepository+'/'+RefRelease+'/'+RefTag+'/'+sample+'/'+'val.'+sample+'.root'

        if (os.path.isfile(NewRelease+'/'+NewTag+'/'+sample+'/val'+sample+'.root')==True):
            # FOR SOME REASON THIS DOES NOT WORK: to be checked...
            print "New file found at: ",NewRelease+'/'+NewTag+'/'+sample+'/val'+sample+'.root'+' -> Use that one'

        elif (GetFilesFromCastor):
# Check the number of events in the harvested samples, needed to retrieve the path on castor
            if (NewFastSim):
                NEVT='27000'
            else:
                NEVT='9000'
                if (sample=='RelValSingleMuPt10'):
                    NEVT='25000'
                elif(sample=='RelValZMM'):
                    NEVT='8995'
                elif((sample=='RelValTTbar')&(NewCondition=='STARTUP')):
                    NEVT='34000'
            os.system('rfcp '+CastorRepository+'/'+NewVersion+'/'+sample+'__'+NewRelease+'-'+NewLabel+'__'+NewFormat+'/run_1/nevents_'+NEVT+'/DQM_V0001_R000000001__'+sample+'__'+NewRelease+'-'+NewLabel+'__'+NewFormat+'_1.root '+NewRelease+'/'+NewTag+'/'+sample+'/'+'val.'+sample+'.root')
#            os.system('rfcp '+CastorRepository+'/'+NewRelease+'/DQM_V0001_R000000001__'+sample+'__'+NewRelease+'-'+NewLabel+'__'+NewFormat+'.root '+NewRelease+'/'+NewTag+'/'+sample+'/'+'val.'+sample+'.root')
#preprod-hegeman              os.system('rfcp '+CastorRepository+'/DQM_V0001_R000000001__'+sample+'__'+NewRelease+'-'+NewLabel+'__'+NewFormat+'_1.root '+NewRelease+'/'+NewTag+'/'+sample+'/'+'val.'+sample+'.root')

        elif (os.path.isfile(newSample)) :
            os.system('cp '+newSample+' '+NewRelease+'/'+NewTag+'/'+sample)

             
        if (os.path.isfile(RefRelease+'/'+RefTag+'/'+sample+'/val'+sample+'.root')!=True and os.path.isfile(refSample)):
            print '*** Getting reference file from '+RefRelease
            os.system('cp '+refSample+' '+RefRelease+'/'+RefTag+'/'+sample)
        elif (GetRefsFromCastor):
            print '*** Getting reference file from castor'
            os.system('rfcp '+CastorRefRepository+'/'+RefRelease+'_'+RefCondition+'_'+sample+'_val.'+sample+'.root '+RefRelease+'/'+RefTag+'/'+sample+'/'+'val.'+sample+'.root')
        else:
            print '*** WARNING: no reference file was found'

        cfgFileName=sample+'_'+NewRelease+'_'+RefRelease
        hltcfgFileName='HLT'+sample+'_'+NewRelease+'_'+RefRelease

        if os.path.isfile(RefRelease+'/'+RefTag+'/'+sample+'/val.'+sample+'.root'):
            replace_map_RECO = { 'DATATYPE': 'RECO', 'NEW_FILE':NewRelease+'/'+NewTag+'/'+sample+'/val.'+sample+'.root', 'REF_FILE':RefRelease+'/'+RefTag+'/'+sample+'/val.'+sample+'.root', 'REF_LABEL':sample, 'NEW_LABEL': sample, 'REF_RELEASE':RefRelease, 'NEW_RELEASE':NewRelease, 'REFSELECTION':RefTag, 'NEWSELECTION':NewTag, 'TrackValHistoPublisher': cfgFileName}
            if (ValidateHLT):
                replace_map_HLT = { 'DATATYPE': 'HLT', 'NEW_FILE':NewRelease+'/'+NewTag+'/'+sample+'/val.'+sample+'.root', 'REF_FILE':RefRelease+'/'+RefTag+'/'+sample+'/val.'+sample+'.root', 'REF_LABEL':sample, 'NEW_LABEL': sample, 'REF_RELEASE':RefRelease, 'NEW_RELEASE':NewRelease, 'REFSELECTION':RefTag, 'NEWSELECTION':NewTag, 'TrackValHistoPublisher': hltcfgFileName}
        else:
            print "No reference file found at: ", RefRelease+'/'+RefTag+'/'+sample
            replace_map_RECO = { 'DATATYPE': 'RECO', 'NEW_FILE':NewRelease+'/'+NewTag+'/'+sample+'/val.'+sample+'.root', 'REF_FILE':NewRelease+'/'+NewTag+'/'+sample+'/val.'+sample+'.root', 'REF_LABEL':sample, 'NEW_LABEL': sample, 'REF_RELEASE':NewRelease, 'NEW_RELEASE':NewRelease, 'REFSELECTION':NewTag, 'NEWSELECTION':NewTag, 'TrackValHistoPublisher': cfgFileName}
            if (ValidateHLT):
                replace_map_HLT = { 'DATATYPE': 'HLT', 'NEW_FILE':NewRelease+'/'+NewTag+'/'+sample+'/val.'+sample+'.root', 'REF_FILE':NewRelease+'/'+NewTag+'/'+sample+'/val.'+sample+'.root', 'REF_LABEL':sample, 'NEW_LABEL': sample, 'REF_RELEASE':NewRelease, 'NEW_RELEASE':NewRelease, 'REFSELECTION':NewTag, 'NEWSELECTION':NewTag, 'TrackValHistoPublisher': hltcfgFileName}

        templatemacroFile = open(macro, 'r')
        macroFile = open(cfgFileName+'.C' , 'w' )
        replace(replace_map_RECO, templatemacroFile, macroFile)

        if (ValidateHLT):
            templatemacroFile = open(macro, 'r')
            hltmacroFile = open(hltcfgFileName+'.C' , 'w' )
            replace(replace_map_HLT, templatemacroFile, hltmacroFile)

        if(Submit):
            os.system('root -b -q -l '+ cfgFileName+'.C'+ '>  macro.'+cfgFileName+'.log')
            if (ValidateHLT):
                os.system('root -b -q -l '+ hltcfgFileName+'.C'+ '>  macro.'+hltcfgFileName+'.log')

        if(Publish):
            if(os.path.exists(newdir)==False):
                os.makedirs(newdir)
            os.system('rm '+NewRelease+'/'+NewTag+'/'+sample+'/val.'+sample+'.root')  
            os.system('scp -r '+NewRelease+'/'+NewTag+'/'+sample+'/* ' + newdir)
