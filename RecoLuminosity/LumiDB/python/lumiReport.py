import os,sys
from RecoLuminosity.LumiDB import tablePrinter, csvReporter,CommonUtil
from RecoLuminosity.LumiDB.wordWrappers import wrap_always, wrap_onspace, wrap_onspace_strict
def toScreenNorm(normdata):
    result=[]
    labels=[('Name','amode','E(GeV)','Norm')]
    print ' ==  = '
    for name,thisnorm in normdata.items():
        amodetag=str(thisnorm[0])
        normval='%.2f'%thisnorm[1]
        egev='%.0f'%thisnorm[2]
        result.append([name,amodetag,egev,normval])
    print tablePrinter.indent (labels+result, hasHeader = True, separateRows = False,prefix = '| ', postfix = ' |', justify = 'left',delim = ' | ', wrapfunc = lambda x: wrap_onspace (x,20) ) 

def toScreenTotDelivered(lumidata,resultlines,scalefactor,isverbose):
    '''
    inputs:
    lumidata {run:[lumilsnum,cmslsnum,timestamp,beamstatus,beamenergy,deliveredlumi,calibratedlumierror,(bxidx,bxvalues,bxerrs),(bxidx,b1intensities,b2intensities)]}
    resultlines [[resultrow1],[resultrow2],...,] existing result row
    '''
    result=[]
    totOldDeliveredLS=0
    totOldDelivered=0.0
    for r in resultlines:
        dl=0.0
        if(r[2]!='n/a'):            
            dl=float(r[2])#in /ub because it comes from file!
            (rr,lumiu)=CommonUtil.guessUnit(dl)
            r[2]='%.3f'%(rr)+' ('+lumiu+')'
        sls=0
        if(r[1]!='n/a'):
            sls=int(r[1])
        totOldDeliveredLS+=sls
        totOldDelivered+=dl
        if(r[4]!='n/a'):
            egv=float(r[4])
            r[4]='%.1f'%egv
        result.append(r)
    totls=0
    totdelivered=0.0
    totaltable=[]
    for run in lumidata.keys():
        lsdata=lumidata[run]
        if lsdata is None:
            result.append([str(run),'n/a','n/a','n/a','n/a'])
            if isverbose:
                result.extend(['n/a'])
            continue
        nls=len(lsdata)
        totls+=nls
        totlumi=sum([x[5] for x in lsdata])
        totdelivered+=totlumi
        (totlumival,lumiunit)=CommonUtil.guessUnit(totlumi)
        beamenergyPerLS=[float(x[4]) for x in lsdata]
        avgbeamenergy=0.0
        if len(beamenergyPerLS):
            avgbeamenergy=sum(beamenergyPerLS)/len(beamenergyPerLS)
        runstarttime=lsdata[0][2]
        if isverbose:
            selectedls=[(x[0],x[1]) for x in lsdata]
            result.append([str(run),str(nls),'%.3f'%(totlumival*scalefactor)+' ('+lumiunit+')',runstarttime.strftime("%m/%d/%y %H:%M:%S"),'%.1f'%(avgbeamenergy), str(selectedls)])
        else:
            result.append([str(run),str(nls),'%.3f'%(totlumival*scalefactor)+' ('+lumiunit+')',runstarttime.strftime("%m/%d/%y %H:%M:%S"),'%.1f'%(avgbeamenergy)])
    sortedresult=sorted(result,key=lambda x : int(x[0]))
    print ' ==  = '
    if isverbose:
        labels = [('Run', 'Total LS', 'Delivered','Start Time','E(GeV)','Selected LS')]
        print tablePrinter.indent (labels+sortedresult, hasHeader = True, separateRows = False,
                               prefix = '| ', postfix = ' |', justify = 'right',
                               delim = ' | ', wrapfunc = lambda x: wrap_onspace (x,20) )
    else:
        labels = [('Run', 'Total LS', 'Delivered','Start Time','E(GeV)')]
        print tablePrinter.indent (labels+result, hasHeader = True, separateRows = False,
                               prefix = '| ', postfix = ' |', justify = 'right',
                               delim = ' | ', wrapfunc = lambda x: wrap_onspace (x,40) )
    print ' ==  =  Total : '
    #if (totdelivered+totOldDelivered)!=0:
    (totalDeliveredVal,totalDeliveredUni)=CommonUtil.guessUnit(totdelivered+totOldDelivered)
    totrowlabels = [('Delivered LS','Delivered('+totalDeliveredUni+')')]
    totaltable.append([str(totls+totOldDeliveredLS),'%.3f'%(totalDeliveredVal*scalefactor)])
    print tablePrinter.indent (totrowlabels+totaltable, hasHeader = True, separateRows = False, prefix = '| ',
                               postfix = ' |', justify = 'right', delim = ' | ',
                               wrapfunc = lambda x: wrap_onspace (x, 20))
    
def toCSVTotDelivered(lumidata,filename,resultlines,scalefactor,isverbose):
    '''
    input:  {run:[lumilsnum,cmslsnum,timestamp,beamstatus,beamenergy,deliveredlumi,calibratedlumierror,(bxidx,bxvalues,bxerrs),(bxidx,b1intensities,b2intensities)]}
    '''
    result=[]
    r=csvReporter.csvReporter(filename)
    fieldnames = ['Run', 'Total LS', 'Delivered(/ub)','UTCTime','E(GeV)']
    if isverbose:
        fieldnames.append('Selected LS')
    for rline in resultlines:
        result.append(rline)
    for run in lumidata.keys():
        lsdata=lumidata[run]
        if lsdata is None:
            result.append([run,'n/a','n/a','n/a','n/a'])
            if isverbose:
                result.extend(['n/a'])
            continue
        nls=len(lsdata)
        totlumival=sum([x[5] for x in lsdata])
        beamenergyPerLS=[float(x[4]) for x in lsdata]
        avgbeamenergy=0.0
        if len(beamenergyPerLS):
            avgbeamenergy=sum(beamenergyPerLS)/len(beamenergyPerLS)
        runstarttime=lsdata[0][2]
        if isverbose:
            selectedls=[(x[0],x[1]) for x in lsdata]
            result.append([run,nls,totlumival*scalefactor,runstarttime.strftime("%m/%d/%y %H:%M:%S"),avgbeamenergy, str(selectedls)])
        else:
            result.append([run,nls,totlumival*scalefactor,runstarttime.strftime("%m/%d/%y %H:%M:%S"),avgbeamenergy])
    sortedresult=sorted(result,key=lambda x : int(x[0]))
    r.writeRow(fieldnames)
    r.writeRows(sortedresult)
#def toScreenLSDelivered(lumidata,resultlines,isverbose):
#    result=[]
#    for run in sorted(lumidata):
#        rundata=lumidata[run]
#        if rundata is None:
#            result.append([str(run),'n/a','n/a','n/a','n/a','n/a','n/a'])
#            continue
#        for lsdata in rundata:
#            if lsdata is None or len(lsdata)==0:
#                result.append([str(run),'n/a','n/a','n/a','n/a','n/a','n/a'])
#                continue
#            else:
#                lumils=lsdata[0]
#                cmsls=lsdata[1]
#                lsts=lsdata[2]
#                beamstatus=lsdata[3]
#                beamenergy=lsdata[4]
#                delivered=lsdata[5]
#                result.append([str(run),str(lumils),str(cmsls),'%.3f'%delivered,lsts.strftime('%m/%d/%y %H:%M:%S'),beamstatus,'%.1f'%beamenergy])
#    labels = [('Run','lumils','cmsls','Delivered(/ub)','UTCTime','Beam Status','E(GeV)')]
#    print tablePrinter.indent (labels+result, hasHeader = True, separateRows = False,
#                               prefix = '| ', postfix = ' |', justify = 'right',
#                               delim = ' | ', wrapfunc = lambda x: wrap_onspace (x,20) )
         
#def toCSVLSDelivered(lumidata,filename,resultlines,isverbose):
#    result=[]
#    fieldnames=['Run','lumils','cmsls','Delivered(/ub)','UTCTime','BeamStatus','E(GeV)']
#    r=csvReporter.csvReporter(filename)
#    for run in sorted(lumidata):
#        rundata=lumidata[run]
#        if rundata is None:
#            result.append([run,'n/a','n/a','n/a','n/a','n/a','n/a'])
#            continue
#        for lsdata in rundata:
#            if lsdata is None:
#                result.append([run,'n/a','n/a','n/a','n/a','n/a','n/a'])
#                continue
#            else:
#                lumils=lsdata[0]
#                cmsls=lsdata[1]
#                lsts=lsdata[2]
#                beamstatus=lsdata[3]
#                beamenergy=lsdata[4]
#                delivered=lsdata[5]
#                result.append([run,lumils,cmsls,delivered,lsts,beamstatus,beamenergy])
#    r.writeRow(fieldnames)
#    r.writeRows(result)

def toScreenOverview(lumidata,resultlines,scalefactor,isverbose):
    '''
    input:
    lumidata {run:[lumilsnum,cmslsnum,timestamp,beamstatus,beamenergy,deliveredlumi,recordedlumi,calibratedlumierror,(bxidx,bxvalues,bxerrs),(bxidx,b1intensities,b2intensities)]}
    resultlines [[resultrow1],[resultrow2],...,] existing result row
    '''
    result=[]
    labels = [('Run', 'Delivered LS', 'Delivered','Selected LS','Recorded')]
    totOldDeliveredLS=0
    totOldSelectedLS=0
    totOldDelivered=0.0
    totOldRecorded=0.0
    
    totaltable=[]
    totalDeliveredLS = 0
    totalSelectedLS = 0
    totalDelivered = 0.0
    totalRecorded = 0.0

    for r in resultlines:
        dl=0.0
        if(r[2]!='n/a'):            
            dl=float(r[2])#delivered in /ub because it comes from file!
            (rr,lumiu)=CommonUtil.guessUnit(dl)
            r[2]='%.3f'%(rr)+' ('+lumiu+')'
        dls=0
        if(r[1]!='n/a'):
            dls=int(r[1])
        totOldDeliveredLS+=dls
        totOldDelivered+=dl
        rls=0
        if(r[3]!='n/a'):
            rlsstr=r[3]
            listcomp=rlsstr.split(', ')
            for lstr in listcomp:
                enddigs=lstr[1:-1].split('-')
                lsmin=int(enddigs[0])
                lsmax=int(enddigs[1])
                rls=lsmax-lsmin+1
                totOldSelectedLS+=rls
        if(r[4]!='n/a'):
            rcd=float(r[4])#recorded in /ub because it comes from file!
            (rrcd,rlumiu)=CommonUtil.guessUnit(rcd)
            r[4]='%.3f'%(rrcd)+' ('+rlumiu+')'
        totOldRecorded+=rcd
        result.append(r)
    for run in lumidata.keys():
        lsdata=lumidata[run]
        if lsdata is None:
            result.append([str(run),'n/a','n/a','n/a','n/a'])
            continue
        nls=len(lsdata)
        deliveredData=[x[5] for x in lsdata]
        totdelivered=sum(deliveredData)
        totalDelivered+=totdelivered
        totalDeliveredLS+=len(deliveredData)
        (totdeliveredlumi,deliveredlumiunit)=CommonUtil.guessUnit(totdelivered)
        recordedData=[x[6] for x in lsdata if x[6] is not None]
        totrecorded=sum(recordedData)
        totalRecorded+=totrecorded
        (totrecordedlumi,recordedlumiunit)=CommonUtil.guessUnit(totrecorded)
        #print 'x[1] ',[x[1] for x in lsdata]
        selectedcmsls=[x[1] for x in lsdata if x[1]!=0]
        #print 'selectedcmsls ',selectedcmsls
        totalSelectedLS+=len(selectedcmsls)
        if len(selectedcmsls)==0:
            selectedlsStr='n/a'
        else:
            selectedlsStr = CommonUtil.splitlistToRangeString(selectedcmsls)
        result.append([str(run),str(nls),'%.3f'%(totdeliveredlumi*scalefactor)+' ('+deliveredlumiunit+')',selectedlsStr,'%.3f'%(totrecordedlumi*scalefactor)+' ('+recordedlumiunit+')'])
    sortedresult=sorted(result,key=lambda x : int(x[0]))    
    print ' ==  = '
    print tablePrinter.indent (labels+sortedresult, hasHeader = True, separateRows = False,
                               prefix = '| ', postfix = ' |', justify = 'right',
                               delim = ' | ', wrapfunc = lambda x: wrap_onspace (x,20) )
    print ' ==  =  Total : '
    (totalDeliveredVal,totalDeliveredUni)=CommonUtil.guessUnit(totalDelivered+totOldDelivered)
    (totalRecordedVal,totalRecordedUni)=CommonUtil.guessUnit(totalRecorded+totOldRecorded)
    totrowlabels = [('Delivered LS','Delivered('+totalDeliveredUni+')','Selected LS','Recorded('+totalRecordedUni+')')]
    totaltable.append([str(totalDeliveredLS+totOldDeliveredLS),'%.3f'%(totalDeliveredVal*scalefactor),str(totalSelectedLS+totOldSelectedLS),'%.3f'%(totalRecordedVal*scalefactor)])
    print tablePrinter.indent (totrowlabels+totaltable, hasHeader = True, separateRows = False, prefix = '| ',
                               postfix = ' |', justify = 'right', delim = ' | ',
                               wrapfunc = lambda x: wrap_onspace (x, 20))
    
def toCSVOverview(lumidata,filename,resultlines,scalefactor,isverbose):
    '''
    input:
    lumidata {run:[lumilsnum,cmslsnum,timestamp,beamstatus,beamenergy,deliveredlumi,recordedlumi,calibratedlumierror,(bxidx,bxvalues,bxerrs),(bxidx,b1intensities,b2intensities)]}
    resultlines [[resultrow1],[resultrow2],...,] existing result row
    '''
    result=[]
    fieldnames = ['Run', 'DeliveredLS', 'Delivered(/ub)','SelectedLS','Recorded(/ub)']
    r=csvReporter.csvReporter(filename)
    for rline in resultlines:
        result.append(rline)
        
    for run in lumidata.keys():
        lsdata=lumidata[run]
        if lsdata is None:
            result.append([run,'n/a','n/a','n/a','n/a'])
            continue
        nls=len(lsdata)
        deliveredData=[x[5] for x in lsdata]
        recordedData=[x[6] for x in lsdata if x[6] is not None]
        totdeliveredlumi=0.0
        totrecordedlumi=0.0
        if len(deliveredData)!=0:
            totdeliveredlumi=sum(deliveredData)
        if len(recordedData)!=0:
            totrecordedlumi=sum(recordedData)
        selectedcmsls=[x[1] for x in lsdata if x[1]!=0]
        if len(selectedcmsls)==0:
            selectedlsStr='n/a'
        else:
            selectedlsStr = CommonUtil.splitlistToRangeString(selectedcmsls)
        result.append([run,nls,totdeliveredlumi*scalefactor,selectedlsStr,totrecordedlumi*scalefactor])
    sortedresult=sorted(result,key=lambda x : int(x[0]))
    r.writeRow(fieldnames)
    r.writeRows(sortedresult)
def toScreenLumiByLS(lumidata,resultlines,scalefactor,isverbose):
    '''
    input:
    lumidata {run:[lumilsnum,cmslsnum,timestamp,beamstatus,beamenergy,deliveredlumi,recordedlumi,calibratedlumierror,(bxidx,bxvalues,bxerrs),(bxidx,b1intensities,b2intensities)]}
    resultlines [[resultrow1],[resultrow2],...,] existing result row
    '''
    result=[]
    labels = [ ('Run','LS','UTCTime','Beam Status','E(GeV)','Delivered(/ub)','Recorded(/ub)')]
    totalrow = []
    
    totalDeliveredLS = 0
    totalSelectedLS = 0
    totalDelivered = 0.0
    totalRecorded = 0.0

    totOldDeliveredLS = 0
    totOldSelectedLS = 0
    totOldDelivered = 0.0
    totOldRecorded = 0.0

    for rline in resultlines:
        myls=rline[1]
        if myls!='n/a':
            [luls,cmls]=myls.split(':')
            totOldDeliveredLS+=1
            if cmls!='0':
                totOldSelectedLS+=1
        dl=0.0
        if rline[5]!='n/a':
            dl=float(rline[5])#delivered in /ub 
            rline[5]='%.2f'%(dl)
            totOldDelivered+=dl
        if rline[6]!='n/a':
           rl=float(rline[6])#recorded in /ub
           rline[6]='%.2f'%(rl)
           totOldRecorded+=rl
        result.append(rline)
        
    for run in lumidata.keys():
        rundata=lumidata[run]
        if rundata is None:
            result.append([str(run),'n/a','n/a','n/a','n/a','n/a'])
            continue
        for lsdata in rundata:
            lumilsnum=lsdata[0]
            cmslsnum=lsdata[1]
            ts=lsdata[2]
            bs=lsdata[3]
            begev=lsdata[4]
            deliveredlumi=lsdata[5]
            recordedlumi=lsdata[6]
            result.append([str(run),str(lumilsnum)+':'+str(cmslsnum),ts.strftime('%m/%d/%y %H:%M:%S'),bs,'%.1f'%begev,'%.2f'%(deliveredlumi*scalefactor),'%.2f'%(recordedlumi*scalefactor)])
            totalDelivered+=deliveredlumi
            totalRecorded+=recordedlumi
            totalDeliveredLS+=1
            if(cmslsnum!=0):
                totalSelectedLS+=1
    totdeliveredlumi=0.0
    deliveredlumiunit='/ub'
    #if (totalDelivered+totOldDelivered)!=0:
    (totdeliveredlumi,deliveredlumiunit)=CommonUtil.guessUnit(totalDelivered+totOldDelivered)
    totrecordedlumi=0.0
    recordedlumiunit='/ub'
    #if (totalRecorded+totOldRecorded)!=0:
    (totrecordedlumi,recordedlumiunit)=CommonUtil.guessUnit(totalRecorded+totOldRecorded)
    lastrowlabels = [ ('Delivered LS','Selected LS', 'Delivered('+deliveredlumiunit+')', 'Recorded('+recordedlumiunit+')')]
    totalrow.append ([str(totalDeliveredLS+totOldDeliveredLS),str(totalSelectedLS+totOldSelectedLS),'%.3f'%(totdeliveredlumi*scalefactor),'%.3f'%(totrecordedlumi*scalefactor)])
    sortedresult=sorted(result,key=lambda x : int(x[0]))
    print ' ==  = '
    print tablePrinter.indent (labels+sortedresult, hasHeader = True, separateRows = False, prefix = '| ',
                               postfix = ' |', justify = 'right', delim = ' | ',
                               wrapfunc = lambda x: wrap_onspace_strict (x, 22))
    print ' ==  =  Total : '
    print tablePrinter.indent (lastrowlabels+totalrow, hasHeader = True, separateRows = False, prefix = '| ',
                               postfix = ' |', justify = 'right', delim = ' | ',
                               wrapfunc = lambda x: wrap_onspace (x, 20))    

                  
def toCSVLumiByLS(lumidata,filename,resultlines,scalefactor,isverbose):
    result=[]
    fieldnames=['Run','LS','UTCTime','Beam Status','E(GeV)','Delivered(/ub)','Recorded(/ub)']
    r=csvReporter.csvReporter(filename)

    for rline in resultlines:
        result.append(rline)
        
    for run in sorted(lumidata):
        rundata=lumidata[run]
        if rundata is None:
            result.append([run,'n/a','n/a','n/a','n/a','n/a'])
            continue
        for lsdata in rundata:
            lumilsnum=lsdata[0]
            cmslsnum=lsdata[1]
            ts=lsdata[2]
            bs=lsdata[3]
            begev=lsdata[4]
            deliveredlumi=lsdata[5]
            recordedlumi=lsdata[6]
            result.append([run,str(lumilsnum)+':'+str(cmslsnum),ts.strftime('%m/%d/%y %H:%M:%S'),bs,begev,deliveredlumi*scalefactor,recordedlumi*scalefactor])
    sortedresult=sorted(result,key=lambda x : int(x[0]))
    r.writeRow(fieldnames)
    r.writeRows(sortedresult)

def toScreenLSEffective(lumidata,resultlines,scalefactor,isverbose):
    '''
    input:  {run:[lumilsnum(0),cmslsnum(1),timestamp(2),beamstatus(3),beamenergy(4),deliveredlumi(5),recordedlumi(6),calibratedlumierror(7),{hltpath:[l1name,l1prescale,hltprescale,efflumi]},bxdata,beamdata]}
    '''
    result=[]#[run,ls,hltpath,l1bitname,hltpresc,l1presc,efflumi]
    for run in sorted(lumidata):#loop over runs
        rundata=lumidata[run]
        if rundata is None:
            result.append([str(run),'n/a','n/a','n/a','n/a','n/a','n/a','n/a'])
            continue
        for lsdata in rundata:
            efflumiDict=lsdata[8]# this ls has no such path?
            if not efflumiDict:
                continue
            cmslsnum=lsdata[1]
            recorded=lsdata[6]
            if not recorded:
                recorded=0.0
            for hltpathname in sorted(efflumiDict):
                pathdata=efflumiDict[hltpathname]
                l1name=pathdata[0]
                if l1name is None:
                    l1name='n/a'
                else:
                    l1name=l1name.replace('"','')
                l1prescale=pathdata[1]
                hltprescale=pathdata[2]
                lumival=pathdata[3]
                if lumival is not None:
                    result.append([str(run),str(cmslsnum),hltpathname,l1name,str(hltprescale),str(l1prescale),'%.2f'%(recorded*scalefactor),'%.2f'%(lumival*scalefactor)])
                else:
                    result.append([str(run),str(cmslsnum),hltpathname,l1name,str(hltprescale),str(l1prescale),'%.2f'%(recorded*scalefactor),'n/a'])
    labels = [('Run','LS','HLTpath','L1bit','HLTpresc','L1presc','Recorded(/ub)','Effective(/ub)')]
    print ' ==  = '
    print tablePrinter.indent (labels+result, hasHeader = True, separateRows = False,
                               prefix = '| ', postfix = ' |', justify = 'right',
                               delim = ' | ', wrapfunc = lambda x: wrap_onspace_strict(x,22) )
def toCSVLSEffective(lumidata,filename,resultlines,scalefactor,isverbose):
    '''
    input:  {run:[lumilsnum(0),cmslsnum(1),timestamp(2),beamstatus(3),beamenergy(4),deliveredlumi(5),recordedlumi(6),calibratedlumierror(7),{hltpath:[l1name,l1prescale,hltprescale,efflumi]},bxdata,beamdata]}
    '''
    result=[]#[run,ls,hltpath,l1bitname,hltpresc,l1presc,efflumi]
    r=csvReporter.csvReporter(filename)
    
    for run in sorted(lumidata):#loop over runs
        rundata=lumidata[run]
        if rundata is None:
            result.append([str(run),'n/a','n/a','n/a','n/a','n/a','n/a','n/a'])
            continue
        for lsdata in rundata:
            efflumiDict=lsdata[8]# this ls has no such path?
            if not efflumiDict:
                continue
            cmslsnum=lsdata[1]
            recorded=lsdata[6]
            if not recorded:
                recorded=0.0
            for hltpathname in sorted(efflumiDict):
                pathdata=efflumiDict[hltpathname]
                l1name=pathdata[0]
                if l1name is None:
                    l1name='n/a'
                else:
                    l1name=l1name.replace('"','')
                l1prescale=pathdata[1]
                hltprescale=pathdata[2]
                lumival=pathdata[3]
                if lumival is not None:
                    result.append([run,cmslsnum,hltpathname,l1name,hltprescale,l1prescale,recorded*scalefactor,lumival*scalefactor])
                else:
                    result.append([run,cmslsnum,hltpathname,l1name,hltprescale,l1prescale,recorded*scalefactor,'n/a'])
    fieldnames = ['Run','LS','HLTpath','L1bit','HLTpresc','L1presc','Recorded(/ub)','Effective(/ub)']
    r.writeRow(fieldnames)
    r.writeRows(result)

def toScreenTotEffective(lumidata,resultlines,scalefactor,isverbose):
    '''
    input:  {run:[lumilsnum(0),triggeredls(1),timestamp(2),beamstatus(3),beamenergy(4),deliveredlumi(5),recordedlumi(6),calibratedlumierror(7),{hltpath:[l1name,l1prescale,hltprescale,efflumi]},bxdata,beamdata](8)}
    screen Run,SelectedLS,Recorded,HLTPath,L1Bit,Effective
    '''
    result=[]#[run,selectedlsStr,recorded,hltpath,l1bit,efflumi]
    totdict={}#{hltpath:[nls,toteff]}
    hprescdict={}
    lprescdict={}
    alltotrecorded=0.0
    selectedcmsls=[]
    for run in sorted(lumidata):#loop over runs
        rundata=lumidata[run]
        if rundata is None:
            result.append([str(run),'n/a','n/a','n/a','n/a'])
            continue
        selectedcmsls=[x[1] for x in rundata if x[1]!=0]
        totrecorded=sum([x[6] for x in rundata if x[6] is not None])
        alltotrecorded+=totrecorded
        totefflumiDict={}
        pathmap={}#{hltpathname:1lname}
        for lsdata in rundata:
            cmslsnum=lsdata[1]
            efflumiDict=lsdata[8]# this ls has no such path?
            if not efflumiDict:
                if cmslsnum in selectedcmsls:
                    selectedcmsls.remove(cmslsnum)
                continue
            for hltpathname,pathdata in efflumiDict.items():
                if not totefflumiDict.has_key(hltpathname):
                    totefflumiDict[hltpathname]=0.0
                    pathmap[hltpathname]='n/a'
                l1name=pathdata[0]
                l1presc=pathdata[1]
                hltpresc=pathdata[2]
                lumival=pathdata[3]
                if not totdict.has_key(hltpathname):
                    totdict[hltpathname]=[0,0.0]
                if l1presc and hltpresc and l1presc*hltpresc!=0:
                    if not hprescdict.has_key(hltpathname):
                        hprescdict[hltpathname]=[]
                    hprescdict[hltpathname].append(hltpresc)
                    if not lprescdict.has_key(l1name):
                        lprescdict[l1name]=[]
                    lprescdict[l1name].append(l1presc)
                    totdict[hltpathname][0]+=1   
                    if lumival:
                        totdict[hltpathname][1]+=lumival
                        totefflumiDict[hltpathname]+=lumival
                        pathmap[hltpathname]=l1name
                else:
                    if cmslsnum in selectedcmsls:
                        selectedcmsls.remove(cmslsnum)
        if len(selectedcmsls)==0:
            selectedlsStr='n/a'
        else:
            selectedlsStr = CommonUtil.splitlistToRangeString(selectedcmsls)
        for name in sorted(totefflumiDict):
            (efflumival,efflumiunit)=CommonUtil.guessUnit(totefflumiDict[name])
            (totrecval,totrecunit)=CommonUtil.guessUnit(totrecorded)
            lname=pathmap[name]
            hprescs=list(set(hprescdict[hltpathname]))
            lprescs=list(set(lprescdict[lname]))
            hprescStr='('+','.join(['%d'%(x) for x in hprescs])+')'
            lprescStr='('+','.join(['%d'%(x) for x in lprescs])+')'
            #print 'efflumival , efflumiunit ',efflumival,efflumiunit
            result.append([str(run),selectedlsStr,'%.3f'%(totrecval*scalefactor)+'('+totrecunit+')',name+hprescStr,lname+lprescStr,'%.3f'%(efflumival*scalefactor)+'('+efflumiunit+')'])
    labels = [('Run','SelectedLS','Recorded','HLTpath','L1bit','Effective')]
    print ' ==  = '
    print tablePrinter.indent (labels+result, hasHeader = True, separateRows = False,
                               prefix = '| ', postfix = ' |', justify = 'right',
                               delim = ' | ', wrapfunc = lambda x: wrap_onspace_strict(x,22) )
    print ' ==  =  Total : '
    lastrowlabels=[('HLTPath','SelectedLS','Recorded','Effective')]
    totresult=[]
    (alltotrecval,alltotrecunit)=CommonUtil.guessUnit(alltotrecorded)
    for hname in sorted(totdict):
        hdata=totdict[hname]
        totnls=hdata[0]
        (toteffval,toteffunit)=CommonUtil.guessUnit(hdata[1])
        totresult.append([hname,str(totnls),'%.2f'%(alltotrecval*scalefactor)+'('+alltotrecunit+')','%.2f'%(toteffval*scalefactor)+'('+toteffunit+')'])
    print tablePrinter.indent (lastrowlabels+totresult, hasHeader = True, separateRows = False,
                               prefix = '| ', postfix = ' |', justify = 'right',
                               delim = ' | ', wrapfunc = lambda x: wrap_onspace (x,20) )
    
def toCSVTotEffective(lumidata,filename,resultlines,scalefactor,isverbose):
    result=[]#[run,hltpath,l1bitname,totefflumi]
    r=csvReporter.csvReporter(filename)
    for run in sorted(lumidata):#loop over runs
        rundata=lumidata[run]
        if rundata is None:
            result.append([str(run),'n/a','n/a','n/a'])
            continue
        selectedcmslsStr=[x[1] for x in rundata if x[1]!=0]
        totrecorded=sum([x[6] for x in rundata if x[6] is not None])
        totefflumiDict={}
        pathmap={}
        for lsdata in rundata:
            efflumiDict=lsdata[8]# this ls has no such path?
            if not efflumiDict:
                continue
            for hltpathname,pathdata in efflumiDict.items():
                if not totefflumiDict.has_key(hltpathname):
                    totefflumiDict[hltpathname]=0.0
                    pathmap[hltpathname]='n/a'                
                lumival=pathdata[3]
                if lumival:
                    totefflumiDict[hltpathname]+=lumival
                    if pathdata[0] is None:
                        pathmap[hltpathname]='n/a'
                    else:
                        pathmap[hltpathname]=pathdata[0].replace('\"','')
        for name in sorted(totefflumiDict):
            result.append([run,name,pathmap[name],totefflumiDict[name]*scalefactor])
    fieldnames=['Run','HLTpath','L1bit','Effective(/ub)']
    r.writeRow(fieldnames)
    r.writeRows(result)
def toCSVLumiByLSXing(lumidata,scalefactor,filename):
    '''
    input:{run:[lumilsnum(0),cmslsnum(1),timestamp(2),beamstatus(3),beamenergy(4),deliveredlumi(5),recordedlumi(6),calibratedlumierror(7),{hltpath:[l1name,l1prescale,hltprescale,efflumi]},bxdata,beamdata]}
    output:
    fieldnames=['Run','CMSLS','Delivered(/ub)','Recorded(/ub)','BX']
    '''
    result=[]
    fieldnames=['run','ls','delivered(/ub)','recorded(/ub)','bx']
    r=csvReporter.csvReporter(filename)
    for run in sorted(lumidata):
        rundata=lumidata[run]
        if rundata is None:
            result.append([run,'n/a','n/a','n/a','n/a'])
            continue
        for lsdata in rundata:
            cmslsnum=lsdata[1]
            deliveredlumi=lsdata[5]
            recordedlumi=lsdata[6]
            (bxidxlist,bxvaluelist,bxerrorlist)=lsdata[8]
            bxresult=[]
            if bxidxlist and bxvaluelist:
                bxinfo=CommonUtil.transposed([bxidxlist,bxvaluelist])
                bxresult=CommonUtil.flatten([run,cmslsnum,deliveredlumi*scalefactor,recordedlumi*scalefactor,bxinfo])
                result.append(bxresult)
            else:
                result.append([run,cmslsnum,deliveredlumi*scalefactor,recordedlumi*scalefactor])
    r.writeRow(fieldnames)
    r.writeRows(result)
    
def toScreenLSTrg(trgdata,iresults=None,isverbose=False):
    '''
    input:{run:[[cmslsnum,deadfrac,deadtimecount,bitzero_count,bitzero_prescale,[(name,count,presc),]],..]
    '''
    result=[]
    for run in trgdata.keys():
        if trgdata[run] is None:
            ll=[str(run),'n/a','n/a','n/a']
            if isverbose:
                ll.append('n/a')
            result.append(ll)
            continue
        perrundata=trgdata[run]
        deadfrac=0.0
        bitdataStr='n/a'
        for lsdata in perrundata:
            cmslsnum=lsdata[0]
            deadfrac=lsdata[1]
            deadcount=lsdata[2]
            bitdata=lsdata[5]# already sorted by name
            flatbitdata=["("+x[0]+',%d'%x[1]+',%d'%x[2]+")" for x in bitdata if x[0]!='False']
            #print 'flatbit ',flatbitdata
            bitdataStr=', '.join(flatbitdata)
            #print 'bitdataStr ',bitdataStr
            if isverbose:
                result.append([str(run),str(cmslsnum),'%.4f'%(deadfrac),'%d'%deadcount,bitdataStr])
            else:
                result.append([str(run),str(cmslsnum),'%.4f'%(deadfrac),'%d'%deadcount])
    print ' ==  = '
    if isverbose:
        labels = [('Run', 'LS', 'dfrac','dcount','(bit,count,presc)')]
    else:
        labels = [('Run', 'LS', 'dfrac','dcount')]
    print tablePrinter.indent (labels+result, hasHeader = True, separateRows = False,
                               prefix = '| ', postfix = ' |', justify = 'left',
                               delim = ' | ', wrapfunc = lambda x: wrap_onspace (x,70) )

    
def toCSVLSTrg(trgdata,ofilename,iresults=None,isverbose=False):
    '''
    input:[[run,cmslsnum,timestamp,deadfrac,name:(prescale,count)]]
    '''
    print trgdata
def toScreenConfTrg(trgconfdata,iresults=None,isverbose=False):
    '''
    input:{run:[datasource,normbitname,[allbits]]}
    '''
    labels=[('Run','source','bit names')]
    result=[]
    for  run in sorted(trgconfdata):
        if trgconfdata[run] is None:
            ll=[str(run),'n/a','n/a','n/a']
            continue
        source=trgconfdata[run][0]
        source=source.split('/')[-1]
        #normbit=trgconfdata[run][1]
        allbits=trgconfdata[run][2]
        bitnames=', '.join(allbits)
        result.append([str(run),source,bitnames])

    print ' ==  = '
    print tablePrinter.indent (labels+result, hasHeader = True, separateRows = False,
                               prefix = '| ', postfix = ' |', justify = 'left',
                               delim = ' | ', wrapfunc = lambda x: wrap_onspace_strict(x,60) )

def toCSVConfTrg(trgconfdata,ofilename,iresults=None,isverbose=False):
    '''
    input:[[run,l1key,bitnames]]
    '''
    print trgconfdata

def toScreenLSHlt(hltdata,iresults=None,isverbose=False):
    '''
    input:{runnumber:[(cmslsnum,[(hltpath,hltprescale,l1pass,hltaccept),...]),(cmslsnum,[])})}
    '''
    result=[]
    for run in hltdata.keys():
        if hltdata[run] is None:            
            ll=[str(run),'n/a','n/a']
            continue
        perrundata=hltdata[run]
        for lsdata in perrundata:
            cmslsnum=lsdata[0]
            allpathinfo=lsdata[1]
            allpathresult=[]
            for thispathinfo in allpathinfo:
                thispathname=thispathinfo[0]
                thispathpresc=thispathinfo[1]
                thisl1pass=None
                thishltaccept=None
                thispathresult=[]
                thispathresult.append(thispathname)
                thispathresult.append('%d'%thispathpresc)
                if isverbose:
                    if thispathinfo[2] :
                        thisl1pass=thispathinfo[2]
                        thispathresult.append('%d'%thisl1pass)
                    else:
                        thispathresult.append('n/a')
                    if thispathinfo[3]:
                        thishltaccept=thispathinfo[3]
                        thispathresult.append('%d'%thishltaccept)
                    else:
                        thispathresult.append('n/a')
                thispathresultStr='('+','.join(thispathresult)+')'
                allpathresult.append(thispathresultStr)
            result.append([str(run),str(cmslsnum),', '.join(allpathresult)])
    print ' ==  = '
    if isverbose:
        labels = [('Run', 'LS', '(hltpath,presc)')]
    else:
        labels = [('Run', 'LS', '(hltpath,presc,l1pass,hltaccept')]
    print tablePrinter.indent (labels+result, hasHeader = True, separateRows = False,
                               prefix = '| ', postfix = ' |', justify = 'left',
                               delim = ' | ', wrapfunc = lambda x: wrap_onspace (x,70) )
    
def toCSVLSHlt(hltdata,ofilename,iresults=None,isverbose=False):
    '''
    input:[[run,cmslsnum,timestamp,deadfrac,name:(prescale,count)]]
    '''
    print hltdata
def toScreenConfHlt(hltconfdata,iresults=None,isverbose=False):
    '''
    input:[[run,hltkey,pathnames]]
    '''
    print hltconfdata

def toCSVConfHlt(hltconfdata,ofilename,iresults=None,isverbose=False):
    '''
    input:[[run,hltkey,pathnames]]
    '''
    print hltconfdata

