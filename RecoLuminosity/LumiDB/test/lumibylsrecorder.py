import os,os.path,csv,coral,commands
from RecoLuminosity.LumiDB import dbUtil,lumiTime,sessionManager,nameDealer
def getrunsInResult(schema,minrun=132440,maxrun=500000):
    '''
    get runs in result tables in specified range
    output:
         [runnum]
         select distinct runnum from hflumiresult where runnum>=:minrun and runnum<=:maxrun;
    '''
    result=[]
    qHandle=schema.newQuery()
    try:
        qHandle.addToTableList( 'HFLUMIRESULT' )
        qHandle.addToOutputList('distinct RUNNUM')
        qCondition=coral.AttributeList()
        qCondition.extend('minrun','unsigned int')
        qCondition.extend('maxrun','unsigned int')
        qCondition['minrun'].setData(minrun)
        qCondition['maxrun'].setData(maxrun)
        qResult=coral.AttributeList()
        qResult.extend('RUNNUM','unsigned int')
        qHandle.defineOutput(qResult)
        qHandle.setCondition('RUNNUM>=:minrun AND RUNNUM<=:maxrun',qCondition)
        cursor=qHandle.execute()
        while cursor.next():
            runnum=cursor.currentRow()['RUNNUM'].data()
            result.append(runnum)
        del qHandle
    except:
        if qHandle:del qHandle
        raise
    return result
    
def getrunsInCurrentData(schema,minrun=132440,maxrun=500000):
    '''
    get runs in data tables in specified range
    output:
         [runnum]
         select runnum,tagid from tagruns where runnum>=:minrun and runnum<=:maxrun;
    '''
    tmpresult={}
    qHandle=schema.newQuery()
    try:
        qHandle.addToTableList( nameDealer.tagRunsTableName() )
        qHandle.addToOutputList('RUNNUM')
        qHandle.addToOutputList('TAGID')
        qCondition=coral.AttributeList()
        qCondition.extend('minrun','unsigned int')
        qCondition.extend('maxrun','unsigned int')
        qCondition['minrun'].setData(minrun)
        qCondition['maxrun'].setData(maxrun)
        qResult=coral.AttributeList()
        qResult.extend('RUNNUM','unsigned int')
        qResult.extend('TAGID','unsigned long long')
        qHandle.defineOutput(qResult)
        qHandle.setCondition('RUNNUM>=:minrun AND RUNNUM<=:maxrun',qCondition)
        cursor=qHandle.execute()
        while cursor.next():
            runnum=cursor.currentRow()['RUNNUM'].data()
            tagid=cursor.currentRow()['TAGID'].data()
            if not tmpresult.has_key(runnum):
                tmpresult[runnum]=tagid
            else:
                if tagid>tmpresult[runnum]:
                    tmpresult[runnum]=tagid
        del qHandle
    except:
        if qHandle:del qHandle
        raise
    if tmpresult:return tmpresult.keys()
    return []
def execCalc(connectStr,authpath,runnum):
    '''
    run lumiCalc2.py lumibyls for the run
    '''
    outdatafile=str(runnum)+'.csv'
    outheaderfile=str(runnum)+'.txt'
    command = 'lumiCalc2.py lumibyls -c ' +connectStr+' -P '+authpath+' -r '+str(runnum)+' -o '+outdatafile+' --headerfile '+outheaderfile
    statusAndOutput = commands.getstatusoutput(command)
    print statusAndOutput
    
class lslumiParser(object):
    def __init__(self,lslumifilename,headerfilename):
        '''
        '''
        self.__filename=lslumifilename
        self.__headername=headerfilename
        self.lumidata=[]#[fill,run,lumils,cmsls,beamstatus,beamenergy,delivered,recorded,avgpu]
        self.datatag=''
        self.normtag=''
    def parse(self):
        '''
        parse ls lumi file
        '''
        hf=open(self.__headername,'rb')
        for line in hf:
            if "lumitype" in line:
                fields=line.strip().split(',')
                for field in fields:
                    a=field.strip().split(':')
                    if a[0]=='datatag':
                        self.datatag=a[1].strip()
                    if a[0]=='normtag':
                        self.normtag=a[1].strip()
                break
        hf.close()
        f=open(self.__filename,'rb')
        freader=csv.reader(f,delimiter=',')
        idx=0
        for row in freader:
           if idx==0:
               idx=1 # skip header
               continue
           [run,fill]=map(lambda i:int(i),row[0].split(':'))
           [lumils,cmsls]=map(lambda i:int(i),row[1].split(':'))
           chartime=row[2]
           beamstatus=row[3]
           beamenergy=float(row[4])
           beamenergy=int(round(beamenergy))
           delivered=float(row[5])
           recorded=float(row[6])
           avgpu=float(row[7])
           self.lumidata.append([fill,run,lumils,cmsls,chartime,beamstatus,beamenergy,delivered,recorded,avgpu])
        f.close()

def addindb(session,datatag,normtag,lumidata,bulksize):
    '''
    input : [fill,run,lumils,cmsls,lstime,beamstauts,beamenergy,delivered,recorded,avgpu]
    '''
    hfresultDefDict=[('RUNNUM','unsigned int'),('LS','unsigned int'),('CMSLS','unsigned int'),('FILLNUM','unsigned int'),('TIME','time stamp'),('BEAM_STATUS','string'),('ENERGY','unsigned int'),('DELIVERED','float'),('RECORDED','float'),('AVG_PU','float'),('DATA_VERSION','string'),('NORM_VERSION','string'),('INSERT_TIME','time stamp')]
    
    committedrows=0
    nrows=0
    bulkvalues=[]
    lute=lumiTime.lumiTime()
    try:
        for datum in lumidata:
            [fillnum,runnum,lumils,cmsls,lstime_char,beamstatus,beamenergy,delivered,recorded,avgpu]=datum
            inserttime=coral.TimeStamp()
            lstime=lute.StrToDatetime(lstime_char,customfm='%m/%d/%y %H:%M:%S')
            corallstime=coral.TimeStamp(lstime.year,lstime.month,lstime.day,lstime.hour,lstime.minute,lstime.second,0)
            bulkvalues.append([('RUNNUM',runnum),('LS',lumils),('CMSLS',cmsls),('FILLNUM',fillnum),('TIME',corallstime),('BEAM_STATUS',beamstatus),('ENERGY',beamenergy),('DELIVERED',delivered),('RECORDED',recorded),('AVG_PU',avgpu),('DATA_VERSION',datatag),('NORM_VERSION',normtag),('INSERT_TIME',inserttime)])
            nrows+=1
            committedrows+=1
            if nrows==bulksize:
                print 'committing trg in LS chunck ',nrows
                db=dbUtil.dbUtil(session.nominalSchema())
                session.transaction().start(False)
                db.bulkInsert('HFLUMIRESULT',hfresultDefDict,bulkvalues)
                session.transaction().commit()
                nrows=0
                bulkvalues=[]
            elif committedrows==len(lumidata):
                print 'committing at the end '
                db=dbUtil.dbUtil(session.nominalSchema())
                session.transaction().start(False)
                db.bulkInsert('HFLUMIRESULT',hfresultDefDict,bulkvalues)
                session.transaction().commit()
                
    except :
        print 'error in addindb'
        raise 
if __name__ == "__main__" :
    sourcestr='oracle://cms_orcon_adg/cms_lumi_prod'
    pth='/afs/cern.ch/cms/lumi/DB'
    sourcesvc=sessionManager.sessionManager(sourcestr,authpath=pth,debugON=False)
    sourcesession=sourcesvc.openSession(isReadOnly=True,cpp2sqltype=[('unsigned int','NUMBER(10)'),('unsigned long long','NUMBER(20)')])
    sourcesession.transaction().start(True)
    sourcerunlist=getrunsInCurrentData(sourcesession.nominalSchema(),minrun=183339,maxrun=209310)
    sourcesession.transaction().commit()
    print sourcerunlist
    deststr='oracle://cms_orcoff_prep/cms_lumi_dev_offline'
    destsvc=sessionManager.sessionManager(deststr,authpath=pth,debugON=False)
    destsession=destsvc.openSession(isReadOnly=True,cpp2sqltype=[('unsigned int','NUMBER(10)'),('unsigned long long','NUMBER(20)')])
    destsession.transaction().start(True)
    destrunlist=getrunsInResult(destsession.nominalSchema(),minrun=183339,maxrun=209310)
    destsession.transaction().commit()
    print destrunlist
    for r in sourcerunlist:
        if r not in destrunlist:
            execCalc(sourcestr,pth,r)
    #p=lslumiParser(lslumifilename,lumiheaderfilename)
    #p.parse()
    #svc=sessionManager.sessionManager(sourcestr,authpath=pth,debugON=False)
    #dbsession=svc.openSession(isReadOnly=False,cpp2sqltype=[('unsigned int','NUMBER(10)'),('unsigned long long','NUMBER(20)')])
    #addindb(dbsession,p.datatag,p.normtag,p.lumidata,bulksize=500)
    #del dbsession
    #del svc
