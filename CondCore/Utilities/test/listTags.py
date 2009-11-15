import DLFCN, sys, os, time
sys.setdlopenflags(DLFCN.RTLD_GLOBAL+DLFCN.RTLD_LAZY)

from pluginCondDBPyInterface import *
a = FWIncantation()
os.putenv("CORAL_AUTH_PATH","/afs/cern.ch/cms/DB/conddb")

import coral
from CondCore.TagCollection import Node,tagInventory,TagTree


svc = coral.ConnectionService()
# session = svc.connect("sqlite_file:CondCore/TagCollection/data/GlobalTag.db",accessMode = coral.access_ReadOnly )
session = svc.connect("oracle://cms_orcoff_prod/CMS_COND_31X_GLOBALTAG",accessMode = coral.access_ReadOnly )
inv=tagInventory.tagInventory(session)
mytree=TagTree.tagTree(session,"GR09_31X_V5P")
result=mytree.getAllLeaves()
tags=[]
for r in result:
    if r.tagid != 0:
        tagcontent=inv.getEntryById(r.tagid)
        tags.append((tagcontent.recordname, tagcontent.labelname, tagcontent.tagname,tagcontent.pfn))


mytree=0
inv=0
del session
del svc

// rdbms = RDBMS()
rdbms = RDBMS("/afs/cern.ch/cms/DB/conddb")
logName = "oracle://cms_orcoff_prod/CMS_COND_31X_POPCONLOG"
gdbName = "oracle://cms_orcoff_prod/CMS_COND_31X_GLOBALTAG"
gName = "GR09_31X_V5P"
rdbms.setLogger(logName)
globalTag = rdbms.globalTag(gdbName,gName)

for tag in tags:
#    dbname = tag[3][tag[3].rfind('/'):]
#    db_o = rdbms.getDB("oracle://cms_orcoff_prod"+dbname)
    db = rdbms.getDB(tag[3])
    log = db.lastLogEntry(tag[2])
    iov = db.iov(tag[2])
    print tag[0],tag[1],tag[3],tag[2]
    print iov.size(), iov.revision(), time.asctime(time.gmtime(unpackTime(iov.timestamp())[0])), iov.comment()
    print log.getState()
    iov=0
    db=0

def iovSize(rdbms,conn,tag) :
    try :
        db = rdbms.getDB(conn)
        log = db.lastLogEntry(tag)
        iov = db.iov(tag)
        size = iov.size()
        for elem in iov.elements :
            if (elem.till()>4294967295) : 
                print tag, elem.since(), elem.till(), elem.payloadToken()
    except RuntimeError :
        print conn, tag," no iov?"
        size=-1
    iov=0
    db=0
    return size




