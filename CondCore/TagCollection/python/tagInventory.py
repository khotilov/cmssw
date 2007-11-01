import coral
import IdGenerator, Node, DBImpl
class  tagInventory(object):
    """Class manages tag inventory 
    """
    def __init__( self , session ):
        """Input: coral session handle
        """
        self.__session = session
        self.__tagInventoryTableName = 'TAGINVENTORY_TABLE'
        self.__tagInventoryTableColumns = {'tagid':'unsigned int', 'tagname':'string', 'pfn':'string','recordname':'string', 'objectname':'string', 'labelname':'string','timetype':'string','comment':'string'}
        self.__tagInventoryTableNotNullColumns = ['tagname','pfn','recordname','objectname','labelname']
        self.__tagInventoryTableUniqueColumns = ['tagname']
        self.__tagInventoryTablePK = ('tagid')
    def existInventoryTable( self ):
        """Check if inventory table exists
        """
        try:
            transaction=self.__session.transaction()
            transaction.start(True)
            schema = self.__session.nominalSchema()
            result=schema.existsTable(self.__tagInventoryTableName)
            transaction.commit()
            #print result
        except Exception, er:
            transaction.rollback()
            raise Exception, str(er)
        return result
    def createInventoryTable( self ):
        """Create tag inventory table. Existing table will be deleted. 
        """
        try:
            transaction=self.__session.transaction()
            transaction.start()
            schema = self.__session.nominalSchema()
            schema.dropIfExistsTable( self.__tagInventoryTableName )
            description = coral.TableDescription();
            description.setName( self.__tagInventoryTableName )
            for columnName, columnType in self.__tagInventoryTableColumns.items():
                description.insertColumn(columnName, columnType)
            for columnName in self.__tagInventoryTableNotNullColumns :
                description.setNotNullConstraint(columnName,True)
            for columnName in self.__tagInventoryTableUniqueColumns :
                description.setUniqueConstraint(columnName)
            combinedunique=('recordname','objectname','labelname')
            description.setUniqueConstraint(combinedunique)
            description.setPrimaryKey(  self.__tagInventoryTablePK )
            self.__tagInventoryTableHandle = schema.createTable( description )
            self.__tagInventoryTableHandle.privilegeManager().grantToPublic( coral.privilege_Select )
            #create also the associated id table
            generator=IdGenerator.IdGenerator(schema)
            generator.createIDTable(self.__tagInventoryTableName,True)
            transaction.commit()
        except Exception, er:
            transaction.rollback()
            raise Exception, str(er)
    def addEntry( self, leafNode ):
        """Add entry into the inventory.\n
        Input: leafNode
        Output: tagid
        """
        tagid=0
        try:
            transaction=self.__session.transaction()
            transaction.start(True)
            schema = self.__session.nominalSchema()
            generator=IdGenerator.IdGenerator(schema)
            tagid=generator.getNewID(generator.getIDTableName(self.__tagInventoryTableName))
            transaction.commit()
            
            transaction.start(False)
            dbop=DBImpl.DBImpl(schema)
            tabrowValueDict={'tagid':tagid,'tagname':leafNode.tagname,'objectname':leafNode.objectname,'pfn':leafNode.pfn,'labelname':leafNode.labelname,'recordname':leafNode.recordname,'timetype':leafNode.timetype,'comment':leafNode.comment}
            dbop.insertOneRow(self.__tagInventoryTableName,
                              self.__tagInventoryTableColumns,
                              tabrowValueDict)
            transaction.commit()
            transaction.start(False)
            generator.incrementNextID(generator.getIDTableName(self.__tagInventoryTableName))
            transaction.commit()
            return tagid
        except Exception, er:
            transaction.rollback()
            raise Exception, str(er)
    def getEntryByName( self, tagName ):
        """Get basic tag from inventory by lable.\n
        Input: tagname
        Output: leafNode
        """
        leafnode = Node.LeafNode()
        transaction=self.__session.transaction()
        try:
            transaction.start(True)
            query = self.__session.nominalSchema().tableHandle(self.__tagInventoryTableName).newQuery()
            for columnName in self.__tagInventoryTableColumns:
                query.addToOutputList(columnName)
            condition = "tagname=:tagname"
            conditionData = coral.AttributeList()
            conditionData.extend( 'tagname','string' )
            conditionData['tagname'].setData(tagName)
            query.setCondition(condition,conditionData)
            cursor = query.execute()
            while ( cursor.next() ):
                leafnode.tagid=cursor.currentRow()['tagid'].data()
                leafnode.tagname=cursor.currentRow()['tagname'].data()
                leafnode.objectname=cursor.currentRow()['objectname'].data()
                leafnode.pfn=cursor.currentRow()['pfn'].data()
                leafnode.labelname=cursor.currentRow()['labelname'].data()
                leafnode.recordname=cursor.currentRow()['recordname'].data()
                leafnode.timetype=cursor.currentRow()['timetype'].data()
                leafnode.comment=cursor.currentRow()['comment'].data()
            transaction.commit()
            del query
            return leafnode
        except Exception, e:
            transaction.rollback()
            raise Exception, str(e)
    def getEntryById( self, tagId ):
        """Get basic tag from inventory by id.\n
        Input: tagid
        Output: leafNode
        """
        leafnode = Node.LeafNode()
        transaction=self.__session.transaction()
        try:
            transaction.start(True)
            query = self.__session.nominalSchema().tableHandle(self.__tagInventoryTableName).newQuery()
            for columnName in self.__tagInventoryTableColumns:
                query.addToOutputList(columnName)
            condition = "tagid=:tagid"
            conditionData = coral.AttributeList()
            conditionData.extend( 'tagid','unsigned long' )
            cursor = query.execute()
            while ( cursor.next() ):
                leafnode.tagid=cursor.currentRow()['tagid'].data()
                leafnode.tagname=cursor.currentRow()['tagname'].data()
                leafnode.objectname=cursor.currentRow()['objectname'].data()
                leafnode.pfn=cursor.currentRow()['pfn'].data()
                leafnode.labelname=cursor.currentRow()['labelname'].data()
                leafnode.recordname=cursor.currentRow()['recordname'].data()
                leafnode.timetype=cursor.currentRow()['timetype'].data()
                leafnode.comment=cursor.currentRow()['comment'].data()
            transaction.commit()
            del query
            return leafnode
        except Exception, e:
            transaction.rollback()
            raise Exception, str(e)
    def getAllEntries( self ):
        """Get all entries in the inventory
        Output: list of leafNode objects
        """
        result=[]
        transaction=self.__session.transaction()
        try:
            transaction.start(True)
            query = self.__session.nominalSchema().tableHandle(self.__tagInventoryTableName).newQuery()
            for columnName in self.__tagInventoryTableColumns:
                query.addToOutputList(columnName)    
            cursor = query.execute()
            while ( cursor.next() ):
                leafnode = Node.LeafNode()
                leafnode.tagid=cursor.currentRow()['tagid'].data()
                leafnode.tagname=cursor.currentRow()['tagname'].data()
                leafnode.objectname=cursor.currentRow()['objectname'].data()
                leafnode.pfn=cursor.currentRow()['pfn'].data()
                leafnode.recordname=cursor.currentRow()['recordname'].data()
                leafnode.labelname=cursor.currentRow()['labelname'].data()
                leafnode.timetype=cursor.currentRow()['timetype'].data()
                leafnode.comment=cursor.currentRow()['comment'].data()
                result.append(leafnode)
            transaction.commit()
            del query
            return result
        except Exception, e:
            transaction.rollback()
            raise Exception, str(e)
        
    def deleteAllEntries( self ):
        """Delete all entries in the inventory
        """
        try:
            transaction=self.__session.transaction()
            transaction.start(False)
            schema = self.__session.nominalSchema()
            dbop=DBImpl.DBImpl(schema)
            inputData = coral.AttributeList()
            dbop.deleteRows(self.__tagInventoryTableName,
                            '',
                            inputData)
            transaction.commit()
        except Exception, e:
            transaction.rollback()
            raise Exception, str(e)
        
    def deleteEntryByName( self, tagname ):
        """Delete entry with given tag name
        """
        try:
            transaction=self.__session.transaction()
            transaction.start(False)
            schema = self.__session.nominalSchema()
            dbop=DBImpl.DBImpl(schema)
            inputData = coral.AttributeList()
            inputData.extend( "tagname","string" )
            inputData[0].setData(tagname)
            dbop.deleteRows(self.__tagInventoryTableName,
                            'tagname=:tagname',
                            inputData)
            transaction.commit()
        except Exception, e:
            transaction.rollback()
            raise Exception, str(e)
    def replaceTagLabel( self, tagname, label ):
        """Replace the run time label of the given tag
        """
        try:
            transaction=self.__session.transaction()
            transaction.start(False)
            schema = self.__session.nominalSchema()
            inputData = coral.AttributeList()
            inputData.extend( "labelname","string" )
            inputData.extend( "tagname", "string" ) 
            inputData[0].setData(label)
            inputData[1].setData(tagname)
            editor = schema.tableHandle(self.__tagInventoryTableName).dataEditor()
            editor.updateRows( "labelname=:labelname", "tagname=:tagname", inputData )
            transaction.commit()
        except Exception, e:
            transaction.rollback()
            raise Exception, str(e)
if __name__ == "__main__":
    context = coral.Context()
    context.setVerbosityLevel( 'ERROR' )
    svc = coral.ConnectionService( context )
    session = svc.connect( 'sqlite_file:testInventory.db',
                           accessMode = coral.access_Update )
    try:
        inv=tagInventory(session)
        inv.createInventoryTable()
        tagentry=Node.LeafNode()
        tagentry.tagname='ecalpedestalsfromonline'
        tagentry.objectname='EcalPedestals'
        tagentry.pfn='oracle://devdb10/CMS_COND_ECAL'
        tagentry.recordname='EcalPedestalsRcd'
        tagentry.labelname=''
        tagentry.timetype='runnumber'
        tagentry.comment='crapcrapcrapandcrackcrackcrak'
        inv.addEntry(tagentry)
        tagentry.tagname='crap'
        tagentry.objectname='MyPedestals'
        tagentry.pfn='oracle://devdb10/CMS_COND_ME'
        tagentry.recordname='MyPedestalsRcd'
        tagentry.labelname='mylabel'
        tagentry.timetype='runnumber'
        tagentry.comment='no'
        inv.addEntry(tagentry)
        result=inv.getAllEntries()
        print 'get all##\t',result
        result=inv.getEntryByName('ecalpedestalsfromonline')
        print 'get ecalpedestalsfromonline##\t',result
        result=inv.getEntryByName('crap')
        print 'get crap##\t',result
        result=inv.getEntryById(0)
        print 'get by id 0##\t',result
        inv.deleteEntryByName('ecalpedestalsfromonline')
        result=inv.getEntryByName('ecalpedestalsfromonline')
        print 'get ecalpedestalsfromonline##\t',result
        result=inv.getEntryByName('crap')
        print 'get crap##\t',result
        del session
    except Exception, e:
        print "Failed in unit test"
        print str(e)
        del session
        
