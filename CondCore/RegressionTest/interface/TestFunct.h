#include <iostream>
#include <unistd.h>




#include "CondCore/ORA/interface/SchemaUtils.h" //v4
#include <stdexcept>
#include "CondCore/ORA/interface/Database.h"
#include "CondCore/ORA/interface/Container.h"
#include "CondCore/ORA/interface/OId.h"
#include "CondCore/ORA/interface/ScopedTransaction.h"
#include "CondCore/ORA/interface/Transaction.h"
#include "CondCore/ORA/interface/Exception.h"
#include "CondCore/ORA/interface/IBlobStreamingService.h"
#include "Reflex/Member.h"
#include "Reflex/Object.h"
#include "CoralBase/Blob.h"

#include "CondCore/DBCommon/interface/DbConnection.h"
#include "CondCore/DBCommon/interface/DbScopedTransaction.h"
#include "CondCore/DBCommon/interface/DbTransaction.h"
#include "CondCore/DBCommon/interface/Exception.h"
#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/PluginManager/interface/standard.h"
#include "CondCore/DBTest/interface/TestPayloadClass.h"
#include "CondCore/MetaDataService/interface/MetaDataSchemaUtility.h"
#include "CondCore/MetaDataService/interface/MetaData.h"
#include "RelationalAccess/SchemaException.h"
#include "RelationalAccess/ISchema.h"
#include "RelationalAccess/ITable.h"
#include "RelationalAccess/TableDescription.h"
#include "RelationalAccess/ITablePrivilegeManager.h"
#include "RelationalAccess/ICursor.h"
#include "RelationalAccess/IQuery.h"
#include "RelationalAccess/ITableDataEditor.h"
#include "CoralBase/AttributeList.h"
#include "CoralBase/AttributeSpecification.h"
#include "CoralBase/Attribute.h"





class TestFunct {
public :
	cond::DbSession s; 
	TestFunct();
	bool Write(std::string mappingName, int payloadID);
	bool Read(std::string mappingName);
	bool ReadAll();
	bool CreateMetaTable();
	bool DropTables(std::string connStr);
	bool DropItem(std::string mappingName);
};