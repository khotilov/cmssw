#include <stdexcept>
#include <string>
#include "OnlineDB/Oracle/interface/Oracle.h"

#include "OnlineDB/EcalCondDB/interface/MonH4TablePositionDat.h"

using namespace std;
using namespace oracle::occi;

MonH4TablePositionDat::MonH4TablePositionDat()
{
  m_env = NULL;
  m_conn = NULL;
  m_writeStmt = NULL;

  m_tableX = 0;
  m_tableY = 0;
}



MonH4TablePositionDat::~MonH4TablePositionDat()
{
}



void MonH4TablePositionDat::prepareWrite()
  throw(runtime_error)
{
  this->checkConnection();

  try {
    m_writeStmt = m_conn->createStatement();
    m_writeStmt->setSQL("INSERT INTO mon_h4_table_position_dat (iov_id, logic_id, "
			"table_x, table_y) "
			"VALUES (:iov_id, :logic_id, "
			":3, :4)");
  } catch (SQLException &e) {
    throw(runtime_error("MonH4TablePositionDat::prepareWrite():  "+e.getMessage()));
  }
}



void MonH4TablePositionDat::writeDB(const EcalLogicID* ecid, const MonH4TablePositionDat* item, MonRunIOV* iov)
  throw(runtime_error)
{
  this->checkConnection();
  this->checkPrepare();

  int iovID = iov->fetchID();
  if (!iovID) { throw(runtime_error("MonH4TablePositionDat::writeDB:  IOV not in DB")); }

  int logicID = ecid->getLogicID();
  if (!logicID) { throw(runtime_error("MonH4TablePositionDat::writeDB:  Bad EcalLogicID")); }
  
  try {
    m_writeStmt->setInt(1, iovID);
    m_writeStmt->setInt(2, logicID);

    m_writeStmt->setFloat(3, item->getTableX() );
    m_writeStmt->setFloat(4, item->getTableY() );

    m_writeStmt->executeUpdate();
  } catch (SQLException &e) {
    throw(runtime_error("MonH4TablePositionDat::writeDB():  "+e.getMessage()));
  }
}



void MonH4TablePositionDat::fetchData(std::map< EcalLogicID, MonH4TablePositionDat >* fillMap, MonRunIOV* iov)
  throw(runtime_error)
{
  this->checkConnection();
  fillMap->clear();

  iov->setConnection(m_env, m_conn);
  int iovID = iov->fetchID();
  if (!iovID) { 
    //  throw(runtime_error("MonH4TablePositionDat::writeDB:  IOV not in DB")); 
    return;
  }

  try {
    Statement* stmt = m_conn->createStatement();
    stmt->setSQL("SELECT cv.name, cv.logic_id, cv.id1, cv.id2, cv.id3, cv.maps_to, "
		 "d.table_x, d.table_y "
		 "FROM channelview cv JOIN mon_h4_table_position_dat d "
		 "ON cv.logic_id = d.logic_id AND cv.name = cv.maps_to "
		 "WHERE d.iov_id = :iov_id");
    stmt->setInt(1, iovID);
    ResultSet* rset = stmt->executeQuery();
    
    std::pair< EcalLogicID, MonH4TablePositionDat > p;
    MonH4TablePositionDat dat;
    while(rset->next()) {
      p.first = EcalLogicID( rset->getString(1),     // name
			     rset->getInt(2),        // logic_id
			     rset->getInt(3),        // id1
			     rset->getInt(4),        // id2
			     rset->getInt(5),        // id3
			     rset->getString(6));    // maps_to

      dat.setTableX( rset->getFloat(7) );
      dat.setTableY( rset->getFloat(8) );

      p.second = dat;
      fillMap->insert(p);
    }
  } catch (SQLException &e) {
    throw(runtime_error("MonH4TablePositionDat::fetchData():  "+e.getMessage()));
  }
}
