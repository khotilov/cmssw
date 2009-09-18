#include <stdexcept>
#include <string>
#include "OnlineDB/Oracle/interface/Oracle.h"

#include "OnlineDB/EcalCondDB/interface/FEConfigFgrEETowerDat.h"
#include "OnlineDB/EcalCondDB/interface/FEConfigFgrInfo.h"

using namespace std;
using namespace oracle::occi;

FEConfigFgrEETowerDat::FEConfigFgrEETowerDat()
{
  m_env = NULL;
  m_conn = NULL;
  m_writeStmt = NULL;
  m_readStmt = NULL;

  m_thresh = 0;
  m_lut_fg = 0;

}



FEConfigFgrEETowerDat::~FEConfigFgrEETowerDat()
{
}



void FEConfigFgrEETowerDat::prepareWrite()
  throw(runtime_error)
{
  this->checkConnection();

  try {
    m_writeStmt = m_conn->createStatement();
    m_writeStmt->setSQL("INSERT INTO "+getTable()+" (fgr_conf_id, logic_id, "
		      "threshold, lut_fg ) "
		      "VALUES (:fgr_conf_id, :logic_id, "
		      ":threshold, :lut_fg )" );
  } catch (SQLException &e) {
    throw(runtime_error("FEConfigFgrEETowerDat::prepareWrite():  "+e.getMessage()));
  }
}



void FEConfigFgrEETowerDat::writeDB(const EcalLogicID* ecid, const FEConfigFgrEETowerDat* item, FEConfigFgrInfo* iconf)
  throw(runtime_error)
{
  this->checkConnection();
  this->checkPrepare();

  int iconfID = iconf->fetchID();
  if (!iconfID) { throw(runtime_error("FEConfigFgrEETowerDat::writeDB:  ICONF not in DB")); }

  int logicID = ecid->getLogicID();
  if (!logicID) { throw(runtime_error("FEConfigFgrEETowerDat::writeDB:  Bad EcalLogicID")); }
 
  try {
    m_writeStmt->setInt(1, iconfID);
    m_writeStmt->setInt(2, logicID);
    m_writeStmt->setInt(3, item->getThreshold());
    m_writeStmt->setInt(4, item->getLutFg());

    m_writeStmt->executeUpdate();
  } catch (SQLException &e) {
    throw(runtime_error("FEConfigFgrEETowerDat::writeDB():  "+e.getMessage()));
  }
}



void FEConfigFgrEETowerDat::fetchData(map< EcalLogicID, FEConfigFgrEETowerDat >* fillMap, FEConfigFgrInfo* iconf)
  throw(runtime_error)
{
  this->checkConnection();
  fillMap->clear();

  iconf->setConnection(m_env, m_conn);
  int iconfID = iconf->fetchID();
  if (!iconfID) { 
    //  throw(runtime_error("FEConfigFgrEETowerDat::writeDB:  ICONF not in DB")); 
    return;
  }
  
  try {

    m_readStmt->setSQL("SELECT cv.name, cv.logic_id, cv.id1, cv.id2, cv.id3, cv.maps_to, "
		 "d.sliding "
		 "FROM channelview cv JOIN "+getTable()+" d "
		 "ON cv.logic_id = d.logic_id AND cv.name = cv.maps_to "
		 "WHERE fgr_conf_id = :fgr_conf_id");
    m_readStmt->setInt(1, iconfID);
    ResultSet* rset = m_readStmt->executeQuery();
    
    std::pair< EcalLogicID, FEConfigFgrEETowerDat > p;
    FEConfigFgrEETowerDat dat;
    while(rset->next()) {
      p.first = EcalLogicID( rset->getString(1),     // name
			     rset->getInt(2),        // logic_id
			     rset->getInt(3),        // id1
			     rset->getInt(4),        // id2
			     rset->getInt(5),        // id3
			     rset->getString(6));    // maps_to

      dat.setThreshold( rset->getInt(7) );  
      dat.setLutFg( rset->getInt(7) );  

      p.second = dat;
      fillMap->insert(p);
    }
  } catch (SQLException &e) {
    throw(runtime_error("FEConfigFgrEETowerDat::fetchData:  "+e.getMessage()));
  }
}

void FEConfigFgrEETowerDat::writeArrayDB(const std::map< EcalLogicID, FEConfigFgrEETowerDat >* data, FEConfigFgrInfo* iconf)
  throw(runtime_error)
{
  this->checkConnection();
  this->checkPrepare();

  int iconfID = iconf->fetchID();
  if (!iconfID) { throw(runtime_error("FEConfigFgrEETowerDat::writeArrayDB:  ICONF not in DB")); }


  int nrows=data->size(); 
  int* ids= new int[nrows];
  int* iconfid_vec= new int[nrows];
  int* xx= new int[nrows];
  int* yy= new int[nrows];


  ub2* ids_len= new ub2[nrows];
  ub2* iconf_len= new ub2[nrows];
  ub2* x_len= new ub2[nrows];
  ub2* y_len= new ub2[nrows];


  const EcalLogicID* channel;
  const FEConfigFgrEETowerDat* dataitem;
  int count=0;
  typedef map< EcalLogicID, FEConfigFgrEETowerDat >::const_iterator CI;
  for (CI p = data->begin(); p != data->end(); ++p) {
        channel = &(p->first);
	int logicID = channel->getLogicID();
	if (!logicID) { throw(runtime_error("FEConfigFgrEETowerDat::writeArrayDB:  Bad EcalLogicID")); }
	ids[count]=logicID;
	iconfid_vec[count]=iconfID;

	dataitem = &(p->second);
	// dataIface.writeDB( channel, dataitem, iconf);
	int x=dataitem->getThreshold();
	int y=dataitem->getLutFg();

	xx[count]=x;
	yy[count]=y;


	ids_len[count]=sizeof(ids[count]);
	iconf_len[count]=sizeof(iconfid_vec[count]);
	
	x_len[count]=sizeof(xx[count]);
	y_len[count]=sizeof(yy[count]);

	count++;
     }


  try {
    m_writeStmt->setDataBuffer(1, (dvoid*)iconfid_vec, OCCIINT, sizeof(iconfid_vec[0]),iconf_len);
    m_writeStmt->setDataBuffer(2, (dvoid*)ids, OCCIINT, sizeof(ids[0]), ids_len );
    m_writeStmt->setDataBuffer(3, (dvoid*)xx, OCCIINT , sizeof(xx[0]), x_len );
    m_writeStmt->setDataBuffer(4, (dvoid*)yy, OCCIINT , sizeof(yy[0]), y_len );

    m_writeStmt->executeArrayUpdate(nrows);

    delete [] ids;
    delete [] iconfid_vec;
    delete [] xx;
    delete [] yy;

    delete [] ids_len;
    delete [] iconf_len;
    delete [] x_len;
    delete [] y_len;

  } catch (SQLException &e) {
    throw(runtime_error("FEConfigFgrEETowerDat::writeArrayDB():  "+e.getMessage()));
  }
}
