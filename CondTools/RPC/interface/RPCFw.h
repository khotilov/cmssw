#ifndef RPC_DB_FW_H
#define RPC_DB_FW_H

/*
 * \class RPCFw
 *  Reads data from OMDS and creates conditioning objects
 *
 *  $Date: 2008/07/31 14:49:53 $
 *  $Revision: 1.1 $
 *  \author D. Pagano - Dip. Fis. Nucl. e Teo. & INFN Pavia
 */



#include "CondTools/RPC/interface/RPCDBCom.h"
#include "CoralBase/TimeStamp.h"
#include "CondTools/RPC/interface/RPCSourceHandler.h"
#include "CondTools/RPC/interface/RPCGasSH.h"


struct dbread{
    float alias;
    float value;
};


class RPCFw : virtual public RPCDBCom
{
public:
  RPCFw( const std::string& connectionString,
         const std::string& userName,
         const std::string& password);
  virtual ~RPCFw();
  void run();

  coral::TimeStamp UTtoT(int utime);


  coral::TimeStamp thr;
  std::vector<RPCdbData::Item> createIMON(int from);
  std::vector<RPCdbData::Item> createVMON(int from); 
  std::vector<RPCdbData::Item> createSTATUS(int from); 
  std::vector<RPCGas::GasItem> createGAS(int from);
  std::vector<RPCGas::TempItem> createT(int from);

  
private:
  std::string m_connectionString;
  std::string m_userName;
  std::string m_password;
};

#endif
