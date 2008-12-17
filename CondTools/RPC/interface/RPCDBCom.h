#ifndef RPCDBCOM_H
#define RPCDBCOM_H

#include <string>
#include "SealKernel/Context.h"
#include "SealKernel/IMessageService.h"

namespace coral {
  class IConnection;
  class ISession;
}

class RPCDBCom
{
public:
  RPCDBCom();
  virtual ~RPCDBCom();
  virtual void run() = 0;
  void setVerbosityLevel( seal::Msg::Level level );

protected:
  coral::ISession* connect( const std::string& connectionString,
                            const std::string& userName,
                            const std::string& password );

private:
  seal::Handle<seal::Context> m_context;
  coral::IConnection*         m_connection;
};

#endif

