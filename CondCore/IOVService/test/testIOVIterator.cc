#include "CondCore/DBCommon/interface/DBSession.h"
#include "CondCore/DBCommon/interface/Connection.h"
#include "CondCore/DBCommon/interface/Exception.h"
#include "CondCore/DBCommon/interface/SessionConfiguration.h"
#include "CondCore/DBCommon/interface/PoolTransaction.h"
#include "CondCore/DBCommon/interface/MessageLevel.h"
#include "CondCore/IOVService/interface/IOVService.h"
#include "CondCore/IOVService/interface/IOVEditor.h"
#include "CondCore/IOVService/interface/IOVIterator.h"
#include <iostream>
int main(){
  try{
    cond::DBSession* session=new cond::DBSession;
    session->open();
    cond::Connection myconnection("sqlite_file:mytest.db",0); 
    myconnection.connect(session);
    cond::PoolTransaction& pooldb=myconnection.poolTransaction();
    pooldb.start(false);
    cond::IOVService iovmanager(pooldb);
    cond::IOVEditor* editor=iovmanager.newIOVEditor();
    editor->insert(20,"pay1tok");
    editor->insert(40,"pay2tok");
    editor->insert(60,"pay3tok");
    pooldb.commit();
    std::string iovtok=editor->token();
    ///test iterator
    cond::IOVIterator* it=iovmanager.newIOVIterator(iovtok);
    std::cout<<"test iterator "<<std::endl;
    pooldb.start(true);
    while( it->next() ){
      std::cout<<"payloadToken "<<it->payloadToken()<<std::endl;
      std::cout<<"since "<<it->validity().first<<std::endl;
      std::cout<<"till "<<it->validity().second<<std::endl;
    }
    std::cout<<"is 30 valid? "<<iovmanager.isValid(iovtok,30)<<std::endl;
    pooldb.commit();
    myconnection.disconnect();
    delete editor;
    delete it;
    delete session;
  }catch(const cond::Exception& er){
    std::cout<<"error "<<er.what()<<std::endl;
  }catch(const std::exception& er){
    std::cout<<"std error "<<er.what()<<std::endl;
  }
}
