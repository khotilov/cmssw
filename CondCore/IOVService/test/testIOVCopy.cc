#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/PluginManager/interface/standard.h"
#include "FWCore/PluginManager/interface/SharedLibrary.h"

#include "CondCore/DBCommon/interface/DbConnection.h"
#include "CondCore/DBCommon/interface/DbTransaction.h"
#include "CondCore/DBCommon/interface/Exception.h"
#include "CondCore/IOVService/interface/IOVService.h"
#include "CondCore/IOVService/interface/IOVEditor.h"
#include "testPayloadObj.h"
#include <iostream>

int main(){
  edmplugin::PluginManager::Config config;
  edmplugin::PluginManager::configure(edmplugin::standard::config());

  std::string sourceConnect("sqlite_file:source.db");
  std::string destConnect("sqlite_file:dest.db");
  try{
    cond::DbConnection connection;
    connection.configuration().setMessageLevel(coral::Debug);
    connection.configuration().setPoolAutomaticCleanUp( false );
    connection.configure();
    //session->configuration().setAuthenticationMethod(cond::XML);
    cond::DbSession sourcedb = connection.createSession();
    sourcedb.open("sqlite_file:source.db");
    cond::DbSession destdb = connection.createSession();
    destdb.open("sqlite_file:dest.db");
    
    cond::IOVService iovmanager(sourcedb);
    cond::IOVEditor* editor=iovmanager.newIOVEditor();
    sourcedb.transaction().start(false);
    editor->create(cond::timestamp,1);
    for(int i=0; i<5; ++i){
      std::cout<<"creating test payload obj"<<i<<std::endl;
      testPayloadObj* myobj=new testPayloadObj;
      for(int j=0; j<10; ++j){
        myobj->data.push_back(i+j);
      }
      
      boost::shared_ptr<testPayloadObj> myobjPtr (myobj );
      std::string tok = sourcedb.storeObject(myobjPtr.get(),"testPayloadObjRcd");
      editor->append(i+10, tok);
    }
    std::string iovtoken=editor->token();
    std::cout<<"iov token "<<iovtoken<<std::endl;
    sourcedb.transaction().commit();
    std::cout<<"source db created "<<std::endl;
    sourcedb.transaction().start(true);
    std::cout<<"source db started "<<std::endl;
    destdb.transaction().start(false);
    std::cout<<"dest db started "<<std::endl;
    iovmanager.exportIOVWithPayload( destdb,
				     iovtoken);
    destdb.transaction().commit();
    std::cout<<"destdb committed"<<std::endl;
    sourcedb.transaction().commit();
    std::cout<<"source db committed"<<std::endl;
    delete editor;
    std::cout<<"editor deleted"<<std::endl;
  }catch(const cond::Exception& er){
    std::cout<<"error "<<er.what()<<std::endl;
  }catch(const std::exception& er){
    std::cout<<"std error "<<er.what()<<std::endl;
  }
}
