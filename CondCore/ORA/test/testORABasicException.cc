#include "CondCore/ORA/interface/Database.h"
#include "CondCore/ORA/interface/Container.h"
#include "CondCore/ORA/interface/Transaction.h"
#include "CondCore/ORA/interface/Exception.h"
#include <iostream>
#include <stdexcept>
namespace {

  void testDbWRFunc( ora::Database& db ){
    ora::Object obj;
    std::string data("BlaBla");
    std::string contName("MyCont");
    try {
      db.create();
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.drop();
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.createContainer<std::string>( contName );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.createContainer<std::string>();
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.dropContainer( contName );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.insertItem( contName, obj );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.insert( contName, data );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.insert( data );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
  }
  
  void testDbUPFunc( ora::Database& db ){
    ora::OId oid;
    ora::Object obj;
    std::string data("BlaBla");
    try {
      db.updateItem( oid, obj );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.update( oid, data );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.erase( oid );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
  }
  
  void testDbRDFunction( ora::Database& db ){
    try {
      db.exists();
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.containers();
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    std::string contName( "NewCont");
    try {
      db.containerHandle( contName );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.containerHandle( 1 );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    ora::OId oid;
    try {
      db.fetchItem( oid );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.fetch<std::string>( oid );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
  }
  
  void testAllFunction( ora::Database& db ){
    testDbWRFunc( db );
    testDbUPFunc( db );
    testDbRDFunction( db );
    try {
      db.transaction();
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      db.flush();
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
  }
}

int main(){

    // writing...  
  ora::Database db;
  try {
    //std::string connStr( "sqlite_file:test.db" );
    std::string connStr( "oracle://devdb10/giacomo" );
    // NOT CONNECTED
    bool connected = db.isConnected();
    std::string okConn("");
    if( !connected ) okConn = "NOT ";
    std::cout << "# Database is "<<okConn<<" connected.Connstr=["<<db.connectionString()<<"]"<<std::endl;
    testAllFunction( db );
    try {
      db.disconnect();
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    try {
      std::string invalidConnection("invalid");
      db.connect( invalidConnection );
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
    db.connect( connStr );
    connected = db.isConnected();
    okConn = std::string("");
    if( !connected ) okConn = "NOT ";
    std::cout << "# Database is now "<<okConn<<" connected.Connstr=["<<db.connectionString()<<"]"<<std::endl;
    testAllFunction( db );
    std::cout << "# Starting read only transaction..."<<std::endl;
    db.transaction().start();
    //testDbWRFunc( db );
    //testDbUPFunc( db );
    try {
      std::cout << "# disconnecting..."<<std::endl;
      db.disconnect();
    } catch ( const ora::Exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    } catch ( const std::exception& e ){
      std::cout << "# Error: "<<e.what()<<std::endl;
    }
  }catch ( const ora::Exception& e ){
    std::cout << "# Test Error: "<<e.what()<<std::endl;
  }catch ( const std::exception& e){
    std::cout << "# Test Error: "<<e.what()<<std::endl;
    }
  
}

