#include "CondCore/ORA/interface/Database.h"
#include "CondCore/ORA/interface/Container.h"
#include "CondCore/ORA/interface/Transaction.h"
#include "CondCore/ORA/interface/ScopedTransaction.h"
#include "CondCore/ORA/interface/Exception.h"
#include "CondCore/ORA/test/Serializer.h"
#include <cstdlib>
#include <iostream>

int main(){
  // writing...
  std::string authpath("/afs/cern.ch/cms/DB/conddb");
  std::string pathenv(std::string("CORAL_AUTH_PATH=")+authpath);
  ::putenv(const_cast<char*>(pathenv.c_str()));
  try {

    std::string connStr( "oracle://cms_orcoff_prep/CMS_COND_WEB" );
    //std::string connStr( "sqlite_file:test.db" );
    ora::Serializer serializer( "ORA_TEST" );
    //serializer.lock( connStr );
    ora::Database db0;
    db0.configuration().setMessageVerbosity( coral::Debug );
    db0.connect( connStr );
    ora::ScopedTransaction trans0( db0.transaction() );
    trans0.start( false );
    if(!db0.exists()){
      db0.create();
    }
    std::set< std::string > conts = db0.containers();
    std::cout << "#### creating containers..."<<std::endl;
    if( conts.find( "Cont0" )== conts.end() ) db0.createContainer<int>("Cont0");
    trans0.commit();
    //**
    trans0.start( false );
    ora::Container contH0 = db0.containerHandle( "Cont0" );
    if( contH0.isLocked() ){
      std::cout <<"### Test ERROR: container should not be locked."<<std::endl;
      return -1;
    }
    std::cout << "#### locking..."<<std::endl;
    contH0.lock();
    if( contH0.isLocked() ){
      std::cout <<"### container has been locked..."<<std::endl;
    }
    std::cout << "#### writing..."<<std::endl;
    int myInt0(999);
    int myInt1(1234567890);
    contH0.insert( myInt0 );
    contH0.insert( myInt1 );
    contH0.flush();
    //::sleep(10);
    //db0.dropContainer( "Cont0" );
    trans0.commit();
    db0.disconnect();
    ::sleep(1);
    db0.connect( connStr );
    ora::ScopedTransaction trans1( db0.transaction() );
    trans1.start( false );
    ora::Container cnt = db0.containerHandle( "Cont0" );
    std::cout <<"### Container \""<<cnt.name()<<"\" has got "<<cnt.size()<<" objects."<<std::endl; 
    db0.drop();
    trans1.commit();
    db0.disconnect();
    serializer.release();


  } catch ( const ora::Exception& exc ){
    std::cout << "### ############# ERROR: "<<exc.what()<<std::endl;
  }
}

