#include "CondCore/DBCommon/interface/DBSession.h"
#include "CondCore/DBCommon/interface/Exception.h"
#include "CondCore/DBCommon/interface/SessionConfiguration.h"
#include "CondCore/DBCommon/interface/MessageLevel.h"
#include "CondCore/DBCommon/interface/Connection.h"
#include "CondCore/DBCommon/interface/CoralTransaction.h"
#include "CondCore/DBCommon/interface/TokenBuilder.h"
#include "CondCore/DBCommon/interface/LogDBEntry.h"
#include "CondCore/DBCommon/interface/Logger.h"
#include "CondCore/DBCommon/interface/UserLogInfo.h"

#include <string>
#include <iostream>
//#include <stdio.h>
//#include <time.h>
#include <unistd.h>

int main(){
  cond::TokenBuilder tk;
  tk.set("3E60FA40-D105-DA11-981C-000E0C4DE431",
	 "CondFormatsCalibration",
	 "Pedestals",
	 "PedestalsRcd",
	 0);
  std::string tok1=tk.tokenAsString();
  tk.set("3E60FA40-D105-DA11-981C-000E0C4DE431",
	 "CondFormatsCalibration",
	 "Pedestals",
	 "PedestalsRcd",
	 1);
  std::string tok2=tk.tokenAsString();
  //std::string constr("sqlite_file:mylog.db");
  std::string constr("oracle://devdb10/cms_xiezhen_dev");
  cond::DBSession* session=new cond::DBSession;
  session->configuration().setMessageLevel( cond::Error );
  session->configuration().setAuthenticationMethod(cond::XML);
  cond::Connection con(constr,-1);
  session->open();
  con.connect(session);
  //cond::CoralTransaction& coralTransaction=con.coralTransaction();
  // coralTransaction.start(false);
  cond::Logger mylogger(&con);
  cond::UserLogInfo a;
  a.provenance="me";
  mylogger.createLogDBIfNonExist();
  if(mylogger.getWriteLock()){
    std::cout<<"1. table locked"<<std::endl;
  }else{
    std::cout<<"1. table lock failed"<<std::endl;
  }
  mylogger.logOperationNow(a,constr,tok1,"mytag1","runnumber",0);
  std::cout<<"1. waiting"<<std::endl;
  sleep(5);
  std::cout<<"1. stop waiting"<<std::endl;
  if(mylogger.releaseWriteLock()){
    std::cout<<"1. table lock released"<<std::endl;
  }else{
    std::cout<<"1. failed to release table lock"<<std::endl;
  }
  if(mylogger.getWriteLock()){
    std::cout<<"1. table locked"<<std::endl;
  }else{
    std::cout<<"1. table lock failed"<<std::endl;
  }
  std::cout<<"1. waiting"<<std::endl;
  sleep(5);
  std::cout<<"1. stop waiting"<<std::endl;
  mylogger.logFailedOperationNow(a,constr,tok1,"mytag1","runnumber",1,"EOOROR");
  std::cout<<"1. waiting"<<std::endl;
  sleep(5);
  std::cout<<"1. stop waiting"<<std::endl;
  if(mylogger.releaseWriteLock()){
    std::cout<<"1. table lock released"<<std::endl;
  }else{
    std::cout<<"1. failed to release table lock"<<std::endl;
  }
  
  if(mylogger.getWriteLock()){
    std::cout<<"1. table locked"<<std::endl;
  }else{
    std::cout<<"1. table lock failed"<<std::endl;
  }
  
  std::cout<<"1. waiting"<<std::endl;
  sleep(5);
  std::cout<<"1. stop waiting"<<std::endl;
  mylogger.logOperationNow(a,constr,tok2,"mytag","runnumber",1);
  std::cout<<"1. waiting"<<std::endl;
  sleep(5);
  std::cout<<"1. stop waiting"<<std::endl;
  if(mylogger.releaseWriteLock()){
    std::cout<<"1. table lock released"<<std::endl;
  }else{
    std::cout<<"1. failed to release table lock"<<std::endl;
  }   
  /*std::cout<<"about to lookup last entry"<<std::endl;
  cond::LogDBEntry result;
  mylogger.LookupLastEntryByProvenance("me",result);
  std::cout<<"result \n";
  std::cout<<"logId "<<result.logId<<"\n";
  std::cout<<"destinationDB "<<result.destinationDB<<"\n";
  std::cout<<"provenance "<<result.provenance<<"\n";
  std::cout<<"usertext "<<result.usertext<<"\n";
  std::cout<<"iovtag "<<result.iovtag<<"\n";
  std::cout<<"iovtimetype "<<result.iovtimetype<<"\n";
  std::cout<<"payloadIdx "<<result.payloadIdx<<"\n";
  std::cout<<"payloadName "<<result.payloadName<<"\n";
  std::cout<<"payloadToken "<<result.payloadToken<<"\n";
  std::cout<<"payloadContainer "<<result.payloadContainer<<"\n";
  std::cout<<"exectime "<<result.exectime<<"\n";
  std::cout<<"execmessage "<<result.execmessage<<std::endl;
  cond::LogDBEntry result2;
  mylogger.LookupLastEntryByTag("mytag1",result2);
  std::cout<<"result2 \n";
  std::cout<<"logId "<<result2.logId<<"\n";
  std::cout<<"destinationDB "<<result2.destinationDB<<"\n";
  std::cout<<"provenance "<<result2.provenance<<"\n";
  std::cout<<"usertext "<<result2.usertext<<"\n";
  std::cout<<"iovtag "<<result2.iovtag<<"\n";
  std::cout<<"iovtimetype "<<result2.iovtimetype<<"\n";
  std::cout<<"payloadIdx "<<result2.payloadIdx<<"\n";
  std::cout<<"payloadName "<<result2.payloadName<<"\n";
  std::cout<<"payloadToken "<<result2.payloadToken<<"\n";
  std::cout<<"payloadContainer "<<result2.payloadContainer<<"\n";
  std::cout<<"exectime "<<result2.exectime<<"\n";
  std::cout<<"execmessage "<<result2.execmessage<<std::endl;
  */
  //coralTransaction.commit();
  con.disconnect();
  delete session;
}
