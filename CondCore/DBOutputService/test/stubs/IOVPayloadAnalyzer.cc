#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondCore/DBCommon/interface/Exception.h"
#include "CondFormats/Calibration/interface/Pedestals.h"
#include "IOVPayloadAnalyzer.h"
IOVPayloadAnalyzer::IOVPayloadAnalyzer(const edm::ParameterSet& iConfig ):
  m_record(iConfig.getParameter< std::string >("record")){
  std::cout<<"IOVPayloadAnalyzer::IOVPayloadAnalyzer"<<std::endl;
}
IOVPayloadAnalyzer::~IOVPayloadAnalyzer(){
  std::cout<<"IOVPayloadAnalyzer::~IOVPayloadAnalyzer"<<std::endl;
}
void IOVPayloadAnalyzer::analyze( const edm::Event& evt, const edm::EventSetup& evtSetup){
  std::cout<<"IOVPayloadAnalyzer::analyze "<<std::endl;
  edm::Service<cond::service::PoolDBOutputService> mydbservice;
  if( !mydbservice.isAvailable() ){
    std::cout<<"Service is unavailable"<<std::endl;
    return;
  }
  unsigned int irun=evt.id().run();
  std::string tag=mydbservice->tag(m_record);
  std::cout<<"tag "<<tag<<std::endl;
  std::cout<<"run "<<irun<<std::endl;
  Pedestals* myped=new Pedestals;
  for(int ichannel=1; ichannel<=5; ++ichannel){
    Pedestals::Item item;
    item.m_mean=1.11*ichannel+irun;
    item.m_variance=1.12*ichannel+irun;
    myped->m_pedestals.push_back(item);
    }
  std::cout<<myped->m_pedestals[1].m_mean<<std::endl;
  if( mydbservice->isNewTagRequest(m_record) ){
    //create 
    cond::Time_t firstTillTime=mydbservice->endOfTime();
    std::cout<<myped->m_pedestals[1].m_mean<<std::endl;
    mydbservice->createNewIOV<Pedestals>(myped,firstTillTime,m_record);
    //std::cout<<"about to throw"<<std::endl;
    //throw cond::Exception("throwme");
    //std::cout<<"thrown"<<std::endl;
  }else{
    //append 
    if(mydbservice->currentTime()%5==0){
      mydbservice->appendSinceTime<Pedestals>(myped,mydbservice->currentTime(),m_record);
    }
  }
}
void IOVPayloadAnalyzer::endJob(){ 
  
}
DEFINE_FWK_MODULE(IOVPayloadAnalyzer);
