/*
 * \file EcalEndcapMonitorDbModule.cc
 * 
 * $Date: 2007/04/05 13:56:47 $
 * $Revision: 1.2 $
 * \author G. Della Ricca
 *
*/

#include <unistd.h>

#include "FWCore/ServiceRegistry/interface/Service.h"

#include <iostream>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DQMServices/Daemon/interface/MonitorDaemon.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SealKernel/Context.h"
#include "SealKernel/ComponentLoader.h"
#include "SealKernel/Exception.h"
#include "SealKernel/IMessageService.h"
#include "PluginManager/PluginManager.h"
#include "RelationalAccess/IConnectionService.h"
#include "RelationalAccess/IConnectionServiceConfiguration.h"

#include "CoralBase/Attribute.h"
#include "CoralBase/AttributeList.h"
#include "CoralBase/AttributeSpecification.h"


#include <DQM/EcalEndcapMonitorDbModule/interface/EcalEndcapMonitorDbModule.h>


EcalEndcapMonitorDbModule::EcalEndcapMonitorDbModule(const edm::ParameterSet& ps){

  // get hold of back-end interface
  dbe_ = edm::Service<DaqMonitorBEInterface>().operator->();

  // MonitorDaemon switch
  enableMonitorDaemon_ = ps.getUntrackedParameter<bool>("enableMonitorDaemon", true);

  if ( enableMonitorDaemon_ ) {
    std::cout << " enableMonitorDaemon switch is ON" << std::endl;
    edm::Service<MonitorDaemon> daemon;
    daemon.operator->();
  } else {
    std::cout << " enableMonitorDaemon switch is OFF" << std::endl;
  }

  xmlFile_ = ps.getUntrackedParameter<std::string>( "xmlFile", "" );
  if ( xmlFile_.size() != 0 ) {
    std::cout << "Monitor Elements from DB xml source file is " << xmlFile_ << std::endl;
  }

  sleepTime_ = ps.getUntrackedParameter<int>( "sleepTime", 0 );
  std::cout << "Sleep time is " << sleepTime_ << " second(s)." << std::endl;
  
  // html output directory
  htmlDir_ = ps.getUntrackedParameter<std::string>("htmlDir", ".");

  if ( htmlDir_.size() != 0 ) {
    std::cout << " HTML output will go to"
	      << " htmlDir = " << htmlDir_ << std::endl;
  } else {
    std::cout << " HTML output is disabled" << std::endl;
  }

  ME_Db_ = new MonitorElementsDb( ps, xmlFile_ );
  
  if ( dbe_ ) dbe_->showDirStructure();
  
}

EcalEndcapMonitorDbModule::~EcalEndcapMonitorDbModule(){

  if ( ME_Db_ ) delete ME_Db_;

}

void EcalEndcapMonitorDbModule::beginJob(const edm::EventSetup& c){

  icycle_ = 0;

  if ( ME_Db_ ) ME_Db_->beginJob(c);

}

void EcalEndcapMonitorDbModule::endJob(void) {

  if ( ME_Db_ ) ME_Db_->endJob();

  std::cout << "EcalEndcapMonitorDbModule: endJob, icycle = " << icycle_ << std::endl;

}

void EcalEndcapMonitorDbModule::analyze(const edm::Event& e, const edm::EventSetup& c){
  
  icycle_++;
  
  std::cout << "EcalEndcapMonitorDbModule: icycle = " << icycle_ << std::endl;

  try {
    seal::Handle<seal::Context> context = new seal::Context;
    seal::PluginManager* pm = seal::PluginManager::get();
    pm->initialise ();
    seal::Handle<seal::ComponentLoader> loader = new seal::ComponentLoader(context.get());

    loader->load("SEAL/Services/MessageService");

    std::vector<seal::Handle<seal::IMessageService> > v_msgSvc;
    context->query(v_msgSvc);
    if ( ! v_msgSvc.empty() ) {
      seal::Handle<seal::IMessageService>& msgSvc = v_msgSvc.front();
      msgSvc->setOutputLevel(seal::Msg::Error);
      //msgSvc->setOutputLevel(seal::Msg::Debug);
    }
    
    loader->load("CORAL/Services/ConnectionService");
    
    loader->load("CORAL/Services/EnvironmentAuthenticationService");

    seal::IHandle<coral::IConnectionService> connectionService = context->query<coral::IConnectionService>("CORAL/Services/ConnectionService");

    loader->load("CORAL/RelationalPlugins/oracle");

    // Set configuration parameters
    coral::IConnectionServiceConfiguration& config = connectionService->configuration();
    config.setConnectionRetrialPeriod(1);
    config.setConnectionRetrialTimeOut(10);

    session_ = connectionService->connect("ECAL CondDB", coral::ReadOnly);    

    if ( ME_Db_ ) ME_Db_->analyze(e, c, session_ );
    
  } catch (coral::Exception& e) {
    std::cerr << "CORAL Exception : " << e.what() << std::endl;
  } catch (std::exception& e) {
    std::cerr << "Standard C++ exception : " << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Exception caught (...)" << std::endl;
  }

  if ( htmlDir_.size() != 0 ) {

    ME_Db_->htmlOutput( htmlDir_ );

  }

  
  delete session_;

  sleep( sleepTime_ );

}

