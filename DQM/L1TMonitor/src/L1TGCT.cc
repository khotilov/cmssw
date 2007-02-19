/*
 * \file L1TGCT.cc
 *
 * $Date: 2007/02/02 06:01:40 $
 * $Revision: 1.00 $
 * \author J. Berryhill
 *
 */

#include "DQM/L1TMonitor/interface/L1TGCT.h"

using namespace std;
using namespace edm;

L1TGCT::L1TGCT(const ParameterSet& ps)
{

  // verbosity switch
  verbose_ = ps.getUntrackedParameter<bool>("verbose", false);

  if(verbose_) cout << "L1TGCT: constructor...." << endl;

  logFile_.open("L1TGCT.log");

  dbe = NULL;
  if ( ps.getUntrackedParameter<bool>("DaqMonitorBEInterface", false) ) 
  {
    dbe = Service<DaqMonitorBEInterface>().operator->();
    dbe->setVerbose(0);
  }

  monitorDaemon_ = false;
  if ( ps.getUntrackedParameter<bool>("MonitorDaemon", false) ) {
    Service<MonitorDaemon> daemon;
    daemon.operator->();
    monitorDaemon_ = true;
  }

  outputFile_ = ps.getUntrackedParameter<string>("outputFile", "");
  if ( outputFile_.size() != 0 ) {
    cout << "L1T Monitoring histograms will be saved to " << outputFile_.c_str() << endl;
  }
  else{
    outputFile_ = "L1TDQM.root";
  }

  bool disable = ps.getUntrackedParameter<bool>("disableROOToutput", false);
  if(disable){
    outputFile_="";
  }


  if ( dbe !=NULL ) {
    dbe->setCurrentFolder("L1TMonitor/L1TGCT");
  }


}

L1TGCT::~L1TGCT()
{
}

void L1TGCT::beginJob(const EventSetup& c)
{

  nev_ = 0;

  // get hold of back-end interface
  DaqMonitorBEInterface* dbe = 0;
  dbe = Service<DaqMonitorBEInterface>().operator->();

  if ( dbe ) {
    dbe->setCurrentFolder("L1TMonitor/L1TGCT");
    dbe->rmdir("L1TMonitor/L1TGCT");
  }


  if ( dbe ) 
  {
    dbe->setCurrentFolder("L1TMonitor/L1TGCT");
    
    gcttest = dbe->book1D("GCT test", 
       "GCT test", 128, -0.5, 127.5 ) ;
  }  
}


void L1TGCT::endJob(void)
{
  if(verbose_) cout << "L1TGCT: end job...." << endl;
  LogInfo("L1TGCT") << "analyzed " << nev_ << " events"; 

 if ( outputFile_.size() != 0  && dbe ) dbe->save(outputFile_);

 return;
}

void L1TGCT::analyze(const Event& e, const EventSetup& c)
{
  nev_++; 
  if(verbose_) cout << "L1TGCT: analyze...." << endl;

  int ntest = 5;
      gcttest->Fill(ntest);
      if (verbose_)
	{     
     std::cout << "\tGCT test " << ntest
   	    << std::endl;
	}
}

