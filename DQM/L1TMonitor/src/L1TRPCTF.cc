/*
 * \file L1TRPCTF.cc
 *
 * $Date: 2007/09/28 08:52:46 $
 * $Revision: 1.6 $
 * \author J. Berryhill
 *
 */

#include "DQM/L1TMonitor/interface/L1TRPCTF.h"

using namespace std;
using namespace edm;

L1TRPCTF::L1TRPCTF(const ParameterSet& ps)
  : rpctfSource_( ps.getParameter< InputTag >("rpctfSource") )
 {

  // verbosity switch
  verbose_ = ps.getUntrackedParameter<bool>("verbose", false);

  if(verbose_) cout << "L1TRPCTF: constructor...." << endl;

  logFile_.open("L1TRPCTF.log");

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
    dbe->setCurrentFolder("L1T/L1TRPCTF");
  }


}

L1TRPCTF::~L1TRPCTF()
{
}

void L1TRPCTF::beginJob(const EventSetup& c)
{

  nev_ = 0;

  // get hold of back-end interface
  DaqMonitorBEInterface* dbe = 0;
  dbe = Service<DaqMonitorBEInterface>().operator->();

  if ( dbe ) {
    dbe->setCurrentFolder("L1T/L1TRPCTF");
    dbe->rmdir("L1T/L1TRPCTF");
  }


  if ( dbe ) 
  {
    dbe->setCurrentFolder("L1T/L1TRPCTF");
    
    rpctfetavalue[1] = dbe->book1D("RPCTF_eta_value", 
       "RPCTF eta value", 100, -2.5, 2.5 ) ;
    rpctfetavalue[2] = dbe->book1D("RPCTF_eta_value_+1", 
       "RPCTF eta value bx +1", 100, -2.5, 2.5 ) ;
    rpctfetavalue[0] = dbe->book1D("RPCTF_eta_value_-1", 
       "RPCTF eta value bx -1", 100, -2.5, 2.5 ) ;
    rpctfphivalue[1] = dbe->book1D("RPCTF_phi_value", 
       "RPCTF phi value", 100, 0.0, 6.2832 ) ;
    rpctfphivalue[2] = dbe->book1D("RPCTF_phi_value_+1", 
       "RPCTF phi value bx +1", 100, 0.0, 6.2832 ) ;
    rpctfphivalue[0] = dbe->book1D("RPCTF_phi_value_-1", 
       "RPCTF phi value bx -1", 100, 0.0, 6.2832 ) ;
    rpctfptvalue[1] = dbe->book1D("RPCTF_pt_value", 
       "RPCTF pt value", 160, -0.5, 159.5 ) ;
    rpctfptvalue[2] = dbe->book1D("RPCTF_pt_value_+1", 
       "RPCTF pt value bx +1", 160, -0.5, 159.5 ) ;
    rpctfptvalue[0] = dbe->book1D("RPCTF_pt_value_-1", 
       "RPCTF pt value bx -1", 160, -0.5, 159.5 ) ;
    rpctfchargevalue[1] = dbe->book1D("RPCTF_charge_value", 
       "RPCTF charge value", 3, -1.5, 1.5 ) ;
    rpctfchargevalue[2] = dbe->book1D("RPCTF_charge_value_+1", 
       "RPCTF charge value bx +1", 3, -1.5, 1.5 ) ;
    rpctfchargevalue[0] = dbe->book1D("RPCTF_charge_value_-1", 
       "RPCTF charge value bx -1", 3, -1.5, 1.5 ) ;
    rpctfquality[1] = dbe->book1D("RPCTF_quality", 
       "RPCTF quality", 20, -0.5, 19.5 ) ;
    rpctfquality[2] = dbe->book1D("RPCTF_quality_+1", 
       "RPCTF quality bx +1", 20, -0.5, 19.5 ) ;
    rpctfquality[0] = dbe->book1D("RPCTF_quality_-1", 
       "RPCTF quality bx -1", 20, -0.5, 19.5 ) ;
    rpctfntrack = dbe->book1D("RPCTF_ntrack", 
       "RPCTF ntrack", 20, -0.5, 19.5 ) ;
    rpctfbx = dbe->book1D("RPCTF_bx", 
       "RPCTF bx", 3, -1.5, 1.5 ) ;
  }  
}


void L1TRPCTF::endJob(void)
{
  if(verbose_) cout << "L1TRPCTF: end job...." << endl;
  LogInfo("L1TRPCTF") << "analyzed " << nev_ << " events"; 

 if ( outputFile_.size() != 0  && dbe ) dbe->save(outputFile_);

 return;
}

void L1TRPCTF::analyze(const Event& e, const EventSetup& c)
{
  nev_++; 
  if(verbose_) cout << "L1TRPCTF: analyze...." << endl;


  edm::Handle<L1MuGMTReadoutCollection> pCollection;


  try {
  e.getByLabel(rpctfSource_,pCollection);
  }
  catch (...) {
    edm::LogInfo("L1TRPCTF") << "can't find L1MuGMTReadoutCollection with label "
			       << rpctfSource_.label() ;
    return;
  }

  L1MuGMTReadoutCollection const* gmtrc = pCollection.product();
  vector<L1MuGMTReadoutRecord> gmt_records = gmtrc->getRecords();
  vector<L1MuGMTReadoutRecord>::const_iterator RRItr;

  int nrpctftrack = 0;
  for( RRItr = gmt_records.begin() ;
       RRItr != gmt_records.end() ;
       RRItr++ ) 
  {

    if (verbose_)
    {
     cout << "Readout Record " << RRItr->getBxInEvent()
   	    << endl;
   }
 
   vector<L1MuRegionalCand> RPCTFCands = RRItr->getBrlRPCCands();
 

   if (verbose_) 
    {
     cout << "RPCTFCands " << RPCTFCands.size()
   	    << endl;
    }

    for( vector<L1MuRegionalCand>::const_iterator 
         ECItr = RPCTFCands.begin() ;
         ECItr != RPCTFCands.end() ;
         ++ECItr ) 
    {

      int bxindex = ECItr->bx() + 1;
      if (ECItr->quality() > 0 ) {
      nrpctftrack++;

      if (verbose_)
	{  
     cout << "RPCTFCand bx " << ECItr->bx()
   	    << endl;
	}
     rpctfbx->Fill(ECItr->bx());

      rpctfetavalue[bxindex]->Fill(ECItr->etaValue());
      if (verbose_)
	{     
     cout << "\tRPCTFCand eta value " << ECItr->etaValue()
   	    << endl;
	}

      rpctfphivalue[bxindex]->Fill(ECItr->phiValue());
      if (verbose_)
	{     
     cout << "\tRPCTFCand phi value " << ECItr->phiValue()
   	    << endl;
	}

      rpctfptvalue[bxindex]->Fill(ECItr->ptValue());
      if (verbose_)
	{     
     cout << "\tRPCTFCand pt value " << ECItr->ptValue()
   	    << endl;
	}


      rpctfchargevalue[bxindex]->Fill(ECItr->chargeValue());
      if (verbose_)
	{     
     cout << "\tRPCTFCand charge value " << ECItr->chargeValue()
   	    << endl;
	}

      rpctfquality[bxindex]->Fill(ECItr->quality());
      if (verbose_)
	{     
     cout << "\tRPCTFCand quality " << ECItr->quality()
   	    << endl;
	}

      }
    }


  }

      rpctfntrack->Fill(nrpctftrack);
      if (verbose_)
	{     
     cout << "\tRPCTFCand ntrack " << nrpctftrack
   	    << endl;
	}
}

