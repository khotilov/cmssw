#include "DQM/HcalMonitorTasks/interface/HcalDeadCellMonitor.h"

#define OUT if(fverbosity_)cout
#define BITSHIFT 5

using namespace std;

HcalDeadCellMonitor::HcalDeadCellMonitor()
{
  ievt_=0;
  // Default initialization
  showTiming   = false;
  fVerbosity   = 0;
  deadmon_makeDiagnostics_ = false;
} //constructor

HcalDeadCellMonitor::~HcalDeadCellMonitor()
{
} //destructor


/* ------------------------------------ */ 

void HcalDeadCellMonitor::setup(const edm::ParameterSet& ps,
				DQMStore* dbe)
{
  HcalBaseMonitor::setup(ps,dbe);
  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }
  baseFolder_ = rootFolder_+"DeadCellMonitor_Hcal";
  if (fVerbosity>0)
    std::cout <<"<HcalDeadCellMonitor::setup>  Setting up histograms"<<std::endl;

  // Assume subdetectors not present until shown otherwise
  HBpresent_ =false;
  HEpresent_ =false;
  HOpresent_ =false;
  HFpresent_ =false;
  ZDCpresent_=false;

  // Dead Cell Monitor - specific cfg variables

  if (fVerbosity>1)
    std::cout <<"<HcalDeadCellMonitor::setup>  Getting variable values from cfg files"<<std::endl;

  // deadmon_makeDiagnostics_ will take on base task value unless otherwise specified
  deadmon_makeDiagnostics_ = ps.getUntrackedParameter<bool>("DeadCellMonitor_makeDiagnosticPlots",makeDiagnostics);
  
  // Set checkNevents values
  deadmon_checkNevents_ = ps.getUntrackedParameter<int>("DeadCellMonitor_checkNevents",checkNevents_);
  // increases rate at which neverpresent tests are run
  deadmon_neverpresent_prescale_     = ps.getUntrackedParameter<int>("DeadCellMonitor_neverpresent_prescale",1);  

  // Set which dead cell checks will be performed
  /* Dead cells can be defined in three ways:
     1)  never present -- digi is never present in run
     2)  occupancy -- digi is absent for (checkNevents_) consecutive events
     3)  energy -- cell is present, but rechit energy is never above threshold value
  */
  deadmon_test_neverpresent_           = ps.getUntrackedParameter<bool>("DeadCellMonitor_test_neverpresent",true);
  deadmon_test_occupancy_         = ps.getUntrackedParameter<bool>("DeadCellMonitor_test_occupancy", true);
  deadmon_test_energy_            = ps.getUntrackedParameter<bool>("DeadCellMonitor_test_energy", true);

  // rechit_occupancy duplicates digi occupancy -- ignore it
  //deadmon_test_rechit_occupancy_  = ps.getUntrackedParameter<bool>("DeadCellMonitor_test_rechit_occupancy", true);

  // Old tests (in version 1.55 and earlier) checked pedestals, compared energies to neighbors
  // Are these ever going to be useful?  If so, we could re-enable them.
  //deadmon_test_pedestal_          = ps.getUntrackedParameter<bool>("DeadCellMonitor_test_pedestal",  true);
  //deadmon_test_neighbor_          = ps.getUntrackedParameter<bool>("DeadCellMonitor_test_neighbor",  true);
  
  deadmon_minErrorFlag_ = ps.getUntrackedParameter<double>("DeadCellMonitor_minErrorFlag",0.0);

  // rechit energy test -- cell must be below threshold value for a number of consecutive events to be considered dead
  energyThreshold_       = ps.getUntrackedParameter<double>("DeadCellMonitor_energyThreshold",                  0);
  HBenergyThreshold_     = ps.getUntrackedParameter<double>("DeadCellMonitor_HB_energyThreshold",energyThreshold_);
  HEenergyThreshold_     = ps.getUntrackedParameter<double>("DeadCellMonitor_HE_energyThreshold",energyThreshold_);
  HOenergyThreshold_     = ps.getUntrackedParameter<double>("DeadCellMonitor_HO_energyThreshold",energyThreshold_);
  HFenergyThreshold_     = ps.getUntrackedParameter<double>("DeadCellMonitor_HF_energyThreshold",energyThreshold_);
  ZDCenergyThreshold_    = ps.getUntrackedParameter<double>("DeadCellMonitor_ZDC_energyThreshold",           -999);

  // neighboring-cell tests -- parameters no longer used

  // Set initial event # to 0
  ievt_=0;

  zeroCounters(true);
  
  // Set up histograms
  if (m_dbe)
    {
      if (fVerbosity>1)
	std::cout <<"<HcalDeadCellMonitor::setup>  Setting up histograms"<<std::endl;

      m_dbe->setCurrentFolder(baseFolder_);
      meEVT_ = m_dbe->bookInt("Dead Cell Task Event Number");
      meEVT_->Fill(ievt_);

      // Create problem cell plots
      // Overall plot gets an initial " " in its name
      ProblemDeadCells=m_dbe->book2D(" ProblemDeadCells",
                                     " Problem Dead Cell Rate for all HCAL",
                                     etaBins_,etaMin_,etaMax_,
                                     phiBins_,phiMin_,phiMax_);
      ProblemDeadCells->setAxisTitle("i#eta",1);
      ProblemDeadCells->setAxisTitle("i#phi",2);

      // 1D plots count number of bad cells
      NumberOfDeadCells=m_dbe->book1D("Problem_TotalDeadCells_HCAL",
				      "Total Number of Dead Hcal Cells",
				      9091,-0.5,9090.5);
      NumberOfDeadCellsHB=m_dbe->book1D("Problem_TotalDeadCells_HB",
					"Total Number of Dead HB Cells",
					2593,-0.5,2592.5);
      NumberOfDeadCellsHE=m_dbe->book1D("Problem_TotalDeadCells_HE",
					"Total Number of Dead HE Cells",
					2593,-0.5,2592.5);
      NumberOfDeadCellsHO=m_dbe->book1D("Problem_TotalDeadCells_HO",
					"Total Number of Dead HO Cells",
					2161,-0.5,2160.5);
      NumberOfDeadCellsHF=m_dbe->book1D("Problem_TotalDeadCells_HF",
					"Total Number of Dead HF Cells",
					1729,-0.5,1728.5);
      NumberOfDeadCellsZDC=m_dbe->book1D("Problem_TotalDeadCells_ZDC",
					"Total Number of Dead ZDC Cells",
					 19,-0.5,18.5);

      // Overall Problem plot appears in main directory; plots by depth appear \in subdirectory
      m_dbe->setCurrentFolder(baseFolder_+"/problem_deadcells");
      setupDepthHists2D(ProblemDeadCellsByDepth, " Problem Dead Cell Rate","");

      // Set up plots for each failure mode of dead cells
      stringstream units; // We'll need to set the titles individually, rather than passing units to setupDepthHists2D (since this also would affect the name of the histograms)
      stringstream name;
      if (deadmon_test_neverpresent_)
	{
	  m_dbe->setCurrentFolder(baseFolder_+"/dead_digi_never_present");
	  setupDepthHists2D(DigisNeverPresentByDepth,
			    "Dead Cells with No Digis Ever","");
	  // 1D plots count number of bad cells
	  NumberOfNeverPresentCells=m_dbe->book1D("Problem_TotalNeverPresentCells_HCAL",
					  "Total Number of Never-Present Hcal Cells",
					  9091,-0.5,9090.5);
	  NumberOfNeverPresentCellsHB=m_dbe->book1D("Problem_NeverPresentCells_HB",
					    "Total Number of Never-Present HB Cells",
					    2593,-0.5,2592.5);
	  NumberOfNeverPresentCellsHE=m_dbe->book1D("Problem_NeverPresentCells_HE",
					    "Total Number of Never-Present HE Cells",
					    2593,-0.5,2592.5);
	  NumberOfNeverPresentCellsHO=m_dbe->book1D("Problem_NeverPresentCells_HO",
					    "Total Number of Never-Present HO Cells",
					    2161,-0.5,2160.5);
	  NumberOfNeverPresentCellsHF=m_dbe->book1D("Problem_NeverPresentCells_HF",
					    "Total Number of Never-Present HF Cells",
					    1729,-0.5,1728.5);
	  NumberOfNeverPresentCellsZDC=m_dbe->book1D("Problem_NeverPresentCells_ZDC",
					     "Total Number of Never-Present ZDC Cells",
					     19,-0.5,18.5);
	}
      if (deadmon_test_occupancy_)
	{
	  m_dbe->setCurrentFolder(baseFolder_+"/dead_digi_often_missing");
	  //units<<"("<<deadmon_checkNevents_<<" consec. events)";
	  name<<"Dead Cells with No Digis";
	  setupDepthHists2D(UnoccupiedDeadCellsByDepth,
			    (char*)(name.str().c_str()),
			    "");
	  name.str("");
	  name<<"HBHF Depth 1 Dead Cells with No Digis for "<<deadmon_checkNevents_<<" Consecutive Events";
	  UnoccupiedDeadCellsByDepth[0]->setTitle(name.str().c_str());
	  name.str("");
	  name<<"HBHF Depth 2 Dead Cells with No Digis for "<<deadmon_checkNevents_<<" Consecutive Events";
	  UnoccupiedDeadCellsByDepth[1]->setTitle(name.str().c_str());
	  name.str("");
	  name<<"HE Depth 3 Dead Cells with No Digis for "<<deadmon_checkNevents_<<" Consecutive Events";
	  UnoccupiedDeadCellsByDepth[2]->setTitle(name.str().c_str());
	  name.str("");
	  name<<"HO/ZDC Dead Cells with No Digis for "<<deadmon_checkNevents_<<" Consecutive Events";
	  name.str("");
	  UnoccupiedDeadCellsByDepth[3]->setTitle(name.str().c_str());
	  name<<"HE Depth 1 Dead Cells with No Digis for "<<deadmon_checkNevents_<<" Consecutive Events";
	  UnoccupiedDeadCellsByDepth[4]->setTitle(name.str().c_str());
	  name.str("");
	  name<<"HE Depth 2 Dead Cells with No Digis for "<<deadmon_checkNevents_<<" Consecutive Events";
	  UnoccupiedDeadCellsByDepth[5]->setTitle(name.str().c_str());
	  name.str("");

	  // 1D plots count number of bad cells
	  name<<"Total Number of Hcal Digis Unoccupied for "<<deadmon_checkNevents_<<" Consecutive Events";
	  NumberOfUnoccupiedCells=m_dbe->book1D("Problem_TotalUnoccupiedCells_HCAL",
						(char*)(name.str().c_str()),
						9091,-0.5,9090.5);
	  name.str("");
	  name<<"Total Number of HB Digis Unoccupied for "<<deadmon_checkNevents_<<" Consecutive Events";
	  NumberOfUnoccupiedCellsHB=m_dbe->book1D("Problem_UnoccupiedCells_HB",
						  (char*)(name.str().c_str()),
						  2593,-0.5,2592.5);
	  name.str("");
	  name<<"Total Number of HE Digis Unoccupied for "<<deadmon_checkNevents_<<" Consecutive Events";
	  NumberOfUnoccupiedCellsHE=m_dbe->book1D("Problem_UnoccupiedCells_HE",
						  (char*)(name.str().c_str()),
						  2593,-0.5,2592.5);
	  name.str("");
	  name<<"Total Number of HO Digis Unoccupied for "<<deadmon_checkNevents_<<" Consecutive Events";
	  NumberOfUnoccupiedCellsHO=m_dbe->book1D("Problem_UnoccupiedCells_HO",
						  (char*)(name.str().c_str()),
						  2161,-0.5,2160.5);
	  name.str("");
	  name<<"Total Number of HF Digis Unoccupied for "<<deadmon_checkNevents_<<" Consecutive Events";
	  NumberOfUnoccupiedCellsHF=m_dbe->book1D("Problem_UnoccupiedCells_HF",
						  (char*)(name.str().c_str()),
						  1729,-0.5,1728.5);
	  name.str("");
	  name<<"Total Number of ZDC Digis Unoccupied for "<<deadmon_checkNevents_<<" Consecutive Events";
	  NumberOfUnoccupiedCellsZDC=m_dbe->book1D("Problem_UnoccupiedCells_ZDC",
						   (char*)(name.str().c_str()),
						   19,-0.5,18.5);
	}
      
      if (deadmon_test_energy_)
	{
	  m_dbe->setCurrentFolder(baseFolder_+"/dead_energytest");
	  setupDepthHists2D(BelowEnergyThresholdCellsByDepth,"Dead Cells Failing Energy Threshold Test","");
	  // set more descriptive titles for threshold plots
	  units.str("");
	  units<<"Dead Cells with Consistent Low Energy Depth 1 -- HB <"<<HBenergyThreshold_<<" GeV, HF <"<<HFenergyThreshold_<<" GeV";
	  BelowEnergyThresholdCellsByDepth[0]->setTitle(units.str().c_str());
	  units.str("");
	  units<<"Dead Cells with Consistent Low Energy Depth 2 -- HB <"<<HBenergyThreshold_<<" GeV, HF <"<<HFenergyThreshold_<<" GeV";
	  BelowEnergyThresholdCellsByDepth[1]->setTitle(units.str().c_str());
	  units.str("");
	  units<<"Dead Cells with Consistent Low Energy Depth 3 -- HE <"<<HEenergyThreshold_<<" GeV";
	  BelowEnergyThresholdCellsByDepth[2]->setTitle(units.str().c_str());
	  units.str("");
	  units<<"Dead Cells with Consistent Low Energy Depth 4 -- HO <"<<HOenergyThreshold_<<" GeV, ZDC <"<< ZDCenergyThreshold_ ;
	  BelowEnergyThresholdCellsByDepth[3]->setTitle(units.str().c_str());
	  units.str("");
	  units<<"Dead Cells with Consistent Low Energy Depth 1 -- HE <"<<HEenergyThreshold_<<" GeV";
	  BelowEnergyThresholdCellsByDepth[4]->setTitle(units.str().c_str());
	  units.str("");
	  units<<"Dead Cells with Consistent Low Energy Depth 2 -- HE <"<<HEenergyThreshold_<<" GeV";
	  BelowEnergyThresholdCellsByDepth[5]->setTitle(units.str().c_str());
	  units.str("");

	  // 1D plots count number of bad cells
	  name<<"Total Number of Hcal RecHits with Consistent Low Energy";
	  NumberOfBelowEnergyCells=m_dbe->book1D("Problem_TotalBelowEnergyCells_HCAL",
						(char*)(name.str().c_str()),
						9091,-0.5,9090.5);
	  name.str("");
	  name<<"Total Number of HB RecHits with Consistent Low Energy < "<<HBenergyThreshold_<<" GeV";
	  NumberOfBelowEnergyCellsHB=m_dbe->book1D("Problem_BelowEnergyCells_HB",
						  (char*)(name.str().c_str()),
						  2593,-0.5,2592.5);
	  name.str("");
	  name<<"Total Number of HE RecHits with Consistent Low Energy < "<<HEenergyThreshold_<<" GeV";
	  NumberOfBelowEnergyCellsHE=m_dbe->book1D("Problem_BelowEnergyCells_HE",
						  (char*)(name.str().c_str()),
						  2593,-0.5,2592.5);
	  name.str("");
	  name<<"Total Number of HO RecHits with Consistent Low Energy < "<<HOenergyThreshold_<<" GeV";
	  NumberOfBelowEnergyCellsHO=m_dbe->book1D("Problem_BelowEnergyCells_HO",
						  (char*)(name.str().c_str()),
						  2161,-0.5,2160.5);
	  name.str("");
	  name<<"Total Number of HF RecHits with Consistent Low Energy < "<<HFenergyThreshold_<<" GeV";
	  NumberOfBelowEnergyCellsHF=m_dbe->book1D("Problem_BelowEnergyCells_HF",
						  (char*)(name.str().c_str()),
						  1729,-0.5,1728.5);
	  name.str("");
	  name<<"Total Number of ZDC RecHits with Consistent Low Energy < "<<ZDCenergyThreshold_<<" GeV";
	  NumberOfBelowEnergyCellsZDC=m_dbe->book1D("Problem_BelowEnergyCells_ZDC",
						   (char*)(name.str().c_str()),
						   19,-0.5,18.5);
	}

    } // if (m_dbe)

  return;
} //void HcalDeadCellMonitor::setup(...)

/* --------------------------- */
void HcalDeadCellMonitor::setupNeighborParams(const edm::ParameterSet& ps,
					      neighborParams& N,
					      char* type)
{
  // This is no longer used -- can remove at some point, after further testing
  // sets up parameters for neighboring-cell algorithm for each subdetector
  /*
  ostringstream myname;
  myname<<"DeadCellMonitor_"<<type<<"_neighbor_deltaIphi";
  N.DeltaIphi = ps.getUntrackedParameter<int>(myname.str().c_str(),
					      defaultNeighborParams_.DeltaIphi);
  myname.str("");
  myname<<"DeadCellMonitor_"<<type<<"_neighbor_deltaIeta";
  N.DeltaIeta = ps.getUntrackedParameter<int>(myname.str().c_str(),
					      defaultNeighborParams_.DeltaIeta);
  myname.str("");
  myname<<"DeadCellMonitor_"<<type<<"_neighbor_deltaDepth";
  N.DeltaDepth = ps.getUntrackedParameter<int>(myname.str().c_str(),
					       defaultNeighborParams_.DeltaDepth);
  myname.str("");
  myname<<"DeadCellMonitor_"<<type<<"_neighbor_maxCellEnergy";
  N.maxCellEnergy = ps.getUntrackedParameter<double>(myname.str().c_str(),
						     defaultNeighborParams_.maxCellEnergy);
  myname.str("");
  myname<<"DeadCellMonitor_"<<type<<"_neighbor_minNeighborEnergy";
  N.minNeighborEnergy = ps.getUntrackedParameter<double>(myname.str().c_str(),
							 defaultNeighborParams_.minNeighborEnergy);
  myname.str("");
  myname<<"DeadCellMonitor_"<<type<<"_neighbor_minGoodNeighborFrac";
  N.minGoodNeighborFrac = ps.getUntrackedParameter<double>(myname.str().c_str(),
							   defaultNeighborParams_.minGoodNeighborFrac);
  myname.str("");
  myname<<"DeadCellMonitor_"<<type<<"_neighbor_maxEnergyFrac";
  N.maxEnergyFrac = ps.getUntrackedParameter<double>(myname.str().c_str(),
						     defaultNeighborParams_.maxEnergyFrac);
  */
  return;
} // void HcalDeadCellMonitor::setupNeighborParams

/* --------------------------- */

void HcalDeadCellMonitor::reset(){}  // reset function is empty for now

/* --------------------------- */



void HcalDeadCellMonitor::done(std::map<HcalDetId, unsigned int>& myqual)
{
  if (dump2database==0) 
    return;

  return;  // this is now done within the client, rather than the task (so that it doesn't get done multiple times when runing offline DQM)
  // Dump to ascii file for database -- now taken care of through ChannelStatus objects
  /*
  char buffer [1024];
  
  ofstream fOutput("hcalDeadCells.txt", ios::out);
  sprintf (buffer, "# %15s %15s %15s %15s %8s %10s\n", "eta", "phi", "dep", "det", "value", "DetId");
  fOutput << buffer;
  */

  int eta,phi;
  float binval;
  int mydepth;

  int subdet;
  char* subdetname;
  if (fVerbosity>1)
    {
      std::cout <<"<HcalDeadCellMonitor>  Summary of Dead Cells in Run: "<<std::endl;
      std::cout <<"(Error rate must be >= "<<deadmon_minErrorFlag_*100.<<"% )"<<std::endl;  
    }
  for (int ieta=1;ieta<=etaBins_;++ieta)
    {
      for (int iphi=1;iphi<=phiBins_;++iphi)
        {
          eta=ieta+int(etaMin_)-1;
          phi=iphi+int(phiMin_)-1;

          for (int d=0;d<6;++d)
            {
	      binval=ProblemDeadCellsByDepth[d]->getBinContent(ieta,iphi);
	     
	      // Set subdetector labels for output
	      if (d<2) // HB/HF
		{
		  if (abs(eta)<29)
		    {
		      subdetname="HB";
		      subdet=1;
		    }
		  else
		    {
		      subdetname="HF";
		      subdet=4;
		    }
		}
	      else if (d==3)
		{
		  if (abs(eta)==43)
		    {
		      subdetname="ZDC";
		      subdet=7; // correct value??
		    }
		  else
		    {
		      subdetname="HO";
		      subdet=3;
		    }
		}
	      else
		{
		  subdetname="HE";
		  subdet=2;
		}
	      // Set correct depth label
	      if (d>3)
		mydepth=d-3;
	      else
		mydepth=d+1;
	      HcalDetId myid((HcalSubdetector)(subdet), eta, phi, mydepth);
	      if (!validDetId((HcalSubdetector)(subdet), eta, phi, mydepth))
		continue;

	      if (fVerbosity>0 && binval>deadmon_minErrorFlag_)
		std::cout <<"Dead Cell "<<subdet<<"("<<eta<<", "<<phi<<", "<<mydepth<<"):  "<<binval*100.<<"%"<<std::endl;
	      int value = 0;
	      if (binval>deadmon_minErrorFlag_)
		value=1;
	      
	      if (value==1)
	      if (myqual.find(myid)==myqual.end())
		{
		  myqual[myid]=(value<<BITSHIFT);  // deadcell shifted to bit 5
		}
	      else
		{
		  int mask=(1<<BITSHIFT);
		  if (value==1)
		    myqual[myid] |=mask;

		  else
		    myqual[myid] &=~mask;
		  if (value==1 && fVerbosity>1) std::cout <<"myqual = "<<std::hex<<myqual[myid]<<std::dec<<"  MASK = "<<std::hex<<mask<<std::dec<<std::endl;
		}
	      /*
	      sprintf(buffer, "  %15i %15i %15i %15s %8X %10X \n",eta,phi,mydepth,subdetname,(value<<BITSHIFT),int(myid.rawId()));
	      fOutput<<buffer;
	      */
	    } // for (int d=0;d<6;++d) // loop over depth histograms
	} // for (int iphi=1;iphi<=phiBins_;++iphi)
    } // for (int ieta=1;ieta<=etaBins_;++ieta)
  //fOutput.close();

  return;

} // void HcalDeadCellMonitor::done()



/* --------------------------------- */

void HcalDeadCellMonitor::clearME()
{
  // I don't think this function gets cleared any more.  
  // And need to add code to clear out subfolders as well?
  if (m_dbe)
    {
      m_dbe->setCurrentFolder(baseFolder_);
      m_dbe->removeContents();
    }
  return;
} // void HcalDeadCellMonitor::clearME()

/* -------------------------------- */


void HcalDeadCellMonitor::processEvent(const HBHERecHitCollection& hbHits,
				       const HORecHitCollection& hoHits,
				       const HFRecHitCollection& hfHits,
				       //const ZDCRecHitCollection& zdcHits,
				       const HBHEDigiCollection& hbhedigi,
				       const HODigiCollection& hodigi,
				       const HFDigiCollection& hfdigi
				       //const ZDCDigiCollection& zdcdigi
				       )
{

  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

  ++ievt_;
  if (m_dbe) meEVT_->Fill(ievt_);

  // HBpresent_, HEpresent need to be determined within loop, since HBHE is a single collection
  HOpresent_ = (hodigi.size()>0||hoHits.size()>0);
  HFpresent_ = (hfdigi.size()>0||hfHits.size()>0);
  //ZDCpresent_ = (zdcdigi.size()>0 || zdcHits.size()>0);

  if (fVerbosity>1) std::cout <<"<HcalDeadCellMonitor::processEvent> Processing event..."<<std::endl;

  // Do Digi-Based dead cell searches 

  if (deadmon_test_neverpresent_ || deadmon_test_occupancy_)
    {

      // Dummy fills
      for (unsigned int i=0;i<UnoccupiedDeadCellsByDepth.size();++i)
	{
	  UnoccupiedDeadCellsByDepth[i]->setBinContent(0,0,ievt_);
	  DigisNeverPresentByDepth[i]->setBinContent(0,0,ievt_);
	}
      NumberOfNeverPresentCells->setBinContent(0,ievt_);
      NumberOfNeverPresentCellsHB->setBinContent(0,ievt_);
      NumberOfNeverPresentCellsHE->setBinContent(0,ievt_);
      NumberOfNeverPresentCellsHO->setBinContent(0,ievt_);
      NumberOfNeverPresentCellsHF->setBinContent(0,ievt_);
      NumberOfNeverPresentCellsZDC->setBinContent(0,ievt_);
      
      NumberOfUnoccupiedCells->setBinContent(0,ievt_);
      NumberOfUnoccupiedCellsHB->setBinContent(0,ievt_);
      NumberOfUnoccupiedCellsHE->setBinContent(0,ievt_);
      NumberOfUnoccupiedCellsHO->setBinContent(0,ievt_);
      NumberOfUnoccupiedCellsHF->setBinContent(0,ievt_);
      NumberOfUnoccupiedCellsZDC->setBinContent(0,ievt_);

      if (showTiming)
	{
	  cpu_timer.reset(); cpu_timer.start();
	}
      for (HBHEDigiCollection::const_iterator j=hbhedigi.begin();
	   j!=hbhedigi.end(); ++j)
	{
	  processEvent_HBHEdigi(j);
	}
      
      for (HODigiCollection::const_iterator j=hodigi.begin();
	   j!=hodigi.end(); ++j)
	{
	  processEvent_HOdigi(j);
	}
      
      for (HFDigiCollection::const_iterator j=hfdigi.begin();
	   j!=hfdigi.end(); ++j)
	{
	  processEvent_HFdigi(j);
	}
      /*
	for (ZDCDigiCollection::const_iterator j=zdcdigi.begin();
	j!=zdcdigi.end(); ++j)
	{
	processEvent_ZDCdigi(j);
	}
      */
      if (showTiming)
	{
	  cpu_timer.stop();  std::cout <<"TIMER:: HcalDeadCellMonitor PROCESSEVENT_DIGI -> "<<cpu_timer.cpuTime()<<std::endl;
	}
    }
  
  // Search for "dead" cells below a certain energy
  if (deadmon_test_energy_) 
    {
      if (showTiming)
	{
	  cpu_timer.reset(); cpu_timer.start();
	}

      // Dummy Fills
      for (unsigned int i=0;i<BelowEnergyThresholdCellsByDepth.size();++i)
	{
	BelowEnergyThresholdCellsByDepth[i]->setBinContent(0,0,ievt_);
	}
      NumberOfBelowEnergyCells->setBinContent(0,ievt_);
      NumberOfBelowEnergyCellsHB->setBinContent(0,ievt_);
      NumberOfBelowEnergyCellsHE->setBinContent(0,ievt_);
      NumberOfBelowEnergyCellsHO->setBinContent(0,ievt_);
      NumberOfBelowEnergyCellsHF->setBinContent(0,ievt_);
      NumberOfBelowEnergyCellsZDC->setBinContent(0,ievt_);

      for (HBHERecHitCollection::const_iterator j=hbHits.begin();
	   j!=hbHits.end(); ++j)
	{
	  processEvent_HBHERecHit(j);
	}
      
      for (HORecHitCollection::const_iterator k=hoHits.begin();
	   k!=hoHits.end(); ++k)
	{
	  processEvent_HORecHit(k);
	}
      
      for (HFRecHitCollection::const_iterator j=hfHits.begin();
	   j!=hfHits.end(); ++j)
	{
	  processEvent_HFRecHit(j);
	}
      /*
	for (ZDCRecHitCollection::const_iterator j=zdcHits.begin();
	j!=zdcHits.end(); ++j)
	{
	processEvent_ZDCRecHit(j);
	}
      */
      if (showTiming)
	{
	  cpu_timer.stop();  std::cout <<"TIMER:: HcalDeadCellMonitor PROCESSEVENT_DIGI -> "<<cpu_timer.cpuTime()<<std::endl;
	}
    }
  
  // Fill problem cells every [checkNevents_]

  int scalefactor=deadmon_checkNevents_/deadmon_neverpresent_prescale_;
  if (ievt_>0 && ievt_%scalefactor==0)
    {
      if (deadmon_test_neverpresent_) fillNevents_neverpresent();
      if (deadmon_test_occupancy_) fillNevents_occupancy();
      if (deadmon_test_energy_) fillNevents_energy();
      fillNevents_problemCells();
      if (ievt_%deadmon_checkNevents_==0) zeroCounters(); // reset for next round of checks
    }

  return;
} // void HcalDeadCellMonitor::processEvent(...)

/* --------------------------------------- */

void HcalDeadCellMonitor::fillDeadHistosAtEndRun()
{
  // Fill histograms one last time at endRun call
  
  /*
    I'm not sure I like this feature.  Suppose checkNevents=500, and the end run occurs at 501?
    Then the occupancy plot would create errors for whichever digis were not found in a single event.
    That's not desired behavior.
    We could just exclude the occupancy test from running here, but I'm not sure that's the best solution either.
    For now (28 Oct. 2008), just disable this functionality.  We'll come back to it if necessary.
  */

  return;

  /*
  if (deadmon_test_occupancy_ && ievt_%deadmon_checkNevents_>0) fillNevents_occupancy();
  if (deadmon_test_pedestal_  && ievt_%deadmon_checkNevents_ >0) fillNevents_pedestal();
  if (deadmon_test_neighbor_  && ievt_%deadmon_checkNevents_ >0) fillNevents_neighbor();
  if ((deadmon_test_energy_ || deadmon_test_rechit_occupancy_)    && ievt_%deadmon_checkNevents_   >0) fillNevents_energy();
  if (deadmon_test_occupancy_ || deadmon_test_pedestal_ || 
      deadmon_test_neighbor_  || deadmon_test_energy_)  
   {
     fillNevents_problemCells();
     FillUnphysicalHEHFBins(ProblemDeadCellsByDepth);
   }
  */
} // fillDeadHistosAtEndOfRun()



/* --------------------------------------- */

// Digi-based dead cell checks

void HcalDeadCellMonitor::processEvent_HBHEdigi(HBHEDigiCollection::const_iterator j)
{
  // Simply check whether a digi is present.  If so, increment occupancy counter.

  int ieta=0;
  int iphi=0;
  int depth=0;
  const HBHEDataFrame digi = (const HBHEDataFrame)(*j);
  ieta=digi.id().ieta();
  iphi=digi.id().iphi();
  depth=digi.id().depth();
  if (!digi.id().validDetId(digi.id().subdet(),ieta,iphi,depth)) return;
  
  if ((HcalSubdetector)(digi.id().subdet())==HcalBarrel)
    {
      HBpresent_=true;
      if (!checkHB_)
	return;
    }
  else 
    {
      HEpresent_=true;
      if (!checkHE_)
	return;
      if (depth<3) depth+=4; // shift HE depths 1&2 up by 4
    }
  
  ++occupancy[static_cast<int>(ieta+(etaBins_-2)/2)][iphi-1][depth-1];
  present[static_cast<int>(ieta+(etaBins_-2)/2)][iphi-1][depth-1]=true;
  return;
} //void HcalDeadCellMonitor::processEvent_HBHEdigi(HBHEDigiCollection::const_iterator j)


void HcalDeadCellMonitor::processEvent_HOdigi(HODigiCollection::const_iterator j)
{
  if (!checkHO_) return;
  int ieta=0;
  int iphi=0;
  int depth=0;
  const HODataFrame digi = (const HODataFrame)(*j);
  ieta=digi.id().ieta();
  iphi=digi.id().iphi();
  depth=digi.id().depth();
  if (!digi.id().validDetId(digi.id().subdet(),ieta,iphi,depth)) return;
  ++occupancy[static_cast<int>(ieta+(etaBins_-2)/2)][iphi-1][depth-1];
  present[static_cast<int>(ieta+(etaBins_-2)/2)][iphi-1][depth-1]=true;
  return;
} //void HcalDeadCellMonitor::processEvent_HOdigi(HODigiCollection::const_iterator j)


void HcalDeadCellMonitor::processEvent_HFdigi(HFDigiCollection::const_iterator j)
{
  if (!checkHF_) return;
  int ieta=0;
  int iphi=0;
  int depth=0;
  const HFDataFrame digi = (const HFDataFrame)(*j);
  ieta=digi.id().ieta();
  iphi=digi.id().iphi();
  depth=digi.id().depth();
  if (!digi.id().validDetId(digi.id().subdet(),ieta,iphi,depth)) return;
  ++occupancy[static_cast<int>(ieta+(etaBins_-2)/2)][iphi-1][depth-1];
  present[static_cast<int>(ieta+(etaBins_-2)/2)][iphi-1][depth-1]=true;
  return;
} //void HcalDeadCellMonitor::processEvent_HFdigi(HFDigiCollection::const_iterator j)

void HcalDeadCellMonitor::processEvent_ZDCdigi(ZDCDigiCollection::const_iterator j)
{
  if (!checkZDC_) return;
  return;
  // need to set up mapping of ZDC into eta/phi/depth space before counting bad entries
} //void HcalDeadCellMonitor::processEvent_ZDCdigi(ZDCDigiCollection::const_iterator j)



//RecHit-based daed cell checks

void HcalDeadCellMonitor::processEvent_HBHERecHit(HBHERecHitCollection::const_iterator HBHEiter)
{
  float en = HBHEiter->energy();
  HcalDetId id(HBHEiter->detid().rawId());
  int ieta = id.ieta();
  int iphi = id.iphi();
  int depth = id.depth();
  
  if (id.subdet()==HcalBarrel)
    {
      HBpresent_=true;
      if (!checkHB_) return;
      if (en>=HBenergyThreshold_)
	++aboveenergy[static_cast<int>(ieta+(etaBins_-2)/2)][iphi-1][depth-1];
    }
  else
    {
      HEpresent_=true;
      if (!checkHE_)return;
      if (depth<3) depth=depth+4; // shift HE depths 1 & 2 up by 4
      if (en>=HEenergyThreshold_)
	++aboveenergy[static_cast<int>(ieta+(etaBins_-2)/2)][iphi-1][depth-1];
    }
  return;
} //void HcalDeadCellMonitor::processEvent_HBHERecHit(HBHERecHitCollection::const_iterator HOiter)


void HcalDeadCellMonitor::processEvent_HORecHit(HORecHitCollection::const_iterator HOiter)
{
  float en = HOiter->energy();
  HcalDetId id(HOiter->detid().rawId());
  int ieta = id.ieta();
  int iphi = id.iphi();
  int depth = id.depth();
  if (en>=HOenergyThreshold_)
	++aboveenergy[static_cast<int>(ieta+(etaBins_-2)/2)][iphi-1][depth-1];
  return;
} //void HcalDeadCellMonitor::processEvent_HORecHit(HORecHitCollection::const_iterator HOiter)

void HcalDeadCellMonitor::processEvent_HFRecHit(HFRecHitCollection::const_iterator HFiter)
{
  float en = HFiter->energy();
  HcalDetId id(HFiter->detid().rawId());
  int ieta = id.ieta();
  int iphi = id.iphi();
  int depth = id.depth();
  if (en>=HFenergyThreshold_)
	++aboveenergy[static_cast<int>(ieta+(etaBins_-2)/2)][iphi-1][depth-1];
  return;
} //void HcalDeadCellMonitor::processEvent_HFRecHit(HFRecHitCollection::const_iterator HFiter)

void HcalDeadCellMonitor::processEvent_ZDCRecHit(ZDCRecHitCollection::const_iterator ZDCiter)
{
  float en = ZDCiter->energy();
  HcalDetId id(ZDCiter->detid().rawId());
  int ieta = id.ieta();
  int iphi = id.iphi();
  int depth = id.depth();
  if (en>=ZDCenergyThreshold_)
	++aboveenergy[static_cast<int>(ieta+(etaBins_-2)/2)][iphi-1][depth-1];
  return;
} //void HcalDeadCellMonitor::processEvent_ZDCRecHit(HFRecHitCollection::const_iterator ZDCiter)



// fill histograms every N events
void HcalDeadCellMonitor::fillNevents_neverpresent(void)
{
  if (!deadmon_test_neverpresent_) return; // extra protection
   if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

   if (fVerbosity>0)
     std::cout <<"<HcalDeadCellMonitor::fillNevents_neverpresent> FILLING OCCUPANCY PLOTS"<<std::endl;
   int mydepth=0;
   int ieta=0;
   int iphi=0;
   // bin used for normalization
   for (unsigned int h=0;h<DigisNeverPresentByDepth.size();++h)
     {
       if (DigisNeverPresentByDepth[h]) DigisNeverPresentByDepth[h]->setBinContent(0,0,ievt_);
     }
   for (int eta=0;eta<(etaBins_-2);++eta)
     {
       ieta=eta-int((etaBins_-2)/2);
       for (int phi=0;phi<72;++phi)
	 {
	   iphi=phi+1;
	   for (int depth=1;depth<=4;++depth) 
	     {
	       for (int subdet=1;subdet<=4;++subdet)
		{
		  if (!validDetId((HcalSubdetector)subdet, ieta, iphi, depth))
		    continue;
		  // Ignore subdetectors that weren't in run
		  if ((subdet==1 && !HBpresent_) || (subdet==2 &&!HEpresent_)||(subdet==3 &&!HOpresent_) || (subdet==4 &&!HFpresent_)) continue;
		  // ignore subdetectors we explicitly mask off 
		  if ((!checkHB_ && subdet==1) ||
		      (!checkHE_ && subdet==2) ||
		      (!checkHO_ && subdet==3) ||
		      (!checkHF_ && subdet==4)) continue;
		  mydepth=depth-1 ; // my depth starts at 0, not 1
		  if (subdet==2 && depth<3) // remember that HE depths 1 and 2 are shifter up by 4 in occupancy array
		    mydepth=mydepth+4;
		  if (present[eta][phi][mydepth]==0)
		    {
		      if (fVerbosity>0) 
			std::cout <<"DEAD CELL; NEVER PRESENT: subdet = "<<subdet<<", eta = "<<ieta<<", phi = "<<iphi<<" depth = "<<depth<<" mydepth = "<<mydepth<<std::endl;
		      
		      // no digi was found for the N events; Fill cell as bad for all N events (N = deadmon_checkNevents_/prescale);
		      if (DigisNeverPresentByDepth[mydepth]) DigisNeverPresentByDepth[mydepth]->Fill(ieta,iphi,deadmon_checkNevents_/deadmon_neverpresent_prescale_);
		    }
		  else  // digi found; this is no longer a dead cell -- erase it
		    if (DigisNeverPresentByDepth[mydepth]) DigisNeverPresentByDepth[mydepth]->setBinContent(eta+2,phi+2,0);
		} // subdet loop
	     } //depth loop
	 } //phi loop
     } // eta loop
   FillUnphysicalHEHFBins(DigisNeverPresentByDepth);
   return;
} // void HcalDeadCellMonitor::fillNevents_neverpresent(void)


void HcalDeadCellMonitor::fillNevents_occupancy(void)
{
  // Fill Histograms showing digi cells with no occupancy for the past checkNevents
  if (!deadmon_test_occupancy_) return; // extra protection here against calling histograms than don't exist
  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

  if (fVerbosity>0)
    std::cout <<"<HcalDeadCellMonitor::fillNevents_occupancy> FILLING OCCUPANCY PLOTS"<<std::endl;

  int mydepth=0;
  int ieta=0;
  int iphi=0;

  // Always fill ievt_ counter with every call, so that we know histogram is being filled
  for (unsigned int h=0;h<UnoccupiedDeadCellsByDepth.size();++h)
    {
      if (UnoccupiedDeadCellsByDepth[h]) UnoccupiedDeadCellsByDepth[h]->setBinContent(0,0,ievt_);
    }

  // Only run fills of histogram when ievt%checkNevents=0
  if ((ievt_%deadmon_checkNevents_)!=0)
    return;


  for (int eta=0;eta<(etaBins_-2);++eta)
    {
      ieta=eta-int((etaBins_-2)/2);
      for (int phi=0;phi<72;++phi)
        {
	  iphi=phi+1;
	  for (int depth=1;depth<=4;++depth) 
            {
	      for (int subdet=1;subdet<=4;++subdet)
		{
		  if (!validDetId((HcalSubdetector)subdet, ieta, iphi, depth))
		    continue;
		  // Ignore subdetectors that weren't in run
		  if ((subdet==1 && !HBpresent_) || (subdet==2 &&!HEpresent_)||(subdet==3 &&!HOpresent_) || (subdet==4 &&!HFpresent_)) continue;
		  // ignore subdetectors we explicitly mask off 
		  if ((!checkHB_ && subdet==1) ||
		      (!checkHE_ && subdet==2) ||
		      (!checkHO_ && subdet==3) ||
		      (!checkHF_ && subdet==4)) continue;
		  mydepth=depth-1 ; // my depth starts at 0, not 1
		  if (subdet==2 && depth<3) // remember that HE depths 1 and 2 are shifter up by 4 in occupancy array
		    mydepth=mydepth+4;

		  if (occupancy[eta][phi][mydepth]==0)
		    {
		      if (fVerbosity>0) 
			std::cout <<"DEAD CELL; NO OCCUPANCY: subdet = "<<subdet<<", eta = "<<ieta<<", phi = "<<iphi<<" depth = "<<depth<<" mydepth = "<<mydepth<<std::endl;

		      
		      // no digi was found for the N events; Fill cell as bad for all N events (N = deadmon_checkNevents_);
		      if (UnoccupiedDeadCellsByDepth[mydepth]) UnoccupiedDeadCellsByDepth[mydepth]->Fill(ieta,iphi,deadmon_checkNevents_);

		    }
		} // for (int subdet=1;subdet<=4;++subdet)
	    } // for (int depth=1;depth<=4;++depth)
	} // for (int phi=0;...)
    } // for (int eta=0;...)
  FillUnphysicalHEHFBins(UnoccupiedDeadCellsByDepth);
  if (showTiming)
    {
      cpu_timer.stop();  std::cout <<"TIMER:: HcalDeadCellMonitor FILLNEVENTS_OCCUPANCY -> "<<cpu_timer.cpuTime()<<std::endl;
    }

  return;

} // void HcalDeadCellMonitor::fillNevents_occupancy(void)



/* ----------------------------------- */

void HcalDeadCellMonitor::fillNevents_energy(void)
{
  // Fill Histograms showing unoccupied rechits, or rec hits with low energy

  if (!deadmon_test_energy_) return;
  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

  if (fVerbosity>0)
    std::cout <<"<HcalDeadCellMonitor::fillNevents_energy> BELOW-ENERGY-THRESHOLD PLOTS"<<std::endl;

  int mydepth=0;
  int ieta=0;
  int iphi=0;
  
  for (unsigned int h=0;h<BelowEnergyThresholdCellsByDepth.size();++h)
    BelowEnergyThresholdCellsByDepth[h]->setBinContent(0,0,ievt_);

  // Only run fills of histogram when ievt%checkNevents=0
  if ((ievt_%deadmon_checkNevents_)!=0)
    return;

  for (int eta=0;eta<(etaBins_-2);++eta)
    {
      ieta=eta-int((etaBins_-2)/2);
      for (int phi=0;phi<72;++phi)
        {
	  iphi=phi+1;
	  for (int depth=1;depth<=4;++depth) 
            {
	      for (int subdet=1;subdet<=4;++subdet)
		{
		  if (!validDetId((HcalSubdetector)subdet, ieta, iphi, depth))
		    continue;
		  // Ignore subdetectors that weren't in run
                  if ((subdet==1 && !HBpresent_) || 
		      (subdet==2 &&!HEpresent_) ||
		      (subdet==3 &&!HOpresent_) ||
		      (subdet==4 &&!HFpresent_)) continue;

		  if ((!checkHB_ && subdet==1) ||
		      (!checkHE_ && subdet==2) ||
		      (!checkHO_ && subdet==3) ||
		      (!checkHF_ && subdet==4)) continue;
		  mydepth=depth-1; // my depth index starts at 0, not 1
		  if (subdet==2 && depth<3) // remember that HE depths 1 & 2 are shifted by 4
		    mydepth=mydepth+4;
		  
		  if (aboveenergy[eta][phi][mydepth]>0) continue; // cell exceeded energy at least once, so it's not dead

		  if (fVerbosity>2) 
		    std::cout <<"DEAD CELL; BELOW ENERGY THRESHOLD = "<<subdet<<" eta = "<<ieta<<", phi = "<<iphi<<" depth = "<<depth+1<<std::endl;
		  // Cell is below energy for all 'checkNevents_' consecutive events; update histogram
		  
		  if (BelowEnergyThresholdCellsByDepth[mydepth]) BelowEnergyThresholdCellsByDepth[mydepth]->Fill(ieta,iphi,deadmon_checkNevents_);
		} // for (int subdet=1;subdet<=4;++subdet)
	    } // for (int depth=1;depth<=4;++depth)
	} // for (int phi=0;...)
    } // for (int eta=0;...)

  FillUnphysicalHEHFBins(BelowEnergyThresholdCellsByDepth);
  if (showTiming)
    {
      cpu_timer.stop();  std::cout <<"TIMER:: HcalDeadCellMonitor FILLNEVENTS_ENERGY -> "<<cpu_timer.cpuTime()<<std::endl;
    }

  return;
} // void HcalDeadCellMonitor::fillNevents_energy(void)



void HcalDeadCellMonitor::fillNevents_problemCells(void)
{
  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

  if (fVerbosity>0)
    std::cout <<"<HcalDeadCellMonitor::fillNevents_problemCells> FILLING PROBLEM CELL PLOTS"<<std::endl;

  int ieta=0;
  int iphi=0;

  double problemvalue=0;
  double sumproblemvalue=0; // summed over all depths


  // Count problem cells in each subdetector
  unsigned int deadHB=0;
  unsigned int deadHE=0;
  unsigned int deadHO=0;
  unsigned int deadHF=0;
  unsigned int deadZDC=0;
  
  unsigned int neverpresentHB=0;
  unsigned int neverpresentHE=0;
  unsigned int neverpresentHO=0;
  unsigned int neverpresentHF=0;
  unsigned int neverpresentZDC=0;

  unsigned int unoccupiedHB=0;
  unsigned int unoccupiedHE=0;
  unsigned int unoccupiedHO=0;
  unsigned int unoccupiedHF=0;
  unsigned int unoccupiedZDC=0;
  
  unsigned int belowenergyHB=0;
  unsigned int belowenergyHE=0;
  unsigned int belowenergyHO=0;
  unsigned int belowenergyHF=0;
  unsigned int belowenergyZDC=0;

  int mydepth=0;

  for (int eta=0;eta<(etaBins_-2);++eta)
    {
      ieta=eta-int((etaBins_-2)/2);
      for (int phi=0;phi<72;++phi)
        {
	  iphi=phi+1;
	  for (int depth=1;depth<=4;++depth) 
            {
	      for (int subdet=1;subdet<=4;++subdet)
		{
		  if (!validDetId((HcalSubdetector)subdet, ieta, iphi, depth))
		    continue;
		  // Ignore subdetectors that weren't in run
                  if ((subdet==1 && !HBpresent_) || 
		      (subdet==2 &&!HEpresent_) ||
		      (subdet==3 &&!HOpresent_) ||
		      (subdet==4 &&!HFpresent_)) continue;

		  if ((!checkHB_ && subdet==1) ||
		      (!checkHE_ && subdet==2) ||
		      (!checkHO_ && subdet==3) ||
		      (!checkHF_ && subdet==4)) continue;
		  mydepth=depth-1; // my depth index starts at 0, not 1
		  if (subdet==2 && depth<3) // remember that HE depths 1 & 2 are shifted by 4
		    mydepth=mydepth+4;
		  // now check which dead cell tests failed; increment counter if any failed
		  if ((deadmon_test_neverpresent_ && present[eta][phi][mydepth]==0) ||
		      (deadmon_test_occupancy_ && occupancy[eta][phi][mydepth]==0 && (ievt_%deadmon_checkNevents_)==0) ||
		      (deadmon_test_energy_ && aboveenergy[eta][phi][mydepth]==0  && (ievt_%deadmon_checkNevents_)==0))
		    {
		      if (subdet==1) ++deadHB;
		      else if (subdet==2) ++deadHE;
		      else if (subdet==3) ++deadHO;
		      else if (subdet==4) ++deadHF;
		      // handle ZDC elsewhere -- in its own loop?
		    }
		  if ((deadmon_test_neverpresent_ && present[eta][phi][mydepth]==0))
		    {
		      if (subdet==1) ++neverpresentHB;
		      else if (subdet==2) ++neverpresentHE;
		      else if (subdet==3) ++neverpresentHO;
		      else if (subdet==4) ++neverpresentHF;
		      // handle ZDC elsewhere -- in its own loop?
		    }
		  if ((deadmon_test_occupancy_ && occupancy[eta][phi][mydepth]==0 && (ievt_%deadmon_checkNevents_)==0))
		    {
		      if (subdet==1) ++unoccupiedHB;
		      else if (subdet==2) ++unoccupiedHE;
		      else if (subdet==3) ++unoccupiedHO;
		      else if (subdet==4) ++unoccupiedHF;
		      // handle ZDC elsewhere -- in its own loop?
		    }
		  if ((deadmon_test_energy_ & aboveenergy[eta][phi][mydepth]==0 && (ievt_%deadmon_checkNevents_)==0))
		    {
		      if (subdet==1) ++belowenergyHB;
		      else if (subdet==2) ++belowenergyHE;
		      else if (subdet==3) ++belowenergyHO;
		      else if (subdet==4) ++belowenergyHF;
		      // handle ZDC elsewhere -- in its own loop?
		    }
		} // subdet loop
	    } //depth loop
	} // phi loop
    } //eta loop
  


 // Fill with number of problem cells found on this pass
  if (ievt_%deadmon_checkNevents_==0)
    {
      NumberOfDeadCellsHB->Fill(deadHB,deadmon_checkNevents_);
      NumberOfDeadCellsHE->Fill(deadHE,deadmon_checkNevents_);
      NumberOfDeadCellsHO->Fill(deadHO,deadmon_checkNevents_);
      NumberOfDeadCellsHF->Fill(deadHF,deadmon_checkNevents_);
      NumberOfDeadCellsZDC->Fill(deadZDC,deadmon_checkNevents_);
      NumberOfDeadCells->Fill(deadHB+deadHE+deadHO+deadHF+deadZDC,deadmon_checkNevents_);

      NumberOfUnoccupiedCellsHE->Fill(unoccupiedHE,deadmon_checkNevents_);
      NumberOfUnoccupiedCellsHO->Fill(unoccupiedHO,deadmon_checkNevents_);
      NumberOfUnoccupiedCellsHF->Fill(unoccupiedHF,deadmon_checkNevents_);
      NumberOfUnoccupiedCellsZDC->Fill(unoccupiedZDC,deadmon_checkNevents_);
      NumberOfUnoccupiedCells->Fill(unoccupiedHB+unoccupiedHE+unoccupiedHO+unoccupiedHF+unoccupiedZDC,deadmon_checkNevents_);
      
      NumberOfBelowEnergyCellsHB->Fill(belowenergyHB,deadmon_checkNevents_);
      NumberOfBelowEnergyCellsHE->Fill(belowenergyHE,deadmon_checkNevents_);
      NumberOfBelowEnergyCellsHO->Fill(belowenergyHO,deadmon_checkNevents_);
      NumberOfBelowEnergyCellsHF->Fill(belowenergyHF,deadmon_checkNevents_);
      NumberOfBelowEnergyCellsZDC->Fill(belowenergyZDC,deadmon_checkNevents_);
      NumberOfBelowEnergyCells->Fill(belowenergyHB+belowenergyHE+belowenergyHO+belowenergyHF+belowenergyZDC,deadmon_checkNevents_);
    }

  // Neverpresent cell algorithm gets called more often; fill with smaller value
  NumberOfNeverPresentCellsHB->Fill(neverpresentHB,deadmon_checkNevents_/deadmon_neverpresent_prescale_);
  NumberOfNeverPresentCellsHE->Fill(neverpresentHE,deadmon_checkNevents_/deadmon_neverpresent_prescale_);
  NumberOfNeverPresentCellsHO->Fill(neverpresentHO,deadmon_checkNevents_/deadmon_neverpresent_prescale_);
  NumberOfNeverPresentCellsHF->Fill(neverpresentHF,deadmon_checkNevents_/deadmon_neverpresent_prescale_);
  NumberOfNeverPresentCellsZDC->Fill(neverpresentZDC,deadmon_checkNevents_/deadmon_neverpresent_prescale_);
  NumberOfNeverPresentCells->Fill(neverpresentHB+neverpresentHE+neverpresentHO+neverpresentHF+neverpresentZDC,deadmon_checkNevents_/deadmon_neverpresent_prescale_);

  for (int eta=0;eta<(etaBins_-2);++eta)
    {
      ieta=eta-int((etaBins_-2)/2);
      for (int phi=0;phi<72;++phi)
        {
	  iphi=phi+1;
	  sumproblemvalue=0;
	  for (int mydepth=0;mydepth<6;++mydepth)
	    {
	      // problem value is sum of problems over all tests for a given depth
	      problemvalue=0;
	      if (deadmon_test_neverpresent_)
		{
		  problemvalue+=DigisNeverPresentByDepth[mydepth]->getBinContent(eta+2,phi+2);
		}
	      if (deadmon_test_occupancy_)
		{
		  problemvalue+=UnoccupiedDeadCellsByDepth[mydepth]->getBinContent(eta+2,phi+2);
		}
	      if (deadmon_test_energy_)
		{
		  problemvalue+=BelowEnergyThresholdCellsByDepth[mydepth]->getBinContent(eta+2,phi+2);
		}
	      problemvalue = min((double)ievt_,problemvalue);
	      sumproblemvalue+=problemvalue;
	      ProblemDeadCellsByDepth[mydepth]->setBinContent(eta+2,phi+2,problemvalue);
	      ProblemDeadCellsByDepth[mydepth]->setBinContent(0,0,ievt_);
	      // verified that maximum = ievt_ here
	    } // for (int mydepth=0;mydepth<6;...)
	  sumproblemvalue = min((double)ievt_,sumproblemvalue);
	  ProblemDeadCells->setBinContent(eta+2,phi+2,sumproblemvalue);
	  ProblemDeadCells->setBinContent(0,0,ievt_);
	} // loop on phi=0;phi<72
    } // loop on eta=0; eta<(etaBins_-2)

  if (showTiming)
    {
      cpu_timer.stop();  std::cout <<"TIMER:: HcalDeadCellMonitor FILLNEVENTS_PROBLEMCELLS -> "<<cpu_timer.cpuTime()<<std::endl;
    }
  return;
} // void HcalDeadCellMonitor::fillNevents_problemCells(void)


void HcalDeadCellMonitor::zeroCounters(bool resetpresent)
{

  // zero all counters

  // 2D histogram counters
  for (unsigned int i=0;i<ETABINS;++i)
    {
      for (unsigned int j=0;j<PHIBINS;++j)
	{
	  for (unsigned int k=0;k<6;++k)
	    {
	      if (resetpresent) present[i][j][k]=false; // keeps track of whether digi was ever present
	      occupancy[i][j][k]=0; // counts occupancy in last (checkNevents) events
	      aboveenergy[i][j][k]=0; // counts instances of cell above threshold energy in last (checkNevents)
	    }
	}
    }

  return;
} // void HcalDeadCellMonitor::zeroCounters(bool resetpresent)
