#include "DQM/HcalMonitorTasks/interface/HcalHotCellMonitor.h"

#define OUT if(fverbosity_)cout
#define BITSHIFT 6

using namespace std;

HcalHotCellMonitor::HcalHotCellMonitor()
{
  ievt_=0;
} //constructor

HcalHotCellMonitor::~HcalHotCellMonitor()
{
} //destructor


/* ------------------------------------ */ 

void HcalHotCellMonitor::setup(const edm::ParameterSet& ps,
				DQMStore* dbe)
{
  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }
  
  if (fVerbosity>0)
    cout <<"<HcalHotCellMonitor::setup>  Setting up histograms"<<endl;

  HcalBaseMonitor::setup(ps,dbe);
  baseFolder_ = rootFolder_+"HotCellMonitor_Hcal";
  
  // Hot Cell Monitor - specific cfg variables

  if (fVerbosity>1)
    cout <<"<HcalHotCellMonitor::setup>  Getting variable values from cfg files"<<endl;
  // determine whether database pedestals are in FC or ADC
  doFCpeds_ = ps.getUntrackedParameter<bool>("HotCellMonitor_pedestalsInFC", true);

  // hotmon_makeDiagnostics_ will take on base task value unless otherwise specified
  hotmon_makeDiagnostics_ = ps.getUntrackedParameter<bool>("HotCellMonitor_makeDiagnosticPlots",makeDiagnostics);
  
  // Set checkNevents values
  hotmon_checkNevents_ = ps.getUntrackedParameter<int>("HotCellMonitor_checkNevents",checkNevents_);
  hotmon_checkNevents_persistent_ = ps.getUntrackedParameter<int>("HotCellMonitor_checkNevents_persistent",hotmon_checkNevents_);
  hotmon_checkNevents_pedestal_   = ps.getUntrackedParameter<int>("HotCellMonitor_checkNevents_pedestal" ,hotmon_checkNevents_);
  hotmon_checkNevents_neighbor_   = ps.getUntrackedParameter<int>("HotCellMonitor_checkNevents_neighbor" ,hotmon_checkNevents_);
  hotmon_checkNevents_energy_     = ps.getUntrackedParameter<int>("HotCellMonitor_checkNevents_energy"   ,hotmon_checkNevents_);
 
  // Set which hot cell checks will be performed
  hotmon_test_persistent_         = ps.getUntrackedParameter<bool>("HotCellMonitor_test_persistent",true);
  hotmon_test_pedestal_           = ps.getUntrackedParameter<bool>("HotCellMonitor_test_pedestal",true);
  hotmon_test_neighbor_           = ps.getUntrackedParameter<bool>("HotCellMonitor_test_neighbor",true);
  hotmon_test_energy_             = ps.getUntrackedParameter<bool>("HotCellMonitor_test_energy",true);

  hotmon_minErrorFlag_ = ps.getUntrackedParameter<double>("HotCellMonitor_minErrorFlag",0.0);

  // pedestal test -- cell must be above pedestal+nsigma for a number of consecutive events to be considered hot
  nsigma_       = ps.getUntrackedParameter<double>("HotCellMonitor_pedestal_Nsigma",       -10);
  HBnsigma_     = ps.getUntrackedParameter<double>("HotCellMonitor_pedestal_HB_Nsigma",nsigma_);
  HEnsigma_     = ps.getUntrackedParameter<double>("HotCellMonitor_pedestal_HE_Nsigma",nsigma_);
  HOnsigma_     = ps.getUntrackedParameter<double>("HotCellMonitor_pedestal_HO_Nsigma",nsigma_);
  HFnsigma_     = ps.getUntrackedParameter<double>("HotCellMonitor_pedestal_HF_Nsigma",nsigma_);
  ZDCnsigma_    = ps.getUntrackedParameter<double>("HotCellMonitor_pedestal_ZDC_Nsigma", nsigma_);

  // rechit energy test -- cell must be above threshold value for a number of consecutive events to be considered hot
  energyThreshold_       = ps.getUntrackedParameter<double>("HotCellMonitor_energyThreshold",   1);
  HBenergyThreshold_     = ps.getUntrackedParameter<double>("HotCellMonitor_HB_energyThreshold",energyThreshold_);
  HEenergyThreshold_     = ps.getUntrackedParameter<double>("HotCellMonitor_HE_energyThreshold",energyThreshold_);
  HOenergyThreshold_     = ps.getUntrackedParameter<double>("HotCellMonitor_HO_energyThreshold",energyThreshold_);
  HFenergyThreshold_     = ps.getUntrackedParameter<double>("HotCellMonitor_HF_energyThreshold",energyThreshold_);
  ZDCenergyThreshold_    = ps.getUntrackedParameter<double>("HotCellMonitor_HF_energyThreshold",-999);

  // rechit event-by-event energy test -- cell must be above threshold to be considered hot
  persistentThreshold_       = ps.getUntrackedParameter<double>("HotCellMonitor_persistentThreshold",   1);
  HBpersistentThreshold_     = ps.getUntrackedParameter<double>("HotCellMonitor_HB_persistentThreshold",persistentThreshold_);
  HEpersistentThreshold_     = ps.getUntrackedParameter<double>("HotCellMonitor_HE_persistentThreshold",persistentThreshold_);
  HOpersistentThreshold_     = ps.getUntrackedParameter<double>("HotCellMonitor_HO_persistentThreshold",persistentThreshold_);
  HFpersistentThreshold_     = ps.getUntrackedParameter<double>("HotCellMonitor_HF_persistentThreshold",persistentThreshold_);
  ZDCpersistentThreshold_    = ps.getUntrackedParameter<double>("HotCellMonitor_HF_persistentThreshold",-999);

  // neighboring-cell tests
  defaultNeighborParams_.DeltaIphi = ps.getUntrackedParameter<int>("HotCellMonitor_neighbor_deltaIphi", 1);
  defaultNeighborParams_.DeltaIeta = ps.getUntrackedParameter<int>("HotCellMonitor_neighbor_deltaIeta", 1);
  defaultNeighborParams_.DeltaDepth = ps.getUntrackedParameter<int>("HotCellMonitor_neighbor_deltaDepth", 0);
  defaultNeighborParams_.minCellEnergy = ps.getUntrackedParameter<double>("HotCellMonitor_neighbor_minCellEnergy",3.);
  defaultNeighborParams_.minNeighborEnergy = ps.getUntrackedParameter<double>("HotCellMonitor_neighbor_minNeighborEnergy",0.);
  defaultNeighborParams_.maxEnergy = ps.getUntrackedParameter<double>("HotCellMonitor_neighbor_maxEnergy",50);
  defaultNeighborParams_.HotEnergyFrac = ps.getUntrackedParameter<double>("HotCellMonitor_neighbor_HotEnergyFrac",0.01);

  setupNeighborParams(ps,HBNeighborParams_ ,"HB");
  setupNeighborParams(ps,HENeighborParams_ ,"HE");
  setupNeighborParams(ps,HONeighborParams_ ,"HO");
  setupNeighborParams(ps,HFNeighborParams_ ,"HF");
  setupNeighborParams(ps,ZDCNeighborParams_,"ZDC");
  HFNeighborParams_.DeltaIphi*=2; // HF cell segmentation is 10 degrees, not 5 (mostly).  Need to multiply by 2 to convert from cell range to degree format

  // Set initial event # to 0
  ievt_=0;

  // zero all counters
  for (int i=0;i<ETABINS;++i)
    {
      for (int j=0;j<PHIBINS;++j)
        {
          for (int k=0;k<4;++k)
            {
              abovepersistent[i][j][k]=0;
	      abovepedestal[i][j][k]=0;
              aboveneighbors[i][j][k]=0;
              aboveenergy[i][j][k]=0;
	      rechit_occupancy_sum[i][j][k]=0;
	      rechit_energy_sum[i][j][k]=0.;
            }
        }
    }

  // Set up histograms
  if (m_dbe)
    {
      if (fVerbosity>1)
	cout <<"<HcalHotCellMonitor::setup>  Setting up histograms"<<endl;

      m_dbe->setCurrentFolder(baseFolder_);
      meEVT_ = m_dbe->bookInt("Hot Cell Task Event Number");
      meEVT_->Fill(ievt_);

      // Create problem cell plots
      // Overall plot gets an initial " " in its name
      ProblemHotCells=m_dbe->book2D(" ProblemHotCells",
                                     " Problem Hot Cell Rate for all HCAL",
                                     etaBins_,etaMin_,etaMax_,
                                     phiBins_,phiMin_,phiMax_);
      ProblemHotCells->setAxisTitle("i#eta",1);
      ProblemHotCells->setAxisTitle("i#phi",2);
      
      (ProblemHotCells->getTH2F())->SetMinimum(hotmon_minErrorFlag_);
      (ProblemHotCells->getTH2F())->SetMaximum(1.);
      
      // Overall Problem plot appears in main directory; plots by depth appear \in subdirectory
      m_dbe->setCurrentFolder(baseFolder_+"/problem_hotcells");
      setupDepthHists2D(ProblemHotCellsByDepth, " Problem Hot Cell Rate","");
      
      setMinMaxHists2D(ProblemHotCellsByDepth,hotmon_minErrorFlag_,1.); // set minimum to hotmon_minErrorFlag_?

      // Set up plots for each failure mode of hot cells
      stringstream units; // We'll need to set the titles individually, rather than passing units to setupDepthHists2D (since this also would affect the name of the histograms)
      
      m_dbe->setCurrentFolder(baseFolder_+"/hot_rechit_above_threshold");
      setupDepthHists2D(AboveEnergyThresholdCellsByDepth,
			"Hot Cells Above Energy Threshold","");
      setMinMaxHists2D(AboveEnergyThresholdCellsByDepth,0.,1.);

      // set more descriptive titles for plots
      units.str("");
      units<<"Hot Cells: Depth 1 -- HB > "<<HBenergyThreshold_<<" GeV, HF > "<<HFenergyThreshold_<<" GeV";
      AboveEnergyThresholdCellsByDepth[0]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells: Depth 2 -- HB > "<<HBenergyThreshold_<<" GeV, HF > "<<HFenergyThreshold_<<" GeV";
      AboveEnergyThresholdCellsByDepth[1]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells: Depth 3 -- HE > "<<HEenergyThreshold_<<" GeV";
      AboveEnergyThresholdCellsByDepth[2]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells: HO > "<<HOenergyThreshold_<<" GeV, ZDC TBD";
      AboveEnergyThresholdCellsByDepth[3]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells: Depth 1 -- HE > "<<HEenergyThreshold_<<" GeV";
      AboveEnergyThresholdCellsByDepth[4]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells: Depth 2 -- HE > "<<HEenergyThreshold_<<" GeV";
      AboveEnergyThresholdCellsByDepth[5]->setTitle(units.str().c_str());

      m_dbe->setCurrentFolder(baseFolder_+"/hot_rechit_always_above_threshold");
      setupDepthHists2D(AbovePersistentThresholdCellsByDepth,
			"Hot Cells Persistently Above Energy Threshold","");
      setMinMaxHists2D(AbovePersistentThresholdCellsByDepth,0.,1.);

      // set more descriptive titles for plots
      units.str("");
      units<<"Hot Cells: Depth 1 -- HB > "<<HBpersistentThreshold_<<" GeV, HF > "<<HFpersistentThreshold_<<" GeV for "<< hotmon_checkNevents_persistent_<<" consec. events";
      AbovePersistentThresholdCellsByDepth[0]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells: Depth 2 -- HB > "<<HBpersistentThreshold_<<" GeV, HF > "<<HFpersistentThreshold_<<" GeV for "<<hotmon_checkNevents_persistent_<<" consec. events";
      AbovePersistentThresholdCellsByDepth[1]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells: Depth 3 -- HE > "<<HEpersistentThreshold_<<" GeV for "<<hotmon_checkNevents_persistent_<<" consec. events";
      AbovePersistentThresholdCellsByDepth[2]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells:  HO > "<<HOpersistentThreshold_<<" GeV, ZDC TBD for "<<hotmon_checkNevents_persistent_<<" consec. events";
      AbovePersistentThresholdCellsByDepth[3]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells: Depth 1 -- HE > "<<HEpersistentThreshold_<<" GeV for "<<hotmon_checkNevents_persistent_<<" consec. events";
      AbovePersistentThresholdCellsByDepth[4]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells: Depth 2 -- HE > "<<HEpersistentThreshold_<<" GeV for "<<hotmon_checkNevents_persistent_<<" consec. events";
      AbovePersistentThresholdCellsByDepth[5]->setTitle(units.str().c_str());
      

      m_dbe->setCurrentFolder(baseFolder_+"/hot_pedestaltest");
      setupDepthHists2D(AbovePedestalHotCellsByDepth,"Hot Cells Above Pedestal","");
      setMinMaxHists2D(AbovePedestalHotCellsByDepth,0.,1.);

      // set more descriptive titles for pedestal plots
      units.str("");
      units<<"Hot Cells Above Pedestal Depth 1 -- HB > ped + "<<HBnsigma_<<" #sigma, HF > ped + "<<HFnsigma_<<" #sigma";
      AbovePedestalHotCellsByDepth[0]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells Above Pedestal Depth 2 -- HB > ped + "<<HBnsigma_<<" #sigma, HF > ped + "<<HFnsigma_<<" #sigma";
      AbovePedestalHotCellsByDepth[1]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells Above Pedestal Depth 3 -- HE > ped + "<<HEnsigma_<<" #sigma";
      AbovePedestalHotCellsByDepth[2]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells Above Pedestal Depth 4 -- HO > ped + "<<HOnsigma_<<" #sigma, ZDC TBD";
      AbovePedestalHotCellsByDepth[3]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells Above Pedestal Depth 1 -- HE > ped + "<<HEnsigma_<<" #sigma";
      AbovePedestalHotCellsByDepth[4]->setTitle(units.str().c_str());
      units.str("");
      units<<"Hot Cells Above Pedestal Depth 2 -- HE > ped + "<<HEnsigma_<<" #sigma";
      AbovePedestalHotCellsByDepth[5]->setTitle(units.str().c_str());
      units.str("");

      m_dbe->setCurrentFolder(baseFolder_+"/hot_neighbortest");
      setupDepthHists2D(AboveNeighborsHotCellsByDepth,"Hot Cells Failing Neighbor Test","");
      setMinMaxHists2D(AboveNeighborsHotCellsByDepth,0.,1.);

      if (hotmon_makeDiagnostics_)
	{
	  if (hotmon_test_pedestal_)
	    {
	      m_dbe->setCurrentFolder(baseFolder_+"/diagnostics/pedestal");
	      d_HBnormped=m_dbe->book1D("HB_normped","HB Hot Cell pedestal diagnostic ",300,-10,20);
	      d_HEnormped=m_dbe->book1D("HE_normped","HE Hot Cell pedestal diagnostic",300,-10,20);
	      d_HOnormped=m_dbe->book1D("HO_normped","HO Hot Cell pedestal diagnostic",300,-10,20);
	      d_HFnormped=m_dbe->book1D("HF_normped","HF Hot Cell pedestal diagnostic",300,-10,20);
	      d_HBnormped->setAxisTitle("(avg ADC-pedestal)/#sigma",1);
	      d_HEnormped->setAxisTitle("(avg ADC-pedestal)/#sigma",1);
	      d_HOnormped->setAxisTitle("(avg ADC-pedestal)/#sigma",1);
	      d_HFnormped->setAxisTitle("(avg ADC-pedestal)/#sigma",1);
	    }
	  if (hotmon_test_energy_ || hotmon_test_persistent_)
	    {
	      m_dbe->setCurrentFolder(baseFolder_+"/diagnostics/rechitenergy");
	      d_HBrechitenergy=m_dbe->book1D("HB_rechitenergy","HB rechit energy",1500,-10,140);
	      d_HErechitenergy=m_dbe->book1D("HE_rechitenergy","HE rechit energy",1500,-10,140);
	      d_HOrechitenergy=m_dbe->book1D("HO_rechitenergy","HO rechit energy",1500,-10,140);
	      d_HFrechitenergy=m_dbe->book1D("HF_rechitenergy","HF rechit energy",1500,-10,140);
	      setupDepthHists2D(d_avgrechitenergymap,
				"Average rec hit energy per cell","");
	    }
	  if (hotmon_test_neighbor_)
	    {
	      m_dbe->setCurrentFolder(baseFolder_+"/diagnostics/neighborcells");
	      d_HBenergyVsNeighbor=m_dbe->book2D("HB_energyVsNeighbor","HB  #Sigma Neighbors vs. rec hit energy",100,-5,15,100,0,25);
	      d_HEenergyVsNeighbor=m_dbe->book2D("HE_energyVsNeighbor","HE  #Sigma Neighbors vs. rec hit energy",100,-5,15,100,0,25);
	      d_HOenergyVsNeighbor=m_dbe->book2D("HO_energyVsNeighbor","HO  #Sigma Neighbors vs. rec hit energy",100,-5,15,100,0,25);
	      d_HFenergyVsNeighbor=m_dbe->book2D("HF_energyVsNeighbor","HF  #Sigma Neighbors vs. rec hit energy",100,-5,15,100,0,25);
	    }
	} // if (hotmon_makeDiagnostics_)
    } // if (m_dbe)

  return;
} //void HcalHotCellMonitor::setup(...)

/* --------------------------- */
void HcalHotCellMonitor::setupNeighborParams(const edm::ParameterSet& ps,
					      hotNeighborParams& N,
					      char* type)
{
  // sets up parameters for neighboring-cell algorithm for each subdetector
  ostringstream myname;
  myname<<"HotCellMonitor_"<<type<<"_neighbor_deltaIphi";
  N.DeltaIphi = ps.getUntrackedParameter<int>(myname.str().c_str(),
					      defaultNeighborParams_.DeltaIphi);
  myname.str("");
  myname<<"HotCellMonitor_"<<type<<"_neighbor_deltaIeta";
  N.DeltaIeta = ps.getUntrackedParameter<int>(myname.str().c_str(),
					      defaultNeighborParams_.DeltaIeta);
  myname.str("");
  myname<<"HotCellMonitor_"<<type<<"_neighbor_deltaDepth";
  N.DeltaDepth = ps.getUntrackedParameter<int>(myname.str().c_str(),
					       defaultNeighborParams_.DeltaDepth);
  myname.str("");
  myname<<"HotCellMonitor_"<<type<<"_neighbor_minCellEnergy";
  N.minCellEnergy = ps.getUntrackedParameter<double>(myname.str().c_str(),
						     defaultNeighborParams_.minCellEnergy);
  myname.str("");
  myname<<"HotCellMonitor_"<<type<<"_neighbor_minNeighborEnergy";
  N.minNeighborEnergy = ps.getUntrackedParameter<double>(myname.str().c_str(),
							 defaultNeighborParams_.minNeighborEnergy);
  myname.str("");
  myname<<"HotCellMonitor_"<<type<<"_neighbor_maxEnergy";
  N.maxEnergy = ps.getUntrackedParameter<double>(myname.str().c_str(),
						 defaultNeighborParams_.maxEnergy);
  myname.str("");
  myname<<"HotCellMonitor_"<<type<<"_HotEnergyFrac";
  N.HotEnergyFrac = ps.getUntrackedParameter<double>(myname.str().c_str(),
						     defaultNeighborParams_.HotEnergyFrac);
  return;
} // void HcalHotCellMonitor::setupNeighborParams

/* --------------------------- */

void HcalHotCellMonitor::reset(){}  // reset function is empty for now

/* --------------------------- */

void HcalHotCellMonitor::createMaps(const HcalDbService& cond)
{

  // Creates maps for pedestals, widths, and pedestals+Nsigma*widths, using HcalDetIds as keys
  
  if (!hotmon_test_pedestal_) return; // no need to create maps if we're not running the pedestal-based hot cell finder

  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }
  
  if (fVerbosity>0)
    cout <<"<HcalHotCellMonitor::createMaps>:  Making pedestal maps"<<endl;
  float ped=0;
  float width=0;
  HcalCalibrations calibs;
  const HcalQIEShape* shape = cond.getHcalShape();

  double myNsigma=0;

  for (int ieta=(int)etaMin_;ieta<=(int)etaMax_;++ieta)
    {
      for (int iphi=(int)phiMin_;iphi<=(int)phiMax_;++iphi)
	{
	  for (int depth=1;depth<=4;++depth)
	    {
	      for (int subdet=1;subdet<=4;++subdet)
		{
		  if (!validDetId((HcalSubdetector)subdet, ieta, iphi, depth))
		    continue;
		  HcalDetId hcal((HcalSubdetector)(subdet), ieta, iphi, depth);
		  
		  if (hcal.subdet()==HcalBarrel)
		    myNsigma=HBnsigma_;
		  else if (hcal.subdet()==HcalEndcap)
		    myNsigma=HEnsigma_;
		  else if (hcal.subdet()==HcalOuter)
		    myNsigma=HOnsigma_;
		  else if (hcal.subdet()==HcalForward)
		    myNsigma=HFnsigma_;
		  
		  calibs=cond.getHcalCalibrations(hcal);
		  const HcalPedestalWidth* pedw = cond.getPedestalWidth(hcal);
		   
		  ped=0.;
		  width=0.;

		  // loop over capids
		  for (int capid=0;capid<4;++capid)
		    {
		      if (doFCpeds_)
			{
			  // pedestals in fC
			  const HcalQIECoder* channelCoder=cond.getHcalCoder(hcal);

			  // Convert pedestals to ADC
			  ped+=channelCoder->adc(*shape,
						 (float)calibs.pedestal(capid),
						 capid);

			  // Okay, this definitely isn't right.  Need to figure out how to convert from fC to ADC properly
			  // Right now, take width as half the difference between (ped+width)- (ped-width), converting each to ADC

			  width+=0.5*(channelCoder->adc(*shape,
							(float)calibs.pedestal(capid)+(float)pow((double)pedw->getWidth(capid),(double)0.5),
							capid)
				      - channelCoder->adc(*shape,
							  (float)calibs.pedestal(capid)-(float)pow((double)pedw->getWidth(capid),(double)0.5),
							  capid));
			} // if (doFCpeds_) // (pedestals in fC)
		      else
			{
			  // pedestals in ADC
			  ped+=calibs.pedestal(capid);
			  width+=pedw->getWidth(capid); // add in quadrature?  Make use of correlations?
			} // else //pedestals in ADC
		    } // for (int capid=0;capid<4;++capid)

		  ped/=4.;  // pedestal value is average over capids
		  if (doFCpeds_)
		    width/=4.;
		  else
		    width=pow((double)width/4.,(double)0.5); // getWidth returns width^2

		  pedestals_[hcal]=ped;
		  widths_[hcal]=width;
		  if (fVerbosity>1) cout <<"<HcalHotCellMonitor::createMaps>  Pedestal Value -- ID = "<<(HcalSubdetector)subdet<<"  ("<<ieta<<", "<<iphi<<", "<<depth<<"): "<<ped<<"; width = "<<width<<endl;
		  pedestal_thresholds_[hcal]=ped+myNsigma*width;
		} // for (int subdet=1,...)
	    } // for (int depth=1;...)
	} // for (int phi ...)
    } // for (int ieta...)
  
  return;
} // void HcalHotCellMonitor::createMaps




/* ------------------------- */

void HcalHotCellMonitor::done(std::map<HcalDetId, unsigned int>& myqual)
{
  if (dump2database==0) // don't do anything special unless specifically asked to dump db file
    return;

  // Dump to ascii file for database -- now handled in Channel Status objects
  /*
  char buffer [1024];
  ofstream fOutput("hcalHotCells.txt", ios::out);
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
      cout <<"<HcalHotCellMonitor>  Summary of Hot Cells in Run: "<<endl;
      cout <<"(Error rate must be >= "<<hotmon_minErrorFlag_*100.<<"% )"<<endl;  
    }
  for (int ieta=1;ieta<=etaBins_;++ieta)
    {
      for (int iphi=1;iphi<=phiBins_;++iphi)
        {
          eta=ieta+int(etaMin_)-1;
          phi=iphi+int(phiMin_)-1;
	  
          for (int d=0;d<6;++d)
            {
	      binval=ProblemHotCellsByDepth[d]->getBinContent(ieta,iphi);
	     
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
	      if (fVerbosity>0 && binval>hotmon_minErrorFlag_)
		cout <<"Hot Cell "<<subdet<<"("<<eta<<", "<<phi<<", "<<mydepth<<"):  "<<binval*100.<<"%"<<endl;
	      int value = 0;
	      if (binval>hotmon_minErrorFlag_)
		value=1;

	      if (myqual.find(myid)==myqual.end())
		{
		  myqual[myid]=(value<<BITSHIFT);  // hotcell shifted to bit 6
		}
	      else
		{
		  int mask=(1<<BITSHIFT);
		  if (value==1)
		    myqual[myid] |=mask;

		  else
		    myqual[myid] &=~mask;
		}
	      /*
	      sprintf(buffer, "  %15i %15i %15i %15s %8X %10X \n",eta,phi,mydepth,subdetname,int(value<<BITSHIFT),int(myid.rawId()));
	      fOutput<<buffer;
	      */
	    } // for (int d=0;d<6;++d) // loop over depth histograms
	} // for (int iphi=1;iphi<=phiBins_;++iphi)
    } // for (int ieta=1;ieta<=etaBins_;++ieta)
  //fOutput.close();

  return;

} // void HcalHotCellMonitor::done()



/* --------------------------------- */

void HcalHotCellMonitor::clearME()
{
  // I don't think this function gets cleared any more.  
  // And need to add code to clear out subfolders as well?
  if (m_dbe)
    {
      m_dbe->setCurrentFolder(baseFolder_);
      m_dbe->removeContents();
    }
  return;
} // void HcalHotCellMonitor::clearME()

/* -------------------------------- */


void HcalHotCellMonitor::processEvent(const HBHERecHitCollection& hbHits,
				       const HORecHitCollection& hoHits,
				       const HFRecHitCollection& hfHits,
				       //const ZDCRecHitCollection& zdcHits,
				       const HBHEDigiCollection& hbhedigi,
				       const HODigiCollection& hodigi,
				       const HFDigiCollection& hfdigi,
				       //const ZDCDigiCollection& zdcdigi,
				       const HcalDbService& cond
				       )
{

  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

  ++ievt_;
  if (m_dbe) meEVT_->Fill(ievt_);

  if (fVerbosity>1) cout <<"<HcalHotCellMonitor::processEvent> Processing event..."<<endl;

  // Search for cells consistently above (pedestal + N sigma)
  if (hotmon_test_pedestal_)
    processEvent_pedestal(hbhedigi,hodigi,hfdigi,cond);

  // Search for hot cells above a certain energy
  if (hotmon_test_energy_ || hotmon_test_persistent_)
    {
      processEvent_rechitenergy(hbHits, hoHits,hfHits);
    }

  // Search for cells that are hot compared to their neighbors
  if (hotmon_test_neighbor_)
    {
      processEvent_rechitneighbors(hbHits, hoHits, hfHits);
    }

  // Fill problem cells
  if (((ievt_%hotmon_checkNevents_persistent_ ==0) && hotmon_test_persistent_ )||
      ((ievt_%hotmon_checkNevents_pedestal_  ==0) && hotmon_test_pedestal_  )||
      ((ievt_%hotmon_checkNevents_neighbor_  ==0) && hotmon_test_neighbor_  )||
      ((ievt_%hotmon_checkNevents_energy_    ==0) && hotmon_test_energy_    ))
    {
      fillNevents_problemCells();
    }

  return;
} // void HcalHotCellMonitor::processEvent(...)

/* --------------------------------------- */

void HcalHotCellMonitor::fillHotHistosAtEndRun()
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

  if (hotmon_test_persistent_ && ievt_%hotmon_checkNevents_persistent_>0) fillNevents_persistentenergy();
  if (hotmon_test_pedestal_  && ievt_%hotmon_checkNevents_pedestal_ >0) fillNevents_pedestal();
  if (hotmon_test_neighbor_  && ievt_%hotmon_checkNevents_neighbor_ >0) fillNevents_neighbor();
  if (hotmon_test_energy_    && ievt_%hotmon_checkNevents_energy_   >0) fillNevents_energy();
  if (hotmon_test_persistent_ || hotmon_test_pedestal_ || 
      hotmon_test_neighbor_  || hotmon_test_energy_)  
    fillNevents_problemCells();
}

/* --------------------------------------- */


void HcalHotCellMonitor::processEvent_rechitenergy( const HBHERecHitCollection& hbheHits,
						     const HORecHitCollection& hoHits,
						     const HFRecHitCollection& hfHits)
						
{
  // Looks at rechits of cells and compares to threshold energies.
  // Cells above thresholds get marked as hot candidates

  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

 if (fVerbosity>1) cout <<"<HcalHotCellMonitor::processEvent_rechitenergy> Processing rechits..."<<endl;
 if (hotmon_test_neighbor_)   rechitEnergies_.clear();

 // loop over HBHE
 for (HBHERecHitCollection::const_iterator HBHEiter=hbheHits.begin(); HBHEiter!=hbheHits.end(); ++HBHEiter) 
   { // loop over all hits
     float en = HBHEiter->energy();
     //float ti = HBHEiter->time();

     HcalDetId id(HBHEiter->detid().rawId());
     int ieta = id.ieta();
     int iphi = id.iphi();
     int depth = id.depth();

     if (hotmon_makeDiagnostics_)
       {
	 ++rechit_occupancy_sum[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1];
	 rechit_energy_sum[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1]+=en;
       }
     if (id.subdet()==HcalBarrel)
       {
	 if (!checkHB_) continue;
	 if (hotmon_makeDiagnostics_) d_HBrechitenergy->Fill(en);
	 if (en>=HBenergyThreshold_)
	   ++aboveenergy[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1];
	 if (en>=HBpersistentThreshold_)
	   ++abovepersistent[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1];
       }
     else if (id.subdet()==HcalEndcap)
       {
	 if (!checkHE_) continue;
	 if (hotmon_makeDiagnostics_) d_HErechitenergy->Fill(en);
	 if (en>=HEenergyThreshold_)
	   ++aboveenergy[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1];
	 if (en>=HEpersistentThreshold_)
	   ++abovepersistent[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1];
       }
     if (hotmon_test_neighbor_) rechitEnergies_[id]=en;
   } //for (HBHERecHitCollection::const_iterator HBHEiter=...)

 // loop over HO
 if (checkHO_)
   {
     for (HORecHitCollection::const_iterator HOiter=hoHits.begin(); HOiter!=hoHits.end(); ++HOiter) 
       { // loop over all hits
	 float en = HOiter->energy();
	 
	 HcalDetId id(HOiter->detid().rawId());
	 int ieta = id.ieta();
	 int iphi = id.iphi();
	 int depth = id.depth();
	 
	 if (hotmon_makeDiagnostics_)
	   {
	     ++rechit_occupancy_sum[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1];
	     rechit_energy_sum[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1]+=en;
	   }
	 if (hotmon_makeDiagnostics_) d_HOrechitenergy->Fill(en);
	 if (en>=HOenergyThreshold_)
	   ++aboveenergy[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1];
	 if (en>=HOpersistentThreshold_)
	   ++abovepersistent[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1];
	 if (hotmon_test_neighbor_) rechitEnergies_[id]=en;
       }
   } // if (checkHO_)
 
 // loop over HF
 if (checkHF_)
   {
     for (HFRecHitCollection::const_iterator HFiter=hfHits.begin(); HFiter!=hfHits.end(); ++HFiter) 
       { // loop over all hits
	 float en = HFiter->energy();
	 
	 HcalDetId id(HFiter->detid().rawId());
	 int ieta = id.ieta();
	 int iphi = id.iphi();
	 int depth = id.depth();
	 
	 if (hotmon_makeDiagnostics_)
	   {
	     ++rechit_occupancy_sum[ieta+(int)((etaBins_-2)/2)][iphi-1][depth+1];
	     rechit_energy_sum[ieta+(int)((etaBins_-2)/2)][iphi-1][depth+1]+=en;
	   }
	 if (hotmon_makeDiagnostics_) d_HFrechitenergy->Fill(en);
	 if (en>=HFenergyThreshold_)
	   ++aboveenergy[ieta+(int)((etaBins_-2)/2)][iphi-1][depth+1]; // HF depths get shifted up by +2
	 if (en>=HBpersistentThreshold_)
	   ++abovepersistent[ieta+(int)((etaBins_-2)/2)][iphi-1][depth+1];
	 if (hotmon_test_neighbor_) rechitEnergies_[id]=en;
       }
   } // if (checkHF_)
 
 
 // Fill histograms 
  if (ievt_%hotmon_checkNevents_energy_==0)
    {
	if (fVerbosity) cout <<"<HcalHotCellMonitor::processEvent_digi> Filling HotCell Energy plots"<<endl;
	fillNevents_energy();
    }
  if (ievt_%hotmon_checkNevents_persistent_==0)
    {
	if (fVerbosity) cout <<"<HcalHotCellMonitor::processEvent_digi> Filling HotCell Persistent Energy plots"<<endl;
	fillNevents_persistentenergy();
    }
  
 if (showTiming)
   {
     cpu_timer.stop();  cout <<"TIMER:: HcalHotCellMonitor PROCESSEVENT_RECHITENERGY -> "<<cpu_timer.cpuTime()<<endl;
   }
 
 return;
} // void HcalHotCellMonitor::processEvent_rechitenergy

/* --------------------------------------- */


void HcalHotCellMonitor::processEvent_rechitneighbors( const HBHERecHitCollection& hbheHits,
							const HORecHitCollection& hoHits,
							const HFRecHitCollection& hfHits
							)
{
  // Compares energy to energy of neighboring cells.
  // This is a slightly simplified version of D0's NADA algorithm

  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

 if (fVerbosity>1) cout <<"<HcalHotCellMonitor::processEvent_rechitneighbors> Processing rechits..."<<endl;

 // if Energy tests weren't run, need to create map of Detid:rechitenergy here
 if (!hotmon_test_energy_ && !hotmon_test_persistent_)
   {
     rechitEnergies_.clear(); // clear old map
     for (HBHERecHitCollection::const_iterator HBHEiter=hbheHits.begin(); HBHEiter!=hbheHits.end(); ++HBHEiter) 
       { // loop over all hits
	 float en = HBHEiter->energy();
	 HcalDetId id(HBHEiter->detid().rawId());
	 if (!checkHB_ && id.subdet()==HcalBarrel)
	   continue;
	 if (!checkHE_ && id.subdet()==HcalEndcap)
	   continue;
	 rechitEnergies_[id]=en;
       }
     // HO
     if (checkHO_)
       {
	 for (HORecHitCollection::const_iterator HOiter=hoHits.begin(); HOiter!=hoHits.end(); ++HOiter) 
	   { // loop over all hits
	     float en = HOiter->energy();
	     HcalDetId id(HOiter->detid().rawId());
	     rechitEnergies_[id]=en;
	   }
       } // if (checkHO_)
     //HF
     if (checkHF_)
       {
	 for (HFRecHitCollection::const_iterator HFiter=hfHits.begin(); HFiter!=hfHits.end(); ++HFiter) 
	   { // loop over all hits
	     float en = HFiter->energy();
	     HcalDetId id(HFiter->detid().rawId());
	     rechitEnergies_[id]=en;
	   }
       } // if (checkHF_)

   } // if (!hotmon_test_energy_ && !hotmon_test_persistent_)   

 // Now do "real" loop, checking against each cell against its neighbors
 
 /* Note:  This works a little differently than the other tests.  The other tests check that a cell consistently
    fails its test condition for N consecutive events.  The neighbor test will flag a cell for every event in which
    it's significantly less than its neighbors, regardless of whether that condition persists for a number of events.
 */

 int ieta, iphi, depth;
 float en;

 int neighborsfound=0;
 float enNeighbor=0;

 // loop over HBHE
 for (HBHERecHitCollection::const_iterator HBHEiter=hbheHits.begin(); 
      HBHEiter!=hbheHits.end(); 
      ++HBHEiter) 
   { // loop over all hits
     
     en = HBHEiter->energy();
     HcalDetId id(HBHEiter->detid().rawId());
     ieta = id.ieta();
     iphi = id.iphi();
     depth = id.depth();

     if (id.subdet()==HcalBarrel)
       {
	 if (!checkHB_) continue;
	 // Case 0:  energy > max value; it's marked as hot regardless of neighbors
	 if (en>HBNeighborParams_.maxEnergy)
	   {
	     aboveneighbors[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1]++;
	     continue;
	   }

	 // Search keys for neighboring cells
	 if (en<HBNeighborParams_.minCellEnergy) // cells below minCellEnergy not considered hot
	   continue;
	 neighborsfound=0;
	 enNeighbor=0;
	 for (int nD=-1*HBNeighborParams_.DeltaDepth;nD<=HBNeighborParams_.DeltaDepth;++nD)
	   {
	     for (int nP =-1*HBNeighborParams_.DeltaIphi;nP<=HBNeighborParams_.DeltaIphi;++nP)
	       {
		 for (int nE =-1*HBNeighborParams_.DeltaIeta;nE<=HBNeighborParams_.DeltaIeta;++nE)
		   {
		     if (nD==0 && nE==0 && nP==0) 
		       continue; // don't count the cell itself
		     int myphi=nP+iphi;
		     if (myphi>72) myphi-=72; // allow for wrapping of cells
		     if (myphi<=0) myphi+=72;
		     if (!validDetId((HcalSubdetector)(1),nE+ieta, myphi, nD+depth)) continue;
		     HcalDetId myid((HcalSubdetector)(1), nE+ieta, myphi, nD+depth); // HB
		     if (rechitEnergies_.find(myid)==rechitEnergies_.end())
		       continue;
		     if (rechitEnergies_[myid]<HBNeighborParams_.minNeighborEnergy)
		       continue;
		     ++neighborsfound;
		     enNeighbor+=rechitEnergies_[myid];
		   } // loop over nE (neighbor eta)
	       } // loop over nP (neighbor phi)
	   } // loop over nD depths

	 if (hotmon_makeDiagnostics_)
	   d_HBenergyVsNeighbor->Fill(en,enNeighbor);
	 
	 // Case 1:  Not enough good neighbors found
	 if (neighborsfound==0)
	   continue;
	 // Case 2:  energy/(avg. neighbor energy) too large for cell to be considered hot
	 if ((1.*enNeighbor/en)>HBNeighborParams_.HotEnergyFrac && en>0 && enNeighbor>0)
	   continue;
	 // Case 3:  Tests passed; cell marked as hot
	 aboveneighbors[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1]++;
       }

     else if (id.subdet()==HcalEndcap)
       {
	 if (!checkHE_) continue;

	 // Case 0:  energy > max value; it's marked as hot regardless of neighbors
	 if (en>HENeighborParams_.maxEnergy)
	   {
	     aboveneighbors[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1]++;
	     continue;
	   }
	 if (en<HENeighborParams_.minCellEnergy)
	   continue; // cells below this value can never be considered hot
	 // Search keys for neighboring cells
	 neighborsfound=0;
	 enNeighbor=0;
	 int HEDeltaIphi = HENeighborParams_.DeltaIphi;
	 // now correct for boundaries
	 if (abs(ieta)>20) HEDeltaIphi*=2; // double iphi boundary range when segmentation switches to 10 degrees
	 // This still needs to be worked on to properly deal with boundaries
	 for (int nD=-1*HENeighborParams_.DeltaDepth;nD<=HENeighborParams_.DeltaDepth;++nD)
	   {
	     for (int nP =-1*HEDeltaIphi;nP<=HEDeltaIphi;++nP)
	       {
		 for (int nE =-1*HENeighborParams_.DeltaIeta;nE<=HENeighborParams_.DeltaIeta;++nE)
		   {
		     if (nD==0 && nE==0 && nP==0) 
		       continue; // don't count the cell itself
		     
		     int myphi=nP+iphi;
                     if (myphi>72) myphi-=72; // allow for wrapping of cells
		     if (myphi<=0) myphi+=72;
		     if (!validDetId((HcalSubdetector)(2),nE+ieta, myphi, nD+depth)) continue;
                     HcalDetId myid((HcalSubdetector)(2), nE+ieta, myphi, nD+depth); // HE
		     if (rechitEnergies_.find(myid)==rechitEnergies_.end())
		       continue;
		     if (rechitEnergies_[myid]<HENeighborParams_.minNeighborEnergy)
		       continue;
		     ++neighborsfound;
		     enNeighbor+=rechitEnergies_[myid];
		   } // loop over nE (neighbor eta)
	       } // loop over nP (neighbor phi)
	   } // loop over nD depths

	 if (hotmon_makeDiagnostics_)
	   d_HEenergyVsNeighbor->Fill(en,enNeighbor);
	 
	 // Case 1:  Not enough good neighbors found
	 if (neighborsfound==0)
	   continue;
	 // Case 2:  energy/(avg. neighbor energy) too large for cell to be considered hot
	 if ((1.*enNeighbor/en)>HENeighborParams_.HotEnergyFrac && en>0 && enNeighbor>0)
	   continue;
	 // Case 3:  Tests passed; cell marked as hot
	 aboveneighbors[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1]++;
       }
} //for (HBHERecHitCollection::const_iterator HBHEiter=...)

 // loop over HO
 if (checkHO_)
   {
     for (HORecHitCollection::const_iterator HOiter=hoHits.begin(); HOiter!=hoHits.end(); ++HOiter) 
       { // loop over all hits
	 float en = HOiter->energy();
	 HcalDetId id(HOiter->detid().rawId());
	 int ieta = id.ieta();
	 int iphi = id.iphi();
	 int depth = id.depth();
	 
	 // Case 0:  energy > max value; it's marked as hot regardless of neighbors
	 if (en>HONeighborParams_.maxEnergy)
	   {
	     aboveneighbors[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1]++;
	     continue;
	   }
	 if (en<HONeighborParams_.minCellEnergy)
	    continue; // cells below this value can never be considered hot

	 // Search keys for neighboring cells
	 neighborsfound=0;
	 enNeighbor=0;
	 for (int nD=-1*HONeighborParams_.DeltaDepth;nD<=HONeighborParams_.DeltaDepth;++nD)
	   {
	     for (int nP =-1*HONeighborParams_.DeltaIphi;nP<=HONeighborParams_.DeltaIphi;++nP)
	       {
		 for (int nE =-1*HONeighborParams_.DeltaIeta;nE<=HONeighborParams_.DeltaIeta;++nE)
		   {
		     if (nD==0 && nE==0 && nP==0) 
		       continue; // don't count the cell itself
		     int myphi=nP+iphi;
		     if (myphi>72) myphi-=72; // allow for wrapping of cells
		     if (myphi<=0) myphi+=72;
                     if (!validDetId((HcalSubdetector)(3),nE+ieta, myphi, nD+depth)) continue;
                     HcalDetId myid((HcalSubdetector)(3), nE+ieta, myphi, nD+depth); // HO
		     if (rechitEnergies_.find(myid)==rechitEnergies_.end())
		       continue;
		     if (rechitEnergies_[myid]<HONeighborParams_.minNeighborEnergy)
		       continue;
		     ++neighborsfound;
		     enNeighbor+=rechitEnergies_[myid];
		   } // loop over nE (neighbor eta)
	       } // loop over nP (neighbor phi)
	   } // loop over nD depths

	 if (hotmon_makeDiagnostics_)
	   d_HOenergyVsNeighbor->Fill(en,enNeighbor);

	 // Case 1:  Not enough good neighbors found
	 if (neighborsfound==0)
	   continue;
	 // Case 2:  energy/(avg. neighbor energy) too large for cell to be considered hot
	 if ((1.*enNeighbor/en)>HONeighborParams_.HotEnergyFrac && en>0 && enNeighbor>0)
	   continue;
	 // Case 3:  Tests passed; cell marked as hot
	 aboveneighbors[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1]++;
       } // loop over hits
   } // if (checkHO_)
 
 // loop over HF
 if (checkHF_)
   {
     for (HFRecHitCollection::const_iterator HFiter=hfHits.begin(); HFiter!=hfHits.end(); ++HFiter) 
       { // loop over all hits
	 float en = HFiter->energy();
	 HcalDetId id(HFiter->detid().rawId());
	 int ieta = id.ieta();
	 int iphi = id.iphi();
	 int depth = id.depth();

	 // Case 0:  energy > max value; it's marked as hot regardless of neighbors
	 if (en>HFNeighborParams_.maxEnergy)
	   {
	     aboveneighbors[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1]++;
	     continue;
	   }
	 if (en<HFNeighborParams_.minCellEnergy)
	   continue; // cells below this value can never be considered hot
	  // Search keys for neighboring cells
	 neighborsfound=0;
	 enNeighbor=0;
	 int HFDeltaIphi = HFNeighborParams_.DeltaIphi;
	 if (abs(ieta)>39) HFDeltaIphi*=2;  // double phi range when segmentation switches to 20 degrees
	 // Still need to create a more robust handling of boundary cases
	 for (int nD=-1*HFNeighborParams_.DeltaDepth;nD<=HFNeighborParams_.DeltaDepth;++nD)
	   {
	     for (int nP =-1*HFDeltaIphi;nP<=HFDeltaIphi;++nP)
	       {
		 for (int nE =-1*HFNeighborParams_.DeltaIeta;nE<=HFNeighborParams_.DeltaIeta;++nE)
		   {
		     if (nD==0 && nE==0 && nP==0) 
		       continue; // don't count the cell itself
		     int myphi=nP+iphi;
                     if (myphi>72) myphi-=72; // allow for wrapping of cells
                     if (myphi<=0) myphi+=72;
		     if (!validDetId((HcalSubdetector)(4),nE+ieta, myphi, nD+depth)) continue;
		     HcalDetId myid((HcalSubdetector)(4), nE+ieta, myphi, nD+depth); // HF
		     if (rechitEnergies_.find(myid)==rechitEnergies_.end())
		       continue;
		     if (rechitEnergies_[myid]<HFNeighborParams_.minNeighborEnergy)
		       continue;
		     ++neighborsfound;
		     enNeighbor+=rechitEnergies_[myid];
		   } // loop over nE (neighbor eta)
	       } // loop over nP (neighbor phi)
	   } // loop over nD depths

	 if (hotmon_makeDiagnostics_)
	   d_HFenergyVsNeighbor->Fill(en,enNeighbor);
	 
	 // Case 1:  Not enough good neighbors found
	 if (neighborsfound==0)
	   continue;
	 // Case 2:  energy/(avg. neighbor energy) too large for cell to be considered hot
	 if ((1.*enNeighbor/en)>HFNeighborParams_.HotEnergyFrac && en>0 && enNeighbor>0)
	   continue;
	 // Case 3:  Tests passed; cell marked as hot
	 aboveneighbors[ieta+(int)((etaBins_-2)/2)][iphi-1][depth+1]++;
       } // loop over all hits
   } // if (checkHF_)
 
 
 // Fill histograms 
  if (ievt_%hotmon_checkNevents_neighbor_==0)
    {
	if (fVerbosity) cout <<"<HcalHotCellMonitor::processEvent_digi> Filling HotCell Neighbor plots"<<endl;
	fillNevents_neighbor();
    }

 if (showTiming)
   {
     cpu_timer.stop();  cout <<"TIMER:: HcalHotCellMonitor PROCESSEVENT_RECHITNEIGHBOR -> "<<cpu_timer.cpuTime()<<endl;
   }
 return;
} // void HcalHotCellMonitor::processEvent_rechitneighbor


/* --------------------------------------- */


void HcalHotCellMonitor::processEvent_pedestal( const HBHEDigiCollection& hbhedigi,
						const HODigiCollection& hodigi,
						const HFDigiCollection& hfdigi,
						//const ZDCDigiCollection& zdcdigi, 
						const HcalDbService& cond
						)
{
  if (!hotmon_test_pedestal_) return;
  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

  if (fVerbosity>1) cout <<"<HcalHotCellMonitor::processEvent_pedestal> Processing digis..."<<endl;

  // Variables used in pedestal check
  float digival=0;
  float maxval=0;
  int maxbin=0;
  float ADCsum=0;

  // Variables used in occupancy check
  int ieta=0;
  int iphi=0;
  int depth=0;

  //HcalCalibrationWidths widths;
  HcalCalibrations calibs;
  const HcalQIEShape* shape=cond.getHcalShape();

  // Loop over HBHE digis

  if (fVerbosity>2) cout <<"<HcalHotCellMonitor::processEvent_pedestal> Processing HBHE..."<<endl;

  for (HBHEDigiCollection::const_iterator j=hbhedigi.begin();
       j!=hbhedigi.end(); ++j)
    {
      digival=0;
      maxval=0;
      maxbin=0;
      ADCsum=0;
      const HBHEDataFrame digi = (const HBHEDataFrame)(*j);
      if (!checkHB_ && (HcalSubdetector)(digi.id().subdet())==HcalBarrel) continue;
      if (!checkHE_ && (HcalSubdetector)(digi.id().subdet())==HcalEndcap) continue;
 
      ieta=digi.id().ieta();
      iphi=digi.id().iphi();
      depth=digi.id().depth();
      
      HcalDetId myid = digi.id();
      //cond.makeHcalCalibrationWidth(digi.id(),&widths);
      const HcalCalibrationWidths widths = cond.getHcalCalibrationWidths(digi.id());

      calibs = cond.getHcalCalibrations(digi.id());

      // Find digi time slice with maximum (pedestal-subtracted) ADC count
      for (int i=0;i<digi.size();++i)
	{
	  int thisCapid = digi.sample(i).capid();
	  if (doFCpeds_)
	    {
	      const HcalQIECoder* coder  = cond.getHcalCoder(digi.id());
	      digival = coder->charge(*shape,digi.sample(i).adc(),digi.sample(i).capid())-calibs.pedestal(thisCapid);
	    }
	  else
	    digival=digi.sample(i).adc()-calibs.pedestal(thisCapid);
	  
	  // Find maximum pedestal-subtracted digi value
	  if (digival>maxval)
	    {
	      maxval=digival;
	      maxbin=i;
	    }
	} // for (int i=0;i<digi.size();++i)
      
      // We'll assume steeply-peaked distribution, so that charge deposit occurs
      // in slices (i-1) -> (i+2) around maximum deposit time i
      
      int bins=0;
      for (int i=max(0,maxbin-1);i<=min(digi.size()-1,maxbin+2);++i)
	{
	  ADCsum+=digi.sample(i).adc();
	  ++bins;
	} // for (int i=max(0,maxbin-1);...)      

      // Compare ADCsum to minimum expected value (pedestal+nsigma)
      // we want to compare the average over the sum to the average (ped+nsigma)
      ADCsum*=1./bins;
      // Search for digi in map of pedestal+threshold values
      if (pedestal_thresholds_.find(myid)!=pedestal_thresholds_.end())
	{
	  if (ADCsum > pedestal_thresholds_[myid])
	    abovepedestal[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1]++;
	  if (hotmon_makeDiagnostics_)
	    {
	      if (widths_[myid]==0) continue;
	      if (myid.subdet()==HcalBarrel)
		d_HBnormped->Fill(1.*(ADCsum-pedestals_[myid])/widths_[myid]);
	      else
		d_HEnormped->Fill(1.*(ADCsum-pedestals_[myid])/widths_[myid]);
	    } // if (hotmon_makeDiagnostics)
	}
    } // for (HBHEDigiCollection...)

  // Loop over HO
  if (checkHO_)
    {
      if (fVerbosity>2) cout <<"<HcalHotCellMonitor::processEvent_pedestal> Processing HO..."<<endl;
      
      for (HODigiCollection::const_iterator j=hodigi.begin();
	   j!=hodigi.end(); ++j)
	{
	  digival=0;
	  maxval=0;
	  maxbin=0;
	  ADCsum=0;
	  const HODataFrame digi = (const HODataFrame)(*j);
	  
	  ieta=digi.id().ieta();
	  iphi=digi.id().iphi();
	  depth=digi.id().depth();
	  
	  HcalDetId myid = digi.id();
	  //cond.makeHcalCalibrationWidth(digi.id(),&widths);
	  const HcalCalibrationWidths widths = cond.getHcalCalibrationWidths(digi.id());

	  calibs = cond.getHcalCalibrations(digi.id());
	  
	  for (int i=0;i<digi.size();++i)
	    {
	      int thisCapid = digi.sample(i).capid();
	      if (doFCpeds_)
		{
		  const HcalQIECoder* coder  = cond.getHcalCoder(digi.id());
		  digival = coder->charge(*shape,digi.sample(i).adc(),digi.sample(i).capid())-calibs.pedestal(thisCapid);
		}
	      else
		digival=digi.sample(i).adc()-calibs.pedestal(thisCapid);
	  
	      // Find maximum pedestal-subtracted digi value
	      if (digival>maxval)
		{
		  maxval=digival;
		  maxbin=i;
		}
	    } // for (int i=0;i<digi.size();++i)
      
	  // We'll assume steeply-peaked distribution, so that charge deposit occurs
	  // in slices (i-1) -> (i+2) around maximum deposit time i
      
	  int bins=0;
	  for (int i=max(0,maxbin-1);i<=min(digi.size()-1,maxbin+2);++i)
	    {
	      ADCsum+=digi.sample(i).adc();
	      ++bins;
	    } // for (int i=max(0,maxbin-1);...)      
	  
	  ADCsum*=1./bins;
	  // Search for digi in map of pedestal+threshold values
	  if (pedestal_thresholds_.find(myid)!=pedestal_thresholds_.end())
	    {
	      if (ADCsum > pedestal_thresholds_[myid])
		abovepedestal[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1]++;
	      if (hotmon_makeDiagnostics_)
		{
		  if (widths_[myid]==0) continue;
		  d_HOnormped->Fill(1.*(ADCsum-pedestals_[myid])/widths_[myid]);
		} // if (hotmon_makeDiagnostics)
	    }
	} // for (HODigiCollection...)
    } // if (checkHO_)

  if (checkHF_)
    {
      // Loop over HF
      if (fVerbosity>2) cout <<"<HcalHotCellMonitor::processEvent_pedestal> Processing HF..."<<endl;

      for (HFDigiCollection::const_iterator j=hfdigi.begin();
	   j!=hfdigi.end(); ++j)
	{
	  digival=0;
	  maxval=0;
	  maxbin=0;
	  ADCsum=0;
	  const HFDataFrame digi = (const HFDataFrame)(*j);

	  ieta=digi.id().ieta();
	  iphi=digi.id().iphi();
	  depth=digi.id().depth()+2; // offset depth by 2 for HF

	  HcalDetId myid = digi.id();
	  //cond.makeHcalCalibrationWidth(digi.id(),&widths);
	  const HcalCalibrationWidths widths = cond.getHcalCalibrationWidths(digi.id());

	  calibs = cond.getHcalCalibrations(digi.id());

	  for (int i=0;i<digi.size();++i)
	    {
	      int thisCapid = digi.sample(i).capid();
	      if (doFCpeds_)
		{
		  const HcalQIECoder* coder  = cond.getHcalCoder(digi.id());
		  digival = coder->charge(*shape,digi.sample(i).adc(),digi.sample(i).capid())-calibs.pedestal(thisCapid);
		}
	      else
		digival=digi.sample(i).adc()-calibs.pedestal(thisCapid);
	  
	      // Find maximum pedestal-subtracted digi value
	      if (digival>maxval)
		{
		  maxval=digival;
		  maxbin=i;
		}
	    } // for (int i=0;i<digi.size();++i)
      
	  // We'll assume steeply-peaked distribution, so that charge deposit occurs
	  // in slices (i-1) -> (i+2) around maximum deposit time i
      
	  int bins=0;
	  for (int i=max(0,maxbin-1);i<=min(digi.size()-1,maxbin+2);++i)
	    {
	      ADCsum+=digi.sample(i).adc();
	      ++bins;
	    } // for (int i=max(0,maxbin-1);...)      

	  ADCsum*=1./bins;
	  // Search for digi in map of pedestal+threshold values
	  if (pedestal_thresholds_.find(myid)!=pedestal_thresholds_.end())
	    {
	      if (ADCsum > pedestal_thresholds_[myid])
		abovepedestal[ieta+(int)((etaBins_-2)/2)][iphi-1][depth-1]++;
	      if (hotmon_makeDiagnostics_)
		{
		  if (widths_[myid]==0) continue;
		  d_HFnormped->Fill(1.*(ADCsum-pedestals_[myid])/widths_[myid]);
		} // if (hotmon_makeDiagnostics)
	    }
	} // for (HFDigiCollection...)
    } // if (checkHF_)

  // Fill histograms 
  if (ievt_%hotmon_checkNevents_pedestal_==0)
    {
      if (fVerbosity) cout <<"<HcalHotCellMonitor::processEvent_pedestal> Filling HotCell Pedestal plots"<<endl;
      fillNevents_pedestal();
    }

   if (showTiming)
    {
      cpu_timer.stop();  cout <<"TIMER:: HcalHotCellMonitor PROCESSEVENT_PEDESTAL -> "<<cpu_timer.cpuTime()<<endl;
    }

  return;
} // void HcalHotCellMonitor::processEvent_pedestal


/* ----------------------------------- */


void HcalHotCellMonitor::fillNevents_persistentenergy(void)
{
  // Fill Histograms showing rechits with energies > some threshold for N consecutive events

  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

  if (fVerbosity>0)
    cout <<"<HcalHotCellMonitor::fillNevents_persistentenergy> FILLING PERSISTENT ENERGY PLOTS"<<endl;
  
  int mydepth=0;
  int ieta=0;
  int iphi=0;
  for (int eta=0;eta<(etaBins_-2);++eta)
    {
      ieta=eta-int((etaBins_-2)/2);
      for (int phi=0;phi<72;++phi)
        {
	  iphi=phi+1;
	  for (int depth=0;depth<4;++depth) // this is one unit less "true" depth (for indexing purposes)
            {
	      for (int subdet=1;subdet<=4;++subdet)
		{
		  if (!validDetId((HcalSubdetector)subdet, ieta, iphi, depth+1))
		    continue;
		  mydepth=depth;
		  if (subdet==4) // remember that HF's elements stored in depths (2,3), not (0,1)
		    mydepth=depth+2;
		  if (hotmon_makeDiagnostics_ && rechit_occupancy_sum[eta][phi][mydepth]>0)
		    {
		      // Fill average energy plots
		      int filldepth=depth;
		      if (subdet==2 && depth<2) filldepth+=4; // HB depths 1&2 get shifted to last two histograms
		      d_avgrechitenergymap[filldepth]->setBinContent(eta+2,phi+2,rechit_energy_sum[eta][phi][mydepth]/rechit_occupancy_sum[eta][phi][mydepth]);
		    }

		  int oldevts=(ievt_/hotmon_checkNevents_persistent_);
		  if (ievt_%hotmon_checkNevents_persistent_==0)
		    oldevts-=1;
		  oldevts*=hotmon_checkNevents_persistent_;
		  int newevts=ievt_-oldevts;
		  if (newevts<=0) continue; // shouldn't happen?
		  // MUST BE ABOVE ENERGY THRESHOLD FOR ALL N EVENTS
		  if (abovepersistent[eta][phi][mydepth]<(unsigned int)newevts)
		    {
		      abovepersistent[eta][phi][mydepth]=0;
		      continue;  // should never be > newevts, either, right?
		    }
		  if (fVerbosity>0) cout <<"HOT CELL; PERSISTENT ENERGY = "<<subdet<<" eta = "<<ieta<<", phi = "<<iphi<<" depth = "<<depth+1<<endl;
		  if (subdet==2 && depth<2) // HE depth positions(0,1) found -- shift up to positions (4,5)
		    mydepth=depth+4;
		  else
		    mydepth=depth; // switches back HF to its correct depth
		  // BinContent starts at 1, not 0 (offset by 0)
		  // Offset by another 1 due to empty bins at edges
		  AbovePersistentThresholdCellsByDepth[mydepth]->setBinContent( eta+2,phi+2,
										(oldevts*AbovePersistentThresholdCellsByDepth[mydepth]->getBinContent(eta+2,phi+2)+newevts)*1./ievt_);
		} // for (int subdet=1;subdet<=4;++subdet)
	      abovepersistent[eta][phi][mydepth]=0; // reset counter
	    } // for (int depth=0;depth<4;++depth)
	} // for (int phi=0;...)
    } // for (int eta=0;...)
  
  if (showTiming)
    {
      cpu_timer.stop();  cout <<"TIMER:: HcalHotCellMonitor FILLNEVENTS_PERSISTENTENERGY -> "<<cpu_timer.cpuTime()<<endl;
    }

  return;
} // void HcalHotCellMonitor::fillNevents_persistentenergy(void)



/* ----------------------------------- */

void HcalHotCellMonitor::fillNevents_pedestal(void)
{
  // Fill Histograms showing digi cells consistently above pedestal values

  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

  if (fVerbosity>0)
    cout <<"<HcalHotCellMonitor::fillNevents_pedestal> FILLING HOT CELL PEDESTAL PLOTS"<<endl;

  int mydepth=0;
  int ieta=0;
  int iphi=0;
  for (int eta=0;eta<(etaBins_-2);++eta)
    {
      ieta=eta-int((etaBins_-2)/2);
      for (int phi=0;phi<72;++phi)
        {
	  iphi=phi+1;
	  for (int depth=0;depth<4;++depth) // this is one unit less "true" depth (for indexing purposes)
            {
	      for (int subdet=1;subdet<=4;++subdet)
		{
		  if (!validDetId((HcalSubdetector)subdet, ieta, iphi, depth+1))
		    continue;
		  
		  int oldevts=(ievt_/hotmon_checkNevents_pedestal_);
		  if (ievt_%hotmon_checkNevents_pedestal_==0)
		    oldevts-=1;
		  oldevts*=hotmon_checkNevents_pedestal_;
		  int newevts=ievt_-oldevts;
		  if (newevts<0) newevts=0; // shouldn't happen

		  mydepth=depth;
		  if (subdet==4) // remember that HF's elements stored in depths (2,3), not (0,1)
		    mydepth=depth+2;

		  // Now that we have a valid cell, check whether it was above pedestal threshold for all N events
		  if (abovepedestal[eta][phi][mydepth]<(unsigned int)newevts)
		    {
		      abovepedestal[eta][phi][mydepth]=0;
		      continue; // cell was above pedestal threshold at least once; ignore it
		    }
		  if (fVerbosity>0) 
		    cout <<"HOT CELL; ABOVE PEDESTAL THRESHOLD = "<<subdet<<" eta = "<<ieta<<", phi = "<<iphi<<" depth = "<<depth+1<<endl;
		  if (subdet==2 && depth<2) // HE depth positions(0,1) found -- shift up to positions (4,5)
		    mydepth=depth+4;
		  else
		    mydepth=depth; // switches back HF to its correct depth
		  
		  if (fVerbosity>0)
		    {
		      cout <<"\t MYDEPTH = "<<mydepth<<endl;
		      cout <<"\t oldevents = "<<oldevts<<"  new = "<<newevts<<endl;
		      cout <<"\t\t"<<(oldevts*AbovePedestalHotCellsByDepth[mydepth]->getBinContent(eta+2,phi+2)+newevts)*1./ievt_<<endl;
		    }
		  // BinContent starts at 1, not 0 (offset by 0)
		  // Offset by another 1 due to empty bins at edges
		  AbovePedestalHotCellsByDepth[mydepth]->setBinContent( eta+2,phi+2,
									 (oldevts*AbovePedestalHotCellsByDepth[mydepth]->getBinContent(eta+2,phi+2)+newevts)*1./ievt_);
		  //reset counter -- shouldn't be necessary (counter should already be 0).
		  abovepedestal[eta][phi][depth]=0;
		} // for (int subdet=1;subdet<=4;++subdet)
	      
	    } // for (int depth=0;depth<4;++depth)
	} // for (int phi=0;...)
    } // for (int eta=0;...)
  
  if (showTiming)
    {
      cpu_timer.stop();  cout <<"TIMER:: HcalHotCellMonitor FILLNEVENTS_ABOVEPEDESTAL -> "<<cpu_timer.cpuTime()<<endl;
    }

  return;


} // void HcalHotCellMonitor::fillNevents_pedestal(void)


/* ----------------------------------- */

void HcalHotCellMonitor::fillNevents_energy(void)
{
  // Fill Histograms showing rec hits that are above some energy value 
  // (Fill for each instance when cell is above energy; don't require it to be hot for a number of consecutive events)

  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

  if (fVerbosity>0)
    cout <<"<HcalHotCellMonitor::fillNevents_energy> ABOVE-ENERGY-THRESHOLD PLOTS"<<endl;

  int mydepth=0;
  int ieta=0;
  int iphi=0;
  for (int eta=0;eta<(etaBins_-2);++eta)
    {
      ieta=eta-int((etaBins_-2)/2);
      for (int phi=0;phi<72;++phi)
        {
	  iphi=phi+1;
	  for (int depth=0;depth<4;++depth) // this is one unit less "true" depth (for indexing purposes)
            {
	      for (int subdet=1;subdet<=4;++subdet)
		{
		  if (!validDetId((HcalSubdetector)subdet, ieta, iphi, depth+1))
		    continue;
		  mydepth=depth;
		  if (subdet==4) // remember that HF's elements stored in depths (2,3), not (0,1)
		    mydepth=depth+2;

		  if (aboveenergy[eta][phi][mydepth]==0)
		    {
		      // Cell never was above threshold energy; skip it
		      continue;
		    }
		  // Allow for case when we perform check at end of run, so that ievt_%checkNevents != 0 ?
		  
		  // Calculate total number of events in prior histogram filled (= ievt_-N in most circumstances).
		  int oldevts=(ievt_/hotmon_checkNevents_pedestal_);
		  if (ievt_%hotmon_checkNevents_pedestal_==0)
		    oldevts-=1;
		  oldevts*=hotmon_checkNevents_pedestal_;
		  int newevts=aboveenergy[eta][phi][mydepth];  // number

		  if (fVerbosity>2) 
		    cout <<"HOT CELL; ABOVE ENERGY THRESHOLD = "<<subdet<<" eta = "<<ieta<<", phi = "<<iphi<<" depth = "<<depth+1<<endl;
		  if (subdet==2 && depth<2) // HE depth positions(0,1) found -- shift up to positions (4,5)
		    mydepth=depth+4;
		  else
		    mydepth=depth; // switches back HF to its correct depth
		  
		  // BinContent starts at 1, not 0 (offset by 0)
		  // Offset by another 1 due to empty bins at edges
		  AboveEnergyThresholdCellsByDepth[mydepth]->setBinContent( eta+2,phi+2,
									    (oldevts*AboveEnergyThresholdCellsByDepth[mydepth]->getBinContent(eta+2,phi+2)+newevts)*1./ievt_);
		} // for (int subdet=1;subdet<=4;++subdet)
	      // reset counter -- shouldn't be necessary (should already be 0)
	      aboveenergy[eta][phi][depth]=0;
	    } // for (int depth=0;depth<4;++depth)
	} // for (int phi=0;...)
    } // for (int eta=0;...)

  if (showTiming)
    {
      cpu_timer.stop();  cout <<"TIMER:: HcalHotCellMonitor FILLNEVENTS_ENERGY -> "<<cpu_timer.cpuTime()<<endl;
    }

  return;


} // void HcalHotCellMonitor::fillNevents_energy(void)



/* ----------------------------------- */

void HcalHotCellMonitor::fillNevents_neighbor(void)
{
  // Fill Histograms showing rec hits with energy much less than neighbors' average

  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

  if (fVerbosity>0)
    cout <<"<HcalHotCellMonitor::fillNevents_neighbor> FILLING ABOVE-NEIGHBOR-ENERGY PLOTS"<<endl;

  int mydepth=0;
  int ieta=0;
  int iphi=0;
  for (int eta=0;eta<(etaBins_-2);++eta)
    {
      ieta=eta-int((etaBins_-2)/2);
      for (int phi=0;phi<72;++phi)
        {
	  iphi=phi+1;
	  for (int depth=0;depth<4;++depth) // this is one unit less "true" depth (for indexing purposes)
            {
	      for (int subdet=1;subdet<=4;++subdet)
		{
		  if (!validDetId((HcalSubdetector)subdet, ieta, iphi, depth+1))
		    continue;
		  mydepth=depth;
		  if (subdet==4) // remember that HF's elements stored in depths (2,3), not (0,1)
		    mydepth=depth+2;
		  if (aboveneighbors[eta][phi][mydepth]>0)
		    {
		      if (fVerbosity>2) cout <<"HOT CELL; ABOVE NEIGHBORS = "<<subdet<<" eta = "<<ieta<<", phi = "<<iphi<<" depth = "<<depth+1<<endl;
		      if (subdet==2 && depth<2) // HE depth positions(0,1) found -- shift up to positions (4,5)
			mydepth=depth+4;
		      else
			mydepth=depth; // switches back HF to its correct depth
		      // no digi was found for the N events; set histogram error rate
		      int oldevts=(ievt_/hotmon_checkNevents_neighbor_);
		      if (ievt_%hotmon_checkNevents_neighbor_==0)
			oldevts-=1;
		      oldevts*=hotmon_checkNevents_neighbor_;
		      // BinContent starts at 1, not 0 (offset by 0)
		      // Offset by another 1 due to empty bins at edges
		      AboveNeighborsHotCellsByDepth[mydepth]->setBinContent( eta+2,phi+2,
									  (oldevts*AboveNeighborsHotCellsByDepth[mydepth]->getBinContent(eta+2,phi+2)+aboveneighbors[eta][phi][mydepth])*1./ievt_);
		      //reset counter
		      aboveneighbors[eta][phi][depth]=0;
		    } // if (aboveneighbors[eta][phi][mydepth]>0)
		} // for (int subdet=1;subdet<=4;++subdet)

	    } // for (int depth=0;depth<4;++depth)
	} // for (int phi=0;...)
    } // for (int eta=0;...)

  if (showTiming)
    {
      cpu_timer.stop();  cout <<"TIMER:: HcalHotCellMonitor FILLNEVENTS_NEIGHBOR -> "<<cpu_timer.cpuTime()<<endl;
    }

  return;


} // void HcalHotCellMonitor::fillNevents_neighbor(void)






void HcalHotCellMonitor::fillNevents_problemCells(void)
{
  if (showTiming)
    {
      cpu_timer.reset(); cpu_timer.start();
    }

  if (fVerbosity>0)
    cout <<"<HcalHotCellMonitor::fillNevents_problemCells> FILLING PROBLEM CELL PLOTS"<<endl;

  int ieta=0;
  int iphi=0;

  double problemvalue=0;
  double sumproblemvalue=0; // summed over all depths
  for (int eta=0;eta<(etaBins_-2);++eta)
    {
      ieta=eta-int((etaBins_-2)/2);
      for (int phi=0;phi<72;++phi)
        {
	  iphi=phi+1;
	  sumproblemvalue=0;
	  for (int mydepth=0;mydepth<6;++mydepth)
	    {
	      // total bad fraction is sum of fractions from individual tests
	      // (eventually, do we want to be more careful about how we handle this, in case checkNevents is
	      //  drastically different for the different tests?)
	      problemvalue=0;
	      if (hotmon_test_persistent_)
		{
		  problemvalue+=AbovePersistentThresholdCellsByDepth[mydepth]->getBinContent(eta+2,phi+2);
		  sumproblemvalue+=AbovePersistentThresholdCellsByDepth[mydepth]->getBinContent(eta+2,phi+2);
		}
	      if (hotmon_test_pedestal_)
		{
		  problemvalue+=AbovePedestalHotCellsByDepth[mydepth]->getBinContent(eta+2,phi+2);
		  sumproblemvalue+=AbovePedestalHotCellsByDepth[mydepth]->getBinContent(eta+2,phi+2);
		}
	      if (hotmon_test_neighbor_)
		{
		  problemvalue+=AboveNeighborsHotCellsByDepth[mydepth]->getBinContent(eta+2,phi+2);
		  sumproblemvalue+=AboveNeighborsHotCellsByDepth[mydepth]->getBinContent(eta+2,phi+2);
		}
	      if (hotmon_test_energy_)
		{
		  problemvalue+=AboveEnergyThresholdCellsByDepth[mydepth]->getBinContent(eta+2,phi+2);
		  sumproblemvalue+=AboveEnergyThresholdCellsByDepth[mydepth]->getBinContent(eta+2,phi+2);
		}
	      problemvalue=min(1.,problemvalue);
	      ProblemHotCellsByDepth[mydepth]->setBinContent(eta+2,phi+2,problemvalue);
	    } // for (int mydepth=0;mydepth<6;...)
	  sumproblemvalue=min(1.,sumproblemvalue);
	  ProblemHotCells->setBinContent(eta+2,phi+2,sumproblemvalue);
	} // loop on phi=0;phi<72
    } // loop on eta=0; eta<(etaBins_-2)
  
  if (showTiming)
    {
      cpu_timer.stop();  cout <<"TIMER:: HcalHotCellMonitor FILLNEVENTS_PROBLEMCELLS -> "<<cpu_timer.cpuTime()<<endl;
    }

} // void HcalHotCellMonitor::fillNevents_problemCells(void)
