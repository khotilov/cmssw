/*
 * \file L1TGCT.cc
 *
 * $Date: 2007/09/26 15:26:23 $
 * $Revision: 1.12 $
 * \author J. Berryhill
 *
 * $Log: L1TGCT.cc,v $
 * Revision 1.12  2007/09/26 15:26:23  berryhil
 *
 *
 * restored L1TGCT.cc
 *
 * Revision 1.10  2007/09/05 22:31:36  wittich
 * - Factorize getByLabels to approximate my understanding of what the
 *   HW can do.
 * - tested (loosely speaking) on GREJ' data.
 *
 * Revision 1.9  2007/09/04 02:54:19  wittich
 * - fix dupe ME in RCT
 * - put in rank>0 req in GCT
 * - various small other fixes
 *
 * Revision 1.8  2007/08/31 18:14:21  wittich
 * update GCT packages to reflect GctRawToDigi, and move to raw plots
 *
 * Revision 1.7  2007/08/31 11:02:56  wittich
 * cerr -> LogInfo
 *
 * Revision 1.6  2007/02/22 19:43:53  berryhil
 *
 *
 *
 * InputTag parameters added for all modules
 *
 *
 *
 */

#include "DQM/L1TMonitor/interface/L1TGCT.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Trigger Headers

// GCT and RCT data formats
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEtSums.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

// L1Extra
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"


using namespace edm;
using namespace l1extra;


// Define statics for bins etc.
const unsigned int ETABINS = 22;
const float ETAMIN = -0.5;
const float ETAMAX = 21.5;

const unsigned int METPHIBINS = 72;
const float METPHIMIN = -0.5;
const float METPHIMAX = 71.5;

const unsigned int PHIBINS = 18;
const float PHIMIN = -0.5;
const float PHIMAX = 17.5;


// const unsigned int L1EETABINS = 22;
// const float L1EETAMIN = -5;
// const float L1EETAMAX = 5;

// const unsigned int L1EPHIBINS = 18;
// const float L1EPHIMIN = -M_PI;
// const float L1EPHIMAX = M_PI;

// Ranks 6, 10 and 12 bits
const unsigned int R6BINS = 64;
const float R6MIN = -0.5;
const float R6MAX = 63.5;
const unsigned int R10BINS = 1024;
const float R10MIN = -0.5;
const float R10MAX = 1023.5;
const unsigned int R12BINS = 4096;
const float R12MIN = -0.5;
const float R12MAX = 4095.5;

// // Physical bins 1 Gev - 1 TeV in 1 GeV steps
// const unsigned int TEVBINS = 1001;
// const float TEVMIN = -0.5;
// const float TEVMAX = 1000.5;

const unsigned int METBINS = 256;
const float METMIN = -0.5;
const float METMAX = 4095.5;

// simple helper functions to take eta to a signed quantity
static int etaBin( unsigned int etaIndex)
{
  return (etaIndex&0x7U)*(etaIndex&0x08U?-1:1);
}


L1TGCT::L1TGCT(const edm::ParameterSet & ps) :
  gctSource_(ps.getParameter<edm::InputTag>("gctSource"))
{

  // verbosity switch
  verbose_ = ps.getUntrackedParameter < bool > ("verbose", false);

  if (verbose_)
    std::cout << "L1TGCT: constructor...." << std::endl;

  logFile_.open("L1TGCT.log");

  dbe = NULL;
  if (ps.getUntrackedParameter < bool > ("DaqMonitorBEInterface", false)) {
    dbe = edm::Service < DaqMonitorBEInterface > ().operator->();
    dbe->setVerbose(0);
  }

  monitorDaemon_ = false;
  if (ps.getUntrackedParameter < bool > ("MonitorDaemon", false)) {
    edm::Service<MonitorDaemon> daemon;
    daemon.operator->();
    monitorDaemon_ = true;
  }

  outputFile_ = ps.getUntrackedParameter < std::string > ("outputFile", "");
  if (outputFile_.size() != 0) {
    std::cout << "L1T Monitoring histograms will be saved to "
	      << outputFile_ << std::endl;
  }
  else {
    outputFile_ = "L1TDQM.root";
  }

  bool disable =
      ps.getUntrackedParameter < bool > ("disableROOToutput", false);
  if (disable) {
    outputFile_ = "";
  }


  if (dbe != NULL) {
    dbe->setCurrentFolder("L1TMonitor/L1TGCT");
  }


}

L1TGCT::~L1TGCT()
{
}

void L1TGCT::beginJob(const edm::EventSetup & c)
{

  nev_ = 0;

  // get hold of back-end interface
  DaqMonitorBEInterface *dbe = 0;
  dbe = edm::Service < DaqMonitorBEInterface > ().operator->();

  if (dbe) {
    dbe->setCurrentFolder("L1TMonitor/L1TGCT");
    dbe->rmdir("L1TMonitor/L1TGCT");
  }


  if (dbe) {
    dbe->setCurrentFolder("L1TMonitor/L1TGCT");


    const int nJetEta = 256;
    const float JetEtaMin = -127.5;
    const float JetEtaMax =  127.5;
    l1GctCenJetsEtEtaPhi_ =
	dbe->book2D("CenJetsEtEtaPhi", "CENTRAL JET E_{T}",
		    PHIBINS, PHIMIN, PHIMAX, 
		    nJetEta, JetEtaMin, JetEtaMax);
    l1GctForJetsEtEtaPhi_ =
	dbe->book2D("ForJetsEtEtaPhi", "FORWARD JET E_{T}",
		    PHIBINS, PHIMIN, PHIMAX, 
		    nJetEta, JetEtaMin, JetEtaMax);
    l1GctTauJetsEtEtaPhi_ =
	dbe->book2D("TauJetsEtEtaPhi", "TAU JET E_{T}", PHIBINS,
		    PHIMIN, PHIMAX, 
		    nJetEta, JetEtaMin, JetEtaMax);
    l1GctIsoEmEtEtaPhi_ =
	dbe->book2D("IsoEmEtEtaPhi", "ISO EM E_{T}", PHIBINS,
		    PHIMIN, PHIMAX, 		    
		    nJetEta, JetEtaMin, JetEtaMax);
    l1GctNonIsoEmEtEtaPhi_ =
	dbe->book2D("NonIsoEmEtEtaPhi", "NON-ISO EM E_{T}",
		    PHIBINS, PHIMIN, PHIMAX, 
		    nJetEta, JetEtaMin, JetEtaMax);

    l1GctCenJetsOccEtaPhi_ =
	dbe->book2D("CenJetsOccEtaPhi", "CENTRAL JET OCCUPANCY",
		    PHIBINS, PHIMIN, PHIMAX, 
		    nJetEta, JetEtaMin, JetEtaMax);
    l1GctForJetsOccEtaPhi_ =
	dbe->book2D("ForJetsOccEtaPhi", "FORWARD JET OCCUPANCY",
		    PHIBINS, PHIMIN, PHIMAX,
		    nJetEta, JetEtaMin, JetEtaMax);
    l1GctTauJetsOccEtaPhi_ =
	dbe->book2D("TauJetsOccEtaPhi", "TAU JET OCCUPANCY",
		    PHIBINS, PHIMIN, PHIMAX, 
		    nJetEta, JetEtaMin, JetEtaMax);
    l1GctIsoEmOccEtaPhi_ =
	dbe->book2D("IsoEmOccEtaPhi", "ISO EM OCCUPANCY",
		    PHIBINS, PHIMIN, PHIMAX, 
		    nJetEta, JetEtaMin, JetEtaMax);
    l1GctNonIsoEmOccEtaPhi_ =
	dbe->book2D("NonIsoEmOccEtaPhi", "NON-ISO EM OCCUPANCY",
		    PHIBINS, PHIMIN, PHIMAX, 
		    nJetEta, JetEtaMin, JetEtaMax);

    l1GctCenJetsRank_ =
	dbe->book1D("CenJetsRank", "CENTRAL JET RANK", R6BINS,
		    R6MIN, R6MAX);
    l1GctForJetsRank_ =
	dbe->book1D("ForJetsRank", "FORWARD JET RANK", R6BINS,
		    R6MIN, R6MAX);
    l1GctTauJetsRank_ =
	dbe->book1D("TauJetsRank", "TAU JET RANK", R6BINS, R6MIN,
		    R6MAX);
    l1GctIsoEmRank_ =
	dbe->book1D("IsoEmRank", "ISO EM RANK", R6BINS, R6MIN,
		    R6MAX);
    l1GctNonIsoEmRank_ =
	dbe->book1D("NonIsoEmRank", "NON-ISO EM RANK", R6BINS,
		    R6MIN, R6MAX);

    l1GctEtMiss_ =
	dbe->book1D("EtMiss", "MISSING E_{T}", METBINS, METMIN,
		    METMAX);
    l1GctEtMissPhi_ =
	dbe->book1D("EtMissPhi", "MISSING E_{T} #phi", METPHIBINS,
		    PHIMIN, PHIMAX);
    l1GctEtTotal_ =
	dbe->book1D("EtTotal", "TOTAL E_{T}", METBINS, METMIN,
		    METMAX);
    l1GctEtHad_ =
	dbe->book1D("EtHad", "TOTAL HAD E_{T}", METBINS, METMIN,
		    METMAX);
  }
}


void L1TGCT::endJob(void)
{
  if (verbose_)
    std::cout << "L1TGCT: end job...." << std::endl;
  edm::LogInfo("L1TGCT") << "analyzed " << nev_ << " events";

  if (outputFile_.size() != 0 && dbe) {
    dbe->save(outputFile_);
  }

  return;
}

void L1TGCT::analyze(const edm::Event & e, const edm::EventSetup & c)
{
  nev_++;
  if (verbose_) {
    std::cout << "L1TGCT: analyze...." << std::endl;
  }
  
  // update to those generated in GctRawToDigi
  edm::Handle < L1GctEmCandCollection > l1IsoEm;
  edm::Handle < L1GctEmCandCollection > l1NonIsoEm;
  edm::Handle < L1GctJetCandCollection > l1CenJets;
  edm::Handle < L1GctJetCandCollection > l1ForJets;
  edm::Handle < L1GctJetCandCollection > l1TauJets;
  edm::Handle < L1GctEtMiss >  l1EtMiss;
  edm::Handle < L1GctEtHad >   l1EtHad;
  edm::Handle < L1GctEtTotal > l1EtTotal;

  // Split this into two parts as this appears to be the way the 
  // HW works.
  bool doJet = true;
  bool doEm = true;
  try {
    // observed in emulator data
      e.getByLabel(gctSource_.label(), "cenJets", l1CenJets);
      e.getByLabel(gctSource_.label(), "forJets", l1ForJets);
      e.getByLabel(gctSource_.label(), "tauJets", l1TauJets);

      e.getByLabel(gctSource_, l1EtMiss);
      e.getByLabel(gctSource_, l1EtHad);
      e.getByLabel(gctSource_, l1EtTotal);
   }
   catch (...) {
     edm::LogInfo("L1TGCT") << " Could not find one of the requested data JET"
       " elements, label was " << gctSource_.label() ;
     doJet = false;
   }
  
  // EM data
  try {
    e.getByLabel(gctSource_.label(), "isoEm", l1IsoEm);
    e.getByLabel(gctSource_.label(), "nonIsoEm", l1NonIsoEm);
  }
  catch (...) {
    edm::LogInfo("L1TGCT") << " Could not find one of the requested EM data "
      " elements, label was " << gctSource_.label() ;
    doEm = false;
  }
  if ( (! doEm) && (! doJet) ) {
    if (  verbose_ )
      std::cout << "L1TGCT: Bailing, didn't find squat."<<std::endl;
    return;
  }
  
  
  // Fill the histograms for the jets
  if ( doJet ) {
    // Central jets
    if ( verbose_ ) {
      std::cout << "L1TGCT: number of central jets = " 
		<< l1CenJets->size() << std::endl;
    }
    for (L1GctJetCandCollection::const_iterator cj = l1CenJets->begin();
	 cj != l1CenJets->end(); cj++) {
      if ( cj->rank() == 0 ) continue;
      l1GctCenJetsEtEtaPhi_->Fill(cj->phiIndex(), etaBin(cj->etaIndex()), 
				  cj->rank());
      l1GctCenJetsOccEtaPhi_->Fill(cj->phiIndex(), etaBin(cj->etaIndex()));
      l1GctCenJetsRank_->Fill(cj->rank());
      if ( verbose_ ) {
	std::cout << "L1TGCT: Central jet " 
		  << cj->phiIndex() << ", " << etaBin(cj->etaIndex())
		  << ", (" << cj->etaIndex() << "), " << cj->rank()
		  << std::endl;
      }
    }

    // Forward jets
    if ( verbose_ ) {
      std::cout << "L1TGCT: number of forward jets = " 
		<< l1ForJets->size() << std::endl;
    }
    for (L1GctJetCandCollection::const_iterator fj = l1ForJets->begin();
	 fj != l1ForJets->end(); fj++) {
      if ( fj->rank() == 0 ) continue;
      l1GctForJetsEtEtaPhi_->Fill(fj->phiIndex(), etaBin(fj->etaIndex()), 
				  fj->rank());
      l1GctForJetsOccEtaPhi_->Fill(fj->phiIndex(), etaBin(fj->etaIndex()));
      l1GctForJetsRank_->Fill(fj->rank());
    }

    // Tau jets
    if ( verbose_ ) {
      std::cout << "L1TGCT: number of tau jets = " 
		<< l1TauJets->size() << std::endl;
    }
    for (L1GctJetCandCollection::const_iterator tj = l1TauJets->begin();
	 tj != l1TauJets->end(); tj++) {
      if ( tj->rank() == 0 ) continue;
      l1GctTauJetsEtEtaPhi_->Fill(tj->phiIndex(), etaBin(tj->etaIndex()), 
				  tj->rank());
      l1GctTauJetsOccEtaPhi_->Fill(tj->phiIndex(), etaBin(tj->etaIndex()));
      l1GctTauJetsRank_->Fill(tj->rank());
    }

    // Energy sums
    l1GctEtMiss_->Fill(l1EtMiss->et());
    l1GctEtMissPhi_->Fill(l1EtMiss->phi());

    // these don't have phi values
    l1GctEtHad_->Fill(l1EtHad->et());
    l1GctEtTotal_->Fill(l1EtTotal->et());

  }


  if ( doEm ) {
    // Isolated EM
    if ( verbose_ ) {
      std::cout << "L1TGCT: number of iso em cands: " 
		<< l1IsoEm->size() << std::endl;
    }
    for (L1GctEmCandCollection::const_iterator ie = l1IsoEm->begin();
	 ie != l1IsoEm->end(); ie++) {
      if ( ie->rank() == 0 ) continue;
      l1GctIsoEmEtEtaPhi_->Fill(ie->phiIndex(), etaBin(ie->etaIndex()), 
				ie->rank());
      l1GctIsoEmOccEtaPhi_->Fill(ie->phiIndex(), etaBin(ie->etaIndex()));
      l1GctIsoEmRank_->Fill(ie->rank());
      if ( verbose_ ) {
	std::cout << "L1TGCT: iso em " 
		  << ie->phiIndex() << ", " << etaBin(ie->etaIndex())
		  << ", (" << ie->etaIndex() << "), " << ie->rank()
		  << std::endl;
      }

    }

    // Non-isolated EM
    if ( verbose_ ) {
      std::cout << "L1TGCT: number of non-iso em cands: " 
		<< l1NonIsoEm->size() << std::endl;
    }
    for (L1GctEmCandCollection::const_iterator ne = l1NonIsoEm->begin();
	 ne != l1NonIsoEm->end(); ne++) {
      if ( ne->rank() == 0 ) continue;
      l1GctNonIsoEmEtEtaPhi_->Fill(ne->phiIndex(), etaBin(ne->etaIndex()), 
				   ne->rank());
      l1GctNonIsoEmOccEtaPhi_->Fill(ne->phiIndex(), etaBin(ne->etaIndex()));
      l1GctNonIsoEmRank_->Fill(ne->rank());
    }
  }




}
