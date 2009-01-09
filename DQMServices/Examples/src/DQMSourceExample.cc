/*
 * \file DQMSourceExample.cc
 * \author C.Leonidopoulos
 * Last Update:
 * $Date: 2008/08/12 10:12:29 $
 * $Revision: 1.16 $
 * $Author: markusm $
 *
 * Description: Simple example showing how to create a DQM source creating and filling
 * monitoring elements
*/

#include "DQMServices/Examples/interface/DQMSourceExample.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TRandom.h"
#include <math.h>

using namespace std;
using namespace edm;


// -----------------------------
//  constructors and destructor
// -----------------------------
DQMSourceExample::DQMSourceExample( const edm::ParameterSet& ps ) {
  parameters_ = ps;
  initialize();
}


DQMSourceExample::~DQMSourceExample() {
}


// ----------------------------------
void DQMSourceExample::initialize() {

  counterLS_  = 0;
  counterEvt_ = 0;

  // get back-end interface
  dbe_ = Service<DQMStore>().operator->();

  // base folder for the contents of this job
  monitorName_ = parameters_.getUntrackedParameter<string>("monitorName","YourSubsystemName");
  cout << "DQMSourceExample: Monitor name = " << monitorName_ << endl;
  if (monitorName_ != "" ) monitorName_ = monitorName_+"/" ;

  prescaleLS_  = parameters_.getUntrackedParameter<int>("prescaleLS",  -1);
  cout << "DQMSourceExample: DQM lumi section prescale = " << prescaleLS_ << " lumi section(s)"<< endl;

  prescaleEvt_ = parameters_.getUntrackedParameter<int>("prescaleEvt", -1);
  cout << "DQMSourceExample: DQM event prescale = " << prescaleEvt_ << " events(s)"<< endl;

}

// ---------------------------------------------------------
void DQMSourceExample::beginJob(const EventSetup& context) {

  // get back-end interface
  dbe_ = Service<DQMStore>().operator->();

  //  create summ and cd into new folder
  dbe_->setCurrentFolder(monitorName_+"DQMsource/Summary");
  summ = dbe_->book1D("summary", "Run Summary", 100, 0, 100); 

  //-------------------------------------
  // testing of Quality Tests 
  //-------------------------------------

   // create and cd into new folder
   dbe_->setCurrentFolder(monitorName_+"DQMsource/QTests");

   // define histogram binning
   NBINS = 40 ; XMIN  =  0.; XMAX  = 40.;

  // book histograms for testsing of quality tests .
  xTrue     = dbe_->book1D("XTrue",       "X Range QTest",                  NBINS, XMIN, XMAX);
  xFalse    = dbe_->book1D("XFalse",      "X Range QTest",                  NBINS, XMIN, XMAX);
  yTrue     = dbe_->book1D("YTrue",       "Y Range QTest",                  NBINS, XMIN, XMAX);
  yFalse    = dbe_->book1D("YFalse",      "Y Range QTest",                  NBINS, XMIN, XMAX);
  wExpTrue  = dbe_->book2D("WExpTrue",    "Contents Within Expected QTest", NBINS, XMIN, XMAX, NBINS, XMIN, XMAX);
  wExpFalse = dbe_->book2D("WExpFalse",   "Contents Within Expected QTest", NBINS, XMIN, XMAX, NBINS, XMIN, XMAX);
  meanTrue  = dbe_->book1D("MeanTrue",    "Mean Within Expected QTest",     NBINS, XMIN, XMAX);
  meanFalse = dbe_->book1D("MeanFalse",   "Mean Within Expected QTest",     NBINS, XMIN, XMAX);
  deadTrue  = dbe_->book1D("DeadTrue",    "Dead Channel QTest",             NBINS, XMIN, XMAX);
  deadFalse = dbe_->book1D("DeadFalse",   "Dead Channel QTest",             NBINS, XMIN, XMAX);
  noisyTrue  = dbe_->book1D("NoisyTrue",  "Noisy Channel QTest",            NBINS, XMIN, XMAX);
  noisyFalse = dbe_->book1D("NoisyFalse", "Noisy Channel QTest",            NBINS, XMIN, XMAX);


  //-------------------------------------
  // book several ME more  
  //-------------------------------------

  //  create and cd into new folder
  dbe_->setCurrentFolder(monitorName_+"DQMsource/C1");
  const int NBINS2 = 10;
 
  i1        = dbe_->bookInt("int1");
  f1        = dbe_->bookFloat("float1");
  s1        = dbe_->bookString("s1", "My string");
  h1        = dbe_->book1D("histo1", "Example 1 1D histogram.", NBINS2, XMIN, XMAX);
  p1        = dbe_->bookProfile(  "prof1", "My profile 1D", NBINS,XMIN,XMAX,NBINS,XMIN,XMAX,"");
  p2        = dbe_->bookProfile2D("prof2", "My profile 2D", NBINS,XMIN,XMAX,NBINS,XMIN,XMAX,NBINS,XMIN,XMAX,"");
 
  // set labels for h1
  char temp[1024];
  for(int i = 1; i <= NBINS2; ++i) {
    sprintf(temp, " bin no. %d", i);
    h1->setBinLabel(i, temp);
  }

  // assign tag to MEs h1
  const unsigned int detector_id = 17;
  dbe_->tag(h1, detector_id);

  // tag full directory
  dbe_->tagContents(monitorName_+"DQMsource/C1", detector_id);

  /*
  // contents of h5 & h6 will be reset at end of monitoring cycle
  h5->setResetMe(true);
  h6->setResetMe(true);
  dbe_->showDirStructure();
  std::vector<std::string> tags;
  dbe_->getAllTags(tags);
  for (size_t i = 0, e = tags.size(); i < e; ++i)
    std::cout << "TAGS [" << i << "] = " << tags[i] << std::endl;
  */

}


// ----------------------------------------------------------------------------
void DQMSourceExample::beginRun(const edm::Run& r, const EventSetup& context) {
}


// ------------------------------------------------------------------------
void DQMSourceExample::beginLuminosityBlock(const LuminosityBlock& lumiSeg,
					    const EventSetup& context) {
}


// ----------------------------------------------------------------------------
void DQMSourceExample::analyze(const Event& iEvent, const EventSetup& iSetup) {
  counterEvt_++;
  if (prescaleEvt_<1)  return;
  if (prescaleEvt_ > 0 && counterEvt_%prescaleEvt_!=0)  return;
  //  cout << " processing conterEvt_: " << counterEvt_ <<endl;

  // fill integer and float
  i1->Fill(4);
  f1->Fill(-3.14);
 
  //----------------------------------------
  // Filling the histograms with random data
  //----------------------------------------

  srand( 0 );
  // fill summ histo
  if(counterEvt_%1000 == 0) {
    cout << " # of events = " << counterEvt_ << endl;
    summ->Fill(counterEvt_/1000., counterEvt_);
  }

  float z  = gRandom->Uniform(XMAX);
  xTrue->Fill(  z, 1./log(z+1.) );
  xFalse->Fill( z+(XMAX/2.),  z );
  yTrue->Fill(  z, 1./log(z+1.) );
  yFalse->Fill( z, z );
  meanTrue->Fill(  gRandom->Gaus(10,  2), 1.);
  meanFalse->Fill( gRandom->Gaus(12,  3), 1.);
  wExpTrue->Fill(  gRandom->Gaus(12,  1), gRandom->Gaus(12, 1), 1.);
  wExpFalse->Fill( gRandom->Gaus(20,  2), gRandom->Gaus(20, 2), 1.);
  deadTrue->Fill(  gRandom->Gaus(20, 10), 2.);
  deadFalse->Fill( gRandom->Gaus(20,  4), 1.);

  for ( int i = 0; i != 10; ++i ) {
    float w = gRandom->Uniform(XMAX);
    noisyTrue->Fill(  w, 1.);
    noisyFalse->Fill( z, 1.);
    float x = gRandom->Gaus(12, 1);
    float y = gRandom->Gaus(20, 2);
    p1->Fill(x, y);
    p2->Fill(x, y, (x+y)/2.);
    h1->Fill(y, 1.);
  }

  // usleep(100);

}


// ----------------------------------------------------------------------

void DQMSourceExample::endLuminosityBlock(const LuminosityBlock& lumiSeg,
					  const EventSetup& context) {

}


// ---------------------------------------------------------------------

void DQMSourceExample::endRun(const Run& r, const EventSetup& context) {

}


// ------------------------------

void DQMSourceExample::endJob() {
 
}
