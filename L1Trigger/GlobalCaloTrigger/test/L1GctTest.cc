#include "L1Trigger/GlobalCaloTrigger/test/L1GctTest.h"

#include "L1Trigger/GlobalCaloTrigger/test/gctTestFunctions.h"

#include "FWCore/Framework/interface/ESHandle.h"

// Trigger configuration includes
#include "CondFormats/L1TObjects/interface/L1GctJetEtCalibrationFunction.h"
#include "CondFormats/DataRecord/interface/L1GctJetCalibFunRcd.h"

// GCT include files
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetEtCalibrationLut.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GlobalCaloTrigger.h"

using std::cout;
using std::endl;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1GctTest::L1GctTest(const edm::ParameterSet& iConfig) :
  theElectronTestIsEnabled   (iConfig.getUntrackedParameter<bool>("doElectrons",   false)),
  theEnergyAlgosTestIsEnabled(iConfig.getUntrackedParameter<bool>("doEnergyAlgos", false)),
  theFirmwareTestIsEnabled   (iConfig.getUntrackedParameter<bool>("doFirmware",    false)),
  theInputDataFileName       (iConfig.getUntrackedParameter<std::string>("inputFile",     "")),
  theReferenceDataFileName   (iConfig.getUntrackedParameter<std::string>("referenceFile", "")),
  m_eventNo(0)
{
  //now do what ever initialization is needed
  // check the files are specified if required
  if (theElectronTestIsEnabled && theInputDataFileName=="") {
    throw cms::Exception ("L1GctTestInitialisationError")
      << "no input filename provided for electron tests.\n"
      << "Specify non-blank parameter inputFile in cmsRun configuration\n"; }
  if (theFirmwareTestIsEnabled && theInputDataFileName=="") {
    throw cms::Exception ("L1GctTestInitialisationError")
      << "no input filename provided for firmware tests.\n" 
      << "Specify non-blank parameter inputFile in cmsRun configuration\n"; }
  if (theFirmwareTestIsEnabled && theReferenceDataFileName=="") {
    throw cms::Exception ("L1GctTestInitialisationError")
      << "no reference filename provided for firmware tests.\n"
      << "Specify non-blank parameter referenceFile in cmsRun configuration\n"; }

  // instantiate the GCT
  m_gct = new L1GlobalCaloTrigger(false,L1GctJetLeafCard::hardwareJetFinder);
  m_tester = new gctTestFunctions();

}


L1GctTest::~L1GctTest()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

   delete m_gct;
   delete m_tester;

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
L1GctTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  bool endOfFile;

   // Initialise the gct
   m_gct->reset();

   // Load data into the gct according to the tests to be carried out
   if (theElectronTestIsEnabled) {
     m_tester->loadNextEvent(m_gct,theInputDataFileName); }

   if (theEnergyAlgosTestIsEnabled) {
     m_tester->loadNextEvent(m_gct, (100>m_eventNo)); }

   if (theFirmwareTestIsEnabled) {
     m_tester->loadNextEvent(m_gct, theInputDataFileName, endOfFile); }

   // Run the gct emulator on the input data
   m_gct->process();

   bool passAllTests = true;

   // Check the results of the emulator
   if (theElectronTestIsEnabled) {
     m_tester->fillElectronData(m_gct);
     passAllTests &= m_tester->checkElectrons(m_gct);
   }

   if (theFirmwareTestIsEnabled) {
     m_tester->fillJetsFromFirmware(theReferenceDataFileName);
     passAllTests &= m_tester->checkJetFinder(m_gct);
   }

   if (theEnergyAlgosTestIsEnabled || theFirmwareTestIsEnabled) {
     m_tester->fillRawJetData(m_gct);
     passAllTests &= m_tester->checkEnergySums(m_gct);
     passAllTests &= m_tester->checkHtSums(m_gct);
     passAllTests &= m_tester->checkJetCounts(m_gct);
   }

   m_eventNo++;
   if (theFirmwareTestIsEnabled && endOfFile) {
     std::cout << "Reached the end of input file after " << m_eventNo << " events\n";
   }
   theFirmwareTestIsEnabled &= !endOfFile;

   // bale out if we fail any test
   if (passAllTests)
   {
      //std::cout << "All tests passed for this event!" << std::endl;
   } else {
      throw cms::Exception("L1GctTestError") << "\ntest failed\n\n";
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
L1GctTest::beginJob(const edm::EventSetup& c)
{
  configureGct(c);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1GctTest::endJob() {
  std::cout << "\n\n=== All tests passed Ok! ===\n\n" << std::endl;
}

void 
L1GctTest::configureGct(const edm::EventSetup& c)
{
  assert(&c!=0);

  // get data from EventSetup
  edm::ESHandle< L1GctJetEtCalibrationFunction > calibFun ;
  c.get< L1GctJetCalibFunRcd >().get( calibFun ) ; // which record?

  if (calibFun.product() == 0) {
    throw cms::Exception("L1GctConfigError")
      << "Failed to find a L1GctJetCalibFunRcd:L1GctJetEtCalibrationFunction in EventSetup!" << endl
      << "Cannot continue without this function" << endl;
  }

  // make a jet Et Lut and tell it about the scales
  m_jetEtCalibLut = L1GctJetEtCalibrationLut::setupLut(calibFun.product());
  m_gct->setJetEtCalibrationLut(m_jetEtCalibLut);

}

