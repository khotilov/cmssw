#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "FastSimDataFormats/PileUpEvents/interface/PUEvent.h"

#include "FastSimulation/PileUpProducer/plugins/PileUpProducer.h"
#include "FastSimulation/Event/interface/GaussianPrimaryVertexGenerator.h"
#include "FastSimulation/Event/interface/FlatPrimaryVertexGenerator.h"
#include "FastSimulation/Event/interface/NoPrimaryVertexGenerator.h"
#include "FastSimulation/Utilities/interface/RandomEngine.h"

#include "HepMC/GenEvent.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"

#include <iostream>
#include <memory>
#include <sys/stat.h>
#include <cmath>

PileUpProducer::PileUpProducer(edm::ParameterSet const & p)  
{    

  // This producer produces a HepMCProduct, with all pileup vertices/particles
  produces<edm::HepMCProduct>("PileUpEvents");
  
  // Initialize the random number generator service
  edm::Service<edm::RandomNumberGenerator> rng;
  if ( ! rng.isAvailable() ) {
    throw cms::Exception("Configuration")
      << "PileUpProducer requires the RandomGeneratorService\n"
         "which is not present in the configuration file.\n"
         "You must add the service in the configuration file\n"
         "or remove the module that requires it";
  }

  // TRandom3 or CLHEP?
  bool m_TRandom = p.getParameter<bool>("UseTRandomEngine");
  if ( !m_TRandom ) { 
    random = new RandomEngine(&(*rng));
  } else {
    TRandom3* anEngine = new TRandom3();
    anEngine->SetSeed(rng->mySeed());
    random = new RandomEngine(anEngine);
  }

  // The pile-up event generation condition
  const edm::ParameterSet& pu = p.getParameter<edm::ParameterSet>("PileUpSimulator");
  averageNumber_ = pu.getParameter<double>("averageNumber");
  theFileNames = pu.getUntrackedParameter<std::vector<std::string> >("fileNames");
  inputFile = pu.getUntrackedParameter<std::string>("inputFile");
  theNumberOfFiles = theFileNames.size();
  theFiles.resize(theNumberOfFiles);
  theTrees.resize(theNumberOfFiles);
  theBranches.resize(theNumberOfFiles);
  thePUEvents.resize(theNumberOfFiles);
  theCurrentEntry.resize(theNumberOfFiles);
  theCurrentMinBiasEvt.resize(theNumberOfFiles);
  theNumberOfEntries.resize(theNumberOfFiles);
  theNumberOfMinBiasEvts.resize(theNumberOfFiles);

  // Initialize the primary vertex generator
  const edm::ParameterSet& vtx = p.getParameter<edm::ParameterSet>("VertexGenerator");
  std::string vtxType = vtx.getParameter<std::string>("type");
  if ( vtxType == "Gaussian" ) 
    theVertexGenerator = new GaussianPrimaryVertexGenerator(vtx,random);
  else if ( vtxType == "Flat" ) 
    theVertexGenerator = new FlatPrimaryVertexGenerator(vtx,random);
  else
    theVertexGenerator = new NoPrimaryVertexGenerator();

}

PileUpProducer::~PileUpProducer() { 

  delete theVertexGenerator;

}

void PileUpProducer::beginJob(const edm::EventSetup & es)
{
  
  std::cout << " PileUpProducer initializing " << std::endl;
  gROOT->cd();
  
  std::string fullPath;
  
  // Read the information from a previous run (to keep reproducibility)
  bool input = this->read(inputFile);
  if ( input ) 
    std::cout << "***WARNING*** You are reading pile-up information from the file "
	      << inputFile << " created in an earlier run."
	      << std::endl;

  // Open the file for saving the information of the current run
  myOutputFile.open ("PileUpOutputFile.txt");
  myOutputBuffer = 0;

  // Open the root files
  std::cout << "Opening minimum-bias event files ... " << std::endl;
  for ( unsigned file=0; file<theNumberOfFiles; ++file ) {

    edm::FileInPath myDataFile("FastSimulation/PileUpProducer/data/"+theFileNames[file]);
    fullPath = myDataFile.fullPath();
    //    theFiles[file] = TFile::Open(theFileNames[file].c_str());
    theFiles[file] = TFile::Open(fullPath.c_str());
    if ( !theFiles[file] ) throw cms::Exception("FastSimulation/PileUpProducer") 
      << "File " << theFileNames[file] << " " << fullPath <<  " not found ";
    //
    theTrees[file] = (TTree*) theFiles[file]->Get("MinBiasEvents"); 
    if ( !theTrees[file] ) throw cms::Exception("FastSimulation/PileUpProducer") 
      << "Tree with name MinBiasEvents not found in " << theFileNames[file];
    //
    theBranches[file] = theTrees[file]->GetBranch("puEvent");
    if ( !theBranches[file] ) throw cms::Exception("FastSimulation/PileUpProducer") 
      << "Branch with name puEvent not found in " << theFileNames[file];
    //
    thePUEvents[file] = new PUEvent();
    theBranches[file]->SetAddress(&thePUEvents[file]);
    //
    theNumberOfEntries[file] = theTrees[file]->GetEntries();
    // Add some randomness (if there was no input file)
    if ( !input ) 
      theCurrentEntry[file] 
	= (unsigned) (theNumberOfEntries[file] * random->flatShoot());

    theTrees[file]->GetEntry(theCurrentEntry[file]);
    unsigned NMinBias = thePUEvents[file]->nMinBias();
    theNumberOfMinBiasEvts[file] = NMinBias;
    // Add some randomness (if there was no input file)
    if ( !input )
	theCurrentMinBiasEvt[file] = 
	  (unsigned) (theNumberOfMinBiasEvts[file] * random->flatShoot());
    
    /*
    std::cout << "File " << theFileNames[file]
	      << " is opened with " << theNumberOfEntries[file] 
	      << " entries and will be read from Entry/Event "
	      << theCurrentEntry[file] << "/" << theCurrentMinBiasEvt[file]
	      << std::endl;
    */
  }
  
  // Return Loot in the same state as it was when entering. 
  gROOT->cd();
  
}
 
void PileUpProducer::endJob()
{ 
    std::cout << " PileUpProducer terminating " << std::endl; 
  // Close all local files
  // Among other things, this allows the TROOT destructor to end up 
  // without crashing, while trying to close these files from outside
  std::cout << "Closing minimum-bias event files... " << std::endl;
  for ( unsigned file=0; file<theFiles.size(); ++file ) {
    
    // std::cout << "Closing " << theFileNames[file] << std::endl;
    theFiles[file]->Close();
    
  }
  
  // Close the output file
  myOutputFile.close();
  
  // And return Loot in the same state as it was when entering. 
  gROOT->cd();
  
}
 
void PileUpProducer::produce(edm::Event & iEvent, const edm::EventSetup & es)
{

  // Create the GenEvent and the HepMCProduct
  std::auto_ptr<edm::HepMCProduct> pu_product(new edm::HepMCProduct());  
  HepMC::GenEvent* evt = new HepMC::GenEvent();
  
  // How many pile-up events?
  int PUevts = (int) random->poissonShoot(averageNumber_);

  // Get N events from random files
  for ( int ievt=0; ievt<PUevts; ++ievt ) { 
    
    // Draw a file in a ramdom manner 
    unsigned file = (unsigned) (theNumberOfFiles * random->flatShoot());
    /*
    if ( debug )  
      std::cout << "The file chosen for event " << ievt 
		<< " is the file number " << file << std::endl; 
    */

    // Smear the primary vertex and express it in mm (stupid GenEvent convention...)
    theVertexGenerator->generate();
    HepMC::FourVector smearedVertex =  
      HepMC::FourVector(theVertexGenerator->X()*10.,
			theVertexGenerator->Y()*10.,
			theVertexGenerator->Z()*10.,
			0.);
    HepMC::GenVertex* aVertex = new HepMC::GenVertex(smearedVertex);
    evt->add_vertex(aVertex);

    // Some rotation around the z axis, for more randomness
    double theAngle = random->flatShoot() * 2. * 3.14159265358979323;
    double cAngle = std::cos(theAngle);
    double sAngle = std::sin(theAngle);

    /*
    if ( debug ) 
      std::cout << "File chosen : " << file 
		<< " Current entry in this file " << theCurrentEntry[file] 
		<< " Current minbias in this chunk= " << theCurrentMinBiasEvt[file] 
		<< " Total number of minbias in this chunk = " << theNumberOfMinBiasEvts[file] << std::endl;
    */

    //      theFiles[file]->cd();
    //      gDirectory->ls();
    // Check we are not either at the end of a minbias bunch 
    // or at the end of a file
    if ( theCurrentMinBiasEvt[file] == theNumberOfMinBiasEvts[file] ) {
      // if ( debug ) std::cout << "End of MinBias bunch ! ";
      ++theCurrentEntry[file];
      // if ( debug) std::cout << "Read the next entry " << theCurrentEntry[file] << std::endl;
      theCurrentMinBiasEvt[file] = 0;
      if ( theCurrentEntry[file] == theNumberOfEntries[file] ) { 
	theCurrentEntry[file] = 0;
	// if ( debug ) std::cout << "End of file - Rewind! " << std::endl;
      }
      // if ( debug ) std::cout << "The PUEvent is reset ... "; 
      thePUEvents[file]->reset();
      unsigned myEntry = theCurrentEntry[file];
      /* 
      if ( debug ) std::cout << "The new entry " << myEntry 
			     << " is read ... in TTree " << theTrees[file] << " "; 
      */
      theTrees[file]->GetEntry(myEntry);
      /*
      if ( debug ) 
	std::cout << "The number of interactions in the new entry is ... "; 	
      */
      theNumberOfMinBiasEvts[file] = thePUEvents[file]->nMinBias();
      // if ( debug ) std::cout << theNumberOfMinBiasEvts[file] << std::endl;
  }
  
    // Read a minbias event chunk
    const PUEvent::PUMinBiasEvt& aMinBiasEvt 
      = thePUEvents[file]->thePUMinBiasEvts()[theCurrentMinBiasEvt[file]];
  
    // Find corresponding particles
    unsigned firstTrack = aMinBiasEvt.first; 
    unsigned trackSize = firstTrack + aMinBiasEvt.size;
    /*
    if ( debug ) std::cout << "First and last+1 tracks are " 
			   << firstTrack << " " << trackSize << std::endl;
    */

    // Loop on particles
    for ( unsigned iTrack=firstTrack; iTrack<trackSize; ++iTrack ) {
      
      const PUEvent::PUParticle& aParticle 
	= thePUEvents[file]->thePUParticles()[iTrack];
      /*
      if ( debug) 
	std::cout << "Track " << iTrack 
		  << " id/px/py/pz/mass "
		  << aParticle.id << " " 
		  << aParticle.px << " " 
		  << aParticle.py << " " 
		  << aParticle.pz << " " 
		  << aParticle.mass << " " << std::endl; 
      */
      
      // Create a FourVector, with rotation 
      double energy = std::sqrt( aParticle.px*aParticle.px
			       + aParticle.py*aParticle.py
			       + aParticle.pz*aParticle.pz
			       + aParticle.mass*aParticle.mass );

      HepMC::FourVector myPart(cAngle * aParticle.px + sAngle * aParticle.py,
			      -sAngle * aParticle.px + cAngle * aParticle.py,
			       aParticle.pz, energy);

      // Add a GenParticle
      HepMC::GenParticle* aGenParticle = new HepMC::GenParticle(myPart,aParticle.id);
      aVertex->add_particle_out(aGenParticle);

    }
    // End of particle loop
    
    // Increment for next time
    ++theCurrentMinBiasEvt[file];
    
  }
  // End of pile-up event loop

  // evt->print();

  // Fill the HepMCProduct from the GenEvent
  if ( evt )  pu_product->addHepMCData( evt );

  // Put the HepMCProduct onto the event
  iEvent.put(pu_product,"PileUpEvents");
  // delete evt;

  // Save the current location in each pile-up event files
  this->save();

}

void
PileUpProducer::save() {

  // Size of buffer
  ++myOutputBuffer;

  // Periodically close the current file and open a new one
  if ( myOutputBuffer/1000*1000 == myOutputBuffer ) { 
    myOutputFile.close();
    myOutputFile.open ("PileUpOutputFile.txt");
    //    myOutputFile.seekp(0); // No need to rewind in that case
  }

  // Save the current position to file
  myOutputFile.write((const char*)(&theCurrentEntry.front()),
		     theCurrentEntry.size()*sizeof(unsigned));
  myOutputFile.write((const char*)&theCurrentMinBiasEvt.front(),
		     theCurrentMinBiasEvt.size()*sizeof(unsigned));
  myOutputFile.flush();

}

bool
PileUpProducer::read(std::string inputFile) {

  ifstream myInputFile;
  struct stat results;
  unsigned size1 = theCurrentEntry.size()*sizeof(unsigned);
  unsigned size2 = theCurrentMinBiasEvt.size()*sizeof(unsigned);
  unsigned size = 0;


  // Open the file (if any)
  myInputFile.open (inputFile.c_str());
  if ( myInputFile.is_open() ) { 

    // Get the size of the file
    if ( stat(inputFile.c_str(), &results) == 0 ) size = results.st_size;
    else return false; // Something is wrong with that file !
  
    // Position the pointer just before the last record
    myInputFile.seekg(size-size1-size2);

    // Read the information
    myInputFile.read((char*)(&theCurrentEntry.front()),size1);
    myInputFile.read((char*)&theCurrentMinBiasEvt.front(),size2);
    myInputFile.close();

    return true;

  } 

  return false;

}

DEFINE_FWK_MODULE(PileUpProducer);
