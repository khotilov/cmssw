#include "JetMETCorrections/METPUSubtraction/plugins/GBRForestWriter.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

#include "CondFormats/EgammaObjects/interface/GBRForest.h"

#include <TFile.h>

GBRForestWriter::GBRForestWriter(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  edm::VParameterSet cfgJobs = cfg.getParameter<edm::VParameterSet>("jobs");
  for ( edm::VParameterSet::const_iterator cfgJob = cfgJobs.begin();
	cfgJob != cfgJobs.end(); ++cfgJob ) {
    jobEntryType* job = new jobEntryType();
    job->inputFileName_ = cfgJob->getParameter<edm::FileInPath>("inputFileName");
    job->gbrForestName_ = cfgJob->getParameter<std::string>("gbrForestName");
    job->outputRecord_ = cfgJob->getParameter<std::string>("outputRecord");
    jobs_.push_back(job);
  }
}

GBRForestWriter::~GBRForestWriter()
{
  for ( std::vector<jobEntryType*>::iterator it = jobs_.begin();
	it != jobs_.end(); ++it ) {
    delete (*it);
  }
}

void GBRForestWriter::analyze(const edm::Event&, const edm::EventSetup&)
{
  std::cout << "<GBRForestWriter::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;

  edm::Service<cond::service::PoolDBOutputService> dbService;
  if ( !dbService.isAvailable() ) 
    throw cms::Exception("GBRForestWriter") 
      << " Failed to access PoolDBOutputService !!\n";
  
  for ( std::vector<jobEntryType*>::iterator job = jobs_.begin();
	job != jobs_.end(); ++job ) {
    if ( !(*job)->inputFileName_.isLocal() ) 
      throw cms::Exception("GBRForestWriter") 
	<< " Failed to find File = " << (*job)->inputFileName_ << " !!\n";
    TFile* inputFile = new TFile((*job)->inputFileName_.fullPath().data());
    //const GBRForest* gbrForest = dynamic_cast<GBRForest*>(inputFile->Get((*job)->gbrForestName_.data())); // CV: dynamic_cast<GBRForest*> fails for some reason ?!
    const GBRForest* gbrForest = (GBRForest*)inputFile->Get((*job)->gbrForestName_.data());
    if ( !gbrForest ) 
      throw cms::Exception("GBRForestWriter") 
	<< " Failed to load GBRForest = " << (*job)->gbrForestName_.data() << " from file = " << (*job)->inputFileName_.fullPath().data() << " !!\n";
    std::cout << " writing " << (*job)->outputRecord_ << "..." << std::endl;
    dbService->writeOne(gbrForest, dbService->beginOfTime(), (*job)->outputRecord_);
    delete inputFile;
  }

  std::cout << "done." << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GBRForestWriter);
