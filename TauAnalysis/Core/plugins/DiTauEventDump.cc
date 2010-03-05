#include "TauAnalysis/Core/plugins/DiTauEventDump.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TauAnalysis/Core/interface/eventDumpAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

#include <TMath.h>

DiTauEventDump::DiTauEventDump(const edm::ParameterSet& cfg)
  : GenericEventDump(cfg)
{
  //std::cout << "<DiTauEventDump::DiTauEventDump>:" << std::endl;

  electronDump_ = makeObjectDump(cfg, "PATElectronDump");
  muonDump_ = makeObjectDump(cfg, "PATMuonDump");
  tauDump_ = makeObjectDump(cfg, "PATTauDump");  
  diTauDump_ = makeObjectDump(cfg, "PADiTauPairDump");
}

DiTauEventDump::~DiTauEventDump()
{
//--- nothing to be done yet...
}

void DiTauEventDump::print(const edm::Event& iEvent, const edm::EventSetup& iSetup, 
			   const std::vector<std::pair<std::string, bool> >& filterResults_cumulative, 
			   const std::vector<std::pair<std::string, bool> >& filterResults_individual, 
			   double eventWeight) const
{  
  //std::cout << "<DiTauEventDump::print>:" << std::endl; 

  printEventHeaderInfo(iEvent, eventWeight);

  printEventTriggerInfo(iEvent);

  printEventSelectionInfo(filterResults_cumulative, filterResults_individual, outputStream_);

  *outputStream_ << ">>GENERATOR LEVEL INFORMATION<<" << std::endl;

  edm::Handle<edm::View<reco::GenParticle> > genParticleCollection;
  iEvent.getByLabel(genParticleSource_, genParticleCollection);
  edm::Handle<edm::View<reco::GenJet> > genTauJetCollection;
  iEvent.getByLabel(genTauJetSource_, genTauJetCollection);
  printGenParticleInfo(genParticleCollection, genTauJetCollection, outputStream_);

  *outputStream_ << ">>RECONSTRUCTION LEVEL INFORMATION<<" << std::endl;

  electronDump_->print(iEvent, iSetup);
  muonDump_->print(iEvent, iSetup);
  tauDump_->print(iEvent, iSetup);

  diTauDump_->print(iEvent, iSetup);

  printMissingEtInfo(iEvent);

  printJetInfo(iEvent);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(EventDumpPluginFactory, DiTauEventDump, "DiTauEventDump");
