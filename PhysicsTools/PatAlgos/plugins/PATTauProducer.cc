//
// $Id: PATTauProducer.cc,v 1.32.16.1 2010/02/22 10:37:18 mbluj Exp $
//

#include "PhysicsTools/PatAlgos/plugins/PATTauProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/TauReco/interface/PFTauDecayModeAssociation.h"
#include "DataFormats/TauReco/interface/CaloTau.h"
#include "DataFormats/TauReco/interface/CaloTauDiscriminator.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <vector>
#include <memory>

using namespace pat;

PATTauProducer::PATTauProducer(const edm::ParameterSet & iConfig):
  isolator_(iConfig.exists("userIsolation") ? iConfig.getParameter<edm::ParameterSet>("userIsolation") : edm::ParameterSet(), false) ,
  useUserData_(iConfig.exists("userData"))
{
  // initialize the configurables
  tauSrc_               = iConfig.getParameter<edm::InputTag>( "tauSource" );

  embedIsolationTracks_ = iConfig.getParameter<bool>         ( "embedIsolationTracks" );
  embedLeadTrack_       = iConfig.getParameter<bool>         ( "embedLeadTrack" );
  embedSignalTracks_    = iConfig.getParameter<bool>         ( "embedSignalTracks" );
  embedLeadPFCand_ = iConfig.getParameter<bool>( "embedLeadPFCand" );
  embedLeadPFChargedHadrCand_ = iConfig.getParameter<bool>( "embedLeadPFChargedHadrCand" );
  embedLeadPFNeutralCand_ = iConfig.getParameter<bool>( "embedLeadPFNeutralCand" );
  embedSignalPFCands_ = iConfig.getParameter<bool>( "embedSignalPFCands" );
  embedSignalPFChargedHadrCands_ = iConfig.getParameter<bool>( "embedSignalPFChargedHadrCands" );
  embedSignalPFNeutralHadrCands_ = iConfig.getParameter<bool>( "embedSignalPFNeutralHadrCands" );
  embedSignalPFGammaCands_ = iConfig.getParameter<bool>( "embedSignalPFGammaCands" );
  embedIsolationPFCands_ = iConfig.getParameter<bool>( "embedIsolationPFCands" );
  embedIsolationPFChargedHadrCands_ = iConfig.getParameter<bool>( "embedIsolationPFChargedHadrCands" );
  embedIsolationPFNeutralHadrCands_ = iConfig.getParameter<bool>( "embedIsolationPFNeutralHadrCands" );
  embedIsolationPFGammaCands_ = iConfig.getParameter<bool>( "embedIsolationPFGammaCands" );

  addGenMatch_    = iConfig.getParameter<bool>         ( "addGenMatch" );

  if (addGenMatch_) {
      embedGenMatch_ = iConfig.getParameter<bool>         ( "embedGenMatch" );
      if (iConfig.existsAs<edm::InputTag>("genParticleMatch")) {
          genMatchSrc_.push_back(iConfig.getParameter<edm::InputTag>( "genParticleMatch" ));
      } else {
          genMatchSrc_ = iConfig.getParameter<std::vector<edm::InputTag> >( "genParticleMatch" );
      }
  }

  addGenJetMatch_    = iConfig.getParameter<bool>         ( "addGenJetMatch" );
  if(addGenJetMatch_) {
    embedGenJetMatch_  = iConfig.getParameter<bool>         ( "embedGenJetMatch" );
    genJetMatchSrc_    = iConfig.getParameter<edm::InputTag>( "genJetMatch" );
  }


  // tau ID configurables
  addTauID_       = iConfig.getParameter<bool>         ( "addTauID" );
  if ( addTauID_ ) {
    // it might be a single tau ID
    if (iConfig.existsAs<edm::InputTag>("tauIDSource")) {
      tauIDSrcs_.push_back(NameTag("", iConfig.getParameter<edm::InputTag>("tauIDSource")));
    }
    // or there might be many of them
    if (iConfig.existsAs<edm::ParameterSet>("tauIDSources")) {
      // please don't configure me twice
      if (!tauIDSrcs_.empty()) throw cms::Exception("Configuration") << 
	"PATTauProducer: you can't specify both 'tauIDSource' and 'tauIDSources'\n";
      // read the different tau ID names
      edm::ParameterSet idps = iConfig.getParameter<edm::ParameterSet>("tauIDSources");
      std::vector<std::string> names = idps.getParameterNamesForType<edm::InputTag>();
      for (std::vector<std::string>::const_iterator it = names.begin(), ed = names.end(); it != ed; ++it) {
	tauIDSrcs_.push_back(NameTag(*it, idps.getParameter<edm::InputTag>(*it)));
      }
    }
    // but in any case at least once
    if (tauIDSrcs_.empty()) throw cms::Exception("Configuration") <<
      "PATTauProducer: id addTauID is true, you must specify either:\n" <<
      "\tInputTag tauIDSource = <someTag>\n" << "or\n" <<
      "\tPSet tauIDSources = { \n" <<
      "\t\tInputTag <someName> = <someTag>   // as many as you want \n " <<
      "\t}\n";
  }
  
  // tau decay mode configurables
  addDecayMode_ = iConfig.getParameter<bool>         ( "addDecayMode" );
  if ( addDecayMode_ ) {
    decayModeSrc_ = iConfig.getParameter<edm::InputTag>( "decayModeSrc" );
  }

  // IsoDeposit configurables
  if (iConfig.exists("isoDeposits")) {
    edm::ParameterSet depconf = iConfig.getParameter<edm::ParameterSet>("isoDeposits");
    if ( depconf.exists("tracker")         ) isoDepositLabels_.push_back(std::make_pair(pat::TrackIso, depconf.getParameter<edm::InputTag>("tracker")));
    if ( depconf.exists("ecal")            ) isoDepositLabels_.push_back(std::make_pair(pat::EcalIso, depconf.getParameter<edm::InputTag>("ecal")));
    if ( depconf.exists("hcal")            ) isoDepositLabels_.push_back(std::make_pair(pat::HcalIso, depconf.getParameter<edm::InputTag>("hcal")));
    if ( depconf.exists("pfAllParticles")  ) isoDepositLabels_.push_back(std::make_pair(pat::PfAllParticleIso, depconf.getParameter<edm::InputTag>("pfAllParticles")));
    if ( depconf.exists("pfChargedHadron") ) isoDepositLabels_.push_back(std::make_pair(pat::PfChargedHadronIso, depconf.getParameter<edm::InputTag>("pfChargedHadron")));
    if ( depconf.exists("pfNeutralHadron") ) isoDepositLabels_.push_back(std::make_pair(pat::PfNeutralHadronIso,depconf.getParameter<edm::InputTag>("pfNeutralHadron")));
    if ( depconf.exists("pfGamma")         ) isoDepositLabels_.push_back(std::make_pair(pat::PfGammaIso, depconf.getParameter<edm::InputTag>("pfGamma")));
    
    if ( depconf.exists("user") ) {
      std::vector<edm::InputTag> userdeps = depconf.getParameter<std::vector<edm::InputTag> >("user");
      std::vector<edm::InputTag>::const_iterator it = userdeps.begin(), ed = userdeps.end();
      int key = UserBaseIso;
      for ( ; it != ed; ++it, ++key) {
	isoDepositLabels_.push_back(std::make_pair(IsolationKeys(key), *it));
      }
    }
  }

  // Efficiency configurables
  addEfficiencies_ = iConfig.getParameter<bool>("addEfficiencies");
  if (addEfficiencies_) {
     efficiencyLoader_ = pat::helper::EfficiencyLoader(iConfig.getParameter<edm::ParameterSet>("efficiencies"));
  }

  // Resolution configurables
  addResolutions_ = iConfig.getParameter<bool>("addResolutions");
  if (addResolutions_) {
     resolutionLoader_ = pat::helper::KinResolutionsLoader(iConfig.getParameter<edm::ParameterSet>("resolutions"));
  }

  // Check to see if the user wants to add user data
  if ( useUserData_ ) {
    userDataHelper_ = PATUserDataHelper<Tau>(iConfig.getParameter<edm::ParameterSet>("userData"));
  }

  // produces vector of taus
  produces<std::vector<Tau> >();
}

PATTauProducer::~PATTauProducer() 
{
}

void PATTauProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) 
{ 
  // Get the collection of taus from the event
  edm::Handle<edm::View<reco::BaseTau> > anyTaus;
  try {
    iEvent.getByLabel(tauSrc_, anyTaus);
  } catch (const edm::Exception &e) {
    edm::LogWarning("DataSource") << "WARNING! No Tau collection found. This missing input will not block the job. Instead, an empty tau collection is being be produced.";
    std::auto_ptr<std::vector<Tau> > patTaus(new std::vector<Tau>()); 
    iEvent.put(patTaus);
    return;
  }

  if (isolator_.enabled()) isolator_.beginEvent(iEvent,iSetup);

  if (efficiencyLoader_.enabled()) efficiencyLoader_.newEvent(iEvent);
  if (resolutionLoader_.enabled()) resolutionLoader_.newEvent(iEvent, iSetup);
   
  std::vector<edm::Handle<edm::ValueMap<IsoDeposit> > > deposits(isoDepositLabels_.size());
  for (size_t j = 0, nd = deposits.size(); j < nd; ++j) {
    iEvent.getByLabel(isoDepositLabels_[j].second, deposits[j]);
  }

  // prepare the MC matching
  std::vector<edm::Handle<edm::Association<reco::GenParticleCollection> > > genMatches(genMatchSrc_.size());
  if (addGenMatch_) {
    for (size_t j = 0, nd = genMatchSrc_.size(); j < nd; ++j) {
      iEvent.getByLabel(genMatchSrc_[j], genMatches[j]);
    }
  }

  edm::Handle<edm::Association<reco::GenJetCollection> > genJetMatch;  
  if (addGenJetMatch_) iEvent.getByLabel(genJetMatchSrc_, genJetMatch); 

  std::auto_ptr<std::vector<Tau> > patTaus(new std::vector<Tau>()); 

  for (size_t idx = 0, ntaus = anyTaus->size(); idx < ntaus; ++idx) {
    edm::RefToBase<reco::BaseTau> tausRef = anyTaus->refAt(idx);
    edm::Ptr<reco::BaseTau> tausPtr = anyTaus->ptrAt(idx);
    
    Tau aTau(tausRef);
    if (embedLeadTrack_)       aTau.embedLeadTrack();
    if (embedSignalTracks_)    aTau.embedSignalTracks();
    if (embedIsolationTracks_) aTau.embedIsolationTracks();    
    if (embedLeadPFCand_) {
      if (aTau.isPFTau() )
	aTau.embedLeadPFCand();
      else 
	edm::LogWarning("Type Error") << "Embedding a PFTau-specific information into a pat::Tau which wasn't made from a reco::PFTau is impossible.\n";
    }
    if (embedLeadPFChargedHadrCand_) {
      if (aTau.isPFTau() )
	aTau.embedLeadPFChargedHadrCand();
      else 
	edm::LogWarning("Type Error") << "Embedding a PFTau-specific information into a pat::Tau which wasn't made from a reco::PFTau is impossible.\n";
    }
    if (embedLeadPFNeutralCand_) {
      if (aTau.isPFTau() )
	aTau.embedLeadPFNeutralCand();
      else 
	edm::LogWarning("Type Error") << "Embedding a PFTau-specific information into a pat::Tau which wasn't made from a reco::PFTau is impossible.\n";
    }
    if (embedSignalPFCands_) {
      if (aTau.isPFTau() )
	aTau.embedSignalPFCands();
      else 
	edm::LogWarning("Type Error") << "Embedding a PFTau-specific information into a pat::Tau which wasn't made from a reco::PFTau is impossible.\n";
    }
    if (embedSignalPFChargedHadrCands_) {
      if (aTau.isPFTau() )
	aTau.embedSignalPFChargedHadrCands();
      else 
	edm::LogWarning("Type Error") << "Embedding a PFTau-specific information into a pat::Tau which wasn't made from a reco::PFTau is impossible.\n";
    }
    if (embedSignalPFNeutralHadrCands_) {
      if (aTau.isPFTau() )
	aTau.embedSignalPFNeutralHadrCands();
      else 
	edm::LogWarning("Type Error") << "Embedding a PFTau-specific information into a pat::Tau which wasn't made from a reco::PFTau is impossible.\n";
    }
    if (embedSignalPFGammaCands_) {
      if (aTau.isPFTau() )
	aTau.embedSignalPFGammaCands();
      else 
	edm::LogWarning("Type Error") << "Embedding a PFTau-specific information into a pat::Tau which wasn't made from a reco::PFTau is impossible.\n";
    }
    if (embedIsolationPFCands_) {
      if (aTau.isPFTau() )
	aTau.embedIsolationPFCands();
      else 
	edm::LogWarning("Type Error") << "Embedding a PFTau-specific information into a pat::Tau which wasn't made from a reco::PFTau is impossible.\n";
    }
    if (embedIsolationPFChargedHadrCands_) {
      if (aTau.isPFTau() )
	aTau.embedIsolationPFChargedHadrCands();
      else 
	edm::LogWarning("Type Error") << "Embedding a PFTau-specific information into a pat::Tau which wasn't made from a reco::PFTau is impossible.\n";
    }
    if (embedIsolationPFNeutralHadrCands_) {
      if (aTau.isPFTau() )
	aTau.embedIsolationPFNeutralHadrCands();
      else 
	edm::LogWarning("Type Error") << "Embedding a PFTau-specific information into a pat::Tau which wasn't made from a reco::PFTau is impossible.\n";
    }
    if (embedIsolationPFGammaCands_) {
      if (aTau.isPFTau() )
	aTau.embedIsolationPFGammaCands();
      else 
	edm::LogWarning("Type Error") << "Embedding a PFTau-specific information into a pat::Tau which wasn't made from a reco::PFTau is impossible.\n";
    }
    
    // store the match to the generated final state muons
    if (addGenMatch_) {
      for(size_t i = 0, n = genMatches.size(); i < n; ++i) {
          reco::GenParticleRef genTau = (*genMatches[i])[tausRef];
          aTau.addGenParticleRef(genTau);
      }
      if (embedGenMatch_) aTau.embedGenParticle();
    }
    
    // store the match to the visible part of the generated tau
    if (addGenJetMatch_) {
      reco::GenJetRef genJetTau = (*genJetMatch)[tausRef];
      if (genJetTau.isNonnull() && genJetTau.isAvailable() ) {
        aTau.setGenJet( genJetTau );
      } // leave empty if no match found
    }

    // prepare ID extraction 
    if ( addTauID_ ) {
      std::vector<pat::Tau::IdPair> ids(tauIDSrcs_.size());
      for ( size_t i = 0; i < tauIDSrcs_.size(); ++i ) {
	edm::Handle<reco::CaloTauDiscriminator> caloTauIdDiscr;
	iEvent.getByLabel(tauIDSrcs_[i].second, caloTauIdDiscr);

	edm::Handle<reco::PFTauDiscriminator> pfTauIdDiscr;
	iEvent.getByLabel(tauIDSrcs_[i].second, pfTauIdDiscr);
	
	if ( typeid(*tausRef) == typeid(reco::PFTau) ) {
	  //std::cout << "filling PFTauDiscriminator '" << tauIDSrcs_[i].first << "' into pat::Tau object..." << std::endl;
	  edm::Handle<reco::PFTauCollection> pfTauCollection; 
	  iEvent.getByLabel(tauSrc_, pfTauCollection);

	  edm::Handle<reco::PFTauDiscriminator> pfTauIdDiscr;
	  iEvent.getByLabel(tauIDSrcs_[i].second, pfTauIdDiscr);

	  ids[i].first = tauIDSrcs_[i].first;
	  ids[i].second = getTauIdDiscriminator(pfTauCollection, idx, pfTauIdDiscr);
	} else if ( typeid(*tausRef) == typeid(reco::CaloTau) ) {
	  //std::cout << "filling CaloTauDiscriminator '" << tauIDSrcs_[i].first << "' into pat::Tau object..." << std::endl;
	  edm::Handle<reco::CaloTauCollection> caloTauCollection; 
	  iEvent.getByLabel(tauSrc_, caloTauCollection);

	  edm::Handle<reco::CaloTauDiscriminator> caloTauIdDiscr;
	  iEvent.getByLabel(tauIDSrcs_[i].second, caloTauIdDiscr);

	  ids[i].first = tauIDSrcs_[i].first;
	  ids[i].second = getTauIdDiscriminator(caloTauCollection, idx, caloTauIdDiscr);
	} else {
	  throw cms::Exception("Type Mismatch") <<
	    "PATTauProducer: unsupported datatype '" << typeid(*tausRef).name() << "' for tauSource\n";
	}
      }

      aTau.setTauIDs(ids);
    }

    // extraction of reconstructed tau decay mode 
    if ( addDecayMode_ ) {
      edm::Handle<reco::PFTauDecayModeAssociation> pfDecayModeAssoc;
      iEvent.getByLabel(decayModeSrc_, pfDecayModeAssoc);
      edm::Handle<reco::PFTauCollection> pfTaus;
      iEvent.getByLabel(tauSrc_, pfTaus);
      reco::PFTauRef pfTauRef(pfTaus, idx);
      // need PFTauRef (edm::RefToBase<reco::BaseTau> does not suffice) 
      // for PFTauDecayMode look-up
      //const reco::PFTauDecayMode& pfDecayMode = (*pfDecayModeAssoc)[tausRef];
      const reco::PFTauDecayMode& pfDecayMode = (*pfDecayModeAssoc)[pfTauRef];
      aTau.setDecayMode(pfDecayMode.getDecayMode());
    }

    // Isolation
    if (isolator_.enabled()) {
      isolator_.fill(*anyTaus, idx, isolatorTmpStorage_);
      typedef pat::helper::MultiIsolator::IsolationValuePairs IsolationValuePairs;
      // better to loop backwards, so the vector is resized less times
      for ( IsolationValuePairs::const_reverse_iterator it = isolatorTmpStorage_.rbegin(), 
	      ed = isolatorTmpStorage_.rend(); it != ed; ++it) {
	aTau.setIsolation(it->first, it->second);
      }
    }
    
    for (size_t j = 0, nd = deposits.size(); j < nd; ++j) {
      aTau.setIsoDeposit(isoDepositLabels_[j].first, (*deposits[j])[tausRef]);
    }

    if (efficiencyLoader_.enabled()) {
      efficiencyLoader_.setEfficiencies( aTau, tausRef );
    }

    if (resolutionLoader_.enabled()) {
      resolutionLoader_.setResolutions(aTau);
    }

    if ( useUserData_ ) {
      userDataHelper_.add( aTau, iEvent, iSetup );
    }

    patTaus->push_back(aTau);
  }

  // sort taus in pT
  std::sort(patTaus->begin(), patTaus->end(), pTTauComparator_);

  // put genEvt object in Event
  iEvent.put(patTaus);

  // clean up
  if (isolator_.enabled()) isolator_.endEvent();
}

template <typename TauCollectionType, typename TauDiscrType>
float PATTauProducer::getTauIdDiscriminator(const edm::Handle<TauCollectionType>& tauCollection, size_t tauIdx, const edm::Handle<TauDiscrType>& tauIdDiscr)
{
  edm::Ref<TauCollectionType> tauRef(tauCollection, tauIdx);
  return (*tauIdDiscr)[tauRef];
}     

// ParameterSet description for module
void PATTauProducer::fillDescriptions(edm::ConfigurationDescriptions & descriptions)
{
  edm::ParameterSetDescription iDesc;
  iDesc.setComment("PAT tau producer module");

  // input source 
  iDesc.add<edm::InputTag>("tauSource", edm::InputTag())->setComment("input collection");

  // embedding
  iDesc.add<bool>("embedIsolationTracks", false)->setComment("embed external isolation tracks");
  iDesc.add<bool>("embedLeadTrack", false)->setComment("embed external leading track");
  iDesc.add<bool>("embedLeadTracks", false)->setComment("embed external signal tracks");

  // MC matching configurables
  iDesc.add<bool>("addGenMatch", true)->setComment("add MC matching");
  iDesc.add<bool>("embedGenMatch", false)->setComment("embed MC matched MC information");
  std::vector<edm::InputTag> emptySourceVector;
  iDesc.addNode( edm::ParameterDescription<edm::InputTag>("genParticleMatch", edm::InputTag(), true) xor 
                 edm::ParameterDescription<std::vector<edm::InputTag> >("genParticleMatch", emptySourceVector, true)
		 )->setComment("input with MC match information");

  // MC jet matching variables
  iDesc.add<bool>("addGenJetMatch", true)->setComment("add MC jet matching");
  iDesc.add<bool>("embedGenJetMatch", false)->setComment("embed MC jet matched jet information");
  iDesc.add<edm::InputTag>("genJetMatch", edm::InputTag("tauGenJetMatch"));


  pat::helper::KinResolutionsLoader::fillDescription(iDesc);

  // tau ID configurables
  iDesc.add<bool>("addTauID", true)->setComment("add tau ID variables");
  edm::ParameterSetDescription tauIDSourcesPSet;
  tauIDSourcesPSet.setAllowAnything(); 
  iDesc.addNode( edm::ParameterDescription<edm::InputTag>("tauIDSource", edm::InputTag(), true) xor
                 edm::ParameterDescription<edm::ParameterSetDescription>("tauIDSources", tauIDSourcesPSet, true)
               )->setComment("input with electron ID variables");

  // tau decay mode configurables
  iDesc.add<bool>("addDecayMode", false)->setComment("add tau ID variables");  
  iDesc.add<edm::InputTag>("decayModeSrc", edm::InputTag("fixedConePFTauDecayModeProducer"));

  // IsoDeposit configurables
  edm::ParameterSetDescription isoDepositsPSet;
  isoDepositsPSet.addOptional<edm::InputTag>("tracker"); 
  isoDepositsPSet.addOptional<edm::InputTag>("ecal");
  isoDepositsPSet.addOptional<edm::InputTag>("hcal");
  isoDepositsPSet.addOptional<edm::InputTag>("pfAllParticles");
  isoDepositsPSet.addOptional<edm::InputTag>("pfChargedHadron");
  isoDepositsPSet.addOptional<edm::InputTag>("pfNeutralHadron");
  isoDepositsPSet.addOptional<edm::InputTag>("pfGamma");
  isoDepositsPSet.addOptional<std::vector<edm::InputTag> >("user");
  iDesc.addOptional("isoDeposits", isoDepositsPSet);

  // Efficiency configurables
  edm::ParameterSetDescription efficienciesPSet;
  efficienciesPSet.setAllowAnything(); // TODO: the pat helper needs to implement a description.
  iDesc.add("efficiencies", efficienciesPSet);
  iDesc.add<bool>("addEfficiencies", false);

  // Check to see if the user wants to add user data
  edm::ParameterSetDescription userDataPSet;
  PATUserDataHelper<Tau>::fillDescription(userDataPSet);
  iDesc.addOptional("userData", userDataPSet);

  edm::ParameterSetDescription isolationPSet;
  isolationPSet.setAllowAnything(); // TODO: the pat helper needs to implement a description.
  iDesc.add("userIsolation", isolationPSet);

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATTauProducer);


