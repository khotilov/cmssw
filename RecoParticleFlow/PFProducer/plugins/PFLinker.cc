#include "RecoParticleFlow/PFProducer/plugins/PFLinker.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/MuonToMuonMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "RecoParticleFlow/PFProducer/interface/GsfElectronEqual.h"
#include "RecoParticleFlow/PFProducer/interface/PhotonEqual.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


PFLinker::PFLinker(const edm::ParameterSet & iConfig) {
  inputTagPFCandidates_ 
    = iConfig.getParameter<edm::InputTag>("PFCandidate");
  inputTagGsfElectrons_
    = iConfig.getParameter<edm::InputTag>("GsfElectrons");
  inputTagPhotons_
    = iConfig.getParameter<edm::InputTag>("Photons");
  inputTagMuons_
    = iConfig.getParameter<edm::InputTag>("Muons");
  
  nameOutputPF_ 
    = iConfig.getParameter<std::string>("OutputPF");
  
  nameOutputElectronsPF_ 
    = iConfig.getParameter<std::string>("ValueMapElectrons");

  nameOutputPhotonsPF_ 
    = iConfig.getParameter<std::string>("ValueMapPhotons");

  producePFCandidates_  
    = iConfig.getParameter<bool>("ProducePFCandidates");

  nameOutputMergedPF_ 
    = iConfig.getParameter<std::string>("ValueMapMerged"); 
  
  fillMuonRefs_
    = iConfig.getParameter<bool>("FillMuonRefs");


  if(producePFCandidates_) {
    produces<reco::PFCandidateCollection>(nameOutputPF_);
  }
  produces<edm::ValueMap<reco::PFCandidatePtr> > (nameOutputElectronsPF_);
  produces<edm::ValueMap<reco::PFCandidatePtr> > (nameOutputPhotonsPF_);
  produces<edm::ValueMap<reco::PFCandidatePtr> > (nameOutputMergedPF_);
  if(fillMuonRefs_)  produces<edm::ValueMap<reco::PFCandidatePtr> > (inputTagMuons_.label());

}

PFLinker::~PFLinker() {;}

void PFLinker::beginRun(edm::Run& run,const edm::EventSetup & es) {;}

void PFLinker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  std::auto_ptr<reco::PFCandidateCollection>
    pfCandidates_p(new reco::PFCandidateCollection);
  
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  bool status=fetchCollection<reco::PFCandidateCollection>(pfCandidates, 
							   inputTagPFCandidates_, 
							   iEvent );
  
  edm::Handle<reco::GsfElectronCollection> gsfElectrons;
  status=fetchCollection<reco::GsfElectronCollection>(gsfElectrons,
						      inputTagGsfElectrons_,
						      iEvent );
  std::map<reco::GsfElectronRef,reco::PFCandidatePtr> electronCandidateMap;  


  edm::Handle<reco::PhotonCollection> photons;
  status=fetchCollection<reco::PhotonCollection>(photons,
						 inputTagPhotons_,
						 iEvent );
  std::map<reco::PhotonRef,reco::PFCandidatePtr> photonCandidateMap;


  edm::Handle<reco::MuonToMuonMap> muonMap;
  if(fillMuonRefs_)
    status=fetchCollection<reco::MuonToMuonMap>(muonMap,
						inputTagMuons_,
						iEvent);
  std::map<reco::MuonRef,reco::PFCandidatePtr> muonCandidateMap;



  unsigned ncand=(status)?pfCandidates->size():0;

  for( unsigned i=0; i<ncand; ++i) {
    edm::Ptr<reco::PFCandidate> candPtr(pfCandidates,i);
    reco::PFCandidate cand(candPtr);
    
    bool isphoton   = cand.particleId() == reco::PFCandidate::gamma && cand.mva_nothing_gamma()>0.;
    bool iselectron = cand.particleId() == reco::PFCandidate::e;
    bool ismuon     = cand.particleId() == reco::PFCandidate::mu && fillMuonRefs_;

    // if not an electron or a photon or a muon just fill the PFCandidate collection
    if ( !(isphoton || iselectron || ismuon)){pfCandidates_p->push_back(cand); continue;}
    
    
    if (ismuon && fillMuonRefs_) {
      reco::MuonRef muRef = (*muonMap)[cand.muonRef()];
      cand.setMuonRef(muRef);
      muonCandidateMap[muRef] = candPtr;
    }
     



    // if it is an electron. Find the GsfElectron with the same GsfTrack
    if (iselectron) {
      const reco::GsfTrackRef & gsfTrackRef(cand.gsfTrackRef());
      GsfElectronEqual myEqual(gsfTrackRef);
      std::vector<reco::GsfElectron>::const_iterator itcheck=find_if(gsfElectrons->begin(),gsfElectrons->end(),myEqual);
      if(itcheck==gsfElectrons->end()) {
	std::ostringstream err;
	err << " Problem in PFLinker: no GsfElectron " << std::endl;
	edm::LogError("PFLinker") << err.str();
	continue; // Watch out ! Continue
      } 
      reco::GsfElectronRef electronRef(gsfElectrons,itcheck-gsfElectrons->begin());
      cand.setGsfElectronRef(electronRef);
      cand.setSuperClusterRef(electronRef->superCluster());
      electronCandidateMap[electronRef] = candPtr;
    }  
  
    // if it is a photon, find the one with the same PF super-cluster
    if (isphoton) {
      const reco::SuperClusterRef & scRef(cand.superClusterRef());
      PhotonEqual myEqual(scRef);
      std::vector<reco::Photon>::const_iterator itcheck=find_if(photons->begin(),photons->end(),myEqual);
      if(itcheck==photons->end()) {
	std::ostringstream err;
	err << " Problem in PFLinker: no Photon " << std::endl;
	edm::LogError("PFLinker") << err.str();
	continue; // Watch out ! Continue
      } 
      reco::PhotonRef photonRef(photons,itcheck-photons->begin());
      cand.setPhotonRef(photonRef);
      cand.setSuperClusterRef(photonRef->superCluster());
      photonCandidateMap[photonRef] = candPtr;
    }

    pfCandidates_p->push_back(cand);
    
  }
  // save the PFCandidates and get a valid handle

  const edm::OrphanHandle<reco::PFCandidateCollection> pfCandidateRefProd = (producePFCandidates_) ? iEvent.put(pfCandidates_p,nameOutputPF_) :
    edm::OrphanHandle<reco::PFCandidateCollection>();

  
  std::cout<<"1"<<std::endl;

  // now make the valuemaps

  edm::ValueMap<reco::PFCandidatePtr> pfMapGsfElectrons = fillValueMap<reco::GsfElectronCollection>(iEvent, 
								nameOutputElectronsPF_, 
								gsfElectrons, 
								electronCandidateMap,
								pfCandidateRefProd);
  
  std::cout<<"2"<<std::endl;

  edm::ValueMap<reco::PFCandidatePtr> pfMapPhotons = fillValueMap<reco::PhotonCollection>(iEvent, 
						      nameOutputPhotonsPF_, 
						      photons, 
						      photonCandidateMap,
						      pfCandidateRefProd);
  
  std::cout<<"3"<<std::endl;

  edm::ValueMap<reco::PFCandidatePtr> pfMapMuons;

  if(fillMuonRefs_){
    edm::Handle<reco::MuonCollection> muons; 
    iEvent.getByLabel(inputTagMuons_.label(), muons);
    
    pfMapMuons = fillValueMap<reco::MuonCollection>(iEvent, 
						    inputTagMuons_.label(), 
						    muons, 
						    muonCandidateMap,
						    pfCandidateRefProd);
    std::cout<<"4"<<std::endl;
  }

  std::cout<<"5"<<std::endl;
  
  std::auto_ptr<edm::ValueMap<reco::PFCandidatePtr> > 
    pfMapMerged(new edm::ValueMap<reco::PFCandidatePtr>());
  edm::ValueMap<reco::PFCandidatePtr>::Filler pfMapMergedFiller(*pfMapMerged);
  
  std::cout<<"6"<<std::endl;

   *pfMapMerged                   += pfMapGsfElectrons;
   std::cout<<"7"<<std::endl;
   *pfMapMerged                   += pfMapPhotons;
   std::cout<<"8"<<std::endl;
   if(fillMuonRefs_) *pfMapMerged += pfMapMuons;
   std::cout<<"9"<<std::endl;
   iEvent.put(pfMapMerged,nameOutputMergedPF_);
   std::cout<<"10"<<std::endl;

}

template<typename T>
bool PFLinker::fetchCollection(edm::Handle<T>& c, 
			       const edm::InputTag& tag, 
			       const edm::Event& iEvent) const {  

  bool found = iEvent.getByLabel(tag, c);
  
  if(!found )
    {
      std::ostringstream  err;
      err<<" cannot get " <<tag<<std::endl;
      edm::LogError("PFLinker")<<err.str();
    }
  return found;
}



template<typename TYPE>
edm::ValueMap<reco::PFCandidatePtr>  PFLinker::fillValueMap(edm::Event & event,
							    std::string label,
							    edm::Handle<TYPE>& inputObjCollection,
							    const std::map<edm::Ref<TYPE>, reco::PFCandidatePtr> & mapToTheCandidate,
							    const edm::OrphanHandle<reco::PFCandidateCollection> & newPFCandColl) const {
  
  std::auto_ptr<edm::ValueMap<reco::PFCandidatePtr> > pfMap_p(new edm::ValueMap<reco::PFCandidatePtr>());
  edm::ValueMap<reco::PFCandidatePtr>::Filler filler(*pfMap_p);
  
  typedef typename std::map<edm::Ref<TYPE>, reco::PFCandidatePtr>::const_iterator MapTYPE_it; 
  
  unsigned nObj=inputObjCollection->size();
  std::vector<reco::PFCandidatePtr> values(nObj);

  for(unsigned iobj=0; iobj < nObj; ++iobj) {

    edm::Ref<TYPE> objRef(inputObjCollection, iobj);
    MapTYPE_it  itcheck = mapToTheCandidate.find(objRef);

    reco::PFCandidatePtr candPtr;

    if(itcheck != mapToTheCandidate.end())
      candPtr = producePFCandidates_ ? reco::PFCandidatePtr(newPFCandColl,itcheck->second.key()) : itcheck->second;
    
    values[iobj] = candPtr;    
  }

  filler.insert(inputObjCollection,values.begin(),values.end());
  filler.fill();
  edm::ValueMap<reco::PFCandidatePtr> returnValue = *pfMap_p;
  event.put(pfMap_p,label);
  return returnValue;
}
