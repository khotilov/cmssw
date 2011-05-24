#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoParticleFlow/PFProducer/plugins/PFPhotonTranslator.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFClusterWidthAlgo.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonCoreFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "RecoEgamma/PhotonIdentification/interface/PhotonIsolationCalculator.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "RecoEgamma/PhotonIdentification/interface/PhotonIsolationCalculator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionFactory.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h" 
#include "CondFormats/EcalObjects/interface/EcalFunctionParameters.h" 
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidatePhotonExtra.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Common/interface/AssociationVector.h"

#include <Math/VectorUtil.h>
#include <vector>
#include "TLorentzVector.h"
#include "TMath.h"

using namespace edm;
using namespace std;
using namespace reco;

using namespace ROOT::Math::VectorUtil;
typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;
typedef math::XYZVector Vector;


PFPhotonTranslator::PFPhotonTranslator(const edm::ParameterSet & iConfig) {

  //std::cout << "PFPhotonTranslator" << std::endl;

  inputTagPFCandidates_ 
    = iConfig.getParameter<edm::InputTag>("PFCandidate");
  
  edm::ParameterSet isoVals  = iConfig.getParameter<edm::ParameterSet> ("isolationValues");
  inputTagIsoVals_.push_back(isoVals.getParameter<edm::InputTag>("pfChargedHadrons"));
  inputTagIsoVals_.push_back(isoVals.getParameter<edm::InputTag>("pfPhotons"));
  inputTagIsoVals_.push_back(isoVals.getParameter<edm::InputTag>("pfNeutralHadrons"));
  

  PFBasicClusterCollection_ = iConfig.getParameter<std::string>("PFBasicClusters");
  PFPreshowerClusterCollection_ = iConfig.getParameter<std::string>("PFPreshowerClusters");
  PFSuperClusterCollection_ = iConfig.getParameter<std::string>("PFSuperClusters");

  PFPhotonCoreCollection_ = iConfig.getParameter<std::string>("PFPhotonCores");
  PFPhotonCollection_ = iConfig.getParameter<std::string>("PFPhotons");

  vertexProducer_   = iConfig.getParameter<std::string>("primaryVertexProducer");

  barrelEcalHits_   = iConfig.getParameter<edm::InputTag>("barrelEcalHits");
  endcapEcalHits_   = iConfig.getParameter<edm::InputTag>("endcapEcalHits");

  hcalTowers_ = iConfig.getParameter<edm::InputTag>("hcalTowers");
  hOverEConeSize_   = iConfig.getParameter<double>("hOverEConeSize");

  if (iConfig.exists("emptyIsOk")) emptyIsOk_ = iConfig.getParameter<bool>("emptyIsOk");
  else emptyIsOk_=false;

  produces<reco::BasicClusterCollection>(PFBasicClusterCollection_); 
  produces<reco::PreshowerClusterCollection>(PFPreshowerClusterCollection_); 
  produces<reco::SuperClusterCollection>(PFSuperClusterCollection_); 
  produces<reco::PhotonCoreCollection>(PFPhotonCoreCollection_);
  produces<reco::PhotonCollection>(PFPhotonCollection_); 
  produces<reco::ConversionCollection>(PFConversionCollection_);
}

PFPhotonTranslator::~PFPhotonTranslator() {}

void PFPhotonTranslator::beginRun(edm::Run& run,const edm::EventSetup & es) {

}

void PFPhotonTranslator::produce(edm::Event& iEvent,  
				    const edm::EventSetup& iSetup) { 

  std::auto_ptr<reco::BasicClusterCollection> 
    basicClusters_p(new reco::BasicClusterCollection);

  std::auto_ptr<reco::PreshowerClusterCollection>
    psClusters_p(new reco::PreshowerClusterCollection);

  std::auto_ptr<reco::ConversionCollection>
    SingleLeg(new reco::ConversionCollection);

  reco::SuperClusterCollection outputSuperClusterCollection;
  reco::PhotonCoreCollection outputPhotonCoreCollection;
  reco::PhotonCollection outputPhotonCollection;

  outputSuperClusterCollection.clear();
  outputPhotonCoreCollection.clear();
  outputPhotonCollection.clear();


  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  bool status=fetchCandidateCollection(pfCandidates, 
				       inputTagPFCandidates_, 
				       iEvent );

  Handle<reco::VertexCollection> vertexHandle;

  
  IsolationValueMaps isolationValues(inputTagIsoVals_.size());
  for (size_t j = 0; j<inputTagIsoVals_.size(); ++j) {
    iEvent.getByLabel(inputTagIsoVals_[j], isolationValues[j]);
  }
  

  // clear the vectors
  photPFCandidateIndex_.clear();
  basicClusters_.clear();
  pfClusters_.clear();
  preshowerClusters_.clear();
  superClusters_.clear();
  basicClusterPtr_.clear();
  preshowerClusterPtr_.clear();
  CandidatePtr_.clear();
  egSCRef_.clear();
  pfConv_.clear();


  // loop on the candidates 
  //CC@@
  // we need first to create AND put the SuperCluster, 
  // basic clusters and presh clusters collection 
  // in order to get a working Handle
  unsigned ncand=(status)?pfCandidates->size():0;

  unsigned iphot=0;
  for( unsigned i=0; i<ncand; ++i ) {

    const reco::PFCandidate& cand = (*pfCandidates)[i];    
    if(cand.particleId()!=reco::PFCandidate::gamma) continue; 
    if(cand. mva_nothing_gamma()>0.001)//Found PFPhoton with PFPhoton Extras saved
      {
	
	//std::cout << "nconv="<<cand.photonExtraRef()->conversionRef().size()<<std::endl;
	pfConv_.push_back(cand.photonExtraRef()->conversionRef());	
	//
	//cand.photonExtraRef()->singleLegConvTrackRef();	
	
      }

    photPFCandidateIndex_.push_back(i);

    basicClusters_.push_back(reco::BasicClusterCollection());
    pfClusters_.push_back(std::vector<const reco::PFCluster *>());
    preshowerClusters_.push_back(reco::PreshowerClusterCollection());
    superClusters_.push_back(reco::SuperClusterCollection());

    reco::PFCandidatePtr ptrToPFPhoton(pfCandidates,i);
    CandidatePtr_.push_back(ptrToPFPhoton);  
    egSCRef_.push_back(cand.superClusterRef());
    //std::cout << "PFPhoton cand " << iphot << std::endl;

    //std::cout << "Cand elements in blocks : " << cand.elementsInBlocks().size() << std::endl;

    for(unsigned iele=0; iele<cand.elementsInBlocks().size(); ++iele) {
      // first get the block 
      reco::PFBlockRef blockRef = cand.elementsInBlocks()[iele].first;
      //
      unsigned elementIndex = cand.elementsInBlocks()[iele].second;
      // check it actually exists 
      if(blockRef.isNull()) continue;
      
      // then get the elements of the block
      const edm::OwnVector< reco::PFBlockElement >&  elements = (*blockRef).elements();
      
      const reco::PFBlockElement & pfbe (elements[elementIndex]); 
      // The first ECAL element should be the cluster associated to the GSF; defined as the seed
      if(pfbe.type()==reco::PFBlockElement::ECAL)
	{	  

	  //std::cout << "BlockElement ECAL" << std::endl;
	  // the Brem photons are saved as daughter PFCandidate; this 
	  // is convenient to access the corrected energy
	  //	  std::cout << " Found candidate "  << correspondingDaughterCandidate(coCandidate,pfbe) << " " << coCandidate << std::endl;
	  createBasicCluster(pfbe,basicClusters_[iphot],pfClusters_[iphot],correspondingDaughterCandidate(cand,pfbe));
	}
      if(pfbe.type()==reco::PFBlockElement::PS1)
	{
	  //std::cout << "BlockElement PS1" << std::endl;
	  createPreshowerCluster(pfbe,preshowerClusters_[iphot],1);
	}
      if(pfbe.type()==reco::PFBlockElement::PS2)
	{
	  //std::cout << "BlockElement PS2" << std::endl;
	  createPreshowerCluster(pfbe,preshowerClusters_[iphot],2);
	}    


    }   // loop on the elements

        // save the basic clusters
    basicClusters_p->insert(basicClusters_p->end(),basicClusters_[iphot].begin(), basicClusters_[iphot].end());
    // save the preshower clusters
    psClusters_p->insert(psClusters_p->end(),preshowerClusters_[iphot].begin(),preshowerClusters_[iphot].end());

    ++iphot;

  } // loop on PFCandidates


   //Save the basic clusters and get an handle as to be able to create valid Refs (thanks to Claude)
  //  std::cout << " Number of basic clusters " << basicClusters_p->size() << std::endl;
  const edm::OrphanHandle<reco::BasicClusterCollection> bcRefProd = 
    iEvent.put(basicClusters_p,PFBasicClusterCollection_);

  //preshower clusters
  const edm::OrphanHandle<reco::PreshowerClusterCollection> psRefProd = 
    iEvent.put(psClusters_p,PFPreshowerClusterCollection_);
  
  // now that the Basic clusters are in the event, the Ref can be created
  createBasicClusterPtrs(bcRefProd);
  // now that the preshower clusters are in the event, the Ref can be created
  createPreshowerClusterPtrs(psRefProd);

  // and now the Super cluster can be created with valid references  
  //if(status) createSuperClusters(*pfCandidates,*superClusters_p);
  if(status) createSuperClusters(*pfCandidates,outputSuperClusterCollection);
  
  //std::cout << "nb superclusters in collection : "<<outputSuperClusterCollection.size()<<std::endl;

  // Let's put the super clusters in the event
  std::auto_ptr<reco::SuperClusterCollection> superClusters_p(new reco::SuperClusterCollection(outputSuperClusterCollection));  
  const edm::OrphanHandle<reco::SuperClusterCollection> scRefProd = iEvent.put(superClusters_p,PFSuperClusterCollection_); 

  const edm::OrphanHandle<reco::ConversionCollection> ConvRefProd = 
    iEvent.put(SingleLeg,PFConversionCollection_);
  /*
  int ipho=0;
  for (reco::SuperClusterCollection::const_iterator gamIter = scRefProd->begin(); gamIter != scRefProd->end(); ++gamIter){
    std::cout << "SC i="<<ipho<<" energy="<<gamIter->energy()<<std::endl;
    ipho++;
  }
  */
  //create photon cores
  //if(status) createPhotonCores(pfCandidates, scRefProd, *photonCores_p);
  if(status) createPhotonCores(scRefProd, outputPhotonCoreCollection);
  
  //std::cout << "nb photoncores in collection : "<<outputPhotonCoreCollection.size()<<std::endl;

  // Put the photon cores in the event
  std::auto_ptr<reco::PhotonCoreCollection> photonCores_p(new reco::PhotonCoreCollection(outputPhotonCoreCollection));  
  //std::cout << "photon core collection put in auto_ptr"<<std::endl;
  const edm::OrphanHandle<reco::PhotonCoreCollection> pcRefProd = iEvent.put(photonCores_p,PFPhotonCoreCollection_); 
  
  //std::cout << "photon core have been put in the event"<<std::endl;
  /*
  int ipho=0;
  for (reco::PhotonCoreCollection::const_iterator gamIter = pcRefProd->begin(); gamIter != pcRefProd->end(); ++gamIter){
    std::cout << "PhotonCore i="<<ipho<<" energy="<<gamIter->pfSuperCluster()->energy()<<std::endl;
    //for (unsigned int i=0; i<)

    std::cout << "PhotonCore i="<<ipho<<" nconv="<<gamIter->conversions().size()<<std::endl;
    ipho++;
  }
  */

  //load vertices
  reco::VertexCollection vertexCollection;
  bool validVertex=true;
  iEvent.getByLabel(vertexProducer_, vertexHandle);
  if (!vertexHandle.isValid()) {
    edm::LogWarning("PhotonProducer") << "Error! Can't get the product primary Vertex Collection "<< "\n";
    validVertex=false;
  }
  if (validVertex) vertexCollection = *(vertexHandle.product());

  //load Ecal rechits
  bool validEcalRecHits=true;
  Handle<EcalRecHitCollection> barrelHitHandle;
  EcalRecHitCollection barrelRecHits;
  iEvent.getByLabel(barrelEcalHits_, barrelHitHandle);
  if (!barrelHitHandle.isValid()) {
    edm::LogError("PhotonProducer") << "Error! Can't get the product "<<barrelEcalHits_.label();
    validEcalRecHits=false; 
  }
  if (  validEcalRecHits)  barrelRecHits = *(barrelHitHandle.product());
  
  Handle<EcalRecHitCollection> endcapHitHandle;
  iEvent.getByLabel(endcapEcalHits_, endcapHitHandle);
  EcalRecHitCollection endcapRecHits;
  if (!endcapHitHandle.isValid()) {
    edm::LogError("PhotonProducer") << "Error! Can't get the product "<<endcapEcalHits_.label();
    validEcalRecHits=false; 
  }
  if( validEcalRecHits) endcapRecHits = *(endcapHitHandle.product());

  //load detector topology & geometry
  iSetup.get<CaloGeometryRecord>().get(theCaloGeom_);

  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopo_);
  const CaloTopology *topology = theCaloTopo_.product();

  // get Hcal towers collection 
  Handle<CaloTowerCollection> hcalTowersHandle;
  iEvent.getByLabel(hcalTowers_, hcalTowersHandle);

  //create photon collection
  if(status) createPhotons(vertexCollection, pcRefProd, topology, &barrelRecHits, &endcapRecHits, hcalTowersHandle, isolationValues, outputPhotonCollection);

  // Put the photons in the event
  std::auto_ptr<reco::PhotonCollection> photons_p(new reco::PhotonCollection(outputPhotonCollection));  
  //std::cout << "photon collection put in auto_ptr"<<std::endl;
  const edm::OrphanHandle<reco::PhotonCollection> photonRefProd = iEvent.put(photons_p,PFPhotonCollection_); 
  //std::cout << "photons have been put in the event"<<std::endl;
  /*
  ipho=0;
  for (reco::PhotonCollection::const_iterator gamIter = photonRefProd->begin(); gamIter != photonRefProd->end(); ++gamIter){
    std::cout << "Photon i="<<ipho<<" energy="<<gamIter->pfSuperCluster()->energy()<<std::endl;
    ipho++;
  }
  */
  
}

bool PFPhotonTranslator::fetchCandidateCollection(edm::Handle<reco::PFCandidateCollection>& c, 
					      const edm::InputTag& tag, 
					      const edm::Event& iEvent) const {  
  bool found = iEvent.getByLabel(tag, c);

  if(!found && !emptyIsOk_)
    {
      std::ostringstream  err;
      err<<" cannot get PFCandidates: "
	 <<tag<<std::endl;
      edm::LogError("PFPhotonTranslator")<<err.str();
    }
  return found;
      
}

// The basic cluster is a copy of the PFCluster -> the energy is not corrected 
// It should be possible to get the corrected energy (including the associated PS energy)
// from the PFCandidate daugthers ; Needs some work 
void PFPhotonTranslator::createBasicCluster(const reco::PFBlockElement & PFBE, 
					      reco::BasicClusterCollection & basicClusters, 
					      std::vector<const reco::PFCluster *> & pfClusters,
					      const reco::PFCandidate & coCandidate) const
{
  reco::PFClusterRef myPFClusterRef = PFBE.clusterRef();
  if(myPFClusterRef.isNull()) return;  

  const reco::PFCluster & myPFCluster (*myPFClusterRef);
  pfClusters.push_back(&myPFCluster);
  //std::cout << " Creating BC " << myPFCluster.energy() << " " << coCandidate.ecalEnergy() <<" "<<  coCandidate.rawEcalEnergy() <<std::endl;
  //std::cout << " # hits " << myPFCluster.hitsAndFractions().size() << std::endl;

//  basicClusters.push_back(reco::CaloCluster(myPFCluster.energy(),
  basicClusters.push_back(reco::CaloCluster(coCandidate.rawEcalEnergy(),
					    myPFCluster.position(),
					    myPFCluster.caloID(),
					    myPFCluster.hitsAndFractions(),
					    myPFCluster.algo(),
					    myPFCluster.seed()));
}

void PFPhotonTranslator::createPreshowerCluster(const reco::PFBlockElement & PFBE, reco::PreshowerClusterCollection& preshowerClusters,unsigned plane) const
{
  reco::PFClusterRef  myPFClusterRef= PFBE.clusterRef();
  preshowerClusters.push_back(reco::PreshowerCluster(myPFClusterRef->energy(),myPFClusterRef->position(),
					       myPFClusterRef->hitsAndFractions(),plane));
}

void PFPhotonTranslator::createBasicClusterPtrs(const edm::OrphanHandle<reco::BasicClusterCollection> & basicClustersHandle )
{
  unsigned size=photPFCandidateIndex_.size();
  unsigned basicClusterCounter=0;
  basicClusterPtr_.resize(size);

  for(unsigned iphot=0;iphot<size;++iphot) // loop on tracks
    {
      unsigned nbc=basicClusters_[iphot].size();
      for(unsigned ibc=0;ibc<nbc;++ibc) // loop on basic clusters
	{
	  //	  std::cout <<  "Track "<< iGSF << " ref " << basicClusterCounter << std::endl;
	  reco::CaloClusterPtr bcPtr(basicClustersHandle,basicClusterCounter);
	  basicClusterPtr_[iphot].push_back(bcPtr);
	  ++basicClusterCounter;
	}
    }
}

void PFPhotonTranslator::createPreshowerClusterPtrs(const edm::OrphanHandle<reco::PreshowerClusterCollection> & preshowerClustersHandle )
{
  unsigned size=photPFCandidateIndex_.size();
  unsigned psClusterCounter=0;
  preshowerClusterPtr_.resize(size);

  for(unsigned iphot=0;iphot<size;++iphot) // loop on tracks
    {
      unsigned nbc=preshowerClusters_[iphot].size();
      for(unsigned ibc=0;ibc<nbc;++ibc) // loop on basic clusters
	{
	  //	  std::cout <<  "Track "<< iGSF << " ref " << basicClusterCounter << std::endl;
	  reco::CaloClusterPtr psPtr(preshowerClustersHandle,psClusterCounter);
	  preshowerClusterPtr_[iphot].push_back(psPtr);
	  ++psClusterCounter;
	}
    }
}

void PFPhotonTranslator::createSuperClusters(const reco::PFCandidateCollection & pfCand,
					       reco::SuperClusterCollection &superClusters) const
{
  unsigned nphot=photPFCandidateIndex_.size();
  for(unsigned iphot=0;iphot<nphot;++iphot)
    {

      //cout << "SC iphot=" << iphot << endl;

      // Computes energy position a la e/gamma 
      double sclusterE=0;
      double posX=0.;
      double posY=0.;
      double posZ=0.;
      
      unsigned nbasics=basicClusters_[iphot].size();
      for(unsigned ibc=0;ibc<nbasics;++ibc)
	{
	  //cout << "BC in SC : iphot="<<iphot<<endl;
	  
	  double e = basicClusters_[iphot][ibc].energy();
	  sclusterE += e;
	  posX += e * basicClusters_[iphot][ibc].position().X();
	  posY += e * basicClusters_[iphot][ibc].position().Y();
	  posZ += e * basicClusters_[iphot][ibc].position().Z();	  
	}
      posX /=sclusterE;
      posY /=sclusterE;
      posZ /=sclusterE;
      
      /*
      if(pfCand[gsfPFCandidateIndex_[iphot]].gsfTrackRef()!=GsfTrackRef_[iphot])
	{
	  edm::LogError("PFElectronTranslator") << " Major problem in PFElectron Translator" << std::endl;
	}
      */      

      // compute the width
      PFClusterWidthAlgo pfwidth(pfClusters_[iphot]);
      
      double correctedEnergy=pfCand[photPFCandidateIndex_[iphot]].ecalEnergy();
      reco::SuperCluster mySuperCluster(correctedEnergy,math::XYZPoint(posX,posY,posZ));
      // protection against empty basic cluster collection ; the value is -2 in this case
      if(nbasics)
	{
//	  std::cout << "SuperCluster creation; energy " << pfCand[gsfPFCandidateIndex_[iphot]].ecalEnergy();
//	  std::cout << " " <<   pfCand[gsfPFCandidateIndex_[iphot]].rawEcalEnergy() << std::endl;
//	  std::cout << "Seed energy from basic " << basicClusters_[iphot][0].energy() << std::endl;
	  mySuperCluster.setSeed(basicClusterPtr_[iphot][0]);
	}
      else
	{
	  //	  std::cout << "SuperCluster creation ; seed energy " << 0 << std::endl;
	  //std::cout << "SuperCluster creation ; energy " << pfCand[photPFCandidateIndex_[iphot]].ecalEnergy();
	  //std::cout << " " <<   pfCand[photPFCandidateIndex_[iphot]].rawEcalEnergy() << std::endl;
//	  std::cout << " No seed found " << 0 << std::endl;	  
//	  std::cout << " MVA " << pfCand[gsfPFCandidateIndex_[iphot]].mva_e_pi() << std::endl;
	  mySuperCluster.setSeed(reco::CaloClusterPtr());
	}
      // the seed should be the first basic cluster

      for(unsigned ibc=0;ibc<nbasics;++ibc)
	{
	  mySuperCluster.addCluster(basicClusterPtr_[iphot][ibc]);
	  //	  std::cout <<"Adding Ref to SC " << basicClusterPtr_[iphot][ibc].index() << std::endl;
	  const std::vector< std::pair<DetId, float> > & v1 = basicClusters_[iphot][ibc].hitsAndFractions();
	  //	  std::cout << " Number of cells " << v1.size() << std::endl;
	  for( std::vector< std::pair<DetId, float> >::const_iterator diIt = v1.begin();
	       diIt != v1.end();
	       ++diIt ) {
	    //	    std::cout << " Adding DetId " << (diIt->first).rawId() << " " << diIt->second << std::endl;
	    mySuperCluster.addHitAndFraction(diIt->first,diIt->second);
	  } // loop over rechits      
	}      

      unsigned nps=preshowerClusterPtr_[iphot].size();
      for(unsigned ips=0;ips<nps;++ips)
	{
	  mySuperCluster.addPreshowerCluster(preshowerClusterPtr_[iphot][ips]);
	}
      

      // Set the preshower energy
      mySuperCluster.setPreshowerEnergy(pfCand[photPFCandidateIndex_[iphot]].pS1Energy()+
					pfCand[photPFCandidateIndex_[iphot]].pS2Energy());

      // Set the cluster width
      mySuperCluster.setEtaWidth(pfwidth.pflowEtaWidth());
      mySuperCluster.setPhiWidth(pfwidth.pflowPhiWidth());
      // Force the computation of rawEnergy_ of the reco::SuperCluster
      mySuperCluster.rawEnergy();

      //cout << "SC energy="<< mySuperCluster.energy()<<endl;

      superClusters.push_back(mySuperCluster);
      //std::cout << "nb super clusters in collection : "<<superClusters.size()<<std::endl;
    }
}

void PFPhotonTranslator::createPhotonCores(const edm::OrphanHandle<reco::SuperClusterCollection> & superClustersHandle, reco::PhotonCoreCollection &photonCores)
{
  
  //std::cout << "createPhotonCores" << std::endl;

  unsigned nphot=photPFCandidateIndex_.size();

  for(unsigned iphot=0;iphot<nphot;++iphot)
    {
      //std::cout << "iphot="<<iphot<<std::endl;

      reco::PhotonCore myPhotonCore;

      reco::SuperClusterRef SCref(reco::SuperClusterRef(superClustersHandle, iphot));
      
      myPhotonCore.setPFlowPhoton(true);
      myPhotonCore.setStandardPhoton(false);
      myPhotonCore.setPflowSuperCluster(SCref);
      myPhotonCore.setSuperCluster(egSCRef_[iphot]);
      reco::ConversionRefVector ConvToAdd=pfConv_[iphot];
      for(unsigned int iConv=0; iConv<ConvToAdd.size(); iConv++)myPhotonCore.addConversion(ConvToAdd[iConv]);

      photonCores.push_back(myPhotonCore);
      
    }

  //std::cout << "end of createPhotonCores"<<std::endl;
}

void PFPhotonTranslator::createPhotons(reco::VertexCollection &vertexCollection, const edm::OrphanHandle<reco::PhotonCoreCollection> & photonCoresHandle, const CaloTopology* topology, const EcalRecHitCollection* barrelRecHits, const EcalRecHitCollection* endcapRecHits, const edm::Handle<CaloTowerCollection> & hcalTowersHandle, const IsolationValueMaps& isolationValues, reco::PhotonCollection &photons)
{

  //cout << "createPhotons" << endl;
  
  unsigned nphot=photPFCandidateIndex_.size();

  for(unsigned iphot=0;iphot<nphot;++iphot)
    {
      //std::cout << "iphot="<<iphot<<std::endl;

      reco::PhotonCoreRef PCref(reco::PhotonCoreRef(photonCoresHandle, iphot));

      math::XYZPoint vtx(0.,0.,0.);
      if (vertexCollection.size()>0) vtx = vertexCollection.begin()->position();
      //std::cout << "vtx made" << std::endl;

      math::XYZVector direction =  PCref->pfSuperCluster()->position() - vtx;

      //It could be that pfSC energy gives not the best resolution : use smaller agregates for some cases ?
      math::XYZVector P3 = direction.unit() * PCref->pfSuperCluster()->energy();
      LorentzVector P4(P3.x(), P3.y(), P3.z(), PCref->pfSuperCluster()->energy());

      reco::Photon myPhoton(P4, PCref->pfSuperCluster()->position(), PCref, vtx);
      //cout << "photon created"<<endl;
      

      reco::Photon::PflowIsolationVariables myPFIso;
      myPFIso.chargedHadronIso=(*isolationValues[0])[CandidatePtr_[iphot]];
      myPFIso.photonIso=(*isolationValues[1])[CandidatePtr_[iphot]];
      myPFIso.neutralHadronIso=(*isolationValues[2])[CandidatePtr_[iphot]];   
      myPhoton.setPflowIsolationVariables(myPFIso);
      
      //cout << "chargedHadronIso="<<myPhoton.chargedHadronIso()<<" photonIso="<<myPhoton.photonIso()<<" neutralHadronIso="<<myPhoton.neutralHadronIso()<<endl;
      

      if (basicClusters_[iphot].size()>0){
      // Cluster shape variables
      //Algorithms from EcalClusterTools could be adapted to PF photons ? (not based on 5x5 BC)
      //It happens that energy computed in eg e5x5 is greater than pfSC energy (EcalClusterTools gathering energies from adjacent crystals even if not belonging to the SC)
      const EcalRecHitCollection* hits = 0 ;
      int subdet = PCref->pfSuperCluster()->seed()->hitsAndFractions()[0].first.subdetId();
      if (subdet==EcalBarrel) hits = barrelRecHits;
      else if  (subdet==EcalEndcap) hits = endcapRecHits;
      const CaloGeometry* geometry = theCaloGeom_.product();

      float maxXtal =   EcalClusterTools::eMax( *(PCref->pfSuperCluster()->seed()), &(*hits) );
      //cout << "maxXtal="<<maxXtal<<endl;
      float e1x5    =   EcalClusterTools::e1x5(  *(PCref->pfSuperCluster()->seed()), &(*hits), &(*topology)); 
      //cout << "e1x5="<<e1x5<<endl;
      float e2x5    =   EcalClusterTools::e2x5Max(  *(PCref->pfSuperCluster()->seed()), &(*hits), &(*topology)); 
      //cout << "e2x5="<<e2x5<<endl;
      float e3x3    =   EcalClusterTools::e3x3(  *(PCref->pfSuperCluster()->seed()), &(*hits), &(*topology)); 
      //cout << "e3x3="<<e3x3<<endl;
      float e5x5    =   EcalClusterTools::e5x5( *(PCref->pfSuperCluster()->seed()), &(*hits), &(*topology)); 
      //cout << "e5x5="<<e5x5<<endl;
      std::vector<float> cov =  EcalClusterTools::covariances( *(PCref->pfSuperCluster()->seed()), &(*hits), &(*topology), geometry); 
      float sigmaEtaEta = sqrt(cov[0]);
      //cout << "sigmaEtaEta="<<sigmaEtaEta<<endl;
      std::vector<float> locCov =  EcalClusterTools::localCovariances( *(PCref->pfSuperCluster()->seed()), &(*hits), &(*topology)); 
      float sigmaIetaIeta = sqrt(locCov[0]);
      //cout << "sigmaIetaIeta="<<sigmaIetaIeta<<endl;
      //float r9 =e3x3/(PCref->pfSuperCluster()->rawEnergy());


      // calculate HoE
      const CaloTowerCollection* hcalTowersColl = hcalTowersHandle.product();
      EgammaTowerIsolation towerIso1(hOverEConeSize_,0.,0.,1,hcalTowersColl) ;  
      EgammaTowerIsolation towerIso2(hOverEConeSize_,0.,0.,2,hcalTowersColl) ;  
      double HoE1=towerIso1.getTowerESum(&(*PCref->pfSuperCluster()))/PCref->pfSuperCluster()->energy();
      double HoE2=towerIso2.getTowerESum(&(*PCref->pfSuperCluster()))/PCref->pfSuperCluster()->energy(); 
      //cout << "HoE1="<<HoE1<<endl;
      //cout << "HoE2="<<HoE2<<endl;  

      reco::Photon::ShowerShape  showerShape;
      showerShape.e1x5= e1x5;
      showerShape.e2x5= e2x5;
      showerShape.e3x3= e3x3;
      showerShape.e5x5= e5x5;
      showerShape.maxEnergyXtal =  maxXtal;
      showerShape.sigmaEtaEta =    sigmaEtaEta;
      showerShape.sigmaIetaIeta =  sigmaIetaIeta;
      showerShape.hcalDepth1OverEcal = HoE1;
      showerShape.hcalDepth2OverEcal = HoE2;
      myPhoton.setShowerShapeVariables ( showerShape ); 
      //cout << "shower shape variables filled"<<endl;
      }

      photons.push_back(myPhoton);

    }

  //std::cout << "end of createPhotons"<<std::endl;
}


const reco::PFCandidate & PFPhotonTranslator::correspondingDaughterCandidate(const reco::PFCandidate & cand, const reco::PFBlockElement & pfbe) const
{
  unsigned refindex=pfbe.index();
  //  std::cout << " N daughters " << cand.numberOfDaughters() << std::endl;
  reco::PFCandidate::const_iterator myDaughterCandidate=cand.begin();
  reco::PFCandidate::const_iterator itend=cand.end();

  for(;myDaughterCandidate!=itend;++myDaughterCandidate)
    {
      const reco::PFCandidate * myPFCandidate = (const reco::PFCandidate*)&*myDaughterCandidate;
      if(myPFCandidate->elementsInBlocks().size()!=1)
	{
	  //	  std::cout << " Daughter with " << myPFCandidate.elementsInBlocks().size()<< " element in block " << std::endl;
	  return cand;
	}
      if(myPFCandidate->elementsInBlocks()[0].second==refindex) 
	{
	  //	  std::cout << " Found it " << cand << std::endl;
	  return *myPFCandidate;
	}      
    }
  return cand;
}

