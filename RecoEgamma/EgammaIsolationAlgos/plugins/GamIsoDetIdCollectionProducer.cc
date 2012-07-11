#include "RecoEgamma/EgammaIsolationAlgos/plugins/GamIsoDetIdCollectionProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CommonTools/Utils/interface/StringToEnumValue.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "RecoCaloTools/Selectors/interface/CaloDualConeSelector.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "DataFormats/DetId/interface/DetIdCollection.h"

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"

GamIsoDetIdCollectionProducer::GamIsoDetIdCollectionProducer(const edm::ParameterSet& iConfig) :
            recHitsLabel_(iConfig.getParameter< edm::InputTag > ("recHitsLabel")),
            emObjectLabel_(iConfig.getParameter< edm::InputTag > ("emObjectLabel")),
            energyCut_(iConfig.getParameter<double>("energyCut")),
            etCut_(iConfig.getParameter<double>("etCut")),
            etCandCut_(iConfig.getParameter<double> ("etCandCut")),
            outerRadius_(iConfig.getParameter<double>("outerRadius")),
            innerRadius_(iConfig.getParameter<double>("innerRadius")),
            interestingDetIdCollection_(iConfig.getParameter<std::string>("interestingDetIdCollection"))
   {

     const std::vector<std::string> flagnamesEB = 
       iConfig.getParameter<std::vector<std::string> >("RecHitFlagToBeExcludedEB");
     
     const std::vector<std::string> flagnamesEE =
       iConfig.getParameter<std::vector<std::string> >("RecHitFlagToBeExcludedEE");
     
     flagsexclEB_= 
       StringToEnumValue<EcalRecHit::Flags>(flagnamesEB);
     
     flagsexclEE_=
       StringToEnumValue<EcalRecHit::Flags>(flagnamesEE);
     
     const std::vector<std::string> severitynamesEB = 
       iConfig.getParameter<std::vector<std::string> >("RecHitSeverityToBeExcludedEB");
     
     severitiesexclEB_= 
       StringToEnumValue<EcalSeverityLevel::SeverityLevel>(severitynamesEB);

     const std::vector<std::string> severitynamesEE = 
       iConfig.getParameter<std::vector<std::string> >("RecHitSeverityToBeExcludedEE");
     
     severitiesexclEE_= 
       StringToEnumValue<EcalSeverityLevel::SeverityLevel>(severitynamesEE);
     
    //register your products
    produces< DetIdCollection > (interestingDetIdCollection_) ;
}

GamIsoDetIdCollectionProducer::~GamIsoDetIdCollectionProducer() 
{}

void GamIsoDetIdCollectionProducer::beginJob () 
{}

// ------------ method called to produce the data  ------------
    void
GamIsoDetIdCollectionProducer::produce (edm::Event& iEvent, 
        const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    //Get EM Object
    Handle<reco::PhotonCollection> emObjectH;
    iEvent.getByLabel(emObjectLabel_,emObjectH);

    // take EcalRecHits
    Handle<EcalRecHitCollection> recHitsH;
    iEvent.getByLabel(recHitsLabel_,recHitsH);
    std::auto_ptr<CaloRecHitMetaCollectionV> recHits_(0); 
    recHits_ = std::auto_ptr<CaloRecHitMetaCollectionV>(new EcalRecHitMetaCollection(*recHitsH));

    edm::ESHandle<CaloGeometry> pG;
    iSetup.get<CaloGeometryRecord>().get(pG);    
    const CaloGeometry* caloGeom = pG.product();

    edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
    iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);
    const EcalSeverityLevelAlgo* sevLevel = sevlv.product();

    CaloDualConeSelector *doubleConeSel_ = 0;
    if(recHitsLabel_.instance() == "EcalRecHitsEB")
        doubleConeSel_= new CaloDualConeSelector(innerRadius_,outerRadius_, &*pG, DetId::Ecal, EcalBarrel);
    else if(recHitsLabel_.instance() == "EcalRecHitsEE")
        doubleConeSel_= new CaloDualConeSelector(innerRadius_,outerRadius_, &*pG, DetId::Ecal, EcalEndcap);

    //Create empty output collections
    std::auto_ptr< DetIdCollection > detIdCollection (new DetIdCollection() ) ;

    reco::PhotonCollection::const_iterator emObj;
    if(doubleConeSel_) { //if cone selector was created
        for (emObj = emObjectH->begin(); emObj != emObjectH->end();  emObj++) { //Loop over candidates

            if(emObj->et() < etCandCut_) continue;
            
            GlobalPoint pclu (emObj->caloPosition().x(),emObj->caloPosition().y(),emObj->caloPosition().z());
            std::auto_ptr<CaloRecHitMetaCollectionV> chosen = doubleConeSel_->select(pclu,*recHits_);

            CaloRecHitMetaCollectionV::const_iterator recIt;
            for (recIt = chosen->begin(); recIt!= chosen->end () ; ++recIt) { // Select RecHits 

                if ( (recIt->energy()) < energyCut_) continue;  //dont fill if below E noise value

                double et = recIt->energy() *
                            caloGeom->getPosition(recIt->detid()).perp() /
                            caloGeom->getPosition(recIt->detid()).mag();
                
                if ( et < etCut_) continue;  //dont fill if below ET noise value

		bool isBarrel = false;
		if (fabs(caloGeom->getPosition(recIt->detid()).eta() < 1.479)) 
		  isBarrel = true;

		int severityFlag = sevLevel->severityLevel(recIt->detid(), *recHitsH);
		std::vector<int>::const_iterator sit;
		if (isBarrel) {
		  sit = std::find(severitiesexclEB_.begin(), severitiesexclEB_.end(), severityFlag);
		  if (sit!= severitiesexclEB_.end())
		    continue;
		} else {
		  sit = std::find(severitiesexclEE_.begin(), severitiesexclEE_.end(), severityFlag);
		  if (sit!= severitiesexclEE_.end())
		    continue;
		}

		std::vector<int>::const_iterator vit;
		if (isBarrel) {
		  vit = std::find(flagsexclEB_.begin(), flagsexclEB_.end(), ((EcalRecHit*)(&*recIt))->recoFlag());
		  if (vit != flagsexclEB_.end())
		    continue;
		} else {
		  vit = std::find(flagsexclEE_.begin(), flagsexclEE_.end(), ((EcalRecHit*)(&*recIt))->recoFlag());
		  if (vit != flagsexclEE_.end())
		    continue;
		}

                if(std::find(detIdCollection->begin(),detIdCollection->end(),recIt->detid()) == detIdCollection->end()) 
                    detIdCollection->push_back(recIt->detid());
            } //end rechits

        } //end candidates

        delete doubleConeSel_;
    } //end if cone selector was created
    
    iEvent.put( detIdCollection, interestingDetIdCollection_ );
}
