#include "RecoEgamma/EgammaIsolationAlgos/plugins/EleIsoDetIdCollectionProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "RecoCaloTools/Selectors/interface/CaloDualConeSelector.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
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

EleIsoDetIdCollectionProducer::EleIsoDetIdCollectionProducer(const edm::ParameterSet& iConfig) :
            recHitsLabel_(iConfig.getParameter< edm::InputTag > ("recHitsLabel")),
            emObjectLabel_(iConfig.getParameter< edm::InputTag > ("emObjectLabel")),
            energyCut_(iConfig.getParameter<double>("energyCut")),
            etCut_(iConfig.getParameter<double>("etCut")),
            etCandCut_(iConfig.getParameter<double> ("etCandCut")),
            outerRadius_(iConfig.getParameter<double>("outerRadius")),
            innerRadius_(iConfig.getParameter<double>("innerRadius")),
            interestingDetIdCollection_(iConfig.getParameter<std::string>("interestingDetIdCollection")),
            severityLevelCut_(iConfig.getParameter<int>("severityLevelCut")),
            //severityRecHitThreshold_(iConfig.getParameter<double>("severityRecHitThreshold")),
            //spIdString_(iConfig.getParameter<std::string>("spikeIdString")),
            //spIdThreshold_(iConfig.getParameter<double>("spikeIdThreshold")),
            v_chstatus_(iConfig.getParameter<std::vector<int> >("recHitFlagsToBeExcluded")) {


  //if     ( !spIdString_.compare("kE1OverE9") )   spId_ = EcalSeverityLevelAlgo::kE1OverE9;
  // else if( !spIdString_.compare("kSwissCross") ) spId_ = EcalSeverityLevelAlgo::kSwissCross;
  //  else                                           spId_ = EcalSeverityLevelAlgo::kSwissCross;
    
    //register your products
    produces< DetIdCollection > (interestingDetIdCollection_) ;
}

EleIsoDetIdCollectionProducer::~EleIsoDetIdCollectionProducer() {}

void EleIsoDetIdCollectionProducer::beginJob () 
{}

// ------------ method called to produce the data  ------------
    void
EleIsoDetIdCollectionProducer::produce (edm::Event& iEvent, 
        const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    //Get EM Object
    Handle<reco::GsfElectronCollection> emObjectH;
    iEvent.getByLabel(emObjectLabel_,emObjectH);

    // take EcalRecHits
    Handle<EcalRecHitCollection> recHitsH;
    iEvent.getByLabel(recHitsLabel_,recHitsH);
    std::auto_ptr<CaloRecHitMetaCollectionV> recHits_(0); 
    recHits_ = std::auto_ptr<CaloRecHitMetaCollectionV>(new EcalRecHitMetaCollection(*recHitsH));

    edm::ESHandle<CaloGeometry> pG;
    iSetup.get<CaloGeometryRecord>().get(pG);    
    const CaloGeometry* caloGeom = pG.product();

    //Get the channel status from the db
    edm::ESHandle<EcalChannelStatus> chStatus;
    iSetup.get<EcalChannelStatusRcd>().get(chStatus);

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

    reco::GsfElectronCollection::const_iterator emObj;
    if(doubleConeSel_) { //if cone selector was created
        for (emObj = emObjectH->begin(); emObj != emObjectH->end();  emObj++) { //Loop over candidates

            if(emObj->et() < etCandCut_) continue; //don't calculate if object hasn't enough energy
            
            GlobalPoint pclu (emObj->caloPosition().x(),emObj->caloPosition().y(),emObj->caloPosition().z());
            std::auto_ptr<CaloRecHitMetaCollectionV> chosen = doubleConeSel_->select(pclu,*recHits_);

            CaloRecHitMetaCollectionV::const_iterator recIt;
            for (recIt = chosen->begin(); recIt!= chosen->end () ; ++recIt) { // Select RecHits 

                if ( fabs(recIt->energy()) < energyCut_) continue;  //dont fill if below E noise value


                double et = recIt->energy() * 
                            caloGeom->getPosition(recIt->detid()).perp() / 
                            caloGeom->getPosition(recIt->detid()).mag();

                if ( fabs(et) < etCut_) continue;  //dont fill if below ET noise value

                //make sure we have a barrel rechit                                     
                //call the severity level method                                        
                //passing the EBDetId                                                   
                //the rechit collection in order to calculate the swiss crss            
                //and the EcalChannelRecHitRcd                                          
                //only consider rechits with ET >                                       
                //the SpikeId method (currently kE1OverE9 or kSwissCross)               
                //cut value for above                                                   
                //then if the severity level is too high, we continue to the next rechit
                
                if(recHitsLabel_.instance() == "EcalRecHitsEB" && 
                   sevLevel->severityLevel(EBDetId(recIt->detid()), *recHitsH)  >= severityLevelCut_) continue;
                //                       *chStatus,                                 
                  //      severityRecHitThreshold_,                 
                //    spId_,                                    
                //      spIdThreshold_                            
                  // ) >= severityLevelCut_) continue;              

                //Check based on flags to protect from recovered channels from non-read towers
                //Assumption is that v_chstatus_ is empty unless doFlagChecks() has been called
                std::vector<int>::const_iterator vit = std::find( v_chstatus_.begin(), v_chstatus_.end(),  ((EcalRecHit*)(&*recIt))->recoFlag() );
                if ( vit != v_chstatus_.end() ) continue; // the recHit has to be excluded from the iso sum


                if(std::find(detIdCollection->begin(),detIdCollection->end(),recIt->detid()) == detIdCollection->end()) 
		            detIdCollection->push_back(recIt->detid()); 
            } //end rechits

        } //end candidates

        delete doubleConeSel_;
    } //end if cone selector was created
    
    iEvent.put( detIdCollection, interestingDetIdCollection_ );
}
