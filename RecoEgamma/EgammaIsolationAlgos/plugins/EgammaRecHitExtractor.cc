//*****************************************************************************
// File:      EgammaRecHitExtractor.cc
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer, hacked by Sam Harper (ie the ugly stuff is mine)
// Institute: IIHE-VUB, RAL
//=============================================================================
//*****************************************************************************
//C++ includes
#include <vector>
#include <functional>

//ROOT includes
#include <Math/VectorUtil.h>

//CMSSW includes
#include "RecoEgamma/EgammaIsolationAlgos/plugins/EgammaRecHitExtractor.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;
using namespace egammaisolation;
using namespace reco::isodeposit;

EgammaRecHitExtractor::EgammaRecHitExtractor(const edm::ParameterSet& par) : 
    etMin_(par.getParameter<double>("etMin")),
    energyMin_(par.getParameter<double>("energyMin")),
    extRadius_(par.getParameter<double>("extRadius")),
    intRadius_(par.getParameter<double>("intRadius")),
    intStrip_(par.getParameter<double>("intStrip")),
    barrelEcalHitsTag_(par.getParameter<edm::InputTag>("barrelEcalHits")), 
    endcapEcalHitsTag_(par.getParameter<edm::InputTag>("endcapEcalHits")),
    fakeNegativeDeposit_(par.getParameter<bool>("subtractSuperClusterEnergy")),
    tryBoth_(par.getParameter<bool>("tryBoth")),
    vetoClustered_(par.getParameter<bool>("vetoClustered")),
    sameTag_(false),
    severityLevelCut_(par.getParameter<int>("severityLevelCut")),
    severityRecHitThreshold_(par.getParameter<double>("severityRecHitThreshold")),
    spIdString_(par.getParameter<std::string>("spikeIdString")),
    spIdThreshold_(par.getParameter<double>("spikeIdThreshold"))
{ 

    if     ( !spIdString_.compare("kE1OverE9") )   spId_ = EcalSeverityLevelAlgo::kE1OverE9;
    else if( !spIdString_.compare("kSwissCross") ) spId_ = EcalSeverityLevelAlgo::kSwissCross;
    else                                           spId_ = EcalSeverityLevelAlgo::kSwissCross;

    if ((intRadius_ != 0.0) && (fakeNegativeDeposit_)) {
        throw cms::Exception("Configuration Error") << "EgammaRecHitExtractor: " << 
            "If you use 'subtractSuperClusterEnergy', you *must* set 'intRadius' to ZERO; it does not make sense, otherwise.";
    }
    std::string isoVariable = par.getParameter<std::string>("isolationVariable");
    if (isoVariable == "et") {
        useEt_ = true;
    } else if (isoVariable == "energy") {
        useEt_ = false;
    } else {
        throw cms::Exception("Configuration Error") << "EgammaRecHitExtractor: isolationVariable '" << isoVariable << "' not known. " 
            << " Supported values are 'et', 'energy'. ";
    }
    if (endcapEcalHitsTag_.encode() ==  barrelEcalHitsTag_.encode()) {
        sameTag_ = true;
        if (tryBoth_) {
            edm::LogWarning("EgammaRecHitExtractor") << "If you have configured 'barrelRecHits' == 'endcapRecHits', so I'm switching 'tryBoth' to FALSE.";
            tryBoth_ = false;
        }
    }

}

EgammaRecHitExtractor::~EgammaRecHitExtractor() { }

reco::IsoDeposit EgammaRecHitExtractor::deposit(const edm::Event & iEvent, 
        const edm::EventSetup & iSetup, const reco::Candidate &emObject ) const {
    edm::ESHandle<CaloGeometry> pG;
    iSetup.get<CaloGeometryRecord>().get(pG);

    //Get the channel status from the db
    edm::ESHandle<EcalChannelStatus> chStatus;
    iSetup.get<EcalChannelStatusRcd>().get(chStatus);

    const CaloGeometry* caloGeom = pG.product(); 
    const CaloSubdetectorGeometry* barrelgeom = caloGeom->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
    const CaloSubdetectorGeometry* endcapgeom = caloGeom->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);

    static std::string metname = "EgammaIsolationAlgos|EgammaRecHitExtractor";

    std::auto_ptr<const CaloRecHitMetaCollectionV> barrelRecHits(0), endcapRecHits(0);

    //Get barrel ECAL RecHits 
    edm::Handle<EcalRecHitCollection> barrelEcalRecHitsH;
    iEvent.getByLabel(barrelEcalHitsTag_, barrelEcalRecHitsH);

    //Get endcap ECAL RecHits 
    edm::Handle<EcalRecHitCollection> endcapEcalRecHitsH;
    iEvent.getByLabel(endcapEcalHitsTag_, endcapEcalRecHitsH);


    //define isodeposit starting from candidate
    reco::SuperClusterRef sc = emObject.get<reco::SuperClusterRef>();
    math::XYZPoint caloPosition = sc->position();

    Direction candDir(caloPosition.eta(), caloPosition.phi());
    reco::IsoDeposit deposit( candDir );
    deposit.setVeto( reco::IsoDeposit::Veto(candDir, intRadius_) ); 
    double sinTheta = sin(2*atan(exp(-sc->eta())));
    deposit.addCandEnergy(sc->energy() * (useEt_ ? sinTheta : 1.0)) ;

    // subtract supercluster if desired
    double fakeEnergy = -sc->rawEnergy();
    if (fakeNegativeDeposit_) {
        deposit.addDeposit(candDir, fakeEnergy * (useEt_ ?  sinTheta : 1.0)); // not exactly clean...
    }

    // fill rechits
    bool inBarrel = sameTag_ || ( abs(sc->eta()) < 1.479 ); //check for barrel. If only one collection is used, use barrel
    if (inBarrel || tryBoth_) {
        collect(deposit, sc, barrelgeom, caloGeom, *barrelEcalRecHitsH, chStatus.product(), true);
    } 
    if ((!inBarrel) || tryBoth_) {
        collect(deposit, sc, endcapgeom, caloGeom, *endcapEcalRecHitsH, chStatus.product(), false);
    }

    return deposit;
}

void EgammaRecHitExtractor::collect(reco::IsoDeposit &deposit, 
        const reco::SuperClusterRef& sc, const CaloSubdetectorGeometry* subdet, 
        const CaloGeometry* caloGeom,
        const EcalRecHitCollection &hits, 
        const EcalChannelStatus* chStatus,
        bool barrel) const 
{

    GlobalPoint caloPosition(sc->position().x(), sc->position().y() , sc->position().z());
    CaloSubdetectorGeometry::DetIdSet chosen = subdet->getCells(caloPosition,extRadius_);
    EcalRecHitCollection::const_iterator j=hits.end();
    double caloeta=caloPosition.eta();
    double calophi=caloPosition.phi();
    double r2 = intRadius_*intRadius_;

    std::vector< std::pair<DetId, float> >::const_iterator rhIt;


    for (CaloSubdetectorGeometry::DetIdSet::const_iterator i = chosen.begin(), end = chosen.end() ; i != end;  ++i)  {  
        j=hits.find(*i);
        if(j != hits.end()){
            const  GlobalPoint & position = caloGeom->getPosition(*i);
            double eta = position.eta();
            double phi = position.phi();
            double energy = j->energy();
            double et = energy*position.perp()/position.mag();
            double phiDiff= reco::deltaPhi(phi,calophi);

            //check if we are supposed to veto clustered and then do so
            if(vetoClustered_) {

                //Loop over basic clusters:
                bool isClustered = false;
                for(    reco::CaloCluster_iterator bcIt = sc->clustersBegin();bcIt != sc->clustersEnd(); ++bcIt) {
                    for(rhIt = (*bcIt)->hitsAndFractions().begin();rhIt != (*bcIt)->hitsAndFractions().end(); ++rhIt) {
                        if( rhIt->first == *i ) isClustered = true;
                        if( isClustered ) break;
                    }
                    if( isClustered ) break;
                } //end loop over basic clusters

                if(isClustered) continue;
            }  //end if removeClustered


            if(barrel &&                                 //make sure we have a barrel rechit
               EcalSeverityLevelAlgo::severityLevel(     //call the severity level method
                   EBDetId(j->id()),                     //passing the EBDetId
                   hits,                                 //the rechit collection in order to calculate the swiss crss
                   *chStatus,                            //and the EcalChannelRecHitRcd
                   severityRecHitThreshold_,             //only consider rechits with ET >
                   spId_,                                //the SpikeId method (currently kE1OverE9 or kSwissCross)
                   spIdThreshold_                        //cut value for above
               ) >= severityLevelCut_) continue;         //then if the severity level is too high, we continue to the next rechit

            if(   fabs(et) > etMin_ 
                    && fabs(energy) > energyMin_  //Changed to fabs
                    && fabs(eta-caloeta) > intStrip_ 
                    && (eta-caloeta)*(eta-caloeta) + phiDiff*phiDiff >r2 ) {

                deposit.addDeposit( Direction(eta, phi), (useEt_ ? et : energy) );

            }
        }
    }
} 


