#include "RecoLocalCalo/EcalRecProducers/plugins/EcalRecHitWorkerSimple.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

EcalRecHitWorkerSimple::EcalRecHitWorkerSimple(const edm::ParameterSet&ps) :
        EcalRecHitWorkerBaseClass(ps)
{
        rechitMaker_ = new EcalRecHitSimpleAlgo();
        v_chstatus_ = ps.getParameter<std::vector<int> >("ChannelStatusToBeExcluded");
}


void EcalRecHitWorkerSimple::set(const edm::EventSetup& es)
{
        es.get<EcalIntercalibConstantsRcd>().get(ical);
        es.get<EcalTimeCalibConstantsRcd>().get(itime);
        es.get<EcalADCToGeVConstantRcd>().get(agc);
        es.get<EcalChannelStatusRcd>().get(chStatus);
        es.get<EcalLaserDbRecord>().get(laser);
}


bool
EcalRecHitWorkerSimple::run( const edm::Event & evt,
                const EcalUncalibratedRecHit& uncalibRH,
                EcalRecHitCollection & result )
{
        DetId detid=uncalibRH.id();

        EcalChannelStatusMap::const_iterator chit = chStatus->find(detid);
        EcalChannelStatusCode chStatusCode = 1;
        if ( chit != chStatus->end() ) {
                chStatusCode = *chit;
        } else {
                edm::LogError("EcalRecHitError") << "No channel status found for xtal " 
                        << detid.rawId() 
                        << "! something wrong with EcalChannelStatus in your DB? ";
        }
        if ( v_chstatus_.size() > 0) {
                std::vector<int>::const_iterator res = std::find( v_chstatus_.begin(), v_chstatus_.end(), chStatusCode.getStatusCode() );
                if ( res != v_chstatus_.end() ) {
                        return false;
                }
        }

        const EcalIntercalibConstantMap& icalMap = ical->getMap();  
        if ( detid.subdetId() == EcalEndcap ) {
                rechitMaker_->setADCToGeVConstant( float(agc->getEEValue()) );
        } else {
                rechitMaker_->setADCToGeVConstant( float(agc->getEBValue()) );
        }

        // first intercalibration constants
        EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(detid);
        EcalIntercalibConstant icalconst = 1;
        if( icalit!=icalMap.end() ) {
                icalconst = (*icalit);
        } else {
                edm::LogError("EcalRecHitError") << "No intercalib const found for xtal "
                        << detid.rawId()
                        << "! something wrong with EcalIntercalibConstants in your DB? ";
        }

        // get laser coefficient
        float lasercalib = laser->getLaserCorrection( detid, evt.time());

        // get time calibration coefficient
        const EcalTimeCalibConstantMap & itimeMap = itime->getMap();  
        EcalTimeCalibConstantMap::const_iterator itime = itimeMap.find(detid);
        EcalTimeCalibConstant itimeconst = 0;
        if( icalit!=icalMap.end() ) {
                itimeconst = (*itime);
        } else {
                edm::LogError("EcalRecHitError") << "No time calib const found for xtal "
                        << detid.rawId()
                        << "! something wrong with EcalTimeCalibConstants in your DB? ";
        }

        // make the rechit and put in the output collection
        result.push_back(EcalRecHit( rechitMaker_->makeRecHit(uncalibRH, icalconst * lasercalib, itimeconst) ));
        return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalRecHitWorkerFactory.h"
DEFINE_EDM_PLUGIN( EcalRecHitWorkerFactory, EcalRecHitWorkerSimple, "EcalRecHitWorkerSimple" );
