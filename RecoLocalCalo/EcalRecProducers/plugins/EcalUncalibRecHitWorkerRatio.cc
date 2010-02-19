#include "RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerRatio.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

EcalUncalibRecHitWorkerRatio::EcalUncalibRecHitWorkerRatio(const edm::ParameterSet&ps) :
        EcalUncalibRecHitWorkerBaseClass(ps)
{
  EBtimeFitParameters_ = ps.getParameter<std::vector<double> >("EBtimeFitParameters"); 
  EEtimeFitParameters_ = ps.getParameter<std::vector<double> >("EEtimeFitParameters"); 
 
  EBamplitudeFitParameters_ = ps.getParameter<std::vector<double> >("EBamplitudeFitParameters"); 
  EEamplitudeFitParameters_ = ps.getParameter<std::vector<double> >("EEamplitudeFitParameters"); 
 
  EBtimeFitLimits_.first  = ps.getParameter<double>("EBtimeFitLimits_Lower"); 
  EBtimeFitLimits_.second = ps.getParameter<double>("EBtimeFitLimits_Upper"); 
 
  EEtimeFitLimits_.first  = ps.getParameter<double>("EEtimeFitLimits_Lower"); 
  EEtimeFitLimits_.second = ps.getParameter<double>("EEtimeFitLimits_Upper"); 
}

void
EcalUncalibRecHitWorkerRatio::set(const edm::EventSetup& es)
{
        es.get<EcalGainRatiosRcd>().get(gains);
        es.get<EcalPedestalsRcd>().get(peds);

}


bool
EcalUncalibRecHitWorkerRatio::run( const edm::Event & evt,
                const EcalDigiCollection::const_iterator & itdg,
                EcalUncalibratedRecHitCollection & result )
{
        DetId detid(itdg->id());

        const EcalPedestals::Item * aped = 0;
        const EcalMGPAGainRatio * aGain = 0;

        if (detid.subdetId()==EcalEndcap) {
                unsigned int hashedIndex = EEDetId(detid).hashedIndex();
                aped  = &peds->endcap(hashedIndex);
                aGain = &gains->endcap(hashedIndex);
        } else {
                unsigned int hashedIndex = EBDetId(detid).hashedIndex();
                aped  = &peds->barrel(hashedIndex);
                aGain = &gains->barrel(hashedIndex);
        }

        pedVec[0] = aped->mean_x12;
        pedVec[1] = aped->mean_x6;
        pedVec[2] = aped->mean_x1;
        pedRMSVec[0] = aped->rms_x12;
        pedRMSVec[1] = aped->rms_x6;
        pedRMSVec[2] = aped->rms_x1;
        gainRatios[0] = 1.;
        gainRatios[1] = aGain->gain12Over6();
        gainRatios[2] = aGain->gain6Over1()*aGain->gain12Over6();

	if (detid.subdetId()==EcalEndcap) {
	  result.push_back(uncalibMaker_endcap_.makeRecHit(*itdg, pedVec, pedRMSVec, gainRatios, EBtimeFitParameters_, EBamplitudeFitParameters_, EBtimeFitLimits_));
        } else {
	  result.push_back(uncalibMaker_barrel_.makeRecHit(*itdg, pedVec, pedRMSVec, gainRatios, EBtimeFitParameters_, EBamplitudeFitParameters_, EBtimeFitLimits_));
        }

        return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerFactory.h"
DEFINE_EDM_PLUGIN( EcalUncalibRecHitWorkerFactory, EcalUncalibRecHitWorkerRatio, "EcalUncalibRecHitWorkerRatio" );
