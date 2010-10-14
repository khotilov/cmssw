#include "RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerGlobal.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalWeightXtalGroupsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTBWeightsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"

EcalUncalibRecHitWorkerGlobal::EcalUncalibRecHitWorkerGlobal(const edm::ParameterSet&ps) :
        EcalUncalibRecHitWorkerBaseClass(ps)
{
        // ratio method parameters
        EBtimeFitParameters_ = ps.getParameter<std::vector<double> >("EBtimeFitParameters"); 
        EEtimeFitParameters_ = ps.getParameter<std::vector<double> >("EEtimeFitParameters"); 
        EBamplitudeFitParameters_ = ps.getParameter<std::vector<double> >("EBamplitudeFitParameters");
        EEamplitudeFitParameters_ = ps.getParameter<std::vector<double> >("EEamplitudeFitParameters");
        EBtimeFitLimits_.first  = ps.getParameter<double>("EBtimeFitLimits_Lower");
        EBtimeFitLimits_.second = ps.getParameter<double>("EBtimeFitLimits_Upper");
        EEtimeFitLimits_.first  = ps.getParameter<double>("EEtimeFitLimits_Lower");
        EEtimeFitLimits_.second = ps.getParameter<double>("EEtimeFitLimits_Upper");
        EBtimeConstantTerm_=ps.getParameter<double>("EBtimeConstantTerm");
        EBtimeNconst_=ps.getParameter<double>("EBtimeNconst");
        EEtimeConstantTerm_=ps.getParameter<double>("EEtimeConstantTerm");
        EEtimeNconst_=ps.getParameter<double>("EEtimeNconst");
        outOfTimeThreshEB_ = ps.getParameter<double>("outOfTimeThresholdEB");
        outOfTimeThreshEE_ = ps.getParameter<double>("outOfTimeThresholdEE");
        amplitudeThreshEB_ = ps.getParameter<double>("amplitudeThresholdEB");
        amplitudeThreshEE_ = ps.getParameter<double>("amplitudeThresholdEE");
	outOfTimeIfGain12OnlyEB_ = ps.getParameter<bool>("outOfTimeIfGain12OnlyEB");
	outOfTimeIfGain12OnlyEE_ = ps.getParameter<bool>("outOfTimeIfGain12OnlyEE");
        ebSpikeThresh_ = ps.getParameter<double>("ebSpikeThreshold");
        // leading edge parameters
        ebPulseShape_ = ps.getParameter<std::vector<double> >("ebPulseShape");
        eePulseShape_ = ps.getParameter<std::vector<double> >("eePulseShape");
}

void
EcalUncalibRecHitWorkerGlobal::set(const edm::EventSetup& es)
{
        // common setup
        es.get<EcalGainRatiosRcd>().get(gains);
        es.get<EcalPedestalsRcd>().get(peds);

        // for the weights method
        es.get<EcalWeightXtalGroupsRcd>().get(grps);
        es.get<EcalTBWeightsRcd>().get(wgts);

        // for the ratio method

        // for the leading edge method
        es.get<EcalTimeCalibConstantsRcd>().get(itime);
}


// check saturation: 5 samples with gainId = 0
template < class C >
int EcalUncalibRecHitWorkerGlobal::isSaturated(const C & dataFrame)
{
        //bool saturated_ = 0;
        int cnt;
        for (int j = 0; j < C::MAXSAMPLES - 5; ++j) {
                cnt = 0;
                for (int i = j; i < (j + 5) && i < C::MAXSAMPLES; ++i) {
                        if ( dataFrame.sample(i).gainId() == 0 ) ++cnt;
                }
                if ( cnt == 5 ) return j-1 ; // the last unsaturated sample
        }
        return -1; // no saturation found
}



bool
EcalUncalibRecHitWorkerGlobal::run( const edm::Event & evt,
                const EcalDigiCollection::const_iterator & itdg,
                EcalUncalibratedRecHitCollection & result )
{
        DetId detid(itdg->id());
        
        // intelligence for recHit computation
        EcalUncalibratedRecHit uncalibRecHit;
        
        
        const EcalPedestals::Item * aped = 0;
        const EcalMGPAGainRatio * aGain = 0;
        const EcalXtalGroupId * gid = 0;

        if (detid.subdetId()==EcalEndcap) {
                unsigned int hashedIndex = EEDetId(detid).hashedIndex();
                aped  = &peds->endcap(hashedIndex);
                aGain = &gains->endcap(hashedIndex);
                gid   = &grps->endcap(hashedIndex);
        } else {
                unsigned int hashedIndex = EBDetId(detid).hashedIndex();
                aped  = &peds->barrel(hashedIndex);
                aGain = &gains->barrel(hashedIndex);
                gid   = &grps->barrel(hashedIndex);
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

	// compute the right bin of the pulse shape using time calibration constants
	EcalTimeCalibConstantMap::const_iterator it = itime->find( detid );
	EcalTimeCalibConstant itimeconst = 0;
	if( it != itime->end() ) {
		  itimeconst = (*it);
	} else {
		  edm::LogError("EcalRecHitError") << "No time intercalib const found for xtal "
		  << detid.rawId()
		  << "! something wrong with EcalTimeCalibConstants in your DB? ";
	}


        // === amplitude computation ===
        int leadingSample = -1;
        if (detid.subdetId()==EcalEndcap) {
                leadingSample = ((EcalDataFrame)(*itdg)).lastUnsaturatedSample();
        } else {
                leadingSample = ((EcalDataFrame)(*itdg)).lastUnsaturatedSample();
        }

        if ( leadingSample >= 0 ) { // saturation
                if ( leadingSample != 4 ) {
                        // all samples different from the fifth are not reliable for the amplitude estimation
                        // put by default the energy at the saturation threshold and flag as saturated
                        float sratio = 1;
                        if ( detid.subdetId()==EcalBarrel) {
                                sratio = ebPulseShape_[5] / ebPulseShape_[4];
                        } else {
                                sratio = eePulseShape_[5] / eePulseShape_[4];
                        }
                        uncalibRecHit = EcalUncalibratedRecHit( (*itdg).id(), 4095*12*sratio, 0, 0, 0);
                        uncalibRecHit.setRecoFlag( EcalUncalibratedRecHit::kSaturated );
                } else {
                        // float clockToNsConstant = 25.;
                        // reconstruct the rechit
                        if (detid.subdetId()==EcalEndcap) {
                                leadingEdgeMethod_endcap_.setPulseShape( eePulseShape_ );
                                // float mult = (float)eePulseShape_.size() / (float)(*itdg).size();
                                // bin (or some analogous mapping) will be used instead of the leadingSample
                                //int bin  = (int)(( (mult * leadingSample + mult/2) * clockToNsConstant + itimeconst ) / clockToNsConstant);
                                // bin is not uset for the moment
                                leadingEdgeMethod_endcap_.setLeadingEdgeSample( leadingSample );
                                uncalibRecHit = leadingEdgeMethod_endcap_.makeRecHit(*itdg, pedVec, gainRatios, 0, 0);
                                uncalibRecHit.setRecoFlag( EcalUncalibratedRecHit::kLeadingEdgeRecovered );
                                leadingEdgeMethod_endcap_.setLeadingEdgeSample( -1 );
                        } else {
                                leadingEdgeMethod_barrel_.setPulseShape( ebPulseShape_ );
                                // float mult = (float)ebPulseShape_.size() / (float)(*itdg).size();
                                // bin (or some analogous mapping) will be used instead of the leadingSample
                                //int bin  = (int)(( (mult * leadingSample + mult/2) * clockToNsConstant + itimeconst ) / clockToNsConstant);
                                // bin is not uset for the moment
                                leadingEdgeMethod_barrel_.setLeadingEdgeSample( leadingSample );
                                uncalibRecHit = leadingEdgeMethod_barrel_.makeRecHit(*itdg, pedVec, gainRatios, 0, 0);
                                uncalibRecHit.setRecoFlag( EcalUncalibratedRecHit::kLeadingEdgeRecovered );
                                leadingEdgeMethod_barrel_.setLeadingEdgeSample( -1 );
                        }
                }
        } else {
                // weights method
                EcalTBWeights::EcalTDCId tdcid(1);
                EcalTBWeights::EcalTBWeightMap const & wgtsMap = wgts->getMap();
                EcalTBWeights::EcalTBWeightMap::const_iterator wit;
                wit = wgtsMap.find( std::make_pair(*gid,tdcid) );
                if( wit == wgtsMap.end() ) {
                        edm::LogError("EcalUncalibRecHitError") << "No weights found for EcalGroupId: " 
                                << gid->id() << " and  EcalTDCId: " << tdcid
                                << "\n  skipping digi with id: " << detid.rawId();

                        return false;
                }
                const EcalWeightSet& wset = wit->second; // this is the EcalWeightSet

                const EcalWeightSet::EcalWeightMatrix& mat1 = wset.getWeightsBeforeGainSwitch();
                const EcalWeightSet::EcalWeightMatrix& mat2 = wset.getWeightsAfterGainSwitch();

                weights[0] = &mat1;
                weights[1] = &mat2;

                // get uncalibrated recHit from weights
		if (detid.subdetId()==EcalEndcap) {
	    	     uncalibRecHit = weightsMethod_endcap_.makeRecHit(*itdg, pedVec, pedRMSVec, gainRatios, weights, testbeamEEShape);
		} else {
		     uncalibRecHit = weightsMethod_barrel_.makeRecHit(*itdg, pedVec, pedRMSVec, gainRatios, weights, testbeamEBShape);
		}

                // === time computation ===
                // ratio method
                float clockToNsConstant = 25.;
                if (detid.subdetId()==EcalEndcap) {
                                ratioMethod_endcap_.init( *itdg, pedVec, pedRMSVec, gainRatios );
                                ratioMethod_endcap_.computeTime( EEtimeFitParameters_, EEtimeFitLimits_, EEamplitudeFitParameters_ );
                                ratioMethod_endcap_.computeAmplitude( EEamplitudeFitParameters_);
                                EcalUncalibRecHitRatioMethodAlgo<EEDataFrame>::CalculatedRecHit crh = ratioMethod_endcap_.getCalculatedRecHit();
                                uncalibRecHit.setJitter( crh.timeMax - 5 );
                                uncalibRecHit.setJitterError( sqrt(pow(crh.timeError,2) + pow(EEtimeConstantTerm_,2)/pow(clockToNsConstant,2)) );
                                uncalibRecHit.setOutOfTimeEnergy( crh.amplitudeMax );
				
                                if (uncalibRecHit.amplitude() > pedRMSVec[1] * amplitudeThreshEE_){ // why pedRMSVec[1] ? 
                                  // determine if gain has switched away from gainId==1 (x12 gain) to possibly veto flagging kOutOfTime
				  bool allSamplesInGain12(true);
                                  for (int iSample = 0; iSample < EEDataFrame::MAXSAMPLES; iSample++) {
                                    int GainId = ((EcalDataFrame)(*itdg)).sample(iSample).gainId();
                                    if (GainId!=1)  allSamplesInGain12 = false;
                                  }
                                  if ( (outOfTimeIfGain12OnlyEE_ && allSamplesInGain12) || (!outOfTimeIfGain12OnlyEE_) ){
                                    float correctedTime = (crh.timeMax-5) * clockToNsConstant + itimeconst;    
                                    float cterm=EEtimeConstantTerm_;
                                    float sigmaped=pedRMSVec[0];
                                    float nterm=EEtimeNconst_*sigmaped/uncalibRecHit.amplitude();
                                    float sigmat=std::sqrt( nterm*nterm  + cterm*cterm   );
                                    
                                    if ( fabs(correctedTime) > sigmat*outOfTimeThreshEE_ ) {
                                      uncalibRecHit.setRecoFlag( EcalUncalibratedRecHit::kOutOfTime );
                                    }
                                  }
                                }

                } else {
                                ratioMethod_barrel_.init( *itdg, pedVec, pedRMSVec, gainRatios );
                                ratioMethod_barrel_.computeTime( EBtimeFitParameters_, EBtimeFitLimits_, EBamplitudeFitParameters_ );
                                ratioMethod_barrel_.computeAmplitude( EBamplitudeFitParameters_);
                                EcalUncalibRecHitRatioMethodAlgo<EBDataFrame>::CalculatedRecHit crh = ratioMethod_barrel_.getCalculatedRecHit();
                                uncalibRecHit.setJitter( crh.timeMax - 5 );
                                uncalibRecHit.setJitterError( sqrt(pow(crh.timeError,2) + pow(EBtimeConstantTerm_,2)/pow(clockToNsConstant,2)) );
                                uncalibRecHit.setOutOfTimeEnergy( crh.amplitudeMax );

                                if (uncalibRecHit.amplitude() > pedRMSVec[1] * amplitudeThreshEB_){ // why pedRMSVec[1] ? 
                                  // determine if gain has switched away from gainId==1 (x12 gain) to possibly veto flagging kOutOfTime
				  bool allSamplesInGain12(true);
                                  for (int iSample = 0; iSample < EBDataFrame::MAXSAMPLES; iSample++) {
                                    int GainId = ((EcalDataFrame)(*itdg)).sample(iSample).gainId();
                                    if (GainId!=1)  allSamplesInGain12 = false;
                                  }
                                  if ( (outOfTimeIfGain12OnlyEB_ && allSamplesInGain12) || (!outOfTimeIfGain12OnlyEB_) ){
                                    float correctedTime = (crh.timeMax-5) * clockToNsConstant + itimeconst;    
                                    float cterm=EBtimeConstantTerm_;
                                    float sigmaped=pedRMSVec[0];
                                    float nterm=EBtimeNconst_*sigmaped/uncalibRecHit.amplitude();
                                    float sigmat=std::sqrt( nterm*nterm  + cterm*cterm   );
                                    
                                    if ( fabs(correctedTime) > sigmat*outOfTimeThreshEB_ ) {
                                      uncalibRecHit.setRecoFlag( EcalUncalibratedRecHit::kOutOfTime );
                                    }
                                  }
				}
                }
		    
		// === chi2express ===
		if (detid.subdetId()==EcalEndcap) {
		      
		    double amplitude = uncalibRecHit.amplitude();
		    double amplitudeOutOfTime = uncalibRecHit.outOfTimeEnergy();
		    double timePulse= uncalibRecHit.jitter()*25.0; // multiply by 25 to translate ADC clocks to ns
		
		    EcalUncalibRecHitRecChi2Algo<EEDataFrame>chi2expressEE_(
				  					    *itdg, 
				  					    amplitude, 
				  					    itimeconst, 
				  					    amplitudeOutOfTime, 
				  					    timePulse, 
				  					    pedVec, 
				  					    pedRMSVec, 
				  					    gainRatios, 
				  					    testbeamEEShape
		    );
		    double chi2 = chi2expressEE_.chi2();
		    uncalibRecHit.setChi2(chi2);
		    double chi2OutOfTime = chi2expressEE_.chi2OutOfTime();
		    uncalibRecHit.setOutOfTimeChi2(chi2OutOfTime);

		} else {
		    double amplitude = uncalibRecHit.amplitude();
		    double amplitudeOutOfTime = uncalibRecHit.outOfTimeEnergy();
		    double timePulse= uncalibRecHit.jitter()*25.0; // multiply by 25 to translate ADC clocks to ns
		  
		    EcalUncalibRecHitRecChi2Algo<EBDataFrame>chi2expressEB_(
		  							    *itdg, 
		  							    amplitude, 
		  							    itimeconst, 
		  							    amplitudeOutOfTime, 
		  							    timePulse, 
		  							    pedVec, 
		 							    pedRMSVec, 
		  							    gainRatios, 
		  							    testbeamEBShape
		    );
		    double chi2 = chi2expressEB_.chi2();
		    uncalibRecHit.setChi2(chi2);
		    double chi2OutOfTime = chi2expressEB_.chi2OutOfTime();
		    uncalibRecHit.setOutOfTimeChi2(chi2OutOfTime);
		}
        }
        // remove setting of kFake, which can be misleading for the time being
        //if ( detid.subdetId()==EcalBarrel ) {
        //        if ( uncalibRecHit.jitter()*25. > -5 ) {
        //                EBDataFrame dt(*itdg);
        //                if ( dt.spikeEstimator() < ebSpikeThresh_ ) uncalibRecHit.setRecoFlag( EcalUncalibratedRecHit::kFake );
        //        }
        //}


        // put the recHit in the collection
        if (detid.subdetId()==EcalEndcap) {
                result.push_back( uncalibRecHit );
        } else {
                result.push_back( uncalibRecHit );
        }
        return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerFactory.h"
DEFINE_EDM_PLUGIN( EcalUncalibRecHitWorkerFactory, EcalUncalibRecHitWorkerGlobal, "EcalUncalibRecHitWorkerGlobal" );
