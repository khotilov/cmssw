#ifndef RecoTBCalo_EcalTBRecProducers_EcalTBWeightUncalibRecHitProducer_HH
#define RecoTBCalo_EcalTBRecProducers_EcalTBWeightUncalibRecHitProducer_HH

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecWeightsAlgo.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "TBDataFormats/EcalTBObjects/interface/EcalTBTDCRecInfo.h"
#include "CondFormats/EcalObjects/interface/EcalWeight.h"

// forward declaration
class EcalTBWeightUncalibRecHitProducer : public edm::EDProducer {

  public:
    typedef std::vector<double> EcalRecoAmplitudes;
    explicit EcalTBWeightUncalibRecHitProducer(const edm::ParameterSet& ps);
    ~EcalTBWeightUncalibRecHitProducer();
    virtual void produce(edm::Event& evt, const edm::EventSetup& es);

  private:

    std::string digiProducer_; // name of module/plugin/producer making digis
    std::string EBdigiCollection_; // secondary name given to collection of digis
    std::string tdcRecInfoCollection_; // secondary name given to collection of digis
    std::string tdcRecInfoProducer_; // secondary name given to collection of digis

    std::string EBhitCollection_; // secondary name to be given to collection of hit

    EcalUncalibRecHitRecWeightsAlgo<EBDataFrame> EBalgo_;

    HepMatrix makeMatrixFromVectors(const std::vector< std::vector<EcalWeight> >& vecvec);
    HepMatrix makeDummySymMatrix(int size);

    int nbTimeBin_;
/*     int nMaxPrintout_; // max # of printouts */
/*     int counter_; // internal verbosity counter */

    //    bool counterExceeded() const { return ( (counter_>nMaxPrintout_) || (counter_<0) ) ; }
};
#endif
