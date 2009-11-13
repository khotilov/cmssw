#ifndef CalibTracker_SiStripESProducers_SiStripGainESProducerTemplate_h
#define CalibTracker_SiStripESProducers_SiStripGainESProducerTemplate_h

// system include files
#include <memory>
#include "boost/shared_ptr.hpp"

// user include files
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CalibFormats/SiStripObjects/interface/SiStripGain.h"
#include "CondFormats/SiStripObjects/interface/SiStripApvGain.h"
#include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"
#include "CalibTracker/Records/interface/SiStripDependentRecords.h"
//
// class declaration
//



template<typename TDependentRecord, typename TInputRecord>
class SiStripGainESProducerTemplate : public edm::ESProducer {
 public:
  SiStripGainESProducerTemplate(const edm::ParameterSet&);
  ~SiStripGainESProducerTemplate(){};
  
  std::auto_ptr<SiStripGain> produce(const TDependentRecord&);

 private:

  SiStripGain* SiStripGainNormalizationFunction(const TDependentRecord& iRecord);
  double getNFactor(const int apvGainIndex);


  std::vector<edm::ParameterSet> apvGainLabels_;
  std::vector<std::string> apvgain_;
  double norm_;
  bool automaticMode_;
  bool  printdebug_;
  SiStripGain * gain_;
  std::vector<edm::ESHandle<SiStripApvGain> > pDD;
  // edm::ESHandle<SiStripApvGain> pDD;

};

template<typename TDependentRecord, typename TInputRecord>
SiStripGainESProducerTemplate<TDependentRecord,TInputRecord>::SiStripGainESProducerTemplate(const edm::ParameterSet& iConfig)
{
  setWhatProduced(this);

  automaticMode_ = iConfig.getParameter<bool>("AutomaticNormalization");
  norm_=iConfig.getParameter<double>("NormalizationFactor");
  printdebug_ = iConfig.getUntrackedParameter<bool>("printDebug", false);
  apvGainLabels_ = iConfig.getParameter<std::vector<edm::ParameterSet> >("APVGain");

  // Fill the vector of apv labels
  std::vector<edm::ParameterSet>::const_iterator gainPSetIt = apvGainLabels_.begin();
  for( ; gainPSetIt != apvGainLabels_.end(); ++gainPSetIt ) {
    apvgain_.push_back(gainPSetIt->getParameter<std::string>("Label"));
  }

  if(!automaticMode_ && norm_<=0){
    edm::LogError("SiStripGainESProducer") << "[SiStripGainESProducer] - ERROR: negative or zero Normalization factor provided. Assuming 1 for such factor" << std::endl;
    norm_=1.;
  }
}

template<typename TDependentRecord, typename TInputRecord>
std::auto_ptr<SiStripGain> SiStripGainESProducerTemplate<TDependentRecord,TInputRecord>::produce(const TDependentRecord& iRecord)
{  
  std::auto_ptr<SiStripGain> ptr(SiStripGainNormalizationFunction(iRecord));
  return ptr;
}

template<typename TDependentRecord, typename TInputRecord>
SiStripGain* SiStripGainESProducerTemplate<TDependentRecord,TInputRecord>::SiStripGainNormalizationFunction(const TDependentRecord& iRecord)
{
  if(typeid(TDependentRecord)==typeid(SiStripGainRcd) && typeid(TInputRecord)==typeid(SiStripApvGainRcd)){

    const SiStripGainRcd& a = dynamic_cast<const SiStripGainRcd&>(iRecord);

    pDD.push_back(edm::ESHandle<SiStripApvGain>());
    a.getRecord<SiStripApvGainRcd>().get(apvgain_[0], pDD[0]);
    SiStripGain * gain = new SiStripGain( *(pDD[0].product()), getNFactor(0));

    if( apvgain_.size() > 1 ) {
      // std::vector<edm::ESHandle<SiStripApvGain> >::const_iterator apvGainIt = apvgain_.begin()+1;
      for( unsigned int i=1; i<apvgain_.size(); ++i ) {
      // for( ; apvGainIt != apvgain_.end(); ++apvGainIt ) {
        pDD.push_back(edm::ESHandle<SiStripApvGain>());
        // a.getRecord<SiStripApvGainRcd>().get(*apvGainIt, *(pDD.back()));
        a.getRecord<SiStripApvGainRcd>().get(apvgain_[i], pDD[i]);
        gain->multiply(*(pDD[i].product()), getNFactor(i));
      }
    }
    return gain;

  }else if(typeid(TDependentRecord)==typeid(SiStripGainSimRcd) && typeid(TInputRecord)==typeid(SiStripApvGainSimRcd)){

    const SiStripGainSimRcd& a = dynamic_cast<const SiStripGainSimRcd&>(iRecord);

    pDD.push_back(edm::ESHandle<SiStripApvGain>());
    a.getRecord<SiStripApvGainSimRcd>().get(apvgain_[0], pDD[0]);
    SiStripGain * gain = new SiStripGain( *(pDD[0].product()), getNFactor(0));

    if( apvgain_.size() > 1 ) {
      // std::vector<edm::ESHandle<SiStripApvGain> >::const_iterator apvGainIt = apvgain_.begin()+1;
      for( unsigned int i=1; i<apvgain_.size(); ++i ) {
      // for( ; apvGainIt != apvgain_.end(); ++apvGainIt ) {
        pDD.push_back(edm::ESHandle<SiStripApvGain>());
        // a.getRecord<SiStripApvGainRcd>().get(*apvGainIt, *(pDD.back()));
        a.getRecord<SiStripApvGainSimRcd>().get(apvgain_[i], pDD[i]);
        gain->multiply(*(pDD[i].product()), getNFactor(i));
      }
    }
    return gain;
  }
    
  edm::LogError("SiStripGainESProducer") << "[SiStripGainNormalizationFunction] - ERROR: asking for a pair of records different from <SiStripGainRcd,SiStripApvGainRcd> and <SiStripGainSimRcd,SiStripApvGainSimRcd>" << std::endl;
  return new SiStripGain();
}

template<typename TDependentRecord, typename TInputRecord>
double  SiStripGainESProducerTemplate<TDependentRecord,TInputRecord>::getNFactor(const int apvGainIndex){
  double NFactor;

  if(automaticMode_ || printdebug_ ){

    std::vector<uint32_t> DetIds;
    pDD[apvGainIndex]->getDetIds(DetIds);

    double SumOfGains=0;
    int NGains=0;
    
    for(std::vector<uint32_t>::const_iterator detit=DetIds.begin(); detit!=DetIds.end(); detit++){
      
      SiStripApvGain::Range detRange = pDD[apvGainIndex]->getRange(*detit);
      
      int iComp=0;
       
      for(std::vector<float>::const_iterator apvit=detRange.first; apvit!=detRange.second; apvit++){
	 
	SumOfGains+=(*apvit);
	NGains++;
	if (printdebug_)
	  edm::LogInfo("SiStripGainESProducer::produce()")<< "detid/component: " << *detit <<"/"<<iComp<< "   gain factor " <<*apvit ;
	iComp++;
      }      
    }
    
    if(automaticMode_){
      if(SumOfGains>0 && NGains>0){
	NFactor=SumOfGains/NGains;
      }
      else{
	edm::LogError("SiStripGainESProducer::produce() - ERROR: empty set of gain values received. Cannot compute normalization factor. Assuming 1 for such factor") << std::endl;
	NFactor=1.;
      }
    }
  }
  
  if(!automaticMode_){
    NFactor=norm_;
  }

  if (printdebug_)  edm::LogInfo("SiStripGainESProducer")<< " putting A SiStrip Gain object in eventSetup with normalization factor " << NFactor ;
  return NFactor;
}  
#endif
