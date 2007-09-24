#ifndef ECAL_FENIXSTRIP_H
#define ECAL_FENIXSTRIP_H

#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixLinearizer.h>
#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixAmplitudeFilter.h>
#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixPeakFinder.h>
#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixEtStrip.h>
#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixLinearizer.h>
#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixEtStrip.h>
#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixStripFgvbEE.h>

#include "CondFormats/DataRecord/interface/EcalTPGLinearizationConstRcd.h"
#include "CondFormats/DataRecord/interface/EcalTPGPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTPGSlidingWindowRcd.h"
#include "CondFormats/DataRecord/interface/EcalTPGWeightIdMapRcd.h"
#include "CondFormats/DataRecord/interface/EcalTPGWeightGroupRcd.h"
#include "CondFormats/DataRecord/interface/EcalTPGFineGrainStripEERcd.h"
#include "CondFormats/EcalObjects/interface/EcalTPGLinearizationConst.h"
#include "CondFormats/EcalObjects/interface/EcalTPGPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalTPGSlidingWindow.h"
#include "CondFormats/EcalObjects/interface/EcalTPGWeightIdMap.h"
#include "CondFormats/EcalObjects/interface/EcalTPGWeightGroup.h"
#include "CondFormats/EcalObjects/interface/EcalTPGFineGrainStripEE.h"

#include "DataFormats/EcalDetId/interface/EcalTriggerElectronicsId.h"
#include <DataFormats/EcalDigi/interface/EBDataFrame.h>
#include <DataFormats/EcalDigi/interface/EEDataFrame.h>
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"

class EBDataFrame;
class EcalTriggerPrimitiveSample;
class EcalTPGPedestals;
class EcalTPGSlidingWindow;
class EcalTPGFineGrainStripEE;
class EcalFenixStripFgvbEE;
class EcalFenixStripFormatEB;
class EcalFenixStripFormatEE;

/** 
    \class EcalFenixStrip
    \brief class representing the Fenix chip, format strip
*/
class EcalFenixStrip {
 public:

  // constructor, destructor
  EcalFenixStrip(const edm::EventSetup& setup, const EcalElectronicsMapping* theMapping,bool debug,bool famos,int maxNrSamples, int nbMaxXtals);
  virtual ~EcalFenixStrip() ;

 private:
 const EcalElectronicsMapping* theMapping_;

  bool debug_;
  bool famos_;
  int nbMaxXtals_;
    
  std::vector <EcalFenixLinearizer *> linearizer_; 

  EcalFenixAmplitudeFilter *amplitude_filter_; 

  EcalFenixPeakFinder *peak_finder_; 
    
  EcalFenixStripFormatEB *fenixFormatterEB_;   
    
  EcalFenixStripFormatEE *fenixFormatterEE_;
    

  EcalFenixEtStrip *adder_;
    
  EcalFenixStripFgvbEE *fgvbEE_;

  // data formats for each event
  std::vector<std::vector<int> > lin_out_;
  std::vector<int> add_out_;
  std::vector<int> filt_out_;
  std::vector<int> peak_out_;
  std::vector<int> format_out_;
  std::vector<int> fgvb_out_;

 public:

  // main methods
  // process method is splitted in 2 parts:
  //   the first one is templated, the same except input
  //   the second part is slightly different for barrel/endcap
  template <class T> 
  void process(const edm::EventSetup&, std::vector<const T> &, int nrxtals, std::vector<int> & out);
  void process_part2_barrel(uint32_t stripid,const EcalTPGSlidingWindow * ecaltpgSlidW);

  void process_part2_endcap(uint32_t stripid, const EcalTPGSlidingWindow * ecaltpgSlidW,const EcalTPGFineGrainStripEE * ecaltpgFgStripEE);


  // getters for the algorithms  ;

  EcalFenixLinearizer *getLinearizer (int i) const { return linearizer_[i];}
  EcalFenixEtStrip *getAdder() const { return  adder_;}
  EcalFenixAmplitudeFilter *getFilter() const { return amplitude_filter_;}
  EcalFenixPeakFinder *getPeakFinder() const { return peak_finder_;}
    
  EcalFenixStripFormatEB *getFormatterEB() const { return fenixFormatterEB_;}
  EcalFenixStripFormatEE *getFormatterEE() const { return fenixFormatterEE_;}
    
  EcalFenixStripFgvbEE *getFGVB()      const { return fgvbEE_;}

  // ========================= implementations ==============================================================
  void process(const edm::EventSetup &setup, std::vector<EBDataFrame> &samples, int nrXtals,std::vector<int> &out){

    // first get parameter records
    //for xtals
    edm::ESHandle<EcalTPGLinearizationConst> theEcalTPGLinearization_handle;
    setup.get<EcalTPGLinearizationConstRcd>().get(theEcalTPGLinearization_handle);
    const EcalTPGLinearizationConst * ecaltpLin = theEcalTPGLinearization_handle.product();
    edm::ESHandle<EcalTPGPedestals> theEcalTPGPedestals_handle;
    setup.get<EcalTPGPedestalsRcd>().get(theEcalTPGPedestals_handle);
    const EcalTPGPedestals * ecaltpPed = theEcalTPGPedestals_handle.product();

    //for strips
    edm::ESHandle<EcalTPGSlidingWindow> theEcalTPGSlidingWindow_handle;
    setup.get<EcalTPGSlidingWindowRcd>().get(theEcalTPGSlidingWindow_handle);
    const EcalTPGSlidingWindow * ecaltpgSlidW = theEcalTPGSlidingWindow_handle.product();
    edm::ESHandle<EcalTPGWeightIdMap> theEcalTPGWEightIdMap_handle;
    setup.get<EcalTPGWeightIdMapRcd>().get(theEcalTPGWEightIdMap_handle);
    const EcalTPGWeightIdMap * ecaltpgWeightMap = theEcalTPGWEightIdMap_handle.product();
    edm::ESHandle<EcalTPGWeightGroup> theEcalTPGWEightGroup_handle;
    setup.get<EcalTPGWeightGroupRcd>().get(theEcalTPGWEightGroup_handle);
    const EcalTPGWeightGroup * ecaltpgWeightGroup = theEcalTPGWEightGroup_handle.product();

    // now call processing
    if (samples.size()==0) {
      std::cout<<" Warning: 0 size vector found in EcalFenixStripProcess!!!!!"<<std::endl;
      return;

    }
    const EcalTriggerElectronicsId elId = theMapping_->getTriggerElectronicsId(samples[0].id());
    uint32_t stripid=elId.rawId() & 0xfffffff8;   //from Pascal
    process_part1(samples,nrXtals,stripid,ecaltpPed,ecaltpLin,ecaltpgWeightMap,ecaltpgWeightGroup);//templated part
    process_part2_barrel(stripid,ecaltpgSlidW);//part different for barrel/endcap
    out=format_out_;
  }

 void  process(const edm::EventSetup &setup, std::vector<EEDataFrame> &samples, int nrXtals, std::vector<int> & out){

  // first get parameter records
  //for xtals
  edm::ESHandle<EcalTPGLinearizationConst> theEcalTPGLinearization_handle;
  setup.get<EcalTPGLinearizationConstRcd>().get(theEcalTPGLinearization_handle);
  const EcalTPGLinearizationConst * ecaltpLin = theEcalTPGLinearization_handle.product();
  edm::ESHandle<EcalTPGPedestals> theEcalTPGPedestals_handle;
  setup.get<EcalTPGPedestalsRcd>().get(theEcalTPGPedestals_handle);
  const EcalTPGPedestals * ecaltpPed = theEcalTPGPedestals_handle.product();

  //for strips
  edm::ESHandle<EcalTPGSlidingWindow> theEcalTPGSlidingWindow_handle;
  setup.get<EcalTPGSlidingWindowRcd>().get(theEcalTPGSlidingWindow_handle);
  const EcalTPGSlidingWindow * ecaltpgSlidW = theEcalTPGSlidingWindow_handle.product();
  edm::ESHandle<EcalTPGWeightIdMap> theEcalTPGWEightIdMap_handle;
  setup.get<EcalTPGWeightIdMapRcd>().get(theEcalTPGWEightIdMap_handle);
  const EcalTPGWeightIdMap * ecaltpgWeightMap = theEcalTPGWEightIdMap_handle.product();
  edm::ESHandle<EcalTPGWeightGroup> theEcalTPGWEightGroup_handle;
  setup.get<EcalTPGWeightGroupRcd>().get(theEcalTPGWEightGroup_handle);
  const EcalTPGWeightGroup * ecaltpgWeightGroup = theEcalTPGWEightGroup_handle.product();
  edm::ESHandle<EcalTPGFineGrainStripEE> theEcalTPGFineGrainStripEE_handle;
  setup.get<EcalTPGFineGrainStripEERcd>().get(theEcalTPGFineGrainStripEE_handle);
  const EcalTPGFineGrainStripEE * ecaltpgFgStripEE = theEcalTPGFineGrainStripEE_handle.product();

// now call processing
   if (samples.size()==0) {
     std::cout<<" Warning: 0 size vector found in EcalFenixStripProcess!!!!!"<<std::endl;
     return;
   }
   const EcalTriggerElectronicsId elId = theMapping_->getTriggerElectronicsId(samples[0].id());
   uint32_t stripid=elId.rawId() & 0xfffffff8;   //from Pascal
   process_part1(samples,nrXtals,stripid,ecaltpPed,ecaltpLin,ecaltpgWeightMap,ecaltpgWeightGroup); //templated part
   process_part2_endcap(stripid,ecaltpgSlidW,ecaltpgFgStripEE);
   out=format_out_;
   return;
 }

 template <class T> 
 void  process_part1(std::vector<T> & df,int nrXtals, uint32_t stripid, const EcalTPGPedestals * ecaltpPed, const EcalTPGLinearizationConst *ecaltpLin,const EcalTPGWeightIdMap * ecaltpgWeightMap,const EcalTPGWeightGroup * ecaltpgWeightGroup)
   {
  
      if(debug_)  std::cout<<"\n\nEcalFenixStrip input is a vector of size: "<<nrXtals<< std::endl;

      //loop over crystals
      for (int cryst=0;cryst<nrXtals;cryst++) {
	if(debug_){
	  std::cout<<std::endl;
	  std::cout <<"cryst= "<<cryst<<" EBDataFrame/EEDataFrame is: "<<std::endl; 
	  for ( int i = 0; i<df[cryst].size();i++){
	    std::cout <<" "<<std::dec<<df[cryst][i].adc();
	  }
	  std::cout<<std::endl;
	}
	// call linearizer
	this->getLinearizer(cryst)->setParameters(df[cryst].id().rawId(),ecaltpPed,ecaltpLin) ; 
	this->getLinearizer(cryst)->process(df[cryst],lin_out_[cryst]);
      }

      if(debug_){
	std::cout<< "output of linearizer is a vector of size: "
              <<std::dec<<lin_out_.size()<<" of which used "<<nrXtals<<std::endl; 
	for (int ix=0;ix<nrXtals;ix++){
	  std::cout<< "cryst: "<<ix<<"  value : "<<std::dec<<std::endl;
	  std::cout<<" lin_out[ix].size()= "<<std::dec<<lin_out_[ix].size()<<std::endl;
	  for (unsigned int i =0; i<lin_out_[ix].size();i++){
	    std::cout <<" "<<std::dec<<(lin_out_[ix])[i];
	  }
	  std::cout<<std::endl;
	}
    
	std::cout<<std::endl;
      }
  
      // call adder
      this->getAdder()->process(lin_out_,nrXtals,add_out_);  //add_out is of size SIZEMAX=maxNrSamples
 
      if(debug_){
	std::cout<< "output of adder is a vector of size: "<<std::dec<<add_out_.size()<<std::endl; 
	for (unsigned int ix=0;ix<add_out_.size();ix++){
	  std::cout<< "cryst: "<<ix<<"  value : "<<std::dec<<add_out_[ix]<<std::endl;
	}
	std::cout<<std::endl;
      }
 

      if (famos_) {
	filt_out_[0]= add_out_[0];
	peak_out_[0]= add_out_[0];
	return;
      }else {
	// call amplitudefilter
	this->getFilter()->setParameters(stripid,ecaltpgWeightMap,ecaltpgWeightGroup); 
	this->getFilter()->process(add_out_,filt_out_); 

	if(debug_){
	  std::cout<< "output of filter is a vector of size: "<<std::dec<<filt_out_.size()<<std::endl; 
	  for (unsigned int ix=0;ix<filt_out_.size();ix++){
	    std::cout<< "cryst: "<<ix<<"  value : "<<std::dec<<filt_out_[ix]<<std::endl;
	  }
	  std::cout<<std::endl;
	}

	// call peakfinder
	this->getPeakFinder()->process(filt_out_,peak_out_);
	if(debug_){
	  std::cout<< "output of peakfinder is a vector of size: "<<peak_out_.size()<<std::endl; 
	  for (unsigned int ix=0;ix<peak_out_.size();ix++){
	    std::cout<< "cryst: "<<ix<<"  value : "<<peak_out_[ix]<<std::endl;
	  }
	  std::cout<<std::endl;
	}
	return;
      }
   }

};
#endif

