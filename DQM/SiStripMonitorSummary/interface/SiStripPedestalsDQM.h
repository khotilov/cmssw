#ifndef SiStripMonitorSummary_SiStripPedestalsDQM_h
#define SiStripMonitorSummary_SiStripPedestalsDQM_h


#include "FWCore/Framework/interface/ESHandle.h"

#include "DQM/SiStripMonitorSummary/interface/SiStripBaseCondObjDQM.h"

#include "CondFormats/SiStripObjects/interface/SiStripPedestals.h"
#include "CondFormats/DataRecord/interface/SiStripPedestalsRcd.h"


class SiStripPedestalsDQM : public SiStripBaseCondObjDQM{
 
  public:
  
  SiStripPedestalsDQM(const edm::EventSetup & eSetup,
                      edm::ParameterSet const& hPSet,
                      edm::ParameterSet const& fPSet);
  
  virtual ~SiStripPedestalsDQM();
  
  void getActiveDetIds(const edm::EventSetup & eSetup);

  void fillModMEs(const std::vector<uint32_t> & selectedDetIds);
  void fillSummaryMEs(const std::vector<uint32_t> & selectedDetIds);
 	       
  void fillMEsForDet(ModMEs selModME_,uint32_t selDetId_);
  void fillMEsForLayer( /*std::map<uint32_t, ModMEs> selModMEsMap_,*/ uint32_t selDetId_);
  
  unsigned long long getCache(const edm::EventSetup & eSetup){ return eSetup.get<SiStripPedestalsRcd>().cacheIdentifier();}
  
  void getConditionObject(const edm::EventSetup & eSetup){
    eSetup.get<SiStripPedestalsRcd>().get(pedestalHandle_);
    cacheID_memory = cacheID_current;
  }

  private:
    edm::ESHandle<SiStripPedestals> pedestalHandle_;

};

#endif
