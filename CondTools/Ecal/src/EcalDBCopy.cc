#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "CondTools/Ecal/interface/EcalDBCopy.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibErrors.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibErrorsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/EcalObjects/interface/EcalWeightXtalGroups.h"
#include "CondFormats/DataRecord/interface/EcalWeightXtalGroupsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalTBWeights.h"
#include "CondFormats/DataRecord/interface/EcalTBWeightsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAlphas.h"
#include "CondFormats/DataRecord/interface/EcalLaserAlphasRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatios.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatiosRef.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRefRcd.h"


#include <vector>



EcalDBCopy::EcalDBCopy(const edm::ParameterSet& iConfig) :
  m_timetype(iConfig.getParameter<std::string>("timetype")),
  m_cacheIDs(),
  m_records()
{

  std::string container;
  std::string tag;
  std::string record;
  typedef std::vector< edm::ParameterSet > Parameters;
  Parameters toCopy = iConfig.getParameter<Parameters>("toCopy");
  for(Parameters::iterator i = toCopy.begin(); i != toCopy.end(); ++i) {
    container = i->getParameter<std::string>("container");
    record = i->getParameter<std::string>("record");
    m_cacheIDs.insert( std::make_pair(container, 0) );
    m_records.insert( std::make_pair(container, record) );
  }
  
}


EcalDBCopy::~EcalDBCopy()
{
  
}

void EcalDBCopy::analyze( const edm::Event& evt, const edm::EventSetup& evtSetup)
{

  std::string container;
  std::string record;
  typedef std::map<std::string, std::string>::const_iterator recordIter;
  for (recordIter i = m_records.begin(); i != m_records.end(); ++i) {
    container = (*i).first;
    record = (*i).second;
    if ( shouldCopy(evtSetup, container) ) {
      copyToDB(evtSetup, container);
    }
  }
  
}



bool EcalDBCopy::shouldCopy(const edm::EventSetup& evtSetup, std::string container)
{

  unsigned long long cacheID = 0;
  if (container == "EcalPedestals") {
    cacheID = evtSetup.get<EcalPedestalsRcd>().cacheIdentifier();
  } else if (container == "EcalADCToGeVConstant") {
    cacheID = evtSetup.get<EcalADCToGeVConstantRcd>().cacheIdentifier();
  } else if (container == "EcalIntercalibConstants") {
    cacheID = evtSetup.get<EcalIntercalibConstantsRcd>().cacheIdentifier();
  } else if (container == "EcalIntercalibErrors") {
    cacheID = evtSetup.get<EcalIntercalibErrorsRcd>().cacheIdentifier();
  } else if (container == "EcalGainRatios") {
    cacheID = evtSetup.get<EcalGainRatiosRcd>().cacheIdentifier();
  } else if (container == "EcalWeightXtalGroups") {
    cacheID = evtSetup.get<EcalWeightXtalGroupsRcd>().cacheIdentifier();
  } else if (container == "EcalTBWeights") {
    cacheID = evtSetup.get<EcalTBWeightsRcd>().cacheIdentifier();
  } else if (container == "EcalLaserAPDPNRatios") {
    cacheID = evtSetup.get<EcalLaserAPDPNRatiosRcd>().cacheIdentifier();
  } else if (container == "EcalLaserAPDPNRatiosRef") {
    cacheID = evtSetup.get<EcalTBWeightsRcd>().cacheIdentifier();
  } else if (container == "EcalLaserAlphas") {
    cacheID = evtSetup.get<EcalTBWeightsRcd>().cacheIdentifier();
  } else if (container == "EcalChannelStatus") {
    cacheID = evtSetup.get<EcalChannelStatusRcd>().cacheIdentifier();
  } 

  else {
    throw cms::Exception("Unknown container");
  }
  
  if (m_cacheIDs[container] == cacheID) {
    return 0;
  } else {
    m_cacheIDs[container] = cacheID;
    return 1;
  }

}



void EcalDBCopy::copyToDB(const edm::EventSetup& evtSetup, std::string container)
{
  edm::Service<cond::service::PoolDBOutputService> dbOutput;
  if ( !dbOutput.isAvailable() ) {
    throw cms::Exception("PoolDBOutputService is not available");
  }

  std::string recordName = m_records[container];

  if (container == "EcalPedestals") {
    edm::ESHandle<EcalPedestals> handle;
    evtSetup.get<EcalPedestalsRcd>().get(handle);
    const EcalPedestals* obj = handle.product();
    cout << "ped pointer is: "<< obj<< endl;
    dbOutput->createNewIOV<const EcalPedestals>( new EcalPedestals(*obj), dbOutput->endOfTime(),recordName);

  }  else if (container == "EcalADCToGeVConstant") {
    edm::ESHandle<EcalADCToGeVConstant> handle;
    evtSetup.get<EcalADCToGeVConstantRcd>().get(handle);
    const EcalADCToGeVConstant* obj = handle.product();
    cout << "adc pointer is: "<< obj<< endl;

   dbOutput->createNewIOV<const EcalADCToGeVConstant>( new EcalADCToGeVConstant(*obj), dbOutput->endOfTime(),recordName);


  }  else if (container == "EcalChannelStatus") {
    edm::ESHandle<EcalChannelStatus> handle;
    evtSetup.get<EcalChannelStatusRcd>().get(handle);
    const EcalChannelStatus* obj = handle.product();
    cout << "channel status pointer is: "<< obj<< endl;

   dbOutput->createNewIOV<const EcalChannelStatus>( new EcalChannelStatus(*obj), dbOutput->endOfTime(),recordName);


  }
else if (container == "EcalIntercalibConstants") {
    edm::ESHandle<EcalIntercalibConstants> handle;
    evtSetup.get<EcalIntercalibConstantsRcd>().get(handle);
    const EcalIntercalibConstants* obj = handle.product();
    cout << "inter pointer is: "<< obj<< endl;
   dbOutput->createNewIOV<const EcalIntercalibConstants>( new EcalIntercalibConstants(*obj), dbOutput->endOfTime(),recordName);

  } else if (container == "EcalIntercalibErrors") {
    edm::ESHandle<EcalIntercalibErrors> handle;
    evtSetup.get<EcalIntercalibErrorsRcd>().get(handle);
    const EcalIntercalibErrors* obj = handle.product();
    cout << "inter pointer is: "<< obj<< endl;
   dbOutput->createNewIOV<const EcalIntercalibErrors>( new EcalIntercalibErrors(*obj), dbOutput->endOfTime(),recordName);

  } else if (container == "EcalGainRatios") {
    edm::ESHandle<EcalGainRatios> handle;
    evtSetup.get<EcalGainRatiosRcd>().get(handle);
    const EcalGainRatios* obj = handle.product();
    cout << "gain pointer is: "<< obj<< endl;
   dbOutput->createNewIOV<const EcalGainRatios>( new EcalGainRatios(*obj), dbOutput->endOfTime(),recordName);

  } else if (container == "EcalWeightXtalGroups") {
    edm::ESHandle<EcalWeightXtalGroups> handle;
    evtSetup.get<EcalWeightXtalGroupsRcd>().get(handle);
    const EcalWeightXtalGroups* obj = handle.product();
    cout << "weight pointer is: "<< obj<< endl;
   dbOutput->createNewIOV<const EcalWeightXtalGroups>( new EcalWeightXtalGroups(*obj), dbOutput->endOfTime(),recordName);

  } else if (container == "EcalTBWeights") {
    edm::ESHandle<EcalTBWeights> handle;
    evtSetup.get<EcalTBWeightsRcd>().get(handle);
    const EcalTBWeights* obj = handle.product();
    cout << "tbweight pointer is: "<< obj<< endl;
   dbOutput->createNewIOV<const EcalTBWeights>( new EcalTBWeights(*obj), dbOutput->endOfTime(),recordName);

  } else if (container == "EcalLaserAlphas") {
    edm::ESHandle<EcalLaserAlphas> handle;
    evtSetup.get<EcalLaserAlphasRcd>().get(handle);
    const EcalLaserAlphas* obj = handle.product();
    cout << "ecalLaserAlpha pointer is: "<< obj<< endl;
   dbOutput->createNewIOV<const EcalLaserAlphas>( new EcalLaserAlphas(*obj), dbOutput->endOfTime(),recordName);

  } else if (container == "EcalLaserAPDPNRatios") {
    edm::ESHandle<EcalLaserAPDPNRatios> handle;
    evtSetup.get<EcalLaserAPDPNRatiosRcd>().get(handle);
    const EcalLaserAPDPNRatios* obj = handle.product();
    cout << "tbweight pointer is: "<< obj<< endl;
   dbOutput->createNewIOV<const EcalLaserAPDPNRatios>( new EcalLaserAPDPNRatios(*obj), dbOutput->endOfTime(),recordName);

  } else if (container == "EcalLaserAPDPNRatiosRef") {
    edm::ESHandle<EcalLaserAPDPNRatiosRef> handle;
    evtSetup.get<EcalLaserAPDPNRatiosRefRcd>().get(handle);
    const EcalLaserAPDPNRatiosRef* obj = handle.product();
    cout << "tbweight pointer is: "<< obj<< endl;
   dbOutput->createNewIOV<const EcalLaserAPDPNRatiosRef>( new EcalLaserAPDPNRatiosRef(*obj), dbOutput->endOfTime(),recordName);

  } else {
    throw cms::Exception("Unknown container");
  }

  cout<< "EcalDBCopy wrote " << recordName << endl;
}
