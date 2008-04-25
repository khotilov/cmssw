#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"

#include "DataFormats/EcalDigi/interface/ESDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "EventFilter/ESDigiToRaw/interface/ESDigiToRaw.h"
#include "EventFilter/ESDigiToRaw/src/ESDataFormatterV1_1.h"
#include "EventFilter/ESDigiToRaw/src/ESDataFormatterV4.h"

using namespace std;
using namespace edm;

ESDigiToRaw::ESDigiToRaw(const edm::ParameterSet& ps) : ESDataFormatter_(0)
{
  
  label_ = ps.getParameter<string>("Label");
  instanceName_ = ps.getParameter<string>("InstanceES");
  debug_ = ps.getUntrackedParameter<bool>("debugMode", false);
  formatMajor_ = ps.getUntrackedParameter<int>("formatMajor",1);
  formatMinor_ = ps.getUntrackedParameter<int>("formatMinor",1);
  lookup_ = ps.getUntrackedParameter<FileInPath>("LookupTable");

  counter_ = 0;
  kchip_ec_ = 0; 
  kchip_bc_ = 0; 

  produces<FEDRawDataCollection>();

  if (formatMajor_==4 && formatMinor_==0) 
    ESDataFormatter_ = new ESDataFormatterV4(ps);
  else 
    ESDataFormatter_ = new ESDataFormatterV1_1(ps);

  // initialize look-up table
  for (int i=0; i<2; ++i)
    for (int j=0; j<2; ++j)
      for (int k=0 ;k<40; ++k)
        for (int m=0; m<40; m++)
          fedId_[i][j][k][m] = -1;

  // read in look-up table
  int iz, ip, ix, iy, fed, kchip, pace, bundle, fiber;
  ifstream file;
  file.open(lookup_.fullPath().c_str());
  if( file.is_open() ) {
    for (int i=0; i<4288; ++i) {
      file>> iz >> ip >> ix >> iy >> fed >> kchip >> pace >> bundle >> fiber;
      fedId_[(3-iz)/2-1][ip-1][ix-1][iy-1] = fed;
    }
  } else {
    cout<<"Look up table file can not be found in "<<lookup_.fullPath().c_str()<<endl;
  }

}

ESDigiToRaw::~ESDigiToRaw() {
  if (ESDataFormatter_) delete ESDataFormatter_;
}

void ESDigiToRaw::beginJob(const edm::EventSetup& es) {
}

void ESDigiToRaw::produce(edm::Event& ev, const edm::EventSetup& es) {

  run_number_ = ev.id().run();
  orbit_number_ = counter_ / LHC_BX_RANGE;
  bx_ = (counter_ % LHC_BX_RANGE);
   
  //lv1_ = counter_;
  lv1_ = ev.id().event();
  kchip_ec_ = (lv1_ % KCHIP_EC_RANGE); 
  kchip_bc_ = (counter_ % KCHIP_BC_RANGE);
  counter_++;

  ESDataFormatter_->setRunNumber(run_number_);
  ESDataFormatter_->setOrbitNumber(orbit_number_);
  ESDataFormatter_->setBX(bx_);
  ESDataFormatter_->setLV1(lv1_);
  ESDataFormatter_->setKchipBC(kchip_bc_);
  ESDataFormatter_->setKchipEC(kchip_ec_);

  pair<int,int> ESFEDIds = FEDNumbering::getPreShowerFEDIds();

  edm::Handle<ESDigiCollection> digis;
  ev.getByLabel(label_, instanceName_, digis);

  int ifed;
  ESDataFormatter::Digis Digis;
  Digis.clear();

  for (ESDigiCollection::const_iterator it=digis->begin(); it!=digis->end(); ++it) {

    const ESDataFrame& df = *it;
    const ESDetId& detId = it->id();

    ifed = fedId_[(3-detId.zside())/2-1][detId.plane()-1][detId.six()-1][detId.siy()-1] - 1;
    if (ifed < 0) continue;

    int fedId = ESFEDIds.first + ifed;

    Digis[fedId].push_back(df);
  }

  auto_ptr<FEDRawDataCollection> productRawData( new FEDRawDataCollection );

  for (int fId=ESFEDIds.first; fId<=ESFEDIds.second; ++fId) {
    FEDRawData *rawData = ESDataFormatter_->DigiToRaw(fId, Digis);
    FEDRawData& fedRawData = productRawData->FEDData(fId); 
    fedRawData = *rawData;
    if (debug_) cout<<"FED : "<<fId<<" Data size : "<<fedRawData.size()<<" (Bytes)"<<endl;
  } 

  ev.put(productRawData);

  return;
}

void ESDigiToRaw::endJob() {
}
