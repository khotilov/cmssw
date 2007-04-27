/** \class EcalTrigPrimFunctionalAlgo
 *
 * EcalTrigPrimFunctionalAlgo is the main algorithm class for TPG
 * It coordinates all the aother algorithms
 * Structure is very close to electronics
 *
 *
 * \author Ursula Berthon, Stephanie Baffioni,  LLR Palaiseau
 *
 * \version   1st Version may 2006
 * \version   2nd Version jul 2006

 *
 ************************************************************/
#include <string>
#include <algorithm>
#include <numeric>
#include <functional>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"

#include "SimCalorimetry/EcalTrigPrimAlgos/interface/EcalTrigPrimFunctionalAlgo.h"
#include "SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixLinearizer.h"
#include "SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixStrip.h"
#include "SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixTcp.h"

//#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "DataFormats/EcalDetId/interface/EcalTriggerElectronicsId.h"

#include "CondFormats/DataRecord/interface/EcalTPParametersRcd.h"

#include <TTree.h>
#include <TMath.h>
//----------------------------------------------------------------------

EcalTrigPrimFunctionalAlgo::EcalTrigPrimFunctionalAlgo(const edm::EventSetup & setup,int binofmax,int nrsamples, bool tcpFormat, bool barrelOnly,bool debug, double ebDccAdcToGeV,double eeDccAdcToGeV):
  valid_(false),valTree_(NULL),binOfMaximum_(binofmax),nrSamplesToWrite_(nrsamples),
  tcpFormat_(tcpFormat), barrelOnly_(barrelOnly), debug_(debug)

{this->init(setup);}

//----------------------------------------------------------------------
EcalTrigPrimFunctionalAlgo::EcalTrigPrimFunctionalAlgo(const edm::EventSetup & setup,TTree *tree,int binofmax, int nrsamples,bool tcpFormat, bool barrelOnly, bool debug, double ebDccAdcToGeV,double eeDccAdcToGeV):
  valid_(true),valTree_(tree),binOfMaximum_(binofmax),nrSamplesToWrite_(nrsamples),
  tcpFormat_(tcpFormat), barrelOnly_(barrelOnly),debug_(debug)

{this->init(setup);}

//----------------------------------------------------------------------
void EcalTrigPrimFunctionalAlgo::init(const edm::EventSetup & setup) {
  //FIXME: check validities
  if (!barrelOnly_) {
    edm::ESHandle<CaloGeometry> theGeometry;
    edm::ESHandle<CaloSubdetectorGeometry> theEndcapGeometry_handle;
    setup.get<IdealGeometryRecord>().get( theGeometry );
    setup.get<IdealGeometryRecord>().get("EcalEndcap",theEndcapGeometry_handle);
    theEndcapGeometry = &(*theEndcapGeometry_handle);
    setup.get<IdealGeometryRecord>().get(eTTmap_);
  }
  edm::ESHandle<EcalTPParameters> theEcalTPParameters_handle;
  setup.get<EcalTPParametersRcd>().get(theEcalTPParameters_handle);
  ecaltpp_=theEcalTPParameters_handle.product();


  // endcap mapping
  edm::ESHandle< EcalElectronicsMapping > ecalmapping;
  setup.get< EcalMappingRcd >().get(ecalmapping);
  theMapping_ = ecalmapping.product();

  //create main sub algos
  estrip_= new EcalFenixStrip(valTree_,ecaltpp_,theMapping_,debug_);
  etcp_ = new EcalFenixTcp(ecaltpp_,tcpFormat_,debug_) ;
}
//----------------------------------------------------------------------

EcalTrigPrimFunctionalAlgo::~EcalTrigPrimFunctionalAlgo() 
{
    delete estrip_;
    delete etcp_;
}
//----------------------------------------------------------------------
void EcalTrigPrimFunctionalAlgo::updateESRecord(double ttfLowEB, double ttfHighEB, double ttfLowEE, double ttfHighEE)
{
  const_cast <EcalTPParameters *> (ecaltpp_)->changeThresholds(ttfLowEB, ttfHighEB, ttfLowEE, ttfHighEE);
}
//----------------------------------------------------------------------
int EcalTrigPrimFunctionalAlgo::findTowerNrInTcc(const EcalTrigTowerDetId &id)
{
  if (id.subDet()== EcalBarrel) { // finds tower nr in TCC   
   const int nrphis=4;
   int ieta=id.ietaAbs();
    int iphi=id.iphi();
    int basenr=(ieta-1)*nrphis +1;
    int towernr=basenr+(iphi-1)%nrphis;
    return  towernr;
  } 
  else if (id.subDet()==EcalEndcap) {
    return theMapping_->iTT(id);
  }
  else {
    LogDebug("EcalTPG")<<"Wrong EcalTrigTowerDetId ";
    return 0;
  }
}
//----------------------------------------------------------------------
int EcalTrigPrimFunctionalAlgo::findTccNr(const EcalTrigTowerDetId &id)
{
// finds Tcc Nr
  if (id.subDet()== EcalBarrel) { 
    return EcalTPParameters::nrMinTccEB_; //FIXME
  }
  else if (id.subDet()==EcalEndcap) {
    return theMapping_->TCCid(id);
  }
  else {
    LogDebug("EcalTPG")<<"Wrong EcalTrigTowerDetId ";
    return 0;
  }     
} 
//----------------------------------------------------------------------
int  EcalTrigPrimFunctionalAlgo::findStripNr(const EBDetId &id){
      int stripnr;
      int n=((id.ic()-1)%100)/20; //20 corresponds to 4 * ecal_barrel_crystals_per_strip FIXME!!
      if (id.ieta()<0) stripnr = n+1;
      //      else stripnr =ecal_barrel_strips_per_trigger_tower - n; 
      else stripnr =EcalTPParameters::nbMaxStrips_ - n; 
      return stripnr;
}
//----------------------------------------------------------------------
int  EcalTrigPrimFunctionalAlgo::findStripNr(const EEDetId &id){
      int stripnr;
      const EcalTriggerElectronicsId elId = theMapping_->getTriggerElectronicsId(id);
      stripnr=elId.pseudoStripId();
      return stripnr;
}
