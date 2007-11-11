#include "L1Trigger/RegionalCaloTrigger/interface/L1RCTInputProducer.h" 

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"

#include "CondFormats/L1TObjects/interface/L1CaloEtScale.h"
#include "CondFormats/DataRecord/interface/L1EmEtScaleRcd.h"

#include "CondFormats/L1TObjects/interface/L1RCTParameters.h"
#include "CondFormats/DataRecord/interface/L1RCTParametersRcd.h"

#include "L1Trigger/RegionalCaloTrigger/interface/L1RCT.h"
#include "L1Trigger/RegionalCaloTrigger/interface/L1RCTLookupTables.h" 

#include <vector>
using std::vector;

#include <iostream>
using std::cout;
using std::endl;

L1RCTInputProducer::L1RCTInputProducer(const edm::ParameterSet& conf) : 
  rctLookupTables(new L1RCTLookupTables),
  rct(new L1RCT(rctLookupTables)),
  useEcal(conf.getParameter<bool>("useEcal")),
  useHcal(conf.getParameter<bool>("useHcal")),
  ecalDigisLabel(conf.getParameter<edm::InputTag>("ecalDigisLabel")),
  hcalDigisLabel(conf.getParameter<edm::InputTag>("hcalDigisLabel"))
{
  produces<vector<unsigned short> >("rctCrate");
  produces<vector<unsigned short> >("rctCard");
  produces<vector<unsigned short> >("rctTower");
  produces<vector<unsigned int> >("rctEGammaET");
  produces<vector<bool> >("rctHoEFGVetoBit");
  produces<vector<unsigned int> >("rctJetMETET");
  produces<vector<bool> >("rctTowerActivityBit");
  produces<vector<bool> >("rctTowerMIPBit");
  produces<vector<unsigned short> >("rctHFCrate");
  produces<vector<unsigned short> >("rctHFRegion");
  produces<vector<unsigned int> >("rctHFET");
  produces<vector<bool> >("rctHFFG");
}

L1RCTInputProducer::~L1RCTInputProducer()
{
  if(rct != 0) delete rct;
  if(rctLookupTables != 0) delete rctLookupTables;
}

void L1RCTInputProducer::beginJob(const edm::EventSetup& eventSetup)
{
}

void L1RCTInputProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{

  // Refresh configuration information every event
  // Hopefully, this does not take too much time
  // There should be a call back function in future to
  // handle changes in configuration

  edm::ESHandle<L1RCTParameters> rctParameters;
  eventSetup.get<L1RCTParametersRcd>().get(rctParameters);
  const L1RCTParameters* r = rctParameters.product();
  edm::ESHandle<CaloTPGTranscoder> transcoder;
  eventSetup.get<CaloTPGRecord>().get(transcoder);
  const CaloTPGTranscoder* t = transcoder.product();
  edm::ESHandle<L1CaloEtScale> emScale;
  eventSetup.get<L1EmEtScaleRcd>().get(emScale);
  const L1CaloEtScale* s = emScale.product();
  rctLookupTables->setRCTParameters(r);
  rctLookupTables->setTranscoder(t);
  rctLookupTables->setL1CaloEtScale(s);

  edm::Handle<EcalTrigPrimDigiCollection> ecal;
  edm::Handle<HcalTrigPrimDigiCollection> hcal;
  
  if (useEcal) { event.getByLabel(ecalDigisLabel, ecal); }
  if (useHcal) { event.getByLabel(hcalDigisLabel, hcal); }

  EcalTrigPrimDigiCollection ecalColl;
  HcalTrigPrimDigiCollection hcalColl;
  if (ecal.isValid()) { ecalColl = *ecal; }
  if (hcal.isValid()) { hcalColl = *hcal; }

  rct->digiInput(ecalColl, hcalColl);

  // Stuff to create

  std::auto_ptr<std::vector<unsigned short> > 
    rctCrate(new vector<unsigned short>);
  std::auto_ptr<std::vector<unsigned short> > 
    rctCard(new vector<unsigned short>);
  std::auto_ptr<std::vector<unsigned short> > 
    rctTower(new vector<unsigned short>);
  std::auto_ptr<std::vector<unsigned int> > 
    rctEGammaET(new vector<unsigned int>);
  std::auto_ptr<std::vector<bool> > rctHoEFGVetoBit(new vector<bool>);
  std::auto_ptr<std::vector<unsigned int> > 
    rctJetMETET(new vector<unsigned int>);
  std::auto_ptr<std::vector<bool> > rctTowerActivityBit(new vector<bool>);
  std::auto_ptr<std::vector<bool> > rctTowerMIPBit(new vector<bool>);
  
  for(int crate = 0; crate < 18; crate++) {
    for(int card = 0; card < 7; card++) {
      for(int tower = 0; tower < 32; tower++) {
	unsigned short ecalCompressedET = 
	  rct->ecalCompressedET(crate, card, tower);
	unsigned short ecalFineGrainBit = 
	  rct->ecalFineGrainBit(crate, card, tower);
	unsigned short hcalCompressedET = 
	  rct->hcalCompressedET(crate, card, tower);
	unsigned int lutBits = 
	  rctLookupTables->lookup(ecalCompressedET, hcalCompressedET, 
				  ecalFineGrainBit, crate, card, tower);
	unsigned int eGammaETCode = lutBits & 0x0000007F;
	bool hOeFGVetoBit = (lutBits >> 7) & 0x00000001;
	unsigned int jetMETETCode = (lutBits >> 8) & 0x000001FF;
	bool activityBit = (lutBits >> 17) & 0x00000001;
	if(eGammaETCode > 0 || jetMETETCode > 0 || 
	   hOeFGVetoBit || activityBit) {
	  rctCrate->push_back(crate);
	  rctCard->push_back(card);
	  rctTower->push_back(tower);
	  rctEGammaET->push_back(eGammaETCode);
	  rctHoEFGVetoBit->push_back(hOeFGVetoBit);
	  rctJetMETET->push_back(jetMETETCode);
	  rctTowerActivityBit->push_back(activityBit);
	  rctTowerMIPBit->push_back(0); // FIXME: MIP bit is not yet defined
	}
      }
    }
  }

  std::auto_ptr<std::vector<unsigned short> > 
    rctHFCrate(new vector<unsigned short>);
  std::auto_ptr<std::vector<unsigned short> > 
    rctHFRegion(new vector<unsigned short>);
  std::auto_ptr<std::vector<unsigned int> > rctHFET(new vector<unsigned int>);
  std::auto_ptr<std::vector<bool> > rctHFFG(new vector<bool>);
  for(int crate = 0; crate < 18; crate++) {
    for(int hfRegion = 0; hfRegion < 8; hfRegion++) {
      unsigned short hfCompressedET = rct->hfCompressedET(crate, hfRegion);
      unsigned int hfETCode = 
	rctLookupTables->lookup(hfCompressedET, crate, 999, hfRegion);
      if(hfETCode > 0) {
	rctHFCrate->push_back(crate);
	rctHFRegion->push_back(hfRegion);
	rctHFET->push_back(hfETCode);
	rctHFFG->push_back(0); // FIXME: HF FG is not yet defined
      }
    }
  }

  //putting stuff back into event
  event.put(rctCrate, "rctCrate");
  event.put(rctCard, "rctCard");
  event.put(rctTower, "rctTower");
  event.put(rctEGammaET, "rctEGammaET");
  event.put(rctHoEFGVetoBit, "rctHoEFGVetoBit");
  event.put(rctJetMETET, "rctJetMETET");
  event.put(rctTowerActivityBit, "rctTowerActivityBit");
  event.put(rctTowerMIPBit, "rctTowerMIPBit");
  event.put(rctHFCrate, "rctHFCrate");
  event.put(rctHFRegion, "rctHFRegion");
  event.put(rctHFET, "rctHFET");
  event.put(rctHFFG, "rctHFFG");

}
