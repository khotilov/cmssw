#include "RecoLocalCalo/CaloTowersCreator/src/CaloTowersCreator.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoLocalCalo/CaloTowersCreator/interface/EScales.h"


const std::vector<double>& 
CaloTowersCreator::getGridValues()
{
  static std::vector<double> retval;
  
  if (retval.size() == 0)
    {
      retval.push_back(0.);
      retval.push_back(10.);
      retval.push_back(20.);
      retval.push_back(30.);
      retval.push_back(40.);
      retval.push_back(50.);
      retval.push_back(100.);
      retval.push_back(1000.); 
    }

  return retval;
}


CaloTowersCreator::CaloTowersCreator(const edm::ParameterSet& conf) : 
  algo_(conf.getParameter<double>("EBThreshold"),
	      conf.getParameter<double>("EEThreshold"),
	      conf.getParameter<double>("HcalThreshold"),
	      conf.getParameter<double>("HBThreshold"),
	      conf.getParameter<double>("HESThreshold"),
	      conf.getParameter<double>("HEDThreshold"),
	      conf.getParameter<double>("HOThreshold"),
	      conf.getParameter<double>("HF1Threshold"),
	      conf.getParameter<double>("HF2Threshold"),
        conf.getUntrackedParameter<std::vector<double> >("EBGrid",getGridValues() ),
        conf.getUntrackedParameter<std::vector<double> >("EBWeights",std::vector<double>(getGridValues().size(),1.) ),
        conf.getUntrackedParameter<std::vector<double> >("EEGrid",getGridValues() ),
        conf.getUntrackedParameter<std::vector<double> >("EEWeights",std::vector<double>(getGridValues().size(),1.) ),
        conf.getUntrackedParameter<std::vector<double> >("HBGrid",getGridValues() ),
        conf.getUntrackedParameter<std::vector<double> >("HBWeights",std::vector<double>(getGridValues().size(),1.) ),
        conf.getUntrackedParameter<std::vector<double> >("HESGrid",getGridValues() ),
        conf.getUntrackedParameter<std::vector<double> >("HESWeights",std::vector<double>(getGridValues().size(),1.) ),
        conf.getUntrackedParameter<std::vector<double> >("HEDGrid",getGridValues() ),
        conf.getUntrackedParameter<std::vector<double> >("HEDWeights",std::vector<double>(getGridValues().size(),1.) ),
        conf.getUntrackedParameter<std::vector<double> >("HOGrid",getGridValues() ),
        conf.getUntrackedParameter<std::vector<double> >("HOWeights",std::vector<double>(getGridValues().size(),1.) ),
        conf.getUntrackedParameter<std::vector<double> >("HF1Grid",getGridValues() ),
        conf.getUntrackedParameter<std::vector<double> >("HF1Weights",std::vector<double>(getGridValues().size(),1.) ),
        conf.getUntrackedParameter<std::vector<double> >("HF2Grid",getGridValues() ),
        conf.getUntrackedParameter<std::vector<double> >("HF2Weights",std::vector<double>(getGridValues().size(),1.) ),
	      conf.getParameter<double>("EBWeight"),
	      conf.getParameter<double>("EEWeight"),
	      conf.getParameter<double>("HBWeight"),
	      conf.getParameter<double>("HESWeight"),
	      conf.getParameter<double>("HEDWeight"),
	      conf.getParameter<double>("HOWeight"),
	      conf.getParameter<double>("HF1Weight"),
	      conf.getParameter<double>("HF2Weight"),
	      conf.getParameter<double>("EcutTower"),
	      conf.getParameter<double>("EBSumThreshold"),
	      conf.getParameter<double>("EESumThreshold"),
	      conf.getParameter<bool>("UseHO"),
         // (for momentum reconstruction algorithm)
        conf.getParameter<int>("MomConstrMethod"),
        conf.getParameter<double>("MomHBDepth"),
        conf.getParameter<double>("MomHEDepth"),
        conf.getParameter<double>("MomEBDepth"),
        conf.getParameter<double>("MomEEDepth")
	),

  hbheLabel_(conf.getParameter<edm::InputTag>("hbheInput")),
  hoLabel_(conf.getParameter<edm::InputTag>("hoInput")),
  hfLabel_(conf.getParameter<edm::InputTag>("hfInput")),
  ecalLabels_(conf.getParameter<std::vector<edm::InputTag> >("ecalInputs")),
  allowMissingInputs_(conf.getUntrackedParameter<bool>("AllowMissingInputs",false))
{
  EBEScale=EScales.EBScale; 
  EEEScale=EScales.EEScale; 
  HBEScale=EScales.HBScale; 
  HESEScale=EScales.HESScale; 
  HEDEScale=EScales.HEDScale; 
  HOEScale=EScales.HOScale; 
  HF1EScale=EScales.HF1Scale; 
  HF2EScale=EScales.HF2Scale; 
  if (EScales.instanceLabel=="") produces<CaloTowerCollection>();
  else produces<CaloTowerCollection>(EScales.instanceLabel);
}

void CaloTowersCreator::produce(edm::Event& e, const edm::EventSetup& c) {
  // get the necessary event setup objects...
  edm::ESHandle<CaloGeometry> pG;
  edm::ESHandle<HcalTopology> htopo;
  edm::ESHandle<CaloTowerConstituentsMap> cttopo;
  c.get<CaloGeometryRecord>().get(pG);
  c.get<IdealGeometryRecord>().get(htopo);
  c.get<IdealGeometryRecord>().get(cttopo);
 
  algo_.setEBEScale(EBEScale);
  algo_.setEEEScale(EEEScale);
  algo_.setHBEScale(HBEScale);
  algo_.setHESEScale(HESEScale);
  algo_.setHEDEScale(HEDEScale);
  algo_.setHOEScale(HOEScale);
  algo_.setHF1EScale(HF1EScale);
  algo_.setHF2EScale(HF2EScale);
  algo_.setGeometry(cttopo.product(),htopo.product(),pG.product());

  algo_.begin(); // clear the internal buffer

  bool present;

  // Step A/C: Get Inputs and process (repeatedly)
  edm::Handle<HBHERecHitCollection> hbhe;
  present=e.getByLabel(hbheLabel_,hbhe);
  if (present || !allowMissingInputs_)  algo_.process(*hbhe);

  edm::Handle<HORecHitCollection> ho;
  present=e.getByLabel(hoLabel_,ho);
  if (present || !allowMissingInputs_) algo_.process(*ho);

  edm::Handle<HFRecHitCollection> hf;
  present=e.getByLabel(hfLabel_,hf);
  if (present || !allowMissingInputs_) algo_.process(*hf);

  std::vector<edm::InputTag>::const_iterator i;
  for (i=ecalLabels_.begin(); i!=ecalLabels_.end(); i++) {
    edm::Handle<EcalRecHitCollection> ec;
    present=e.getByLabel(*i,ec);
    if (present || !allowMissingInputs_) algo_.process(*ec);
  }

  // Step B: Create empty output
  std::auto_ptr<CaloTowerCollection> prod(new CaloTowerCollection());

  // Step C: Process
  algo_.finish(*prod);

  // Step D: Put into the event
  if (EScales.instanceLabel=="") e.put(prod);
  else e.put(prod,EScales.instanceLabel);
}

