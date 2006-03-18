#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalDigi/interface/ESDataFrame.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloHitResponse.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloTDigitizer.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EcalSimParameterMap.h"
#include "SimCalorimetry/EcalSimAlgos/interface/ESShape.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EcalDigitizerTraits.h"
#include "SimCalorimetry/EcalSimAlgos/interface/ESElectronicsSim.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include <vector>
#include <iostream>
#include <iterator>
using namespace std;
using namespace cms;

int main() {
  // make a silly little hit in each subdetector, which should
  // correspond to a 300 keV particle
  ESDetId ESDetId(1, 1, 1, 1, 1);
  PCaloHit ESHit(ESDetId.rawId(), 0.0003, 0.);

  vector<DetId> ESDetIds;
  ESDetIds.push_back(ESDetId);

  edm::PCaloHitContainer ESHits;
  ESHits.push_back(ESHit);

  string ESName = "EcalHitsES";
  vector<string> caloDets, muonDets, trackingDets;
  caloDets.push_back(ESName);

  CrossingFrame crossingFrame(-5, 5, 25, muonDets, trackingDets, caloDets);
  crossingFrame.addSignalCaloHits(ESName, &ESHits);

  EcalSimParameterMap parameterMap;
  ESShape shape(1);

  CaloHitResponse ESResponse(&parameterMap, &shape);

  ESElectronicsSim electronicsSim(true, 3, 1, 1000, 9, 78.47);

  bool addNoise = false;
  CaloTDigitizer<ESDigitizerTraits> ESDigitizer(&ESResponse, &electronicsSim, addNoise);
  ESDigitizer.setDetIds(ESDetIds);

  auto_ptr<ESDigiCollection> ESResult(new ESDigiCollection);

  MixCollection<PCaloHit> ESHitCollection(&crossingFrame, ESName);

  ESDigitizer.run(ESHitCollection, *ESResult);

  // print out all the digis
  cout << "ES Frames" << endl;
  copy(ESResult->begin(), ESResult->end(), std::ostream_iterator<ESDataFrame>(std::cout, "\n"));

  return 0;
}


