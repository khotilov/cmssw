#include "SimCalorimetry/HcalSimAlgos/interface/HPDNoiseMaker.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HPDNoiseData.h"

int main () {
  HPDNoiseMaker maker ("hpdNoiseLibrary.root");
  maker.addHpd ("HPD01", 0.01);
  maker.addHpd ("HPD02", 0.02);
  maker.addHpd ("HPD03", 0.03);

  HcalDetId id;
  float data[10];

  HPDNoiseData event;
  for (size_t i = 0; i < 10; i++) data[i] = i;
  id = HcalDetId (HcalBarrel, 1, 1, 1);
  event.addChannel (id, data);
  id = HcalDetId (HcalBarrel, 1, 2, 1);
  event.addChannel (id, data);
  id = HcalDetId (HcalBarrel, 1, 3, 1);
  event.addChannel (id, data);
  
  maker.newHpdEvent ("HPD01", event);

  for (size_t i = 0; i < 10; i++) data[i] = i*10;
  id = HcalDetId (HcalBarrel, 2, 1, 1);
  event.addChannel (id, data);
  id = HcalDetId (HcalBarrel, 2, 2, 1);
  event.addChannel (id, data);
  id = HcalDetId (HcalBarrel, 2, 3, 1);
  event.addChannel (id, data);
  
  maker.newHpdEvent ("HPD02", event);
  maker.newHpdEvent ("HPD02", event);

  for (size_t i = 0; i < 10; i++) data[i] = i*100;
  id = HcalDetId (HcalBarrel, 3, 1, 1);
  event.addChannel (id, data);
  id = HcalDetId (HcalBarrel, 3, 2, 1);
  event.addChannel (id, data);
  id = HcalDetId (HcalBarrel, 3, 3, 1);
  event.addChannel (id, data);
  
  maker.newHpdEvent ("HPD03", event);
  maker.newHpdEvent ("HPD03", event);
  maker.newHpdEvent ("HPD03", event);

  return 0;
}
