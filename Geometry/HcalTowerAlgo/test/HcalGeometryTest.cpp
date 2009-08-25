#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalHardcodeGeometryLoader.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include <iostream>

void testTriggerGeometry() {

  HcalTrigTowerGeometry trigTowers;
  std::cout << "HCAL trigger tower eta bounds " << std::endl;
  for(int ieta = 1; ieta <= 32; ++ieta) {
    double eta1, eta2;
    trigTowers.towerEtaBounds(ieta, eta1, eta2);
    std::cout << ieta << " "  << eta1 << " " << eta2 << std::endl;
  }

  // now test some cell mappings
  HcalDetId barrelDet(HcalBarrel, 1, 1, 1);
  HcalDetId endcapDet(HcalEndcap, 29, 1, 1);
  HcalDetId forwardDet1(HcalForward, 29, 71, 1);
  HcalDetId forwardDet2(HcalForward, 29, 71, 2);
  HcalDetId forwardDet3(HcalForward, 40, 71, 1);

  typedef std::vector<HcalTrigTowerDetId> TowerDets;
  TowerDets barrelTowers = trigTowers.towerIds(barrelDet);
  TowerDets endcapTowers = trigTowers.towerIds(endcapDet);
  TowerDets forwardTowers1 = trigTowers.towerIds(forwardDet1);
  TowerDets forwardTowers2 = trigTowers.towerIds(forwardDet2);
  TowerDets forwardTowers3 = trigTowers.towerIds(forwardDet3);

  assert(barrelTowers.size() ==1);
  assert(endcapTowers.size() ==2);
  assert(forwardTowers1.size() ==1);
  assert(forwardTowers2.size() ==1);
  assert(forwardTowers3.size() ==1);

  std::cout << barrelTowers[0] << std::endl;
  std::cout << endcapTowers[0] << std::endl;
  std::cout << endcapTowers[1] << std::endl;
  std::cout << forwardTowers1[0] << std::endl;
  std::cout << forwardTowers3[0] << std::endl;

}


void testClosestCell(const HcalDetId & detId, const CaloSubdetectorGeometry * geom)
{
  const CaloCellGeometry* cell = geom->getGeometry(detId);
  HcalDetId closest = geom->getClosestCell( cell->getPosition() );


  if(closest != detId)
  {
    std::cout << "ERROR mismatch.  Original HCAL cell is "
              << detId << " while nearest is " << closest << std::endl;
  }
}

void testClosestCells() 
{
   HcalHardcodeGeometryLoader l;
   HcalHardcodeGeometryLoader::ReturnType g = l .load();
   // make sure each cel is its own closest cell
   HcalDetId barrelDet(HcalBarrel, 1, 1, 1);
   HcalDetId barrelDet2(HcalBarrel, 16, 50, 1);
   HcalDetId endcapDet1(HcalEndcap, -17, 72, 1);
   HcalDetId endcapDet2(HcalEndcap, 29, 35, 1);
   HcalDetId forwardDet1(HcalForward, 30, 71, 1);
   HcalDetId forwardDet3(HcalForward, -40, 71, 1);

   testClosestCell( barrelDet  , g );
   testClosestCell( barrelDet2 , g );
   testClosestCell( endcapDet1 , g );
   testClosestCell( endcapDet2 , g );
   testClosestCell( forwardDet1, g );
   testClosestCell( forwardDet3, g );

   std::vector<DetId> ids=g->getValidDetIds(DetId::Hcal,HcalBarrel);
   for (std::vector<DetId>::iterator i=ids.begin(); i!=ids.end(); i++) 
   {
      testClosestCell( HcalDetId(*i), g );
   }
}



int main() {

  HcalHardcodeGeometryLoader l;
  HcalHardcodeGeometryLoader::ReturnType b=l.load(DetId::Hcal,HcalBarrel);
  HcalHardcodeGeometryLoader::ReturnType e=l.load(DetId::Hcal,HcalEndcap);
  HcalHardcodeGeometryLoader::ReturnType o=l.load(DetId::Hcal,HcalOuter);
  HcalHardcodeGeometryLoader::ReturnType f=l.load(DetId::Hcal,HcalForward);

  std::cout << std::endl << " BARREL : " << std::endl;
  std::vector<DetId> ids=b->getValidDetIds(DetId::Hcal,HcalBarrel);
  for (std::vector<DetId>::iterator i=ids.begin(); i!=ids.end(); i++) {
    HcalDetId hid=(*i);
    if (hid.iphi()!=1) continue;
    const CaloCellGeometry* geom=b->getGeometry(hid);
    const CaloCellGeometry::CornersVec& corners=geom->getCorners();
    std::cout << hid << std::endl;
    for (CaloCellGeometry::CornersVec::const_iterator j=corners.begin(); j!=corners.end(); j++) {
      std::cout << "  " << *j << std::endl;
    }
  }

  std::cout << std::endl << " FORWARD : " << std::endl;
  ids=f->getValidDetIds(DetId::Hcal,HcalForward);
  for (std::vector<DetId>::iterator i=ids.begin(); i!=ids.end(); i++) {
    HcalDetId hid=(*i);
    //  if (hid.iphi()!=1 && hid.iphi()!=2 && hid.iphi()!=3) continue;
    std::cout << hid << std::endl;
    
    const CaloCellGeometry* geom=f->getGeometry(hid);
    const CaloCellGeometry::CornersVec& corners=geom->getCorners();
    for (CaloCellGeometry::CornersVec::const_iterator j=corners.begin(); j!=corners.end(); j++) {
      std::cout << "  " << *j << std::endl;
    }
  }

  std::cout << std::endl << " ENDCAP : " << std::endl;
  ids=e->getValidDetIds(DetId::Hcal,HcalEndcap);
  for (std::vector<DetId>::iterator i=ids.begin(); i!=ids.end(); i++) {
    HcalDetId hid=(*i);
    if (hid.iphi()!=1 && hid.iphi()!=2 && hid.iphi()!=3) continue;
    std::cout << hid << std::endl;
    
    const CaloCellGeometry* geom=e->getGeometry(hid);
    const CaloCellGeometry::CornersVec& corners=geom->getCorners();
    for (CaloCellGeometry::CornersVec::const_iterator j=corners.begin(); j!=corners.end(); j++) {
      std::cout << "  " << *j << std::endl;
    }
  }

  std::cout << std::endl << " OUTER : " << std::endl;
  ids=o->getValidDetIds(DetId::Hcal,HcalEndcap);
  for (std::vector<DetId>::iterator i=ids.begin(); i!=ids.end(); i++) {
    HcalDetId hid=(*i);
    if (hid.iphi()!=1 && hid.iphi()!=2 && hid.iphi()!=3) continue;
    std::cout << hid << std::endl;
    
    const CaloCellGeometry* geom=o->getGeometry(hid);
    const CaloCellGeometry::CornersVec& corners=geom->getCorners();
    for (CaloCellGeometry::CornersVec::const_iterator j=corners.begin(); j!=corners.end(); j++) {
      std::cout << "  " << *j << std::endl;
    }
  }

  testTriggerGeometry();

  testClosestCells();
  return 0;
}
