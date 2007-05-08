#include "FWCore/Framework/interface/EventSetupProvider.h"
#include "FWCore/Framework/interface/recordGetImplementation.icc"
 
#include "SimG4Core/Geometry/interface/DDDWorld.h"
#include "SimG4Core/Geometry/interface/DDG4Builder.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"

#include "G4RunManagerKernel.hh"
#include "G4PVPlacement.hh"
 
using namespace edm;

DDDWorld::DDDWorld(const DDCompactView* cpv, 
		   G4LogicalVolumeToDDLogicalPartMap & map,
		   SensitiveDetectorCatalog & catalog,
		   bool check) {

  std::auto_ptr<DDG4Builder> theBuilder(new DDG4Builder(cpv, check));

  DDGeometryReturnType ret = theBuilder->BuildGeometry();
  G4LogicalVolume *    world = ret.logicalVolume();
  G4VPhysicalVolume *  pv = 
    new G4PVPlacement(0,G4ThreeVector(),world,"DDDWorld",0,false,0);
  SetAsWorld(pv);
  map     = ret.lvToDDLPMap();
  catalog = ret.sdCatalog();
}

DDDWorld::~DDDWorld() {}

void DDDWorld::SetAsWorld(G4VPhysicalVolume * pv) {
  G4RunManagerKernel * kernel = G4RunManagerKernel::GetRunManagerKernel();
  if (kernel != 0) kernel->DefineWorldVolume(pv);
  std::cout << " World volume defined " << std::endl;
}

