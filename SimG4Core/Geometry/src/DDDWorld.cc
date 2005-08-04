#include "FWCore/Framework/interface/EventSetupProvider.h"
#include "FWCore/Framework/interface/recordGetImplementation.icc"
 
#include "SimG4Core/Geometry/interface/DDDWorld.h"
#include "SimG4Core/Geometry/interface/DDG4Builder.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"

#include "G4RunManagerKernel.hh"
#include "G4PVPlacement.hh"
 
using namespace edm;

DDDWorld::DDDWorld(const DDCompactView* cpv) 
{
    std::auto_ptr<DDG4Builder> theBuilder(new DDG4Builder(cpv));

    G4LogicalVolume * world = theBuilder->BuildGeometry();
    G4VPhysicalVolume * pv = 
	new G4PVPlacement(0,G4ThreeVector(),world,"DDDWorld",0,false,0);
    SetAsWorld(pv);
}

DDDWorld::~DDDWorld() {}

void DDDWorld::SetAsWorld(G4VPhysicalVolume * pv)
{
    G4RunManagerKernel * kernel = G4RunManagerKernel::GetRunManagerKernel();
    if (kernel != 0) kernel->DefineWorldVolume(pv);
    std::cout << " World volume defined " << std::endl;
}

