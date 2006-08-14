#include "SimG4Core/Physics81/interface/PhysicsList.h"
#include "SimG4Core/Physics81/interface/DDG4ProductionCuts.h"

PhysicsList::PhysicsList(const edm::ParameterSet & p) 
  : G4VModularPhysicsList(), m_pPhysics(p),  prodCuts(0)
{
    //if (m_pPhysics.getParameter<bool>("CutsPerRegion")) 
      // prodCuts = new DDG4ProductionCuts();	
}
 
PhysicsList::~PhysicsList() 
{
    if (prodCuts!=0) delete prodCuts;
}

void PhysicsList::SetCuts() 
{ 

    SetDefaultCutValue(m_pPhysics.getParameter<double>("DefaultCutValue")*cm);
    SetCutsWithDefault();
    
    if ( m_pPhysics.getParameter<bool>("CutsPerRegion") )
    {
       std::cout << " Setting Production Cuts per Region" << std::endl ;
       DDG4ProductionCuts prodCuts;
       prodCuts.SetVerbosity( m_pPhysics.getUntrackedParameter<int>("Verbosity",0)) ;
       prodCuts.update();
    }

    if ( m_pPhysics.getUntrackedParameter<int>("Verbosity") > 1 ) 
	G4VUserPhysicsList::DumpCutValuesTable();

    return ;

}

