#include "SimG4Core/Physics/interface/DDG4ProductionCuts.h"
#include "SimG4Core/Notification/interface/SimG4Exception.h"

#include "G4ProductionCuts.hh"
#include "G4RegionStore.hh"

#include <algorithm>

using std::cout;
using std::endl;

DDG4ProductionCuts::DDG4ProductionCuts(const G4LogicalVolumeToDDLogicalPartMap& map) : map_(map), m_Verbosity(0) {
    m_KeywordRegion =  "CMSCutsRegion";
}

DDG4ProductionCuts::~DDG4ProductionCuts(){
}

/** helper function to compare parts through their name instead of comparing them
    by their pointers. 
    It's guaranteed to produce the same order in subsequent application runs,
    while pointers usually can't guarantee this
*/
bool dd_is_greater(const std::pair<G4LogicalVolume*, DDLogicalPart> & p1,
                   const std::pair<G4LogicalVolume*, DDLogicalPart> & p2)
{
  bool result = false;
  if (p1.second.name().ns() > p2.second.name().ns()) {
    result = true;
  }
  if (p1.second.name().ns() == p2.second.name().ns()) {
    if (p1.second.name().name() > p2.second.name().name()) {
      result = true;
    }
    if (p1.second.name().name() == p2.second.name().name()) {
      if (p1.first->GetName() > p2.first->GetName()) {
	result = true;
      }
    }
  }
  return result;
}		   

void DDG4ProductionCuts::update() {
  //
  // Now we can operate, loop over the Map by Martin ...
  //
  G4LogicalVolumeToDDLogicalPartMap::Vector vec = map_.all(m_KeywordRegion);
  // sort all root volumes - to get the same sequence at every run of the application.
  // (otherwise, the sequence will depend on the pointer (memory address) of the 
  // involved objects, because 'new' does no guarantee that you allways get a
  // higher (or lower) address when allocating an object of the same type ...
  sort(vec.begin(),vec.end(),&dd_is_greater);
  if ( m_Verbosity > 0 ) {
    std::cout <<" DDG4ProductionCuts (New) : starting\n"
	      <<" DDG4ProductionCuts : Got "<<vec.size()
	      <<" region roots.\n"
	      <<" DDG4ProductionCuts : List of all roots:\n";
    for ( size_t jj=0; jj<vec.size(); ++jj)
      std::cout << "   DDG4ProductionCuts : root=" << vec[jj].second.name() << std::endl;
  }

  //cout <<" DDG4ProductionCuts : returning, doing nothing! " << endl;
  //return;

  //
  // In this way I have all the DDLP which have a Region info attached
  //
  for (G4LogicalVolumeToDDLogicalPartMap::Vector::iterator tit = vec.begin();
       tit != vec.end(); tit++){
    SetProdCuts((*tit).second,(*tit).first);
  }
}

void DDG4ProductionCuts::SetProdCuts(const DDLogicalPart lpart, 
				     G4LogicalVolume* lvol ) {  
  
  if ( m_Verbosity > 0 ) 
    std::cout <<" DDG4ProductionCuts: inside SetProdCuts\n";
  
  G4Region * region = 0;
  
  std::string  regionName;
  unsigned int num= map_.toString(m_KeywordRegion,lpart,regionName);
  
  if (num != 1){
    throw SimG4Exception("DDG4ProductionCuts: Problem with Region tags.");
  }
  
  if ( m_Verbosity > 0 ) std::cout << "Using region " << regionName << "\n";

  region = GetRegion(regionName);
  region->AddRootLogicalVolume(lvol);
  
  if ( m_Verbosity > 0 )
    std::cout << " MakeRegions: added " << lvol->GetName()
	      << " to region " << region->GetName() << std::endl;
  //
  // search for production cuts
  // you must have three of them: e+ e- gamma
  //
  double gammacut;
  double electroncut;
  double positroncut;
  int temp =  map_.toDouble("ProdCutsForGamma",lpart,gammacut);
  if (temp != 1){
    throw SimG4Exception("DDG4ProductionCuts: Problem with Region tags: no/more than one ProdCutsForGamma.");
  }
  temp =  map_.toDouble("ProdCutsForElectrons",lpart,electroncut);
  if (temp != 1){
    throw SimG4Exception("DDG4ProductionCuts: Problem with Region tags: no/more than one ProdCutsForElectrons.");
  }
  temp =  map_.toDouble("ProdCutsForPositrons",lpart,positroncut);
  if (temp != 1) {
    throw SimG4Exception("DDG4ProductionCuts: Problem with Region tags: no/more than one ProdCutsForPositrons.");
  }
  //
  // For the moment I assume all of the three are set
  //
  G4ProductionCuts * prodCuts = GetProductionCuts(region);
  prodCuts->SetProductionCut( gammacut, idxG4GammaCut );
  prodCuts->SetProductionCut( electroncut, idxG4ElectronCut );
  prodCuts->SetProductionCut( positroncut, idxG4PositronCut );
  if ( m_Verbosity > 0 ) {
    std::cout << "DDG4ProductionCuts : Setting cuts for " << regionName
	      << "\n    Electrons: " << electroncut
	      << "\n    Positrons: " << positroncut
	      << "\n    Gamma    : " << gammacut << std::endl;
  }
}

G4Region * DDG4ProductionCuts::GetRegion(const std::string & regName) 
{
  G4Region * reg =  G4RegionStore::GetInstance()->FindOrCreateRegion (regName);
  return reg;
}

 G4ProductionCuts * DDG4ProductionCuts::GetProductionCuts( G4Region* reg )
{
  G4ProductionCuts * prodCuts = reg->GetProductionCuts();
  if( !prodCuts ) {
    prodCuts = new G4ProductionCuts();
    reg->SetProductionCuts(prodCuts);
  }
  return prodCuts;
}

