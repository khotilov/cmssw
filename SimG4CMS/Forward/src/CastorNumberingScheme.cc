///////////////////////////////////////////////////////////////////////////////
// File: CastorNumberingScheme.cc
// Description: Numbering scheme for Castor
///////////////////////////////////////////////////////////////////////////////
#include "SimG4CMS/Forward/interface/CastorNumberingScheme.h"
#include "DataFormats/HcalDetId/interface/HcalCastorDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include <iostream>

#define debug

CastorNumberingScheme::CastorNumberingScheme() {
  edm::LogInfo("ForwardSim") << "Creating CastorNumberingScheme";
}

CastorNumberingScheme::~CastorNumberingScheme() {
  edm::LogInfo("ForwardSim") << "Deleting CastorNumberingScheme";
}

uint32_t CastorNumberingScheme::getUnitID(const G4Step* aStep) const {

  uint32_t intindex = 0;

  uint32_t index=0;
  int      level = detectorLevel(aStep);

#ifdef debug
  LogDebug("ForwardSim") << "CastorNumberingScheme number of levels= " <<level;
#endif

  if (level > 0) {
    int*      copyno = new int[level];
    G4String* name   = new G4String[level];

    detectorLevel(aStep, level, copyno, name);

    int zside   = 0;
    int sector  = 0;
    int module = 0;

    //    HcalCastorDetId::Section section;

    for (int ich=0; ich  <  level; ich++) {
      if(name[ich] == "CAST") {
	// Z index +Z = 1 ; -Z = 2
	zside   = copyno[ich];
      } else if(name[ich] == "CAES" || name[ich] == "CEDS") {
	// sector number for dead material 1 - 8
	//section = HcalCastorDetId::EM;
	sector = copyno[ich];
      } else if(name[ich] == "CAHS" || name[ich] == "CHDS") {
	// sector number for dead material 1 - 8
	//section = HcalCastorDetId::HAD;
	sector = copyno[ich];
      } else if(name[ich] == "CAER" || name[ich] == "CEDR") {
	// zmodule number 1-2 for EM section (2 copies)
	module = copyno[ich];
      } else if(name[ich] == "CAHR" || name[ich] == "CHDR") {
	//zmodule number 3-14 for HAD section (12 copies)
	module = copyno[ich] + 2;  
      } else if(name[ich] == "C3EF" || name[ich] == "C3HF") {
	// sector number for sensitive material 1 - 16
	sector = sector*2 - 1 ;
      } else if(name[ich] == "C4EF" || name[ich] == "C4HF") {
	// sector number for sensitive material 1 - 16
	sector = sector*2 ;
      }
    
    
#ifdef debug
      LogDebug("ForwardSim") << "CastorNumberingScheme  " << "ich = " << ich  
			     << "copyno" << copyno[ich] << "name = " 
			     << name[ich];
#endif
    }
    // use for Castor number 9 
    // 
    // Z index +Z = 1 ; -Z = 2
    // sector number 1 - 16
    // zmodule number  1 - 18


    //    int det = 9; 
    //    intindex = packIndex (det, zside, sector, zmodule);

    //intindex = packIndex (section, zside, sector, module);

    intindex = packIndex(zside, sector, module);

    bool true_for_positive_eta;
    if(zside == 1)true_for_positive_eta = true;
    if(zside == -1)true_for_positive_eta = false;

    HcalCastorDetId castorId = HcalCastorDetId(true_for_positive_eta, sector, module);
    index = castorId.rawId();

#ifdef debug
    LogDebug("ForwardSim") << "CastorNumberingScheme :" <<" zside "
			   << zside << " sector " << sector << " module " 
			   << module << " UnitID 0x" << std::hex << intindex 
			   << std::dec;
#endif

    delete[] copyno;
    delete[] name;
  }
  return index;
  
}

//uint32_t CastorNumberingScheme::packIndex(int section, int z, int sector,  int module ) {

uint32_t CastorNumberingScheme::packIndex(int z, int sector, int module) {
  /*
  uint32_t idx=(section&31)<<28;     //bits 28-31   (21-27 are free for now)
  idx+=((z-1)&1)<<20;                //bits  20  (1...2)
  idx+=(sector&15)<<6;               //bits  6-9 (1...16)
  idx+=(module&63);                 //bits  0-5 (1...18)
  return idx;
  */

  uint32_t idx=((z-1)&1)<<8;       //bit 8
  idx+=(sector&15)<<4;             //bits  4-7 (1...16)
  idx+=(module&15);                 //bits  0-3 (1...14)
  return idx;

}

//void CastorNumberingScheme::unpackIndex(const uint32_t& idx, int& section, int& z, int& sector, int& module) {

void CastorNumberingScheme::unpackIndex(const uint32_t& idx, int& z, int& sector, int& module) {
  /*
  section = (idx>>28)&31;
  z   = (idx>>20)&1;
  z  += 1;
  sector = (idx>>6)&15;
  module= (idx&63);
  */
  z   = (idx>>8)&1;
  z  += 1;
  sector = (idx>>4)&15;
  module= (idx&15);
}

int CastorNumberingScheme::detectorLevel(const G4Step* aStep) const {
  
  //Find number of levels
  const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
  int level = 0;
  if (touch) level = ((touch->GetHistoryDepth())+1);
  return level;
}
  
void CastorNumberingScheme::detectorLevel(const G4Step* aStep, int& level,
					  int* copyno, G4String* name) const {
 
  //Get name and copy numbers
  if (level > 0) {
    const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
    for (int ii = 0; ii < level; ii++) {
      int i      = level - ii - 1;
      name[ii]   = touch->GetVolume(i)->GetName();
      copyno[ii] = touch->GetReplicaNumber(i);
    }
  }
}
