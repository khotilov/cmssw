#ifndef Fireworks_Calo_FWHFTowerSliceSelector_h
#define Fireworks_Calo_FWHFTowerSliceSelector_h
// -*- C++ -*-
//
// Package:     Calo
// Class  :     FWHFTowerSliceSelector
// 
/**\class FWHFTowerSliceSelector FWHFTowerSliceSelector.h Fireworks/Calo/interface/FWHFTowerSliceSelector.h

 Description: [one line class summary]

 Usage:
    <usage>

*/
//
// Original Author:  Alja Mrak-Tadel
//         Created:  Wed Jun  2 19:21:13 CEST 2010
// $Id: FWHFTowerSliceSelector.h,v 1.1 2010/06/02 17:34:03 amraktad Exp $
//

// system include files

// user include files
class HcalDetId;

#include "Fireworks/Calo/src/FWFromSliceSelector.h"

// forward declarations

class FWHFTowerSliceSelector : public FWFromSliceSelector
{
public:
   FWHFTowerSliceSelector(TH2F* h, const FWEventItem* i) : FWFromSliceSelector(h, i) {}
   virtual ~FWHFTowerSliceSelector() {}

   virtual void doSelect(const TEveCaloData::CellId_t&);
   virtual void doUnselect(const TEveCaloData::CellId_t&);
   
private:
   bool findBinFromId(HcalDetId& id, int tower) const;
};


#endif
