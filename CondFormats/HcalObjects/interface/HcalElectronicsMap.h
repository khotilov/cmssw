#ifndef HcalElectronicsMap_h
#define HcalElectronicsMap_h

/** 
\class HcalElectronicsMap
\author Fedor Ratnikov (UMd)
POOL object to store map between detector ID, electronics ID and trigger ID
$Author: ratnikov
$Date: 2006/02/15 19:48:02 $
$Revision: 1.7 $
*/

#include <vector>
#include <algorithm>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalCalibDetId.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "CalibFormats/HcalObjects/interface/HcalText2DetIdConverter.h"

// 
class HcalElectronicsMap {
 public:
  HcalElectronicsMap();
  ~HcalElectronicsMap();

  /// lookup the logical detid associated with the given electronics id
  //return Null item if no such mapping
  const DetId lookup(HcalElectronicsId fId, bool fWarning = true) const;

  /// brief lookup the electronics detid associated with the given logical id
  //return Null item if no such mapping
  const HcalElectronicsId lookup(DetId fId, bool fWarning = true) const;

  /// brief lookup the trigger logical detid associated with the given electronics id
  //return Null item if no such mapping
  const DetId lookupTrigger(HcalElectronicsId fId, bool fWarning = true) const;

  /// brief lookup the electronics detid associated with the given trigger logical id
  //return Null item if no such mapping
  const HcalElectronicsId lookupTrigger(DetId fId, bool fWarning = true) const;


  std::vector <HcalElectronicsId> allElectronicsId () const;
  std::vector <DetId> allDetectorId () const;
  std::vector <HcalDetId> allHcalDetectorId () const;
  std::vector <HcalCalibDetId> allCalibrationId () const;
  std::vector <HcalTrigTowerDetId> allTriggerId () const;

  // map channels
  bool setMapping (DetId fId, HcalElectronicsId fElectronicsId, HcalTrigTowerDetId fTriggerId);
  bool mapEId2tId (HcalElectronicsId fElectronicsId, HcalTrigTowerDetId fTriggerId);
  bool mapEId2chId (HcalElectronicsId fElectronicsId, DetId fId);
  // sorting
  void sortById ();
  void sortByElectronicsId ();
  void sortByTriggerId ();
  void sort () {sortByElectronicsId ();}
  class Item { 
  public:
    Item () {mId = mElId = mTrigId = 0;}
    Item (unsigned long fId, unsigned long fElId, unsigned long fTrigId) 
      : mId (fId), mElId (fElId), mTrigId (fTrigId) {}
    class LessById {public: bool operator () (const Item& a, const Item& b) {return a.mId < b.mId;}};
    class LessByElId {public: bool operator () (const Item& a, const Item& b) {return a.mElId < b.mElId;}};
    class LessByTrigId {public: bool operator () (const Item& a, const Item& b) {return a.mTrigId < b.mTrigId;}};
    unsigned long mId;
    unsigned long mElId;
    unsigned long mTrigId;
  };
 protected:
  const Item* findById (unsigned long fId, bool fWarning) const;
  const Item* findByElId (unsigned long fElId, bool fWarning) const;
  const Item* findByTrigId (unsigned long fTrigId, bool fWarning) const;
  
  std::vector<Item> mItems;
  bool mSortedById;
  bool mSortedByElId;
  bool mSortedByTrigId;
};

#endif
