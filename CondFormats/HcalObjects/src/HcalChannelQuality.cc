/** 
\class HcalChannelQuality
\author Fedor Ratnikov (UMd)
POOL object to store pedestal values 4xCapId
$Author: ratnikov
$Date: 2005/11/02 21:49:19 $
$Revision: 1.2 $
*/

#include <iostream>

#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"

HcalChannelQuality::HcalChannelQuality() 
  : mSorted (false) {}

HcalChannelQuality::~HcalChannelQuality(){}

HcalChannelQuality::Quality HcalChannelQuality::quality (unsigned long fId, bool fWarning) const {
  Item target;
  target.mId = fId;
  std::vector<Item>::const_iterator cell;
  if (sorted ()) {
    cell = std::lower_bound (mItems.begin(), mItems.end(), target);
  }
  else {
    if (fWarning) std::cerr << "HcalChannelQuality::quality-> container is not sorted. Sort it to search effectively" << std::endl;
    cell = std::find (mItems.begin(), mItems.end(), target);
  }
  if (cell == mItems.end() || cell->mId != fId) return UNKNOWN;
  return (HcalChannelQuality::Quality) cell->mQuality;
}

std::vector<unsigned long> HcalChannelQuality::getAllChannels () const {
  std::vector<unsigned long> result;
  for (std::vector<Item>::const_iterator item = mItems.begin (); item != mItems.end (); item++) {
    result.push_back (item->mId);
  }
  return result;
}

bool HcalChannelQuality::setChannel (unsigned long fId, Quality fQuality) {
  Item item;
  item.mId = fId;
  item.mQuality = unsigned (fQuality);
  mSorted = false;
  return true;
}

void HcalChannelQuality::sort () {
  if (!mSorted) {
    std::sort (mItems.begin(), mItems.end());
    mSorted = true;
  }
}

const char* HcalChannelQuality::str (Quality fQuality) {
  switch (fQuality) {
  case BAD: return "BAD";
  case GOOD: return "GOOD";
  case HOT: return "HOT";
  case DEAD: return "DEAD";
  case END: return "END";
  default: return "UNKNOWN";
  }
}
