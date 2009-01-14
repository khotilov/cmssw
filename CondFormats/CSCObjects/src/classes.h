#include "CondFormats/CSCObjects/interface/CSCPedestals.h"
#include "CondFormats/CSCObjects/interface/CSCDBPedestals.h"
#include "CondFormats/CSCObjects/interface/CSCGains.h"
#include "CondFormats/CSCObjects/interface/CSCDBGains.h"
#include "CondFormats/CSCObjects/interface/CSCNoiseMatrix.h"
#include "CondFormats/CSCObjects/interface/CSCDBNoiseMatrix.h"
#include "CondFormats/CSCObjects/interface/CSCcrosstalk.h"
#include "CondFormats/CSCObjects/interface/CSCDBCrosstalk.h"
#include "CondFormats/CSCObjects/interface/CSCBadStrips.h"
#include "CondFormats/CSCObjects/interface/CSCBadWires.h"
#include "CondFormats/CSCObjects/interface/CSCChamberMap.h"
#include "CondFormats/CSCObjects/interface/CSCCrateMap.h"
#include "CondFormats/CSCObjects/interface/CSCDDUMap.h"
#include "CondFormats/CSCObjects/interface/CSCChamberIndex.h"

namespace {
  struct dictionary {

    std::vector<CSCPedestals::Item>   pedcontainer1;
    std::map< int, std::vector<CSCPedestals::Item> > pedmap;
    std::vector<CSCDBPedestals::Item> pedcontainer2;

    std::vector<CSCGains::Item>   gcontainer1;
    std::map< int, std::vector<CSCGains::Item> > gmap;
    std::vector<CSCDBGains::Item> gcontainer2;

    std::vector<CSCNoiseMatrix::Item> mcontainer1;
    std::map< int, std::vector< CSCNoiseMatrix::Item> > mmap;
    std::vector<CSCDBNoiseMatrix::Item> mcontainer2;

    std::vector<CSCcrosstalk::Item> ccontainer1;
    std::map< int, std::vector< CSCcrosstalk::Item> > cmap;
    std::vector<CSCDBCrosstalk::Item> ccontainer2;

    std::vector<CSCBadStrips::BadChamber> bschmcontainer;
    std::vector<CSCBadStrips::BadChannel> bschncontainer;

    std::vector<CSCBadWires::BadChamber> bwchmcontainer;
    std::vector<CSCBadWires::BadChannel> bwchncontainer;

  };
}

#include "CondFormats/CSCObjects/interface/CSCChamberMap.h"
namespace{
  std::map< int, CSCMapItem::MapItem > chmap;
}
#include "CondFormats/CSCObjects/interface/CSCCrateMap.h"
namespace{
  std::map< int, CSCMapItem::MapItem > cratemap;
}
#include "CondFormats/CSCObjects/interface/CSCDDUMap.h"
namespace{
  std::map< int, CSCMapItem::MapItem > ddumap;
}
#include "CondFormats/CSCObjects/interface/CSCChamberIndex.h"
namespace{
  std::map< int, CSCMapItem::MapItem > indexmap;
}

#include "CondFormats/CSCObjects/interface/CSCBadChambers.h"
#include "CondFormats/CSCObjects/interface/CSCIdentifier.h"
#include "CondFormats/CSCObjects/interface/CSCReadoutMapping.h"
#include "CondFormats/CSCObjects/interface/CSCTriggerMapping.h"
#include "CondFormats/CSCObjects/interface/CSCL1TPParameters.h"
