#ifndef PixelPortcardMap_h
#define PixelPortcardMap_h
//
// This class provides the maping between
// portcards and the modules controlled by
// the card
//
//
//
 
#include <string>
#include <vector>
#include <map>
#include <set>
#include "CalibFormats/SiPixelObjects/interface/PixelConfigBase.h"
#include "CalibFormats/SiPixelObjects/interface/PixelModuleName.h"

namespace pos{
  class PixelPortcardMap: public PixelConfigBase
  {
  public:

    PixelPortcardMap(std::string filename);

    PixelPortcardMap(std::vector< std::vector < std::string> > &tableMat);

    // Get the port card and AOH associated with this module.  If the module has one(two) channels, this vector contains one(two) element(s).
    const std::set< std::pair< std::string, int > > PortCardAndAOHs(const PixelModuleName& aModule) const;
    //                            portcardname, aoh #

    const std::set< std::string > portcards(const PixelModuleName& aModule) const;

    int numChannels(const PixelModuleName& aModule) {return PortCardAndAOHs(aModule).size();}

    const std::pair< std::string, int > PortCardAndAOH(const PixelModuleName& aModule, const std::string& TBMChannel) const;

    // set of all modules attached to a port card
    std::set< PixelModuleName > modules(std::string portCardName) const;

    // all port cards in the map
    std::set< std::string > portcards();

  private:
    //    (modulename, TBM channel A or B) <--> (portcardname, aoh #)
    std::map< std::pair<PixelModuleName, std::string>, std::pair<std::string, int> > map_;
  };
}
#endif
