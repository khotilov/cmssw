#ifndef PixelTKPCIFECConfig_h
#define PixelTKPCIFECConfig_h
//
// This class specifies the settings on the TKPCIFEC
// and the settings on the portcard
//
//
//
//
#include <vector>
#include <string>
#include <map>
#include "CalibFormats/SiPixelObjects/interface/PixelConfigBase.h"

namespace pos{
  using namespace std;

  class PixelPortCardConfig: public PixelConfigBase{

  public:
  
    PixelPortCardConfig(std::vector < std::vector< std::string> >  &tableMat);
    PixelPortCardConfig(std::string);

    void writeASCII(std::string dir="") const;
  
    unsigned int getdevicesize() const;
    std::string  getTKFECID() const;
    unsigned int getringAddress() const;
    unsigned int getccuAddress() const;
    unsigned int getchannelAddress() const;
    unsigned int geti2cSpeed() const;
    unsigned int getdeviceAddress(unsigned int i) const;
    unsigned int getdeviceValues(unsigned int i) const;
    unsigned int getdeviceAddressForSetting(std::string settingName) const;
    unsigned int getdeviceValuesForSetting(std::string settingName) const;
    void setdeviceValues(unsigned int address, unsigned int value);
    void setdeviceValues(std::string settingName, unsigned int value);
  
  private:
    void fillNameToAddress();

    std::string portcardname_;
 
    std::string  TKFECID_;//FEC ID string, as defined in tkfecconfig.dat
    unsigned int ringAddress_;//ring #
    unsigned int ccuAddress_;//CCU #

    unsigned int channelAddress_;//there are 8? channels on a CCU board
    vector < pair<unsigned int, unsigned int> > device_;//the address on the portcard, and the value of it
    unsigned int i2cSpeed_;//for the portcard, the slow i2c speed is 100kHz
  
    std::map<std::string, unsigned int> nameToAddress_; // translation from name to address, filled in by fillNameToAddress();
  };
}
#endif
