#ifndef PixelTKFECConfig_h
#define PixelTKFECConfig_h
//
// This class specifies which TKFEC boards
// are used and how they are addressed
// 
// 
// 
//
//
//
#include <iostream>
#include <vector>
#include <string>
#include "CalibFormats/SiPixelObjects/interface/PixelConfigBase.h"
#include "CalibFormats/SiPixelObjects/interface/PixelTKFECParameters.h"

class PixelTKFECConfig: public PixelConfigBase {

 public:

 PixelTKFECConfig(std::string filename);  //  <---- Modified for the conversion from parallel vectors to object that contain the configuration
   
 PixelTKFECConfig(std::vector<std::vector<std::string> >& tableMat ); 

    unsigned int getNTKFECBoards() const;

    std::string  getTKFECID(unsigned int i) const;
    unsigned int getCrate(unsigned int i) const;
    std::string  getType(unsigned int i) const;
    unsigned int getAddress(unsigned int i) const;
    unsigned int crateFromTKFECID(std::string TKFECID) const;
    std::string  typeFromTKFECID(std::string TKFECID) const;
    unsigned int addressFromTKFECID(std::string TKFECID) const;
    
 private:
    std::vector< PixelTKFECParameters > TKFECconfig_;
};

#endif
