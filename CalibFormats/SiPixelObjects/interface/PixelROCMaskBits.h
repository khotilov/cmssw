#ifndef PixelROCMaskBits_h
#define PixelROCMaskBits_h
//
//
//

#include <fstream>
#include <string>
#include "CalibFormats/SiPixelObjects/interface/PixelROCName.h"

class PixelROCMaskBits {

 public:

    PixelROCMaskBits();
    
    void setROCMaskBits(PixelROCName& rocid ,std::string bits);

    int read(const PixelROCName& rocid, std::ifstream& in);

    int readBinary(const PixelROCName& rocid, std::ifstream& in);

    unsigned int mask(unsigned int col, unsigned int row) const;

    void setMask(unsigned int col, unsigned int row, unsigned int mask);

    void writeBinary(std::ofstream& out) const;

    void writeASCII(std::ofstream& out) const;

    PixelROCName name() const {return rocid_;}

    friend std::ostream& operator<<(std::ostream& s, const PixelROCMaskBits& maskbits);

 private:

    PixelROCName rocid_;
    unsigned char bits_[520];


};


#endif
