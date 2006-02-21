#ifndef CSCALCTTrailer_h
#define CSCALCTTrailer_h

/** documented in  flags
  http://www.phys.ufl.edu/~madorsky/alctv/alct2000_spec.PDF
*/

#include <string.h> // memcpy

class CSCALCTTrailer
{
public:
  CSCALCTTrailer() { bzero(this, sizeInWords()*2);  e0dLine = 0xe0d; reserved_4=0xd;}
  explicit CSCALCTTrailer(const unsigned short * buf) {
    memcpy(this, buf, sizeInWords()*2);
  }
  unsigned short * data() {return (unsigned short *) this;}
  /// in 16-bit frames
  static int sizeInWords() {return 4;}

  int getCRC() { 
    //printf("crc1 %x crc0 %x \n",crc1,crc0);
    return ((crc1&0x7ff)<<11) | (crc0&0x7ff) ; 
  }

  bool check() const {return (e0dLine & 0xfff) == 0xe0d;}
  int wordCount() { return frameCount; }
  unsigned alctCRCCheck() const { return reserved_3; }

private:

  unsigned crc0:11, reserved_0:5;
  unsigned crc1:11, reserved_1:5;
  unsigned e0dLine:12, reserved_2:4;
  unsigned frameCount:10, reserved_3:2, reserved_4:4;
};

#endif

