#ifndef DataFormats_SiStripDigi_SiStripRawDigi_H
#define DataFormats_SiStripDigi_SiStripRawDigi_H

#include "boost/cstdint.hpp"

/** 
    @brief A Digi for the silicon strip detector, containing only adc
    information, and suitable for storing raw hit information.
*/
class SiStripRawDigi {

 public:

/*   SiStripRawDigi() : adc_(0) {;} */
/*   SiStripRawDigi( const uint16_t& adc ) : adc_(adc) {;} */
/*   ~SiStripRawDigi() {;} */

/*   inline const uint16_t& adc() const { return adc_; } */
  
/*   inline bool operator< ( const SiStripRawDigi& other ) const { return false; }  */
  
  
  SiStripRawDigi() : strip_(0), adc_(0) {;} //@@ temp
  SiStripRawDigi( uint16_t strip, uint16_t adc ) : strip_(strip), adc_(adc) {;} //@@ temp
  ~SiStripRawDigi() {;}
  
  uint16_t strip()   const { return strip_; } //@@ temp
  uint16_t channel() const { return strip(); } //@@ temp
  uint16_t adc()     const { return adc_; } //@@ temp
  
  inline bool operator< ( const SiStripRawDigi& other ) const { return ( this->strip() < other.strip() ); } //@@ temp
  
 private:
  
  uint16_t strip_; //@@ temp
  uint16_t adc_;
  
};

#endif // DataFormats_SiStripDigi_SiStripRawDigi_H



