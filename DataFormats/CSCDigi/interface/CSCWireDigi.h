#ifndef CSCWireDigi_CSCWireDigi_h
#define CSCWireDigi_CSCWireDigi_h

/**\class CSCWireDigi
 *
 * Digi for CSC anode wires. 
 *
 */

#include <vector>
#include <iosfwd>

class CSCWireDigi{

public:

  /// Constructors
  
  CSCWireDigi (int wire, unsigned int tbinb);  /// wiregroup#, tbin bit word
  CSCWireDigi ();                     /// default

  /// return wiregroup number
  int getWireGroup() const {return wire_;}
  /// return the word with time bins bits
  unsigned int getTimeBinWord() const {return tbinb_;}
  /// return tbin number, (obsolete, use getTimeBin() instead)
  int getBeamCrossingTag() const;
  /// return first tbin ON number
  int getTimeBin()         const;
  /// return vector of time bins ON
  std::vector<int> getTimeBinsOn() const;

  /// Print content of digi
  void print() const;

  /// set wiregroup number
  void setWireGroup(unsigned int wiregroup) {wire_= wiregroup;}


private:

  uint16_t wire_;
  uint32_t tbinb_;

};

std::ostream & operator<<(std::ostream & o, const CSCWireDigi& digi);

#endif
