#ifndef L1GCTJETCOUNTS_H
#define L1GCTJETCOUNTS_H

#include <boost/cstdint.hpp>
#include <vector>
#include <ostream>

#include "DataFormats/L1GlobalTrigger/interface/L1TriggerObject.h"

///
/// \class L1GctJetCounts
/// 
/// \author: Jim Brooke
///
/// Class to store the GCT jet count output
/// 


class L1GctJetCounts : public L1TriggerObject {

 public:

  /// default constructor
  L1GctJetCounts();

  /// constructor for unpacking
  L1GctJetCounts(uint32_t data0, uint32_t data1);

  /// constructor for emulator
  L1GctJetCounts(std::vector<unsigned> counts);

  /// destructor
  ~L1GctJetCounts();

  /// name method
  std::string name() const { return "JetCounts"; }

  /// empty method
  bool empty() const { return false; }

  /// get raw word 0
  uint32_t raw0() const { return m_data0; }

  /// get raw word 1
  uint32_t raw1() const { return m_data1; }

  /// get count by index
  unsigned count(unsigned i) const;

 private:

  /// static maximum number of jet counts
  static const unsigned MAX_COUNTS=12;

  uint32_t m_data0;
  uint32_t m_data1;

};

std::ostream& operator<<(std::ostream& s, const L1GctJetCounts& c);

#endif
