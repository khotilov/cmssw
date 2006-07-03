#ifndef L1GCTJETCAND_H
#define L1GCTJETCAND_H

#include <boost/cstdint.hpp>
#include <ostream>
#include <string>

#include "DataFormats/L1GlobalTrigger/interface/L1TriggerObject.h"

/*! \class L1GctJetCand
 * \brief Level-1 Trigger jet candidate
 *
 */

/*! \author Jim Brooke
 *  \date June 2006
 */


class L1GctJetCand : public L1TriggerObject {
public:
  /// default constructor (for vector initialisation etc.)
  L1GctJetCand();

  /// construct from raw data
  L1GctJetCand(uint16_t data, bool isTau, bool isFor);

  /// construct from rank, eta, phi, isolation
  /// NB - eta = -6 to -0, +0 to +6. Sign is bit 3, 1 means -ve Z, 0 means +ve Z
  L1GctJetCand(unsigned rank, unsigned phi, unsigned eta, bool isTau, bool isFor);

  /// destructor
  ~L1GctJetCand();

  /// name of object - inherited from L1TriggerObject
  std::string name() const;

  /// was an object really found? - inherited from L1TriggerObject
  bool empty() const;

  /// get the raw data
  uint16_t raw() const { return m_data; }
  
  /// get rank bits
  unsigned rank() const { return m_data & 0x3f; }

  /// get eta index (bit 3 is sign, 1 for -ve Z, 0 for +ve Z)
  unsigned etaIndex() const { return (m_data>>6) & 0xf; }

  /// get eta sign bit (1 for -ve Z, 0 for +ve Z)
  unsigned etaSign() const { return (m_data>>9) & 0x1; }
  
  /// get phi index (0-17)
  unsigned phiIndex() const { return (m_data>>10) & 0x1f; }

  /// check if this is a central jet
  bool isCentral() const { return (!m_isTau) && (!m_isFor); }

  /// check if this is a tau
  bool isTau() const { return m_isTau; }

  /// check if this is a forward jet
  bool isForward() const { return m_isFor; }

 private:

  uint16_t m_data;
  bool m_isTau;
  bool m_isFor;

 };

std::ostream& operator<<(std::ostream& s, const L1GctJetCand& cand);

#endif
