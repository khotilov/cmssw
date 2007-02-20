#ifndef L1GCTINTERNEMCAND_H
#define L1GCTINTERNEMCAND_H

#include <boost/cstdint.hpp>
#include <ostream>
#include <string>

#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCand.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"

/*! \class L1GctEmCand
 * \brief Level-1 Trigger EM candidate within GCT processing
 *
 */

/*! \author Jim Brooke
 *  \date June 2006
 */


class L1GctInternEmCand : public L1GctCand {
public:

  /// default constructor (for vector initialisation etc.)
  L1GctInternEmCand();

  /// construct from raw data
  L1GctInternEmCand(uint16_t data, bool iso, unsigned block, unsigned index);

  /// destructor (virtual to prevent compiler warnings)
  virtual ~L1GctInternEmCand();
  
  /// region associated with the candidate
  L1CaloRegionDetId regionId() const;

  /// name of object
  std::string name() const;

  /// was an object really found?
  bool empty() const;
  
  /// get the raw data
  uint16_t raw() const { return m_data; }
  
  /// get rank bits
  unsigned rank() const { return m_data & 0x3f; }

  /// get eta index -  Barrel 1:4, Endcap  5:7, HF  8:11
  unsigned etaIndex() const { return (m_data>>6) & 0xf; }

  /// get eta sign (1 for -ve Z, 0 for +ve Z)
  unsigned etaSign() const { return (m_data>>9) & 0x1; }

  /// get phi index (0-17)
  unsigned phiIndex() const { return (m_data>>10) & 0x1f; }

  /// which stream did this come from
  bool isolated() const { return m_iso; }

  /// which capture block did this come from
  unsigned capBlock() const { return (m_source>>9) & 0x7f; }

  /// what index within capture block
  unsigned capIndex() const { return m_source&0x1ff; }

  /// equality operator
  int operator==(const L1GctInternEmCand& c) const { return ((m_data==c.raw() && m_iso==c.isolated())
                                                      || (this->empty() && c.empty())); }

  /// inequality operator
  int operator!=(const L1GctInternEmCand& c) const { return ((m_data!=c.raw() || m_iso!=c.isolated())
                                                     && (!this->empty() || !c.empty())); }

 private:

  // set internal data from rank and region ieta, iphi
  void construct(unsigned rank, unsigned eta, unsigned etaSgn, unsigned phi);

 private:

  uint16_t m_data;
  uint16_t m_source;
  bool m_iso;

 };


std::ostream& operator<<(std::ostream& s, const L1GctInternEmCand& cand);



#endif 
