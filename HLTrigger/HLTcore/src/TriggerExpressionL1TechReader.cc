#include <iostream>
#include <bitset>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMask.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionL1TechReader.h"

namespace triggerExpression {

// define the result of the module from the L1 reults
bool L1TechReader::operator()(const Data & data) const {
  if (not data.hasL1T())
    return false;

  typedef std::pair<std::string, unsigned int> value_type;
  if (data.ignoreL1TechPrescales()) {
    // select PSB#9 and bunch crossing 0
    const L1GtPsbWord & psb = data.l1tResults().gtPsbWord(0xbb09, 0);
    // the four 16-bit words psb.bData(1), psb.aData(1), psb.bData(0) and psb.aData(0) yield
    // (in this sequence) the 64 technical trigger bits from most significant to least significant bit
    std::bitset<64> psbTriggerWord( ((uint64_t) psb.bData(1) << 48) | 
                                    ((uint64_t) psb.aData(1) << 32) | 
                                    ((uint64_t) psb.bData(0) << 16) | 
                                    ((uint64_t) psb.aData(0)) );
    BOOST_FOREACH(const value_type & trigger, m_triggers)
      if (psbTriggerWord[trigger.second])
        return true;
  } else {
    BOOST_FOREACH(const value_type & trigger, m_triggers)
      if (data.l1tResults().technicalTriggerWord()[trigger.second])
        return true;
  }

  return false;
}

void L1TechReader::dump(std::ostream & out) const {
  if (m_triggers.size() == 0) {
    out << "FALSE";
  } else if (m_triggers.size() == 1) {
    out << m_triggers[0].first;
  } else {
    out << "(" << m_triggers[0].first;
    for (unsigned int i = 1; i < m_triggers.size(); ++i)
      out << " OR " << m_triggers[i].first;
    out << ")";
  }
}

void L1TechReader::init(const Data & data) {
  const L1GtTriggerMenu & menu = data.l1tMenu();
  const L1GtTriggerMask & mask = data.l1tTechMask();

  // clear the previous configuration
  m_triggers.clear();

  // check if the pattern has is a glob expression, or a single trigger name 
  if (not edm::is_glob(m_pattern)) {
    // no wildcard expression
    const AlgorithmMap & triggerMap = menu.gtTechnicalTriggerMap();
    AlgorithmMap::const_iterator entry = triggerMap.find(m_pattern);
    if (entry != triggerMap.end()) {
      // single L1 bit
      m_triggers.push_back( std::make_pair(m_pattern, entry->second.algoBitNumber()) );
    } else
      // trigger not found in the current menu
      if (data.shouldThrow())
        throw cms::Exception("Configuration") << "requested L1 trigger \"" << m_pattern << "\" does not exist in the current L1 menu";
      else
        edm::LogWarning("Configuration") << "requested L1 trigger \"" << m_pattern << "\" does not exist in the current L1 menu";
  } else {
    // expand wildcards in the pattern 
    bool match = false;
    boost::regex re(edm::glob2reg(m_pattern));
    const AlgorithmMap & triggerMap = menu.gtTechnicalTriggerMap();
    BOOST_FOREACH(const AlgorithmMap::value_type & entry, triggerMap)
      if (boost::regex_match(entry.first, re)) {
        match = true;
        if (data.ignoreL1Mask() or (mask.gtTriggerMask()[entry.second.algoBitNumber()] & data.daqPartitions()) != data.daqPartitions()) // unmasked in one or more partitions
          m_triggers.push_back( std::make_pair(entry.first, entry.second.algoBitNumber()) );
      }

    if (not match) {
      // m_pattern does not match any L1 bits
      if (data.shouldThrow())
        throw cms::Exception("Configuration") << "requested pattern \"" << m_pattern <<  "\" does not match any L1 trigger in the current menu";
      else
        edm::LogWarning("Configuration") << "requested pattern \"" << m_pattern <<  "\" does not match any L1 trigger in the current menu";
    }
  }

}

} // namespace triggerExpression
