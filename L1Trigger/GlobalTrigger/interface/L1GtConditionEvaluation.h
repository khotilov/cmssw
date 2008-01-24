#ifndef GlobalTrigger_L1GtConditionEvaluation_h
#define GlobalTrigger_L1GtConditionEvaluation_h

/**
 * \class L1GtConditionEvaluation
 * 
 * 
 * Description: Base class for evaluation of the L1 Global Trigger object templates.
 * 
 * Implementation:
 *    <TODO: enter implementation details>
 *   
 * \author: Vasile Mihai Ghete   - HEPHY Vienna 
 * 
 * $Date$
 * $Revision$
 *
 */

// system include files
#include <iostream>

#include <string>
#include <vector>

#include <boost/cstdint.hpp>

// user include files

//   base class

//
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// forward declarations

// class interface
class L1GtConditionEvaluation
{

public:

    /// constructor
    L1GtConditionEvaluation();

    /// copy constructor
    L1GtConditionEvaluation(L1GtConditionEvaluation&);

    /// destructor
    virtual ~L1GtConditionEvaluation();

public:

    /// get / set the maximum number of objects received for the
    /// evaluation of the condition
    inline int condMaxNumberObjects() const {
        return m_condMaxNumberObjects;
    }

    inline void setCondMaxNumberObjects(int condMaxNumberObjectsValue) {
        m_condMaxNumberObjects = condMaxNumberObjectsValue;
    }

    /// get the latest result for the condition
    inline bool condLastResult() const {
        return m_condLastResult;
    }

    /// call evaluateCondition and save last result
    inline void evaluateConditionStoreResult() {
        m_condLastResult = evaluateCondition();
    }

    /// the core function to check if the condition matches
    virtual const bool evaluateCondition() const = 0;

    /// get numeric expression
    virtual std::string getNumericExpression() {
        if (m_condLastResult) {
            return "1";
        }
        else {
            return "0";
        }
    }

    /// get all the object combinations evaluated to true in the condition
    inline CombinationsInCond* getCombinationsInCond() const {
        return m_combinationsInCond;
    }

    /// print condition
    virtual void print(std::ostream& myCout) const;

protected:

    /// check if a value is greater than a threshold or
    /// greater-or-equal depending on the value of the condGEqValue flag
    template<class Type1, class Type2> const bool checkThreshold(const Type1& threshold,
        const Type2& value, const bool condGEqValue) const;

    ///check if a bit with a given number is set in a mask
    template<class Type1> const bool checkBit(const Type1& mask, const unsigned int bitNumber) const;

protected:

    /// maximum number of objects received for the evaluation of the condition
    /// usually retrieved from event setup
    int m_condMaxNumberObjects;

    /// the last result of evaluateCondition()
    bool m_condLastResult;

    /// store all the object combinations evaluated to true in the condition
    CombinationsInCond* m_combinationsInCond;

};

// define templated methods

// check if a value is greater than a threshold or
// greater-or-equal depending on the value of the condGEqValue flag
template<class Type1, class Type2> const bool L1GtConditionEvaluation::checkThreshold(
    const Type1& threshold, const Type2& value, const bool condGEqValue) const {

    //if (value > 0) {
    //    LogTrace("L1GtConditionEvaluation") << "  threshold check for condGEqValue = "
    //        << condGEqValue << "\n    hex: " << std::hex << "threshold = " << threshold
    //        << " value = " << value << "\n    dec: " << std::dec << "threshold = " << threshold
    //        << " value = " << value << std::endl;
    //}

    if (condGEqValue) {
        if (value >= threshold) {

            //LogTrace("L1GtConditionEvaluation") << "    condGEqValue: value >= threshold"
            //    << std::endl;

            return true;
        }

        return false;

    }
    else {

        if (value == threshold) {

            //LogTrace("L1GtConditionEvaluation") << "    condGEqValue: value = threshold"
            //    << std::endl;

            return true;
        }

        return false;
    }
}

// check if a bit with a given number is set in a mask
template<class Type1> const bool L1GtConditionEvaluation::checkBit(const Type1& mask,
    const unsigned int bitNumber) const {

    boost::uint64_t oneBit = 1ULL;

    if (bitNumber >= (sizeof(oneBit)*8)) {

        LogTrace("L1GtConditionEvaluation") << "    checkBit " << "\n     Bit number = "
            << bitNumber << " larger than maximum allowed " << sizeof(oneBit)*8 << std::endl;

        return false;
    }

    oneBit <<= bitNumber;

    //LogTrace("L1GtConditionEvaluation") << "    checkBit " << "\n     mask address = " << &mask
    //    << std::dec << "\n     dec: " << "mask = " << mask << " oneBit = " << oneBit
    //    << " bitNumber = " << bitNumber << std::hex << "\n     hex: " << "mask = " << mask
    //    << " oneBit = " << oneBit << " bitNumber = " << bitNumber << std::dec
    //    << "\n     mask & oneBit result = " << bool ( mask & oneBit ) << std::endl;

    return (mask & oneBit);
}

#endif
