//-------------------------------------------------
//
/** \class L1MuGMTLFDisableHotLUT
 *
 *   LFDisableHot look-up table
 *          
 *   this class was automatically generated by 
 *     L1MuGMTLUT::MakeSubClass()  
*/ 
//   $Date: 2007/03/23 18:51:35 $
//   $Revision: 1.2 $
//
//   Author :
//   H. Sakulin            HEPHY Vienna
//
//   Migrated to CMSSW:
//   I. Mikulec
//
//--------------------------------------------------
#ifndef L1TriggerGlobalMuonTrigger_L1MuGMTLFDisableHotLUT_h
#define L1TriggerGlobalMuonTrigger_L1MuGMTLFDisableHotLUT_h

//---------------
// C++ Headers --
//---------------


//----------------------
// Base Class Headers --
//----------------------
#include "L1Trigger/GlobalMuonTrigger/src/L1MuGMTLUT.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------


//              ---------------------
//              -- Class Interface --
//              ---------------------


class L1MuGMTLFDisableHotLUT : public L1MuGMTLUT {
  
 public:
  enum {DT, CSC};

  /// constuctor using function-lookup
  L1MuGMTLFDisableHotLUT() : L1MuGMTLUT("LFDisableHot", 
				       "DT CSC",
				       "eta(6) phi(8)",
				       "disable_hot(1)", 10, false) {
    InitParameters();
  } ;

  /// destructor
  virtual ~L1MuGMTLFDisableHotLUT() {};

  /// specific lookup function for disable_hot
  unsigned SpecificLookup_disable_hot (int idx, unsigned eta, unsigned phi) const {
    std::vector<unsigned> addr(2);
    addr[0] = eta;
    addr[1] = phi;
    return Lookup(idx, addr) [0];
  };

  /// specific lookup function for entire output field
  unsigned SpecificLookup (int idx, unsigned eta, unsigned phi) const {
    std::vector<unsigned> addr(2);
    addr[0] = eta;
    addr[1] = phi;
    return LookupPacked(idx, addr);
  };



  /// access to lookup function with packed input and output

  virtual unsigned LookupFunctionPacked (int idx, unsigned address) const {
    std::vector<unsigned> addr = u2vec(address, m_Inputs);
    return TheLookupFunction(idx ,addr[0] ,addr[1]);

  };

 private:
  /// Initialize scales, configuration parameters, alignment constants, ...
  void InitParameters();

  /// The lookup function - here the functionality of the LUT is implemented
  unsigned TheLookupFunction (int idx, unsigned eta, unsigned phi) const;

};
#endif



















