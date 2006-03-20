#ifndef MuonNumbering_RPCNumberingScheme_h
#define MuonNumbering_RPCNumberingScheme_h

/** \class RPCNumberingScheme
 *
 * implementation of MuonNumberingScheme for muon rpc,
 * converts the MuonBaseNumber to a unit id
 *  
 *  $Date: 2006/02/15 13:21:24 $
 *  $Revision: 1.1 $
 * \author Arno Straessner, CERN <arno.straessner@cern.ch>
 *
 */

#include "Geometry/MuonNumbering/interface/MuonNumberingScheme.h"

class MuonBaseNumber;

class RPCNumberingScheme : public MuonNumberingScheme {
 public:

  RPCNumberingScheme();
  virtual ~RPCNumberingScheme(){};
  
  virtual int baseNumberToUnitNumber(const MuonBaseNumber&) const;
  
 private:
  int theRegionLevel;
  int theBWheelLevel;
  int theBStationLevel;
  int theBPlaneLevel;
  int theBChamberLevel;
  int theEPlaneLevel;
  int theESectorLevel;
  int theERollLevel;

};

#endif
