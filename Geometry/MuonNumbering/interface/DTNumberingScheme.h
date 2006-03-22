#ifndef MuonNumbering_DTNumberingScheme_h
#define MuonNumbering_DTNumberingScheme_h

/** \class DTNumberingScheme
 *
 * implementation of MuonNumberingScheme for muon barrel,
 * converts the MuonBaseNumber to a unit id
 *  
 *  $Date: 2006/02/15 13:21:24 $
 *  $Revision: 1.1 $
 * \author Arno Straessner, CERN <arno.straessner@cern.ch>
 *
 */

#include "Geometry/MuonNumbering/interface/MuonNumberingScheme.h"

class MuonBaseNumber;
class MuBarDetBuilder;

class DTNumberingScheme : public MuonNumberingScheme {
 public:

  DTNumberingScheme();
  ~DTNumberingScheme(){}
  
  virtual int baseNumberToUnitNumber(const MuonBaseNumber num);
  
 private:
  // Decode MuonBaseNumber to id: no checking
  void decode(const MuonBaseNumber& num,
              int& wire_id,
              int& layer_id,
              int& superlayer_id,
              int& sector_id,
              int& station_id,
              int& wheel_id
             ) const;

  int theRegionLevel;
  int theWheelLevel;
  int theStationLevel;
  int theSuperLayerLevel;
  int theLayerLevel;
  int theWireLevel;

  /** Same as BaseNumberToUnitNumber but w/o check: used by MuBarDetBuilder
   * class to build the geometry from DDD */
  int getDetId(const MuonBaseNumber num) const;
  friend class DTGeometryBuilderFromDDD;
  
};

#endif
