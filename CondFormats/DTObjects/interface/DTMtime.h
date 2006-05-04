#ifndef DTMtime_H
#define DTMtime_H
/** \class DTMtime
 *
 *  Description:
 *       Class to hold drift tubes mean-times
 *             ( SL by SL mean-time calculation )
 *
 *  $Date: 2006/02/28 18:06:29 $
 *  $Revision: 1.3 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "CondFormats/DTObjects/interface/DTTimeUnits.h"
#include "DataFormats/MuonDetId/interface/DTSuperLayerId.h"

//---------------
// C++ Headers --
//---------------
#include <string>
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class DTSLMtimeData {

 public:

  DTSLMtimeData();
  ~DTSLMtimeData();

  int   wheelId;
  int stationId;
  int  sectorId;
  int      slId;
  float mTime;
  float mTrms;

};


class DTMtime {

 public:

  /** Constructor
   */
  DTMtime();
  DTMtime( const std::string& version );

  /** Destructor
   */
  ~DTMtime();

  /** Operations
   */
  /// get content
  int slMtime( int   wheelId,
               int stationId,
               int  sectorId,
               int      slId,
               float&  mTime,
               float&  mTrms,
               DTTimeUnits::type unit = DTTimeUnits::counts ) const;
  int slMtime( const DTSuperLayerId& id,
               float&  mTime,
               float&  mTrms,
               DTTimeUnits::type unit = DTTimeUnits::counts ) const;
  float unit() const;

  /// access version
  const
  std::string& version() const;
  std::string& version();

  /// reset content
  void clear();

  int setSLMtime( int   wheelId,
                  int stationId,
                  int  sectorId,
                  int      slId,
                  float   mTime,
                  float   mTrms,
                  DTTimeUnits::type unit = DTTimeUnits::counts );
  int setSLMtime( const DTSuperLayerId& id,
                  float   mTime,
                  float   mTrms,
                  DTTimeUnits::type unit = DTTimeUnits::counts );
  void setUnit( float unit );

  /// Access methods to data
  typedef std::vector<DTSLMtimeData>::const_iterator const_iterator;
  const_iterator begin() const;
  const_iterator end() const;

 private:

  /// read and store full content
  void initSetup() const;

  std::string dataVersion;
  float nsPerCount;

  std::vector<DTSLMtimeData> slData;

//  static int rmsFactor;

};


#endif // DTMtime_H

