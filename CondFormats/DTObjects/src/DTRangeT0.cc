/*
 *  See header file for a description of this class.
 *
 *  $Date: 2007/11/24 12:29:12 $
 *  $Revision: 1.7.6.3 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// This Class' Header --
//----------------------
#include "CondFormats/DTObjects/interface/DTRangeT0.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CondFormats/DTObjects/interface/DTDataBuffer.h"

//---------------
// C++ Headers --
//---------------
#include <iostream>

//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------
DTRangeT0::DTRangeT0():
  dataVersion( " " ) {
}


DTRangeT0::DTRangeT0( const std::string& version ):
  dataVersion( version ) {
}


DTRangeT0Id::DTRangeT0Id() :
    wheelId( 0 ),
  stationId( 0 ),
   sectorId( 0 ),
       slId( 0 ) {
}


DTRangeT0Data::DTRangeT0Data() :
  t0min( 0 ),
  t0max( 0 ) {
}


//--------------
// Destructor --
//--------------
DTRangeT0::~DTRangeT0() {
  DTDataBuffer<int,int>::dropBuffer( mapName() );
}


DTRangeT0Id::~DTRangeT0Id() {
}


DTRangeT0Data::~DTRangeT0Data() {
}


//--------------
// Operations --
//--------------
int DTRangeT0::get( int   wheelId,
                     int stationId,
                    int  sectorId,
                    int      slId,
                    int&    t0min,
                    int&    t0max ) const {

  t0min =
  t0max = 0;

  std::string mName = mapName();
  DTBufferTree<int,int>* dBuf =
  DTDataBuffer<int,int>::findBuffer( mName );
  if ( dBuf == 0 ) {
    cacheMap();
    dBuf =
    DTDataBuffer<int,int>::findBuffer( mName );
  }

  std::vector<int> chanKey;
  chanKey.push_back(   wheelId );
  chanKey.push_back( stationId );
  chanKey.push_back(  sectorId );
  chanKey.push_back(      slId );
  int ientry;
  int searchStatus = dBuf->find( chanKey.begin(), chanKey.end(), ientry );
  if ( !searchStatus ) {
    const DTRangeT0Data& data( dataList[ientry].second );
    t0min = data.t0min;
    t0max = data.t0max;
  }

  return searchStatus;

}


int DTRangeT0::get( const DTSuperLayerId& id,
                    int&    t0min,
                    int&    t0max ) const {
  return get( id.wheel(),
              id.station(),
              id.sector(),
              id.superLayer(),
              t0min, t0max );
}


const
std::string& DTRangeT0::version() const {
  return dataVersion;
}


std::string& DTRangeT0::version() {
  return dataVersion;
}


void DTRangeT0::clear() {
  DTDataBuffer<int,int>::dropBuffer( mapName() );
  dataList.clear();
  return;
}


int DTRangeT0::set( int   wheelId,
                    int stationId,
                    int  sectorId,
                    int      slId,
                    int     t0min,
                    int     t0max ) {

  std::string mName = mapName();
  DTBufferTree<int,int>* dBuf =
  DTDataBuffer<int,int>::findBuffer( mName );
  if ( dBuf == 0 ) {
    cacheMap();
    dBuf =
    DTDataBuffer<int,int>::findBuffer( mName );
  }
  std::vector<int> chanKey;
  chanKey.push_back(   wheelId );
  chanKey.push_back( stationId );
  chanKey.push_back(  sectorId );
  chanKey.push_back(      slId );
  int ientry;
  int searchStatus = dBuf->find( chanKey.begin(), chanKey.end(), ientry );

  if ( !searchStatus ) {
    DTRangeT0Data& data( dataList[ientry].second );
    data.t0min = t0min;
    data.t0max = t0max;
    return -1;
  }
  else {
    DTRangeT0Id key;
    key.  wheelId =   wheelId;
    key.stationId = stationId;
    key. sectorId =  sectorId;
    key.     slId =      slId;
    DTRangeT0Data data;
    data.t0min = t0min;
    data.t0max = t0max;
    ientry = dataList.size();
    dataList.push_back( std::pair<DTRangeT0Id,DTRangeT0Data>( key, data ) );
    dBuf->insert( chanKey.begin(), chanKey.end(), ientry );
    return 0;
  }

  return 99;

}


int DTRangeT0::set( const DTSuperLayerId& id,
                    int t0min,
                    int t0max ) {
  return set( id.wheel(),
              id.station(),
              id.sector(),
              id.superLayer(),
              t0min, t0max );
}


DTRangeT0::const_iterator DTRangeT0::begin() const {
  return dataList.begin();
}


DTRangeT0::const_iterator DTRangeT0::end() const {
  return dataList.end();
}


std::string DTRangeT0::mapName() const {
  std::string name = dataVersion + "_map_RangeT0";
  char nptr[100];
  sprintf( nptr, "%x", reinterpret_cast<unsigned int>( this ) );
  name += nptr;
  return name;
}


void DTRangeT0::cacheMap() const {

  std::string mName = mapName();
  DTBufferTree<int,int>* dBuf =
  DTDataBuffer<int,int>::openBuffer( mName );

  int entryNum = 0;
  int entryMax = dataList.size();
  while ( entryNum < entryMax ) {

    const DTRangeT0Id& chan = dataList[entryNum].first;

    std::vector<int> chanKey;
    chanKey.push_back( chan.  wheelId );
    chanKey.push_back( chan.stationId );
    chanKey.push_back( chan. sectorId );
    chanKey.push_back( chan.     slId );
    dBuf->insert( chanKey.begin(), chanKey.end(), entryNum++ );

  }

  return;

}

