/*
 *  See header file for a description of this class.
 *
 *  $Date: 2009/09/25 12:03:21 $
 *  $Revision: 1.2 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "CondTools/DT/interface/DTHVStatusHandler.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CondTools/DT/interface/DTDBSession.h"
#include "CondFormats/DTObjects/interface/DTHVStatus.h"
#include "CondCore/DBCommon/interface/AuthenticationMethod.h"
#include "CondCore/DBCommon/interface/DBSession.h"
#include "CondCore/DBCommon/interface/Connection.h"
#include "CondCore/DBCommon/interface/CoralTransaction.h"
#include "CondCore/DBCommon/interface/SessionConfiguration.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RelationalAccess/ISessionProxy.h"
#include "RelationalAccess/ISchema.h"
#include "RelationalAccess/ITable.h"
#include "RelationalAccess/ICursor.h"
#include "RelationalAccess/IQuery.h"
#include "RelationalAccess/TableDescription.h"
#include "RelationalAccess/ITableDataEditor.h"
#include "CoralBase/AttributeList.h"
#include "CoralBase/AttributeSpecification.h"
#include "CoralBase/Attribute.h"
#include "CoralBase/TimeStamp.h"

//---------------
// C++ Headers --
//---------------
#include <map>
#include <sys/time.h>

//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------
DTHVStatusHandler::DTHVStatusHandler( const edm::ParameterSet& ps ) :
 dataTag(               ps.getParameter<std::string> ( "tag" ) ),
 onlineConnect(         ps.getParameter<std::string> ( "onlineDB" ) ),
 onlineAuthentication(  ps.getParameter<std::string> ( 
                        "onlineAuthentication" ) ),
 bufferConnect(         ps.getParameter<std::string> ( "bufferDB" ) ),
 ySince(                ps.getParameter<int> ( "sinceYear"   ) ),
 mSince(                ps.getParameter<int> ( "sinceMonth"  ) ),
 dSince(                ps.getParameter<int> ( "sinceDay"    ) ),
 hSince(                ps.getParameter<int> ( "sinceHour"   ) ),
 pSince(                ps.getParameter<int> ( "sinceMinute" ) ),
 sSince(                ps.getParameter<int> ( "sinceSecond" ) ),
 yUntil(                ps.getParameter<int> ( "untilYear"   ) ),
 mUntil(                ps.getParameter<int> ( "untilMonth"  ) ),
 dUntil(                ps.getParameter<int> ( "untilDay"    ) ),
 hUntil(                ps.getParameter<int> ( "untilHour"   ) ),
 pUntil(                ps.getParameter<int> ( "untilMinute" ) ),
 sUntil(                ps.getParameter<int> ( "untilSecond" ) ),
 mapVersion(            ps.getParameter<std::string> ( "mapVersion"   ) ),
 splitVersion(          ps.getParameter<std::string> ( "splitVersion" ) ) {
  std::cout << " PopCon application for DT HV data export "
            << onlineAuthentication
            << std::endl;
  minHV = new float[4];
  maxHV = new float[4];
  minHV[0] = 3500.0;
  minHV[1] = 3500.0;
  minHV[2] = 1700.0;
  minHV[3] = 1100.0;
  maxHV[0] = 4000.0;
  maxHV[1] = 4000.0;
  maxHV[2] = 2200.0;
  maxHV[3] = 1600.0;
  maxCurrent = 3.0;
  maxPayload = 1000;
}

//--------------
// Destructor --
//--------------
DTHVStatusHandler::~DTHVStatusHandler() {
}

//--------------
// Operations --
//--------------
void DTHVStatusHandler::getNewObjects() {

  std::cout << "get new objects..." << std::endl;

// online DB connection

  cond::DBSession* omds_session;
  cond::Connection* omds_connect;
  cond::CoralTransaction* omds_transaction;

  std::cout << "open omds session... " << onlineAuthentication << std::endl;

  omds_session = new cond::DBSession();
  // to get the username/passwd from $CORAL_AUTH_PATH/authentication.xml
  omds_session->configuration().setAuthenticationMethod( cond::XML );
  omds_session->configuration().setAuthenticationPath( onlineAuthentication );
  // set message level to Error or Debug
  omds_session->configuration().setMessageLevel( cond::Error );
//  omds_session->connectionConfiguration().setConnectionRetrialTimeOut( 60 );
  omds_session->open();

  std::cout << "omds session open, start transaction" << std::endl;

  omds_connect = new cond::Connection( onlineConnect );
  omds_connect->connect( omds_session );
  omds_transaction = &( omds_connect->coralTransaction() );
  omds_transaction->start( true );

  std::cout << "omds transaction started" << std::endl;

  std::cout << "get omds session proxy... " << std::endl;
  omds_s_proxy = &( omds_transaction->coralSessionProxy() );
  std::cout << "omds session proxy got" << std::endl;
  omds_transaction->start( true );

// buffer DB connection

  cond::DBSession* buff_session;
  cond::Connection* buff_connect;
  cond::CoralTransaction* buff_transaction;

  std::cout << "open buffer session..." << std::endl;

  buff_session = new cond::DBSession();
  // to get the username/passwd from $CORAL_AUTH_PATH/authentication.xml
  buff_session->configuration().setAuthenticationMethod( cond::XML );
  buff_session->configuration().setAuthenticationPath( onlineAuthentication );
  // set message level to Error or Debug
  buff_session->configuration().setMessageLevel( cond::Error );
//  buff_session->connectionConfiguration().setConnectionRetrialTimeOut( 60 );
  buff_session->open();

  std::cout << "buffer session open, start transaction" << std::endl;

  buff_connect = new cond::Connection( bufferConnect );
  buff_connect->connect( buff_session );
  buff_transaction = &( buff_connect->coralTransaction() );
  buff_transaction->start( false );

  std::cout << "buffer transaction started" << std::endl;

  std::cout << "get buffer session proxy... " << std::endl;
  buff_s_proxy = &( buff_transaction->coralSessionProxy() );
  std::cout << "buffer session proxy got" << std::endl;
  buff_transaction->start( false );

// offline info

  //to access the information on the tag inside the offline database:
  cond::TagInfo const & ti = tagInfo();
  cond::Time_t last = ti.lastInterval.first;
  std::cout << "latest DCS data (HV) already copied until: "
            << last << std::endl;

  coral::TimeStamp coralSince( ySince, mSince, dSince,
                               hSince, pSince, sSince, 0 );
  procSince = condTime( coralSince );
  coral::TimeStamp coralUntil( yUntil, mUntil, dUntil,
                               hUntil, pUntil, sUntil, 0 );
  procUntil = condTime( coralUntil );
  lastFound = 0;
  nextFound = 0;
  timeLimit = 0;
  lastStamp = 0;

  if ( last == 0 ) {
    DTHVStatus* dummyStatus = new DTHVStatus( dataTag );
    cond::Time_t snc = 1;
    m_to_transfer.push_back( std::make_pair( dummyStatus, snc ) );
    last = procSince + 1;
    std::cout << "no old data... " << last << std::endl;
  }
  coral::TimeStamp coralLast = coralTime( last );
  coral::TimeStamp coralProc = coral::TimeStamp::now();
  cond::Time_t condProc = condTime( coralProc );

  if ( procSince > condProc ) {
      std::cout << "Required time interval in the future: "
                                                       << std::endl
                << " copy since " << ySince << " "
                                  << mSince << " "
                                  << dSince
                << " ( "          << procSince << " )" << std::endl
                << " current time " << coralProc.year( ) << " "
                                    << coralProc.month() << " "
                                    << coralProc.day(  ) << std::endl;
  }
  if ( procUntil > condProc ) procUntil = condProc;
  if ( last > procSince ) {
    if ( last < procUntil ) {
      procSince = last;
      checkNewData();
    }
    else {
      std::cout << "Required time interval already copied: "
                                                       << std::endl
                << " copy until " << yUntil << " "
                                  << mUntil << " "
                                  << dUntil
                << " ( "          << procUntil << " )" << std::endl
                << " data until " << coralLast.year( ) << " "
                                  << coralLast.month() << " "
                                  << coralLast.day(  ) << std::endl;
    }
  }
  else {
    std::cout << "Required time interval not contiguous with copied data: "
                                                     << std::endl
              << " data until " << coralLast.year( ) << " "
                                << coralLast.month() << " "
                                << coralLast.day(  ) << std::endl
              << " copy since " << ySince << " "
                                << mSince << " "
                                << dSince
              << " ( "          << procSince << " )" << std::endl;
  }

  delete omds_connect;
  delete omds_session;
  buff_transaction->commit();
  delete buff_connect;
  delete buff_session;
  return;

}


void DTHVStatusHandler::checkNewData() {

  //to access the information on last successful log entry for this tag:
//  cond::LogDBEntry const & lde = logDBEntry();     

  //to access the lastest payload (Ref is a smart pointer)
  Ref payload = lastPayload();

//  unsigned lastIOV = last;
//  cond::Time_t lastIOV = last;
//  std::cout << "check for new data since " << lastIOV << std::endl;
  std::cout << "check for new data since " << procSince << std::endl;

  std::set<std::string> lt( omds_s_proxy->nominalSchema().listTables() );
  std::set<std::string>::const_iterator iter = lt.begin();
  std::set<std::string>::const_iterator iend = lt.end();
  while ( iter != iend ) {
    const std::string& istr = *iter++;
    std::cout << "TABLE: " << istr << std::endl;
  }

  getLayerSplit();
  getChannelMap();

  std::cout << "open buffer db..." << std::endl;

  if ( !( buff_s_proxy->nominalSchema().existsTable( "HVSNAPSHOT" ) ) )
      createSnapshot();
  updateHVStatus();

  return;

}

std::string DTHVStatusHandler::id() const {
  return "DTHVStatusHandler";
}


void DTHVStatusHandler::getChannelMap() {

  if ( !( buff_s_proxy->nominalSchema().existsTable( "HVALIASES" ) ) ) {
    dumpHVAliases();
  }
  else {
    std::cout << "retrieve aliases table..." << std::endl;
    coral::ITable& hvalTable =
      buff_s_proxy->nominalSchema().tableHandle( "HVALIASES" );
    std::auto_ptr<coral::IQuery> hvalQuery( hvalTable.newQuery() );
    hvalQuery->addToOutputList( "DETID" );
    hvalQuery->addToOutputList(  "DPID" );
    coral::ICursor& hvalCursor = hvalQuery->execute();
    int chId;
    int dpId;
    while ( hvalCursor.next() ) {
      chId = hvalCursor.currentRow()["DETID"].data<int>();
      dpId = hvalCursor.currentRow()[ "DPID"].data<int>();
      aliasMap.insert( std::pair<int,int>( dpId, chId ) );
      layerMap.insert( std::pair<int,int>( chId, dpId ) );
    }
  }

  return;

}


void DTHVStatusHandler::getLayerSplit() {
  std::cout << "retrieve layer split table..." << std::endl;
  int whe;
  int sec;
  int sta;
  int qua;
  int lay;
  int l_p;
  int f_c;
  int l_c;
  coral::ITable& lsplTable =
    omds_s_proxy->nominalSchema().tableHandle( "DT_HV_LAYER_SPLIT" );
  std::auto_ptr<coral::IQuery> lsplQuery( lsplTable.newQuery() );
  coral::AttributeList versionBindVariableList;
  versionBindVariableList.extend( "version", typeid(std::string) );
  versionBindVariableList["version"].data<std::string>() = mapVersion;
  lsplQuery->setCondition( "VERSION=:version", versionBindVariableList );
  lsplQuery->addToOutputList( "WHEEL" );
  lsplQuery->addToOutputList( "SECTOR" );
  lsplQuery->addToOutputList( "STATION" );
  lsplQuery->addToOutputList( "SUPERLAYER" );
  lsplQuery->addToOutputList( "LAYER" );
  lsplQuery->addToOutputList( "PART" );
  lsplQuery->addToOutputList( "FIRST_CELL" );
  lsplQuery->addToOutputList( "LAST_CELL" );
  coral::ICursor& lsplCursor = lsplQuery->execute();
  while ( lsplCursor.next() ) {
    whe = lsplCursor.currentRow()["WHEEL"     ].data<int>();
    sec = lsplCursor.currentRow()["SECTOR"    ].data<int>();
    sta = lsplCursor.currentRow()["STATION"   ].data<int>();
    qua = lsplCursor.currentRow()["SUPERLAYER"].data<int>();
    lay = lsplCursor.currentRow()["LAYER"     ].data<int>();
    l_p = lsplCursor.currentRow()["PART"      ].data<int>();
    f_c = lsplCursor.currentRow()["FIRST_CELL"].data<int>();
    l_c = lsplCursor.currentRow()[ "LAST_CELL"].data<int>();
    DTWireId wireId( whe, sta, sec, qua, lay, 10 + l_p );
    laySplit.insert( std::pair<int,int>( wireId.rawId(), 
                                         ( f_c * 10000 ) + l_c ) );
  }
}


void DTHVStatusHandler::getChannelSplit() {
  std::cout << "retrieve channel split table..." << std::endl;
  int pwhe;
  int psec;
  int psta;
  int pqua;
  int play;
  int pl_p;
  int snum;
  int swhe;
  int ssec;
  int ssta;
  int squa;
  int slay;
  int sl_p;
  coral::ITable& csplTable =
    omds_s_proxy->nominalSchema().tableHandle( "DT_HV_CHANNEL_SPLIT" );
  std::auto_ptr<coral::IQuery> csplQuery( csplTable.newQuery() );
  coral::AttributeList versionBindVariableList;
  versionBindVariableList.extend( "version", typeid(std::string) );
  versionBindVariableList["version"].data<std::string>() = splitVersion;
  csplQuery->setCondition( "VERSION=:version", versionBindVariableList );
  csplQuery->addToOutputList( "P_WHEEL" );
  csplQuery->addToOutputList( "P_SECTOR" );
  csplQuery->addToOutputList( "P_STATION" );
  csplQuery->addToOutputList( "P_SUPERLAYER" );
  csplQuery->addToOutputList( "P_LAYER" );
  csplQuery->addToOutputList( "P_PART" );
  csplQuery->addToOutputList( "S_NUMBER" );
  csplQuery->addToOutputList( "S_WHEEL" );
  csplQuery->addToOutputList( "S_SECTOR" );
  csplQuery->addToOutputList( "S_STATION" );
  csplQuery->addToOutputList( "S_SUPERLAYER" );
  csplQuery->addToOutputList( "S_LAYER" );
  csplQuery->addToOutputList( "S_PART" );
  coral::ICursor& csplCursor = csplQuery->execute();
  while ( csplCursor.next() ) {
    pwhe = csplCursor.currentRow()["P_WHEEL"     ].data<int>();
    psec = csplCursor.currentRow()["P_SECTOR"    ].data<int>();
    psta = csplCursor.currentRow()["P_STATION"   ].data<int>();
    pqua = csplCursor.currentRow()["P_SUPERLAYER"].data<int>();
    play = csplCursor.currentRow()["P_LAYER"     ].data<int>();
    pl_p = csplCursor.currentRow()["P_PART"      ].data<int>();
    snum = csplCursor.currentRow()["S_NUMBER"    ].data<int>();
    swhe = csplCursor.currentRow()["S_WHEEL"     ].data<int>();
    ssec = csplCursor.currentRow()["S_SECTOR"    ].data<int>();
    ssta = csplCursor.currentRow()["S_STATION"   ].data<int>();
    squa = csplCursor.currentRow()["S_SUPERLAYER"].data<int>();
    slay = csplCursor.currentRow()["S_LAYER"     ].data<int>();
    sl_p = csplCursor.currentRow()["S_PART"      ].data<int>();
    DTWireId pId( pwhe, psta, psec, pqua, play, 10 + pl_p );
    DTWireId sId( swhe, ssta, ssec, squa, slay, 10 + sl_p );
    int pRaw = pId.rawId();
    int sRaw = sId.rawId();
    std::vector<int>* splitList = 0;
    std::map< int,std::vector<int>* >::iterator iter =
                                                channelSplit.find( pRaw );
    std::map< int,std::vector<int>* >::iterator iend =
                                                channelSplit.end();
    if ( iter == iend ) {
      channelSplit.insert( std::pair< int,
                                      std::vector<int>* >( pRaw, splitList =
                                      new std::vector<int> ) );
    }
    else {
      splitList = iter->second;
    }
    splitList->push_back( sRaw );
  }
  return;
}


void DTHVStatusHandler::dumpHVAliases() {

  std::cout << "DTHVStatusHandler::dumpHVAliases - begin" << std::endl;

  std::cout << "create aliases description..." << std::endl;
  coral::TableDescription hvalDesc;
  hvalDesc.setName( "HVALIASES" );
  hvalDesc.insertColumn( "DETID",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(int) ) );
  hvalDesc.insertColumn(  "DPID",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(int) ) );
  std::cout << "create aliases table..." << std::endl;
  coral::ITable& hvalTable = 
  buff_s_proxy->nominalSchema().createTable( hvalDesc );

  std::cout << "open DPNAME table..." << std::endl;
  std::map<int,std::string> idMap;
  coral::ITable& dpidTable =
    omds_s_proxy->nominalSchema().tableHandle( "DP_NAME2ID" );
  std::auto_ptr<coral::IQuery> dpidQuery( dpidTable.newQuery() );
  dpidQuery->addToOutputList( "ID" );
  dpidQuery->addToOutputList( "DPNAME" );
  coral::ICursor& dpidCursor = dpidQuery->execute();
  while( dpidCursor.next() ) {
    const coral::AttributeList& row = dpidCursor.currentRow();
    int id         = static_cast<int>( 0.01 +
                     row["ID"    ].data<float>() );
    std::string dp = row["DPNAME"].data<std::string>();
    idMap.insert( std::pair<int,std::string>( id, dp ) );
  }
  std::cout << "DPNAME table read... " << idMap.size() << std::endl;

  std::cout << "open ALIASES table..." << std::endl;
  std::map<std::string,std::string> cnMap;
  coral::ITable& nameTable =
    omds_s_proxy->nominalSchema().tableHandle( "ALIASES" );
  std::auto_ptr<coral::IQuery> nameQuery( nameTable.newQuery() );
  nameQuery->addToOutputList( "DPE_NAME" );
  nameQuery->addToOutputList( "ALIAS" );
  coral::ICursor& nameCursor = nameQuery->execute();
  while( nameCursor.next() ) {
    const coral::AttributeList& row = nameCursor.currentRow();
    std::string dp = row["DPE_NAME"].data<std::string>();
    std::string an = row["ALIAS"   ].data<std::string>();
    if ( an.length() < 20 ) continue;
    cnMap.insert( std::pair<std::string,std::string>( dp, an ) );
  }
  std::cout << "ALIASES table read... " << cnMap.size() << std::endl;

  std::map<int,std::string>::const_iterator idIter = idMap.begin();
  std::map<int,std::string>::const_iterator idIend = idMap.end();
  std::string outChk( "/outputChannel" );
  while ( idIter != idIend ) {
    const std::pair<int,std::string>& ientry = *idIter++;
    int dpId       = ientry.first;
    std::string dp = ientry.second;
    int ldp = dp.length();
    if ( ldp < 20 ) continue;
    std::string subOut( dp.substr( ldp - 17, 17 ) );
    std::string subChk( subOut.substr( 0, 14 ) );
    if ( subChk != outChk ) continue;
    std::string chName( dp.substr( 0, ldp - 17 ) );
    chName += ".actual.OvC";
    int chCode = subOut.c_str()[16] - '0';
    std::map<std::string,std::string>::const_iterator jter =
                                                      cnMap.find( chName );
    if ( jter == cnMap.end() ) continue;
    const std::pair<std::string,std::string>& jentry = *jter;
    std::cout << dp << std::endl << chName << " " << chCode << std::endl;
    std::string an( jentry.second );
    int al = an.length();
    int iofw = 7 + an.find( "DT_HV_W", 0 );
    int iofc = 3 + an.find( "_MB", 0 );
    int iofs = 2 + an.find( "_S" , 0 );
    int iofq = 3 + an.find( "_SL", 0 );
    int iofl = 2 + an.find( "_L" , 0 );
    if ( ( iofw == al ) ||
         ( iofc == al ) ||
         ( iofs == al ) ||
         ( iofq == al ) ||
         ( iofl == al ) ) {
      break;
    }
    int ioew = an.find( "_", iofw );
    int ioec = an.find( "_", iofc );
    int ioes = an.find( "_", iofs );
    int ioeq = an.find( "_", iofq );
    int ioel = an.find( "_", iofl );
    std::string swhe( an.substr( iofw, ioew - iofw ) );
    const char* cwhe = swhe.c_str();
    int whe = cwhe[1] - '0';
    if ( *cwhe != 'P' ) whe = -whe;

    std::string scha( an.substr( iofc, ioec - iofc ) );
    const char* ccha = scha.c_str();
    int cha = *ccha - '0';

    std::string ssec( an.substr( iofs, ioes - iofs ) );
    const char* csec = ssec.c_str();
    int sec = ( ( *csec - '0' ) * 10 ) + ( csec[1] - '0' );
    if ( ( csec[2] == 'R' ) && ( sec == 10 ) ) sec = 14;
    if ( ( csec[2] == 'L' ) && ( sec ==  4 ) ) sec = 13;

    std::string squa( an.substr( iofq, ioeq - iofq ) );
    const char* cqua = squa.c_str();
    int qua = *cqua - '0';

    std::string slay( an.substr( iofl, ioel - iofl ) );
    const char* clay = slay.c_str();
    int lay = *clay - '0';

    DTWireId wireId( whe, cha, sec, qua, lay, 10 + chCode );
    int chId = wireId.rawId();
    coral::AttributeList newChan;
    newChan.extend( "DETID", typeid(int) );
    newChan.extend(  "DPID", typeid(int) );
    newChan["DETID"].data<int>() = chId;
    newChan[ "DPID"].data<int>() = dpId;
    hvalTable.dataEditor().insertRow( newChan );
    aliasMap.insert( std::pair<int,int>( dpId, chId ) );
    layerMap.insert( std::pair<int,int>( chId, dpId ) );
  }

  std::cout << "DTHVStatusHandler::dumpHVAliases - end" << std::endl;
  return;
}


void DTHVStatusHandler::createSnapshot() {
  std::cout << "create snapshot description..." << std::endl;
  coral::TableDescription hvssDesc;
  hvssDesc.setName( "HVSNAPSHOT" );
  hvssDesc.insertColumn( "TIME",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(coral::TimeStamp) ) );
  hvssDesc.insertColumn( "WHEEL",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(int) ) );
  hvssDesc.insertColumn( "STATION",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(int) ) );
  hvssDesc.insertColumn( "SECTOR",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(int) ) );
  hvssDesc.insertColumn( "SUPERLAYER",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(int) ) );
  hvssDesc.insertColumn( "LAYER",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(int) ) );
  hvssDesc.insertColumn( "CHAN",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(int) ) );
  hvssDesc.insertColumn( "TYPE",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(int) ) );
  hvssDesc.insertColumn( "VALUE",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(float) ) );
  std::cout << "create snapshot table..." << std::endl;
  buff_s_proxy->nominalSchema().createTable( hvssDesc );
  coral::ITable& bufferTable = 
    buff_s_proxy->nominalSchema().tableHandle( "HVSNAPSHOT" );
  coral::AttributeList newMeas;
  newMeas.extend( "TIME",       typeid(coral::TimeStamp) );
  newMeas.extend( "WHEEL",      typeid(int) );
  newMeas.extend( "STATION",    typeid(int) );
  newMeas.extend( "SECTOR",     typeid(int) );
  newMeas.extend( "SUPERLAYER", typeid(int) );
  newMeas.extend( "LAYER",      typeid(int) );
  newMeas.extend( "CHAN",       typeid(int) );
  newMeas.extend( "TYPE",       typeid(int) );
  newMeas.extend( "VALUE",      typeid(float) );

  long long int zeroTime = 0LL;
  newMeas["TIME"     ].data<coral::TimeStamp>() =
                            coral::TimeStamp( zeroTime );
  newMeas["VALUE"    ].data<float>() = -999999.0;

  std::map<int,int>::const_iterator iter = aliasMap.begin();
  std::map<int,int>::const_iterator iend = aliasMap.end();
  while ( iter != iend ) {
    const std::pair<int,int>& entry= *iter++;
    int detId = entry.second;
    DTWireId chlId( detId );
    newMeas["WHEEL"     ].data<int>() = chlId.wheel     ();
    newMeas["STATION"   ].data<int>() = chlId.station   ();
    newMeas["SECTOR"    ].data<int>() = chlId.sector    ();
    newMeas["SUPERLAYER"].data<int>() = chlId.superLayer();
    newMeas["LAYER"     ].data<int>() = chlId.layer     ();
    newMeas["CHAN"      ].data<int>() = chlId.wire      () - 10;
    int itype;
    for ( itype = 1; itype <= 2; itype++ ) {
      newMeas["TYPE"      ].data<int>() = itype;
      bufferTable.dataEditor().insertRow( newMeas );
    }
  }

  std::cout << "create logging info..." << std::endl;
  if ( buff_s_proxy->nominalSchema().existsTable( "LOG" ) )
       buff_s_proxy->nominalSchema().  dropTable( "LOG" );
  coral::TableDescription infoDesc;
  infoDesc.setName( "LOG" );
  infoDesc.insertColumn( "EXECTIME",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(coral::TimeStamp) ) );
  infoDesc.insertColumn( "SNAPSHOT",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(coral::TimeStamp) ) );
  buff_s_proxy->nominalSchema().createTable( infoDesc );
  coral::AttributeList newInfo;
  newInfo.extend( "EXECTIME", typeid(coral::TimeStamp) );
  newInfo.extend( "SNAPSHOT", typeid(coral::TimeStamp) );
  newInfo["EXECTIME"].data<coral::TimeStamp>() =
                           coral::TimeStamp( zeroTime );
  newInfo["SNAPSHOT"].data<coral::TimeStamp>() =
                           coral::TimeStamp( zeroTime );
  coral::ITable& infoTable = 
    buff_s_proxy->nominalSchema().tableHandle( "LOG" );
  infoTable.dataEditor().insertRow( newInfo );

  return;

}


void DTHVStatusHandler::updateHVStatus() {
  std::map<int,timedMeasurement> snapshotValues;
  int missingChannels = recoverSnapshot( snapshotValues );
  cond::Time_t snapshotTime = recoverLastTime();
  std::cout << " snapshot at " << snapshotTime << " ( "
                               << coralTime( snapshotTime )
                                  .total_nanoseconds() << " ) "
                               << std::endl;
  if ( snapshotTime > procSince ) {
    coral::TimeStamp coralSnap = coralTime( snapshotTime );
    std::cout << "too recent snapshot: " << std::endl
              << " snapshot at " << coralSnap.year( ) << " "
                                 << coralSnap.month() << " "
                                 << coralSnap.day(  ) << std::endl
              << " copy since " << ySince << " "
                                << mSince << " "
                                << dSince
              << " ( "          << procSince << " )" << std::endl;
    return;
  }
  long long int dTime = 3600;//43200;
  dTime <<= 32;
  cond::Time_t condUntil = procSince;
  cond::Time_t condSince = condUntil - dTime;

  while ( missingChannels ) {
    std::cout << "back iteration: "
              << condSince << " ( " << coralTime( condSince )
                                       .total_nanoseconds() << " ) -> "
              << condUntil << " ( " << coralTime( condUntil )
                                       .total_nanoseconds() << " ) "
              << std::endl;
    if ( condSince <= snapshotTime ) condSince = snapshotTime;
    std::cout << "corrected since: "
              << condSince << " ( " << coralTime( condSince )
                                       .total_nanoseconds() << " ) "
              << std::endl;
    if ( condSince >= condUntil    ) break;
    std::cout << "missing... " << missingChannels << std::endl;
    checkForPeriod( condSince, condUntil,
                    snapshotValues, missingChannels, false );
    condUntil = condSince;
    condSince = condUntil - dTime;
  }

  dumpSnapshot( coralTime( procSince ), snapshotValues );

  copyHVData( snapshotValues );

  dumpSnapshot( coral::TimeStamp( lastStamp ), snapshotValues );

  return;
}

int DTHVStatusHandler::recoverSnapshot( std::map<int,timedMeasurement>&
                                         snapshotValues ) {
  int missingChannels = 0;
  std::map<int,int>::const_iterator layIter = layerMap.begin();
  std::map<int,int>::const_iterator layIend = layerMap.end();
  std::cout << "retrieve snapshot table..." << std::endl;
  coral::ITable& hvssTable =
         buff_s_proxy->nominalSchema().tableHandle( "HVSNAPSHOT" );
  std::auto_ptr<coral::IQuery> hvssQuery( hvssTable.newQuery() );
  hvssQuery->addToOutputList( "TIME" );
  hvssQuery->addToOutputList( "WHEEL" );
  hvssQuery->addToOutputList( "STATION" );
  hvssQuery->addToOutputList( "SECTOR" );
  hvssQuery->addToOutputList( "SUPERLAYER" );
  hvssQuery->addToOutputList( "LAYER" );
  hvssQuery->addToOutputList( "CHAN" );
  hvssQuery->addToOutputList( "TYPE" );
  hvssQuery->addToOutputList( "VALUE" );
  coral::ICursor& hvssCursor = hvssQuery->execute();
  while ( hvssCursor.next() ) {
    coral::TimeStamp time =
           hvssCursor.currentRow()["TIME"].data<coral::TimeStamp>();
    int     whe = hvssCursor.currentRow()["WHEEL"     ].data<int>();
    int     sta = hvssCursor.currentRow()["STATION"   ].data<int>();
    int     sec = hvssCursor.currentRow()["SECTOR"    ].data<int>();
    int     qua = hvssCursor.currentRow()["SUPERLAYER"].data<int>();
    int     lay = hvssCursor.currentRow()["LAYER"     ].data<int>();
    int     l_p = hvssCursor.currentRow()["CHAN"      ].data<int>();
    int     mty = hvssCursor.currentRow()["TYPE"      ].data<int>();
    float value = hvssCursor.currentRow()["VALUE"     ].data<float>();
    if ( mty > 2 ) continue;
    DTWireId wireId( whe, sta, sec, qua, lay, 10 + l_p );
    layIter = layerMap.find( wireId.rawId() );
    if ( layIter == layIend ) {
      continue;
    }
    int dpId = ( layIter->second * 10 ) + mty;
    snapshotValues.insert( std::pair<int,timedMeasurement>( 
                                   dpId, timedMeasurement( 
                                   time.total_nanoseconds(), value ) ) );
    missingChannels++;
  }
  return missingChannels;
}


cond::Time_t DTHVStatusHandler::recoverLastTime() {
  coral::ITable& infoTable =
         buff_s_proxy->nominalSchema().tableHandle( "LOG" );
  std::auto_ptr<coral::IQuery> infoQuery( infoTable.newQuery() );
  infoQuery->addToOutputList( "SNAPSHOT" );
  coral::ICursor& infoCursor = infoQuery->execute();
  coral::TimeStamp time;
  while ( infoCursor.next() ) {
    time = infoCursor.currentRow()["SNAPSHOT"].data<coral::TimeStamp>();
  }
  return condTime( time );
}


void DTHVStatusHandler::dumpSnapshot( const coral::TimeStamp& time,
                                      std::map<int,timedMeasurement>&
                                      snapshotValues ) {

  std::cout << "dump snapshot to buffer db..." << std::endl;
  std::string emptyCondition( "" );
  coral::AttributeList emptyBindVariableList;
  std::map<int,int>::const_iterator mapIter = aliasMap.begin();
  std::map<int,int>::const_iterator mapIend = aliasMap.end();
  coral::ITable& hvssTable =
         buff_s_proxy->nominalSchema().tableHandle( "HVSNAPSHOT" );
  coral::ITableDataEditor& hvssEditor( hvssTable.dataEditor() );
  long nRows = hvssEditor.deleteRows( emptyCondition, emptyBindVariableList );
  std::cout << nRows << " rows deleted" << std::endl;

  coral::AttributeList newMeas;
  newMeas.extend( "TIME",       typeid(coral::TimeStamp) );
  newMeas.extend( "WHEEL",      typeid(int) );
  newMeas.extend( "STATION",    typeid(int) );
  newMeas.extend( "SECTOR",     typeid(int) );
  newMeas.extend( "SUPERLAYER", typeid(int) );
  newMeas.extend( "LAYER",      typeid(int) );
  newMeas.extend( "CHAN",       typeid(int) );
  newMeas.extend( "TYPE",       typeid(int) );
  newMeas.extend( "VALUE",      typeid(float) );

  nRows = 0;
  std::map<int,timedMeasurement>::const_iterator ssvIter =
                                  snapshotValues.begin();
  std::map<int,timedMeasurement>::const_iterator ssvIend =
                                  snapshotValues.end();
  while ( ssvIter != ssvIend ) {
    const std::pair<int,timedMeasurement>& entry = *ssvIter++;
    int dpty = entry.first;
    int dpId = dpty / 10;
    int type = dpty % 10;
    mapIter = aliasMap.find( dpId );
    if ( mapIter == mapIend ) continue;
    DTWireId chlId( mapIter->second );
    const timedMeasurement& tMeas = entry.second;
    long long int newTime = tMeas.first;
    newMeas["TIME"      ].data<coral::TimeStamp>() =
                               coral::TimeStamp( newTime );
    newMeas["WHEEL"     ].data<int>() = chlId.wheel     ();
    newMeas["STATION"   ].data<int>() = chlId.station   ();
    newMeas["SECTOR"    ].data<int>() = chlId.sector    ();
    newMeas["SUPERLAYER"].data<int>() = chlId.superLayer();
    newMeas["LAYER"     ].data<int>() = chlId.layer     ();
    newMeas["CHAN"      ].data<int>() = chlId.wire      () - 10;
    newMeas["TYPE"      ].data<int>() = type;
    newMeas["VALUE"     ].data<float>() = tMeas.second;
    hvssEditor.insertRow( newMeas );
    nRows++;
  }
  std::cout << nRows << " rows updated" << std::endl;

  std::cout << "create logging info..." << std::endl;
  if ( buff_s_proxy->nominalSchema().existsTable( "LOG" ) )
       buff_s_proxy->nominalSchema().  dropTable( "LOG" );
  coral::TableDescription infoDesc;
  infoDesc.setName( "LOG" );
  infoDesc.insertColumn( "EXECTIME",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(coral::TimeStamp) ) );
  infoDesc.insertColumn( "SNAPSHOT",
                         coral::AttributeSpecification::typeNameForId( 
                         typeid(coral::TimeStamp) ) );
  buff_s_proxy->nominalSchema().createTable( infoDesc );
  coral::AttributeList newInfo;
  newInfo.extend( "EXECTIME", typeid(coral::TimeStamp) );
  newInfo.extend( "SNAPSHOT", typeid(coral::TimeStamp) );
  newInfo["EXECTIME"].data<coral::TimeStamp>() = coral::TimeStamp::now();
  newInfo["SNAPSHOT"].data<coral::TimeStamp>() = time;
  coral::ITable& infoTable = 
    buff_s_proxy->nominalSchema().tableHandle( "LOG" );
  infoTable.dataEditor().insertRow( newInfo );

  return;

}


int DTHVStatusHandler::checkForPeriod( cond::Time_t condSince,
                                       cond::Time_t condUntil,
                                       std::map<int,timedMeasurement>&
                                                snapshotValues,
                                       int& missingChannels,
                                       bool copyOffline ) {

  std::map<int,timedMeasurement>::iterator mapIter = snapshotValues.begin();
  std::map<int,timedMeasurement>::iterator mapIend = snapshotValues.end();

  std::map<long long int,channelValue> periodBuffer;

  coral::ITable& fwccTable =
    omds_s_proxy->nominalSchema().tableHandle( "FWCAENCHANNEL" );
  std::auto_ptr<coral::IQuery> fwccQuery( fwccTable.newQuery() );
  fwccQuery->addToOutputList( "DPID"          );
  fwccQuery->addToOutputList( "CHANGE_DATE"   );
  fwccQuery->addToOutputList( "ACTUAL_VMON"   );
  fwccQuery->addToOutputList( "ACTUAL_IMON"   );
  fwccQuery->addToOutputList( "ACTUAL_ISON"   );
  fwccQuery->addToOutputList( "ACTUAL_STATUS" );
  fwccQuery->addToOutputList( "ACTUAL_OVC"    );
  coral::AttributeList timeBindVariableList;
  timeBindVariableList.extend( "since", typeid(coral::TimeStamp) );
  timeBindVariableList.extend( "until", typeid(coral::TimeStamp) );
  coral::TimeStamp coralSince = coralTime( condSince );
  coral::TimeStamp coralUntil = coralTime( condUntil );
  std::cout << "look for data since "            
            << coralSince.year( ) << " "
            << coralSince.month() << " "
            << coralSince.day(  ) << " "
            << coralSince.hour( ) << " until "
            << coralUntil.year( ) << " "
            << coralUntil.month() << " "
            << coralUntil.day(  ) << " "
            << coralUntil.hour( ) << std::endl;
  timeBindVariableList["since"].data<coral::TimeStamp>() =
                                     coralTime( condSince );
  timeBindVariableList["until"].data<coral::TimeStamp>() =
                                     coralTime( condUntil );
  fwccQuery->setCondition( "CHANGE_DATE>:since and CHANGE_DATE<:until",
                           timeBindVariableList );
  fwccQuery->addToOrderList( "CHANGE_DATE" );
  coral::ICursor& fwccCursor = fwccQuery->execute();
  int nrows = 0;
  while ( fwccCursor.next() ) {
    nrows++;
    const coral::Attribute& dp     = fwccCursor.currentRow()["DPID"         ];
    const coral::Attribute& vmon   = fwccCursor.currentRow()["ACTUAL_VMON"  ];
    const coral::Attribute& imon   = fwccCursor.currentRow()["ACTUAL_IMON"  ];
    coral::TimeStamp changeTime =
           fwccCursor.currentRow()["CHANGE_DATE"].data<coral::TimeStamp>();
    long long int cTimeValue = changeTime.total_nanoseconds();
    if ( !copyOffline ) cTimeValue = -cTimeValue;
    if ( dp.isNull() ) {
      std::cout << "------- " << nrows << std::endl;
      continue;
    }
    int dpId = 10 * static_cast<int>( 0.01 +
           fwccCursor.currentRow()["DPID"].data<float>() );
    if ( !( vmon.isNull() ) ) {
      while ( periodBuffer.find( cTimeValue ) !=
              periodBuffer.end() ) cTimeValue++;
      int chan = dpId + 1;
      periodBuffer.insert( std::pair<long long int,channelValue> (
                                     cTimeValue, channelValue( chan,
                                                 vmon.data<float>() ) ) );
    }
    if ( !( imon.isNull() ) ) {
      while ( periodBuffer.find( cTimeValue ) !=
              periodBuffer.end() ) cTimeValue++;
      int chan = dpId + 2;
      periodBuffer.insert( std::pair<long long int,channelValue> (
                                     cTimeValue, channelValue( chan,
                                                 imon.data<float>() ) ) );
    }
  }

  long long int dTime = 1;
  dTime <<= 32;
  std::cout << "data found in period: " << periodBuffer.size() << std::endl;
  std::map<long long int,channelValue>::const_iterator bufIter =
                                                       periodBuffer.begin();
  std::map<long long int,channelValue>::const_iterator bufIend =
                                                       periodBuffer.end();

  bool changedStatus = false;
  while ( bufIter != bufIend ) {
    const std::pair<long long int,channelValue>& entry = *bufIter++;
    long long int mTime = entry.first;
    if ( !copyOffline ) mTime = -mTime;
    channelValue cValue = entry.second;
    int   chan = cValue.first;
    float cont = cValue.second;
    mapIter = snapshotValues.find( chan );
    if ( ( mapIter != mapIend ) &&
         ( mapIter->second.first < mTime ) ) {
      nextFound = condTime( mTime );
      if ( changedStatus ) {
        if ( nextFound > timeLimit ) {
          DTHVStatus* hvStatus = offlineList( snapshotValues );
          m_to_transfer.push_back( std::make_pair( hvStatus, lastFound ) );
          changedStatus = false;
          if ( !( --maxPayload ) ) {
            procUntil = lastFound;
            break;
          }
        }
      }
      if ( copyOffline && !changedStatus &&
           checkStatusChange( chan, mapIter->second.second, cont ) ) {
        timeLimit = nextFound + dTime;
        changedStatus = true;
      }
      mapIter->second = timedMeasurement( lastStamp = mTime, cont );
      lastFound = nextFound;
      missingChannels--;
    }
  }

  std::cout << nrows << std::endl;
  return nrows;

}


void DTHVStatusHandler::copyHVData( std::map<int,timedMeasurement>&
                                         snapshotValues ) {
  long long int dTime = 3600;//43200;
  dTime <<= 32;

  cond::Time_t condSince = procSince;
  cond::Time_t condUntil = condSince + dTime;
  if ( condUntil > procUntil ) condUntil = procUntil;

  int dum = 0;
  lastStatus = 0;
  while ( condSince < condUntil ) {
    checkForPeriod( condSince, condUntil, snapshotValues, dum, true );
    condSince = condUntil;
    condUntil = condSince + dTime;
    if ( condUntil > procUntil ) condUntil = procUntil;
  }
  if ( ( lastFound != 0 ) && ( maxPayload > 0 ) ) {
    DTHVStatus* hvStatus = offlineList( snapshotValues );
    m_to_transfer.push_back( std::make_pair( hvStatus, lastFound ) );
  }

  return;
}


DTHVStatus* DTHVStatusHandler::offlineList(
            std::map<int,timedMeasurement>& snapshotValues ) {
  DTHVStatus* hv = new DTHVStatus( dataTag );
  std::map<int,timedMeasurement>::const_iterator mapIter =
                                                 snapshotValues.begin();
  std::map<int,timedMeasurement>::const_iterator mapIend =
                                                 snapshotValues.end();
  std::map<int,int>::const_iterator aliasIter = aliasMap.begin();
  std::map<int,int>::const_iterator aliasIend = aliasMap.end();
  while ( mapIter != mapIend ) {
    const std::pair<int,timedMeasurement>& entry = *mapIter++;
    int chan = entry.first;
    const timedMeasurement& tMeas = entry.second;
    float value = tMeas.second;
    int dpId = chan / 10;
    int type = chan % 10;
    if ( type > 2 ) continue;
    aliasIter = aliasMap.find( dpId );
    if ( aliasIter == aliasIend ) continue;
    int rawId = aliasIter->second;
    setFlags( hv, rawId, type, value );
    std::map< int,std::vector<int>* >::const_iterator m_iter =
                                       channelSplit.find( rawId );
    std::map< int,std::vector<int>* >::const_iterator m_iend =
                                       channelSplit.end();
    if ( m_iter != m_iend ) {
      std::vector<int>* cList = m_iter->second;
      std::vector<int>::const_iterator l_iter = cList->begin();
      std::vector<int>::const_iterator l_iend = cList->end();
      while ( l_iter != l_iend ) {
        setFlags( hv, *l_iter++, type, value );
      }
    }
  }
  return hv;
}


void DTHVStatusHandler::setFlags( DTHVStatus* hv,
                                  int rawId, int type, float value ) {
  DTWireId chlId( rawId );
  int whe = chlId.wheel     ();
  int sta = chlId.station   ();
  int sec = chlId.sector    ();
  int qua = chlId.superLayer();
  int lay = chlId.layer     ();
  int l_p = chlId.wire      () - 10;
  int chanError = checkCurrentStatus( l_p, type, value );
  if ( !chanError ) return;
  switch ( l_p ) {
  case 0:
    setChannelFlag( hv, whe, sta, sec, qua, lay, 0, 'A', chanError );
    break;
  case 1:
    setChannelFlag( hv, whe, sta, sec, qua, lay, 1, 'A', chanError );
    break;
  case 2:
    setChannelFlag( hv, whe, sta, sec, qua, lay, 0, 'C', chanError );
    setChannelFlag( hv, whe, sta, sec, qua, lay, 1, 'C', chanError );
    break;
  case 3:
    setChannelFlag( hv, whe, sta, sec, qua, lay, 0, 'S', chanError );
    setChannelFlag( hv, whe, sta, sec, qua, lay, 1, 'S', chanError );
    break;
  }
  return;
}


void DTHVStatusHandler::setChannelFlag( DTHVStatus* hv,
                                        int whe, int sta, int sec,
                                        int qua, int lay, int l_p,
                                        char cht, int err ) {
  int fCell = 0;
  int lCell = 99;
  int flagA = 0;
  int flagC = 0;
  int flagS = 0;
  int searchStatus = hv->get( whe, sta, sec, qua, lay, l_p,
                              fCell, lCell, flagA, flagC, flagS );
  if ( searchStatus ) {
    DTWireId wireId( whe, sta, sec, qua, lay, 10 + l_p );
    std::map<int,int>::const_iterator splitIter =
                                      laySplit.find( wireId.rawId() );
    std::map<int,int>::const_iterator splitIend =
                                      laySplit.end();
    if ( splitIter != splitIend ) {
      int code = splitIter->second;
      fCell = code / 10000;
      lCell = code % 10000;
    }
  }
  switch ( cht ) {
  case 'A':
    flagA |= err;
    break;
  case 'C':
    flagC |= err;
    break;
  case 'S':
    flagS |= err;
    break;
  default:
    break;
  }
  hv->set( whe, sta, sec, qua, lay, l_p,
           fCell, lCell, flagA, flagC, flagS );
  return;
}


int DTHVStatusHandler::checkCurrentStatus( int part, int type, float value ) {
  if ( part < 0 ) return 0;
  if ( part > 3 ) return 0;
  int status = 0;
  if ( type == 1 ) {
    if ( value < minHV[part] ) status += 1;
    if ( value > maxHV[part] ) status += 2;
  }
  if ( type == 2 ) {
    float maxCurrent = 10.0;
    if ( value > maxCurrent ) status += 4;
  }
  return status;
}


int DTHVStatusHandler::checkStatusChange( int chan,
                                          float oldValue, float newValue ) {
  int dpId = chan / 10;
  int type = chan % 10;
  std::map<int,int>::const_iterator aliasIter = aliasMap.find( dpId );
  std::map<int,int>::const_iterator aliasIend = aliasMap.end();
  if ( aliasIter == aliasIend ) return false;
  DTWireId chlId( aliasIter->second );
  int l_p = chlId.wire() - 10;
  int oldStatus = checkCurrentStatus( l_p, type, oldValue );
  int newStatus = checkCurrentStatus( l_p, type, newValue );
  if ( newStatus == oldStatus ) return 0;
  if ( newStatus ) return +1;
  return -1;
}


coral::TimeStamp DTHVStatusHandler::coralTime( const  cond::Time_t&    time ) {
  long long int iTime = ( ( ( ( time >> 32 ) & 0xFFFFFFFF ) * 1000000000 ) +
                          ( (   time         & 0xFFFFFFFF ) * 1000       ) );
  coral::TimeStamp cTime( iTime );
  return cTime;
}


cond::Time_t     DTHVStatusHandler::condTime(  const coral::TimeStamp& time ) {
  cond::Time_t cTime = ( ( time.total_nanoseconds() / 1000000000 )  << 32 ) + 
                       ( ( time.total_nanoseconds() % 1000000000 ) / 1000 );
  return cTime;
}


cond::Time_t     DTHVStatusHandler::condTime(  long long int           time ) {
  cond::Time_t cTime = ( ( time                     / 1000000000 )  << 32 ) + 
                       ( ( time                     % 1000000000 ) / 1000 );
  return cTime;
}


