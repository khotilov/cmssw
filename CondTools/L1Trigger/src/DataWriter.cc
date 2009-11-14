#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondTools/L1Trigger/interface/DataWriter.h"
#include "CondTools/L1Trigger/interface/Exception.h"
#include "CondCore/MetaDataService/interface/MetaData.h"
#include "CondCore/IOVService/interface/IOVService.h"
#include "CondCore/DBCommon/interface/ClassInfoLoader.h"
#include "CondCore/DBCommon/interface/Exception.h"

#include <utility>

namespace l1t
{
  DataWriter::DataWriter(){}
  DataWriter::~DataWriter(){}



std::string
DataWriter::writePayload( const edm::EventSetup& setup,
			  const std::string& recordType )
{
  WriterFactory* factory = WriterFactory::get();
  std::auto_ptr<WriterProxy> writer(factory->create( recordType + "@Writer" )) ;
  if( writer.get() == 0 )
    {
      throw cond::Exception( "DataWriter: could not create WriterProxy with name "
			     + recordType + "@Writer" ) ;
    }

  edm::Service<cond::service::PoolDBOutputService> poolDb;
  if (!poolDb.isAvailable())
    {
      throw cond::Exception( "DataWriter: PoolDBOutputService not available."
			     ) ;
    }
  cond::DbSession session = poolDb->session();
  cond::DbScopedTransaction tr(session);
  // if throw transaction will unroll
  tr.start(false);

  // update key to have new payload registered for record-type pair.
  std::string payloadToken = writer->save( setup, session ) ;

  edm::LogVerbatim( "L1-O2O" ) << recordType << " PAYLOAD TOKEN "
			       << payloadToken ;

  tr.commit ();

  return payloadToken ;
}

void
DataWriter::writeKeyList( L1TriggerKeyList* keyList,
			  edm::RunNumber_t sinceRun,
			  bool logTransactions )
{
  edm::Service<cond::service::PoolDBOutputService> poolDb;
  if( !poolDb.isAvailable() )
    {
      throw cond::Exception( "DataWriter: PoolDBOutputService not available."
			     ) ;
    }

  cond::DbSession session = poolDb->session();
  cond::DbScopedTransaction tr(session);
  tr.start(false);

  // Write L1TriggerKeyList payload
  pool::Ref<L1TriggerKeyList> ref = 
    session.storeObject(keyList,
			cond::classNameForTypeId(typeid(L1TriggerKeyList))
			);
			
  // Save payload token before committing.
  std::string payloadToken = ref.toString();
  
  // Commit before calling updateIOV(), otherwise PoolDBOutputService gets
  // confused.
  tr.commit ();
  
  // Set L1TriggerKeyList IOV
  updateIOV( "L1TriggerKeyListRcd",
	     payloadToken,
	     sinceRun,
	     logTransactions ) ;
}

bool
DataWriter::updateIOV( const std::string& esRecordName,
		       const std::string& payloadToken,
		       edm::RunNumber_t sinceRun,
		       bool logTransactions )
{
  edm::LogVerbatim( "L1-O2O" ) << esRecordName
			       << " PAYLOAD TOKEN " << payloadToken ;

  edm::Service<cond::service::PoolDBOutputService> poolDb;
  if (!poolDb.isAvailable())
    {
      throw cond::Exception( "DataWriter: PoolDBOutputService not available."
			     ) ;
    }

  bool iovUpdated = true ;

  if( poolDb->isNewTagRequest( esRecordName ) )
    {
      sinceRun = poolDb->beginOfTime() ;
      poolDb->createNewIOV( payloadToken,
			    sinceRun,
			    poolDb->endOfTime(),
			    esRecordName,
			    logTransactions ) ;
    }
  else
    {	
      cond::TagInfo tagInfo ;
      poolDb->tagInfo( esRecordName, tagInfo ) ;

      if( sinceRun == 0 ) // find last since and add 1
	{
	  sinceRun = tagInfo.lastInterval.first ;
	  ++sinceRun ;
	}

      if( tagInfo.lastPayloadToken != payloadToken )
	{
	  poolDb->appendSinceTime( payloadToken,
				   sinceRun,
				   esRecordName,
				   logTransactions ) ;
	}
      else
	{
	  iovUpdated = false ;
	  edm::LogVerbatim( "L1-O2O" ) << "IOV already up to date." ;
	}
    }

  if( iovUpdated )
    {
      edm::LogVerbatim( "L1-O2O" ) << esRecordName
				   << " SINCE " << sinceRun ;
    }

  return iovUpdated ;
}

std::string
DataWriter::payloadToken( const std::string& recordName,
			  edm::RunNumber_t runNumber )
{
  edm::Service<cond::service::PoolDBOutputService> poolDb;
  if( !poolDb.isAvailable() )
    {
      throw cond::Exception( "DataWriter: PoolDBOutputService not available."
			     ) ;
    }

  // Get tag corresponding to EventSetup record name.
  std::string iovTag = poolDb->tag( recordName ) ;

  // Get IOV token for tag.
  cond::DbSession session = poolDb->session();
  cond::DbScopedTransaction tr(session);
  tr.start(true);
  cond::MetaData metadata(session ) ;
  std::string iovToken ;
  if( metadata.hasTag( iovTag ) )
    {
      iovToken = metadata.getToken( iovTag ) ;
    }
  if( iovToken.empty() )
    {
      return std::string() ;
    }

  // Get payload token for run number.

  cond::IOVService iovService( session ) ;
  std::string payloadToken = iovService.payloadToken( iovToken, runNumber ) ;

  tr.commit() ;
  return payloadToken ;
}

} // ns
