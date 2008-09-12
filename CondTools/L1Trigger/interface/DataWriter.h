#ifndef CondTools_L1Trigger_DataWriter_h
#define CondTools_L1Trigger_DataWriter_h

// Framework
#include "FWCore/Framework/interface/IOVSyncValue.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/DataKey.h"

#include "CondCore/DBCommon/interface/TypedRef.h"
#include "CondCore/MetaDataService/interface/MetaData.h"

#include "DataFormats/Provenance/interface/RunID.h"

// L1T includes
#include "CondFormats/L1TObjects/interface/L1TriggerKey.h"
#include "CondFormats/L1TObjects/interface/L1TriggerKeyList.h"

#include "CondTools/L1Trigger/interface/DataManager.h"
#include "CondTools/L1Trigger/interface/WriterProxy.h"

#include <string>
#include <map>

namespace l1t
{

/* This class is used to write L1 Trigger configuration data to Pool DB. Later this data can be
 * read from pool DB with class l1t::DataReader.
 *
 * L1TriggerKey is always writen to database by metod writeKey (...). Payload can be writen to database
 * in two ways. If user knows data type that she or he is interested to write to database, should use
 * templated version of writePayload. If that is not the case, then non-template version of write payload should
 * be used. Both of these methods will write payload to database and updates passed L1TriggerKey to reflect this
 * new change. This means that L1TriggerKey should be writen to database last, after all payload has been writen.
 *
 * Non-template version of writePayload was writen in such way that it could take data from EventSetup. Other options
 * also can be implemented.
 *
 * In order to use non-template version of this class, user has to make sure to register datatypes that she or he is
 * interested to write to the framework. This should be done with macro REGISTER_L1_WRITER(record, type) found in
 * WriterProxy.h file. Also, one should take care to register these data types to CondDB framework with macro
 * REGISTER_PLUGIN(record, type) from registration_macros.h found in PluginSystem.
 *
 * Validity interval of all data is controled via L1TriggerKey.
 */
class DataWriter : public DataManager
{
    public:
        /* Constructors. This system will use pool db that is configured via provided parameters.
         * connect - connection string to pool db, e.g. sqlite_file:test.db
         * catalog - catalog that should be used for this connection. e.g. file:test.xml
         */
  //        explicit DataWriter (const std::string & connect, const std::string & catalog)
  //            : DataManager (connect, catalog) {};
  explicit DataWriter (const std::string& connectString,
		       const std::string& authenticationPath )
    : DataManager( connectString, authenticationPath ) {};
        virtual ~DataWriter () {};

        /* Data writting functions */

        /* Writes payload to DB and assignes it to provided L1TriggerKey. Key is updated
         * to reflect changes.
         */
        template<typename T>
        void writePayload (L1TriggerKey & key, T * payload, const std::string & recordName);

        void writePayload (L1TriggerKey & key, const edm::EventSetup & setup,
                const std::string & record, const std::string & type);

        /* Writes given key to DB and starts its IOV from provided run number.
         * From here pointer ownership is managed by POOL
         */
        void writeKey (L1TriggerKey * key,
		       const std::string & tag,
		       const edm::RunNumber_t sinceRun);

      // Added by wsun 03/2008

      // Get payload from EventSetup and write to DB with no IOV
      // recordType = "record@type", return value is payload token
      std::string writePayload( const edm::EventSetup& setup,
				const std::string& recordType ) ;

      // Write L1TriggerKeyList to DB
      void writeKeyList( L1TriggerKeyList* keyList,
			 const std::string& tag, // tag for IOV sequence
			 edm::RunNumber_t sinceRun = 0 ) ;

      // Append IOV with sinceRun to IOV sequence with given tag
      void updateIOV( const std::string& tag,
		      const std::string& payloadToken,
		      const edm::RunNumber_t sinceRun ) ;

    protected:
        /* Helper method that maps tag with iovToken */
        void addMappings (const std::string tag, const std::string iovToken);
        std::string findTokenForTag (const std::string & tag);

        // map to speed up detection if we already have mapping from given tag to token
        typedef std::map<std::string, std::string> TagToToken;
        TagToToken tagToToken;
};

/* Part of the implementation: Template functions that can't go to cc failes */

template<typename T>
void DataWriter::writePayload (L1TriggerKey & key, T * payload, const std::string & recordName)
{
/*     pool->connect (); */
  //connection->connect ( session );
/*     pool->startTransaction (false); */
  cond::PoolTransaction& pool = connection->poolTransaction() ;
    pool.start (false);

    cond::TypedRef<T> ref (pool, payload);
    ref.markWrite (recordName);

    std::string typeName = 
        edm::eventsetup::heterocontainer::HCTypeTagTemplate<T, edm::eventsetup::DataKey>::className ();

    key.add (recordName, typeName, ref.token ());

    pool.commit ();
/*     pool->disconnect (); */
    //connection->disconnect ();
}

} // ns

#endif
