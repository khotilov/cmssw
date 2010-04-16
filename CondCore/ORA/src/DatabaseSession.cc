#include "CondCore/ORA/interface/Configuration.h"
#include "CondCore/ORA/interface/Exception.h"
#include "DatabaseSession.h"
#include "IDatabaseSchema.h"
#include "Sequences.h"
#include "MappingDatabase.h"
#include "TransactionCache.h"
#include "DatabaseContainer.h"
#include "ClassUtils.h"
#include "MappingRules.h"
// externals
#include "RelationalAccess/ISessionProxy.h"
#include "RelationalAccess/ITransaction.h"

ora::ContainerUpdateTable::ContainerUpdateTable():
  m_table(){
}

ora::ContainerUpdateTable::~ContainerUpdateTable(){
}

void ora::ContainerUpdateTable::takeNote( int contId,
                                          unsigned int size ){
  std::map<int, unsigned int>::iterator iC = m_table.find( contId );
  if( iC == m_table.end() ){
    iC = m_table.insert( std::make_pair( contId, 0 ) ).first;
  }
  iC->second = size;
}

const std::map<int, unsigned int>& ora::ContainerUpdateTable::table(){
  return m_table;
}

void ora::ContainerUpdateTable::clear(){
  m_table.clear();
}

ora::DatabaseSession::DatabaseSession( Configuration& configuration ):
  m_connectionPool( new ConnectionPool ),
  m_dbSession(),
  m_connectionString( "" ),
  m_schema(),
  m_contIdSequence(),
  m_mappingDb(),
  m_transactionCache(),
  m_containerUpdateTable(),
  m_configuration( configuration ){
}

ora::DatabaseSession::DatabaseSession(boost::shared_ptr<ConnectionPool>& connectionPool,
                                      Configuration& configuration ):
  m_connectionPool( connectionPool ),
  m_dbSession(),
  m_connectionString( "" ),
  m_schema(),
  m_contIdSequence(),
  m_mappingDb(),
  m_transactionCache(),
  m_containerUpdateTable(),
  m_configuration( configuration ){
}

ora::DatabaseSession::~DatabaseSession(){
  disconnect();
}

bool ora::DatabaseSession::connect( const std::string& connectionString,
                                    bool readOnly ){
  m_dbSession = m_connectionPool->connect( connectionString, readOnly?coral::ReadOnly:coral::Update );
  if(m_dbSession.isValid()) {
    m_connectionString = connectionString;
  }
  return isConnected();
}

void ora::DatabaseSession::clearTransaction(){
  m_transactionCache.reset();
  m_mappingDb.reset();
  m_contIdSequence.reset();
  m_schema.reset();
  m_containerUpdateTable.clear();
}

void ora::DatabaseSession::disconnect(){
  if( isConnected() ){
    if( isTransactionActive()) rollbackTransaction();
  }
  clearTransaction();
  m_dbSession.close();
  m_connectionString.clear();
}

bool ora::DatabaseSession::isConnected(){
  return m_dbSession.isValid();
}

const std::string& ora::DatabaseSession::connectionString(){
  return m_connectionString;
}

void ora::DatabaseSession::startTransaction( bool readOnly ){
  if( !m_transactionCache.get() ){
    m_dbSession.get().transaction().start( readOnly );
    m_schema.reset( IDatabaseSchema::createSchemaHandle( m_dbSession.get().nominalSchema() ));
    m_contIdSequence.reset( new NamedSequence( MappingRules::sequenceNameForContainerId(), *m_schema ));
    m_mappingDb.reset( new MappingDatabase( *m_schema ));
    m_transactionCache.reset( new TransactionCache );
  }
}

void ora::DatabaseSession::commitTransaction(){
  if( m_transactionCache.get() ){
    m_schema->containerHeaderTable().updateNumberOfObjects( m_containerUpdateTable.table() );
    m_dbSession.get().transaction().commit();
    clearTransaction();
  }
}

void ora::DatabaseSession::rollbackTransaction(){
  if( m_transactionCache.get() ){
    m_dbSession.get().transaction().rollback();
    clearTransaction();
  }
}

bool ora::DatabaseSession::isTransactionActive(){
  return m_transactionCache.get() != 0;
}

bool ora::DatabaseSession::exists(){
  if(!m_transactionCache->dbExistsLoaded()){
    m_transactionCache->setDbExists( m_schema->exists() );
  }
  return m_transactionCache->dbExists();
}

void ora::DatabaseSession::create(){
  m_schema->create();
  m_transactionCache->setDbExists( true );
}

void ora::DatabaseSession::drop(){
  m_schema->drop();
  m_transactionCache->setDbExists( false );
}

void ora::DatabaseSession::open(){
  if(!m_transactionCache->isLoaded()){
    std::map<std::string, ContainerHeaderData> containersData;
    m_schema->containerHeaderTable().getContainerData( containersData );
    for(std::map<std::string, ContainerHeaderData>::iterator iC = containersData.begin();
        iC != containersData.end(); ++iC){
      Handle<DatabaseContainer> container( new DatabaseContainer( iC->second.id, iC->first,
                                                                  iC->second.className,
                                                                  iC->second.numberOfObjects, *this ) );
      m_transactionCache->addContainer( iC->second.id, iC->first, container );
    }
    m_transactionCache->setLoaded();
  }
}

ora::Handle<ora::DatabaseContainer> ora::DatabaseSession::createContainer( const std::string& containerName,
                                                                           const Reflex::Type& type ){
  // create the container
  int newContId = m_contIdSequence->getNextId( true );
  Handle<DatabaseContainer> newCont ( new DatabaseContainer( newContId, containerName, type, *this ) );
  m_transactionCache->addContainer( newContId, containerName, newCont );
  m_schema->containerHeaderTable().addContainer( newContId, containerName, newCont->className() );
  newCont->create();
  return newCont;
}

void ora::DatabaseSession::dropContainer( const std::string& name ){
  Handle<DatabaseContainer> cont = m_transactionCache->getContainer( name );
  cont->drop();
  m_transactionCache->eraseContainer( cont->id(), name );
  m_schema->containerHeaderTable().removeContainer( cont->id() );
}

ora::Handle<ora::DatabaseContainer> ora::DatabaseSession::containerHandle( const std::string& name ){
  return m_transactionCache->getContainer( name );
}

ora::Handle<ora::DatabaseContainer> ora::DatabaseSession::containerHandle( int contId ){
  return  m_transactionCache->getContainer( contId );
}

const std::map<int, ora::Handle<ora::DatabaseContainer> >& ora::DatabaseSession::containers(){
  return m_transactionCache->containers();
}

ora::IDatabaseSchema& ora::DatabaseSession::schema(){
  return *m_schema;  
}

ora::NamedSequence& ora::DatabaseSession::containerIdSequence(){
  return *m_contIdSequence;
}

ora::MappingDatabase& ora::DatabaseSession::mappingDatabase(){
  return *m_mappingDb;  
}

ora::Configuration& ora::DatabaseSession::configuration(){
  return m_configuration;
}

ora::ContainerUpdateTable& ora::DatabaseSession::containerUpdateTable(){
  return m_containerUpdateTable;
}
