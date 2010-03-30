#include "CondCore/ORA/interface/Configuration.h"
#include "CondCore/ORA/interface/Exception.h"
#include "ContainerSchema.h"
#include "DatabaseSession.h"
#include "IDatabaseSchema.h"
#include "MappingToSchema.h"
#include "MappingDatabase.h"
#include "MappingGenerator.h"
#include "MappingRules.h"
#include "ClassUtils.h"
// externals
#include "RelationalAccess/ISchema.h"

namespace ora {

  void getTableHierarchyFromMappingElement( const MappingElement& source,
                                            std::map<std::string, std::set<std::string> >& tableList ){
    const std::string& tableName = source.tableName();
    std::map<std::string, std::set<std::string> >::iterator iTab = tableList.find( tableName );
    if( iTab ==tableList.end() ){
      std::set<std::string> dependencies;
      tableList.insert(std::make_pair( tableName, dependencies ) );
    }
    for( MappingElement::const_iterator iElem = source.begin();
         iElem != source.end(); iElem++ ){
      std::map<std::string, std::set<std::string> >::iterator iT = tableList.find( tableName );
      const std::string& innerTable = iElem->second.tableName();
      if( innerTable != tableName ){
        iT->second.insert( innerTable );
      }
      getTableHierarchyFromMappingElement( iElem->second, tableList );
    }
  }
  
  //void getTableHierarchyFromMapping( const MappingTree& source,
  //                                   std::map<std::string,
  //                                   std::set<std::string> >& tableList ){
  //  getTableHierarchyFromMappingElement( source.element(), tableList );
  //  for( std::map<std::string, MappingElement >::const_iterator iDepEl = source.dependentElements().begin();
  //       iDepEl != source.dependentElements().end(); iDepEl++ ){
  //    getTableHierarchyFromMappingElement( iDepEl->second, tableList );
  //  }
  //}

  void addFromTableHierarchy( const std::string& tableName,
                              std::map<std::string, std::set<std::string> >& tableList,
                              std::vector<std::string>& orderedList ){
    orderedList.push_back( tableName );
    std::map<std::string, std::set<std::string> >::const_iterator iDeps = tableList.find( tableName );
    if(iDeps != tableList.end() ){
      for( std::set<std::string>::const_iterator iDt = iDeps->second.begin();
           iDt != iDeps->second.end(); iDt++ ){
        addFromTableHierarchy( *iDt, tableList, orderedList );
      }
    }
  }

}

ora::ContainerSchema::ContainerSchema( int containerId,
                                       const std::string& containerName,
                                       const Reflex::Type& containerType,
                                       DatabaseSession& session ):
  m_containerId( containerId ),
  m_containerName( containerName ),
  m_className( containerType.Name( Reflex::SCOPED ) ),
  m_classDict( containerType ),
  m_session( session ),
  m_loaded( false ),
  m_containerSchemaSequences( session.schema() ),
  m_mapping(),
  m_dependentMappings(){
}

ora::ContainerSchema::ContainerSchema( int containerId,
                                       const std::string& containerName,
                                       const std::string& className,
                                       DatabaseSession& session ):
  m_containerId( containerId ),
  m_containerName( containerName ),
  m_className( className ),
  m_classDict(),
  m_session( session ),
  m_loaded( false ),
  m_containerSchemaSequences( session.schema() ),
  m_mapping(),
  m_dependentMappings(){
  m_classDict = ClassUtils::lookupDictionary( className, false );
}

ora::ContainerSchema::~ContainerSchema(){
  for( std::map<std::string,MappingTree*>::iterator iDep = m_dependentMappings.begin();
       iDep != m_dependentMappings.end(); ++iDep ){
    delete iDep->second;
  }
}

void ora::ContainerSchema::checkClassDict(){
  if( !m_classDict ) throwException("Container class \""+m_className+"\" has not been found in the dictionary.",
                                    "ContainerSchema::checkClassDict");
}

void ora::ContainerSchema::create(){

  checkClassDict();
  std::string newMappingVersion = m_session.mappingDatabase().newMappingVersionForContainer( m_containerName );
  MappingGenerator mapGen( m_session.schema().storageSchema() );
  if(!mapGen.createNewMapping( m_containerName, m_classDict, m_mapping )){
    throwException("Mapping generation failed.",
                   "ContainerSchema::create");
  }
  m_mapping.setVersion( newMappingVersion );
  m_session.mappingDatabase().storeMapping( m_mapping );
  m_session.mappingDatabase().insertClassVersion( m_classDict, 0, m_containerId, newMappingVersion, true );
  MappingToSchema mapping2Schema( m_session.schema().storageSchema() );
  mapping2Schema.createOrAlter(  m_mapping, false );
  m_loaded = true;
}

void ora::ContainerSchema::drop(){

  std::set<std::string> containerMappingVersions;
  if( m_session.mappingDatabase().getMappingVersionsForContainer( m_containerId, containerMappingVersions )){
    std::map< std::string, std::set<std::string> > tableHierarchy;
    std::set<std::string> topLevelTables; // should be strictly only one!
    for( std::set<std::string>::const_iterator iV = containerMappingVersions.begin();
         iV!= containerMappingVersions.end(); ++iV ){
      MappingTree mapping;
      if( m_session.mappingDatabase().getMappingByVersion( *iV, mapping ) ){
        topLevelTables.insert( mapping.topElement().tableName() );
        getTableHierarchyFromMappingElement( mapping.topElement(), tableHierarchy );
        m_containerSchemaSequences.erase( MappingRules::sequenceNameForDependentClass( m_containerName, mapping.className()));
      }
    }

    std::vector<std::string> orderedTableList;
    for(std::set<std::string>::const_iterator iMainT = topLevelTables.begin();
        iMainT != topLevelTables.end(); ++iMainT ){
      addFromTableHierarchy( *iMainT, tableHierarchy, orderedTableList );
    }
    for(std::vector<std::string>::reverse_iterator iTable = orderedTableList.rbegin();
        iTable != orderedTableList.rend(); iTable++ ){
      m_session.schema().storageSchema().dropIfExistsTable( *iTable );
    }
  }
  for( std::set<std::string>::const_iterator iM = containerMappingVersions.begin();
       iM != containerMappingVersions.end(); ++iM ){
    m_session.mappingDatabase().removeMapping( *iM );
  }
  m_containerSchemaSequences.erase( MappingRules::sequenceNameForContainer( m_containerName ));

}

void ora::ContainerSchema::evolve(){
  // will be always enabled?
  //bool schemaEvolutionEnabled = true;
  //if(!schemaEvolutionEnabled){
  //  throwException("No mapping available for the provided class version","ContainerSchema::create");
  //}
  checkClassDict();
  MappingGenerator mapGen( m_session.schema().storageSchema() );
  // retrieve the base mapping
  MappingTree baseMapping;
  if( !m_session.mappingDatabase().getBaseMappingForContainer( m_classDict.Name(Reflex::SCOPED), m_containerId, baseMapping )){
    throwException("Base mapping has not been found in the database.",
                   "ContainerSchema::evolve");
  }
  if(!mapGen.createNewMapping( m_containerName, m_classDict, baseMapping,  m_mapping )){
    throwException("Mapping generation failed..",
                   "ContainerSchema::evolve");
  }
  std::string newMappingVersion = m_session.mappingDatabase().newMappingVersionForContainer( m_containerName );
  m_mapping.setVersion( newMappingVersion );
  m_session.mappingDatabase().storeMapping( m_mapping );
  m_session.mappingDatabase().insertClassVersion( m_classDict, 0, m_containerId, newMappingVersion );
  MappingToSchema mapping2Schema( m_session.schema().storageSchema() );
  mapping2Schema.createOrAlter(  m_mapping, true  );
  m_loaded = true;
}

const Reflex::Type& ora::ContainerSchema::type(){
  return m_classDict;
}

ora::MappingTree& ora::ContainerSchema::mapping(){
  checkClassDict();
  if(!m_loaded ){
    if( !m_session.mappingDatabase().getMappingForContainer( m_classDict, m_containerId, m_mapping ) ){
      evolve();
    }
    m_loaded = true;
  }
  return m_mapping;
}

bool ora::ContainerSchema::loadMappingForDependentClass( const Reflex::Type& dependentClassDict ){
  checkClassDict();
  std::string className = dependentClassDict.Name(Reflex::SCOPED);
  std::map<std::string,MappingTree*>::iterator iDep = m_dependentMappings.find( className );
  if( iDep ==  m_dependentMappings.end() ){
    // not in cache, search the database...
    iDep = m_dependentMappings.insert( std::make_pair( className, new MappingTree ) ).first;
    if( ! m_session.mappingDatabase().getMappingForContainer( dependentClassDict, m_containerId, *iDep->second ) ){
      m_dependentMappings.erase( className );
      return false;
    }
  }
  return true;  
}

void ora::ContainerSchema::extend( const Reflex::Type& dependentClassDict ){
  bool okGen = false;
  std::string className = dependentClassDict.Name(Reflex::SCOPED);
  std::map<std::string,MappingTree*>::iterator iDep =
    m_dependentMappings.insert( std::make_pair( className, new MappingTree ) ).first;
  MappingGenerator mapGen( m_session.schema().storageSchema() );
  MappingTree baseMapping;
  bool foundBase = m_session.mappingDatabase().getBaseMappingForContainer( className, m_containerId, baseMapping );
  if( foundBase ){
    // evolution required...
    okGen = mapGen.createNewDependentMapping( dependentClassDict, m_mapping, baseMapping, *iDep->second );
  } else {
    // new mapping from scratch...
    okGen = mapGen.createNewDependentMapping( dependentClassDict, m_mapping, *iDep->second );
  }
  if(!okGen) throwException("Mapping generation failed.",
                            "ContainerSchema::extend");
  std::string newMappingVersion = m_session.mappingDatabase().newMappingVersionForContainer( m_containerName );
  iDep->second->setVersion( newMappingVersion );
  m_session.mappingDatabase().storeMapping( *iDep->second );
  //m_session.mappingDatabase().insertClassVersion( m_classDict, 1, m_containerId, newMappingVersion, !foundBase );
  m_session.mappingDatabase().insertClassVersion( dependentClassDict, 1, m_containerId, newMappingVersion, !foundBase );
  MappingToSchema mapping2Schema( m_session.schema().storageSchema() );
  mapping2Schema.createOrAlter(  *iDep->second, foundBase  );
}

bool ora::ContainerSchema::extendIfRequired( const Reflex::Type& dependentClassDict ){
  bool ret = false;
  if( ! loadMappingForDependentClass( dependentClassDict ) ){
    extend( dependentClassDict );
    ret = true;
  }
  return ret;
}

ora::MappingElement& ora::ContainerSchema::mappingForDependentClass( const Reflex::Type& dependentClassDict ){
  if( ! loadMappingForDependentClass( dependentClassDict ) ){
    // new mapping genereation required...
    if( m_session.configuration().properties().getFlag( Configuration::automaticDatabaseCreation()) ||
        m_session.configuration().properties().getFlag( Configuration::automaticContainerCreation() ) ){
      extend( dependentClassDict );
    } 
  }
  std::string className = dependentClassDict.Name(Reflex::SCOPED);
  std::map<std::string,MappingTree*>::iterator iDep = m_dependentMappings.find( className );
  if( iDep ==  m_dependentMappings.end() ){
    throwException( "Mapping for class \""+ className + "\" is not available in the database.",
                    "ContainerSchema::mappingForDependentClass");
  }
  return iDep->second->topElement();
}

bool ora::ContainerSchema::mappingForDependentClasses( std::vector<ora::MappingElement>& destination ){
  return m_session.mappingDatabase().getDependentMappingsForContainer( m_containerId, destination );
}

ora::Sequences& ora::ContainerSchema::containerSequences(){
  return m_containerSchemaSequences;
}

ora::IBlobStreamingService* ora::ContainerSchema::blobStreamingService(){
  return m_session.configuration().blobStreamingService();
}

ora::IReferenceHandler* ora::ContainerSchema::referenceHandler(){
  return m_session.configuration().referenceHandler();
}
  
const std::string& ora::ContainerSchema::mappingVersion(){
  return m_mapping.version();
}

int ora::ContainerSchema::containerId(){
  return m_containerId;
}

const std::string&  ora::ContainerSchema::containerName(){
  return m_containerName;
}

const std::string&  ora::ContainerSchema::className(){
  return m_className;
}

coral::ISchema& ora::ContainerSchema::storageSchema(){
  return m_session.schema().storageSchema();
}


