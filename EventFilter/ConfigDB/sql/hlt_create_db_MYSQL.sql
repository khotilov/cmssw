--
-- CMS High Level Trigger Configuration Database Schema: MySQL
-- -----------------------------------------------------------
--
-- CREATED:
-- 12/12/2006 Philipp Schieferdecker <philipp.schieferdecker@cern.ch>
--

-- create the database
DROP DATABASE IF EXISTS hltdb;
CREATE DATABASE hltdb;
USE hltdb;

-- TABLE 'SoftwareReleases'
CREATE TABLE SoftwareReleases
(
	releaseId  	BIGINT UNSIGNED   NOT NULL AUTO_INCREMENT UNIQUE,
	releaseTag     	VARCHAR(32)       NOT NULL UNIQUE,
	PRIMARY KEY(releaseId)
) ENGINE=INNODB;

-- TABLE 'Directories'
CREATE TABLE Directories
(
	dirId		BIGINT UNSIGNED   NOT NULL AUTO_INCREMENT UNIQUE,
	parentDirId     BIGINT UNSIGNED,
        dirName         VARCHAR(512)      NOT NULL UNIQUE,
	created		TIMESTAMP         NOT NULL,
	PRIMARY KEY(dirId)
) ENGINE=INNODB;

-- TABLE 'Configurations'
CREATE TABLE Configurations
(
	configId   	BIGINT UNSIGNED   NOT NULL AUTO_INCREMENT UNIQUE,
	configDescriptor VARCHAR(256)     NOT NULL UNIQUE,
	parentDirId     BIGINT UNSIGNED   NOT NULL,
	config     	VARCHAR(64)       NOT NULL,
	version         SMALLINT UNSIGNED NOT NULL,
	created         TIMESTAMP         NOT NULL,
	UNIQUE (parentDirId,config,version),
	PRIMARY KEY(configId),
	FOREIGN KEY(parentDirId) REFERENCES Directories(dirId)
) ENGINE=INNODB;

-- TABLE 'ConfigurationReleaseAssoc'
CREATE TABLE ConfigurationReleaseAssoc
(
	configId   	BIGINT UNSIGNED   NOT NULL,
	releaseId  	BIGINT UNSIGNED   NOT NULL,
	FOREIGN KEY(configId) REFERENCES Configurations(configId),
	FOREIGN KEY(releaseId) REFERENCES SoftwareReleases(releaseId)
) ENGINE=INNODB;

-- TABLE 'SuperIds'
CREATE TABLE SuperIds
(
	superId    	BIGINT UNSIGNED   NOT NULL AUTO_INCREMENT UNIQUE,
	PRIMARY KEY(superId)
) ENGINE=INNODB;

-- TABLE 'SuperIdReleaseAssoc'
CREATE TABLE SuperIdReleaseAssoc
(
	superId    	BIGINT UNSIGNED   NOT NULL,
	releaseId  	BIGINT UNSIGNED   NOT NULL,
	FOREIGN KEY(superId) REFERENCES SuperIds(superId),
	FOREIGN KEY(releaseId) REFERENCES SoftwareReleases(releaseId)
) ENGINE=INNODB;

-- TABLE 'Paths'
CREATE TABLE Paths
(
	pathId     	BIGINT UNSIGNED   NOT NULL AUTO_INCREMENT UNIQUE,
	configId   	BIGINT UNSIGNED   NOT NULL,
	name       	VARCHAR(64)       NOT NULL,
	sequenceNb 	SMALLINT UNSIGNED NOT NULL,
	isEndPath       BOOL              NOT NULL DEFAULT false,
	UNIQUE(configId,name),
	PRIMARY KEY(pathId),
	FOREIGN KEY(configId) REFERENCES Configurations(configId)
) ENGINE=INNODB;

-- TABLE 'PathInPathAssoc'
CREATE TABLE PathInPathAssoc
(
	parentPathId	BIGINT UNSIGNED   NOT NULL,
	childPathId	BIGINT UNSIGNED   NOT NULL,
	sequenceNb	SMALLINT UNSIGNED NOT NULL,
	FOREIGN KEY (parentPathId) REFERENCES Paths(pathId),
	FOREIGN KEY (childPathId)  REFERENCES Paths(pathId)
) ENGINE=INNODB;

-- TABLE 'Sequences'
CREATE TABLE Sequences
(
	sequenceId	BIGINT UNSIGNED	  NOT NULL AUTO_INCREMENT UNIQUE,
	configId        BIGINT UNSIGNED   NOT NULL,
	name		VARCHAR(64)	  NOT NULL,
	PRIMARY KEY(sequenceId),
	FOREIGN KEY(configId) REFERENCES Configurations(configId)
) ENGINE=INNODB;

-- TABLE 'PathSequenceAssoc'
CREATE TABLE PathSequenceAssoc
(
	pathId		BIGINT UNSIGNED   NOT NULL,
	sequenceId	BIGINT UNSIGNED   NOT NULL,
	sequenceNb      SMALLINT UNSIGNED NOT NULL,
	UNIQUE(pathId,sequenceId),
	FOREIGN KEY(pathId)     REFERENCES Paths(pathId),
	FOREIGN KEY(sequenceId) REFERENCES Sequences(sequenceId)
) ENGINE=INNODB;


--
-- SERVICES
--

-- TABLE 'ServiceTemplates'
CREATE TABLE ServiceTemplates
(
	superId  	BIGINT UNSIGNED   NOT NULL UNIQUE,
	name       	VARCHAR(64)       NOT NULL,
	cvstag       	VARCHAR(64)       NOT NULL,
	PRIMARY KEY(superId),
	FOREIGN KEY(superId) REFERENCES SuperIds(superId)
) ENGINE=INNODB;

-- TABLE 'Services'
CREATE TABLE Services
(
	superId      	BIGINT UNSIGNED   NOT NULL UNIQUE,
	templateId     	BIGINT UNSIGNED   NOT NULL,
	configId   	BIGINT UNSIGNED   NOT NULL,
	sequenceNb	SMALLINT UNSIGNED NOT NULL,
	PRIMARY KEY(superId),
	FOREIGN KEY(superId)    REFERENCES SuperIds(superId),
	FOREIGN KEY(templateId) REFERENCES ServiceTemplates(superId),
	FOREIGN KEY(configId)   REFERENCES Configurations(configId)
) ENGINE=INNODB;


--
-- EDSOURCES
--

-- TABLE 'EDSourceTemplates'
CREATE TABLE EDSourceTemplates
(
	superId  	BIGINT UNSIGNED   NOT NULL UNIQUE,
	name       	VARCHAR(64)       NOT NULL,
	cvstag       	VARCHAR(64)       NOT NULL,
	PRIMARY KEY(superId),
	FOREIGN KEY(superId) REFERENCES SuperIds(superId)
) ENGINE=INNODB;

-- TABLE 'EDSources'
CREATE TABLE EDSources
(
	superId      	BIGINT UNSIGNED   NOT NULL UNIQUE,
	templateId     	BIGINT UNSIGNED   NOT NULL,
	configId   	BIGINT UNSIGNED   NOT NULL,
	sequenceNb	SMALLINT UNSIGNED NOT NULL,
	PRIMARY KEY(superId),
	FOREIGN KEY(superId)    REFERENCES SuperIds(superId),
	FOREIGN KEY(templateId) REFERENCES EDSourceTemplates(superId),
	FOREIGN KEY(configId)   REFERENCES Configurations(configId)
) ENGINE=INNODB;


--
-- ESSOURCES
--

-- TABLE 'ESSourceTemplates'
CREATE TABLE ESSourceTemplates
(
	superId  	BIGINT UNSIGNED   NOT NULL UNIQUE,
	name       	VARCHAR(64)       NOT NULL,
	cvstag       	VARCHAR(64)       NOT NULL,
	PRIMARY KEY(superId),
	FOREIGN KEY(superId) REFERENCES SuperIds(superId)
) ENGINE=INNODB;

-- TABLE 'ESSources'
CREATE TABLE ESSources
(
	superId      	BIGINT UNSIGNED   NOT NULL UNIQUE,
	templateId     	BIGINT UNSIGNED   NOT NULL,
	configId   	BIGINT UNSIGNED   NOT NULL,
	name       	VARCHAR(64)	  NOT NULL,
	sequenceNb      SMALLINT UNSIGNED NOT NULL,
	PRIMARY KEY(superId),
	FOREIGN KEY(superId)    REFERENCES SuperIds(superId),
	FOREIGN KEY(templateId) REFERENCES ESSourceTemplates(superId),
	FOREIGN KEY(configId)   REFERENCES Configurations(configId)
) ENGINE=INNODB;


--
-- MODULES
--

-- TABLE 'ModuleTypes'
CREATE TABLE ModuleTypes
(
	typeId 		BIGINT UNSIGNED    NOT NULL AUTO_INCREMENT UNIQUE,
	type   		VARCHAR(32)        NOT NULL UNIQUE,
	PRIMARY KEY(typeId)
) ENGINE=INNODB;

-- TABLE 'ModuleTemplates'
CREATE TABLE ModuleTemplates
(
	superId  	BIGINT UNSIGNED   NOT NULL UNIQUE,
	typeId  	BIGINT UNSIGNED   NOT NULL,
	name       	VARCHAR(64)       NOT NULL,
	cvstag       	VARCHAR(64)       NOT NULL,
	PRIMARY KEY(superId),
	FOREIGN KEY(superId) REFERENCES SuperIds(superId),
	FOREIGN KEY(typeId)  REFERENCES ModuleTypes(typeId)
) ENGINE=INNODB;

-- TABLE 'Modules'
CREATE TABLE Modules
(
	superId   	BIGINT UNSIGNED   NOT NULL UNIQUE,
	templateId  	BIGINT UNSIGNED   NOT NULL,
	name       	VARCHAR(64)       NOT NULL,
	PRIMARY KEY(superId),
	FOREIGN KEY(superId) REFERENCES SuperIds(superId),
	FOREIGN KEY(templateId) REFERENCES ModuleTemplates(superId)
) ENGINE=INNODB;

-- TABLE 'PathModuleAssoc'
CREATE TABLE PathModuleAssoc
(
	pathId     	BIGINT UNSIGNED   NOT NULL,
        moduleId   	BIGINT UNSIGNED   NOT NULL,
	sequenceNb	SMALLINT UNSIGNED NOT NULL,
	FOREIGN KEY(pathId)   REFERENCES Paths(pathId),
	FOREIGN KEY(moduleId) REFERENCES Modules(superId)
) ENGINE=INNODB;

-- TABLE 'SequenceModuleAssoc'
CREATE TABLE SequenceModuleAssoc
(
	sequenceId     	BIGINT UNSIGNED   NOT NULL,
        moduleId   	BIGINT UNSIGNED   NOT NULL,
	sequenceNb	SMALLINT UNSIGNED NOT NULL,
	FOREIGN KEY(sequenceId) REFERENCES Sequences(sequenceId),
	FOREIGN KEY(moduleId)   REFERENCES Modules(superId)
) ENGINE=INNODB;


--
-- PARAMETER SETS
--

-- TABLE 'ParameterSets'
CREATE TABLE ParameterSets
(
	superId		BIGINT UNSIGNED	  NOT NULL UNIQUE,
	name		VARCHAR(64)	  NOT NULL,
	tracked         BOOLEAN           NOT NULL,
	PRIMARY KEY(superId),
	FOREIGN KEY(superId) REFERENCES SuperIds(superId)
) ENGINE=INNODB;

--TABLE 'VecParameterSets'
CREATE TABLE VecParameterSets
(
	superId		BIGINT UNSIGNED   NOT NULL UNIQUE,
	name		VARCHAR(64)	  NOT NULL,
	tracked         BOOLEAN           NOT NULL,
	PRIMARY KEY(superId),
	FOREIGN KEY(superId) REFERENCES SuperIds(superId)
) ENGINE=INNODB;

-- TABLE 'ConfigurationParamSetAssoc'
CREATE TABLE ConfigurationParamSetAssoc
(
	configId	BIGINT UNSIGNED	  NOT NULL,
	paramSetId	BIGINT UNSIGNED	  NOT NULL,
	sequenceNb	SMALLINT UNSIGNED NOT NULL,
	FOREIGN KEY(configId)    REFERENCES Configurations(configId),
	FOREIGN KEY(paramSetId) REFERENCES ParameterSets(superId)
) ENGINE=INNODB;

-- TABLE 'SuperIdParamSetAssoc'
CREATE TABLE SuperIdParamSetAssoc
(
	superId		BIGINT UNSIGNED	  NOT NULL,
	paramSetId	BIGINT UNSIGNED	  NOT NULL,
	sequenceNb	SMALLINT UNSIGNED NOT NULL,
	FOREIGN KEY(superId)    REFERENCES SuperIds(superId),
	FOREIGN KEY(paramSetId) REFERENCES ParameterSets(superId)
) ENGINE=INNODB;

-- TABLE 'SuperIdVecParamSetAssoc'
CREATE TABLE SuperIdVecParamSetAssoc
(
	superId		BIGINT UNSIGNED	  NOT NULL,
	vecParamSetId	BIGINT UNSIGNED	  NOT NULL,
	sequenceNb	SMALLINT UNSIGNED NOT NULL,
	FOREIGN KEY(superId)       REFERENCES SuperIds(superId),
	FOREIGN KEY(vecParamSetId) REFERENCES VecParameterSets(superId)
) ENGINE=INNODB;


--
-- PARAMETERS
--

-- TABLE 'ParameterTypes'
CREATE TABLE ParameterTypes
(
	paramTypeId	BIGINT UNSIGNED   NOT NULL AUTO_INCREMENT UNIQUE,
	paramType       VARCHAR(32)       NOT NULL UNIQUE,
	PRIMARY KEY(paramTypeId)
) ENGINE=INNODB;

-- TABLE 'Parameters'
CREATE TABLE Parameters
(
	paramId    	BIGINT UNSIGNED   NOT NULL AUTO_INCREMENT UNIQUE,
	paramTypeId    	BIGINT UNSIGNED   NOT NULL,
	name       	VARCHAR(64)       NOT NULL,
	tracked         BOOLEAN           NOT NULL,
	PRIMARY KEY(paramId),
	FOREIGN KEY(paramTypeId) REFERENCES ParameterTypes(paramTypeId)
) ENGINE=INNODB;

-- TABLE 'SuperIdParameterAssoc'
CREATE TABLE SuperIdParameterAssoc
(
	superId		BIGINT UNSIGNED	  NOT NULL,
	paramId		BIGINT UNSIGNED	  NOT NULL,
	sequenceNb	SMALLINT UNSIGNED NOT NULL,
	FOREIGN KEY(superId) REFERENCES SuperIds(superId),
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'Int32ParamValues'
CREATE TABLE Int32ParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL UNIQUE,
	value      	BIGINT            NOT NULL,
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'VInt32ParamValues'
CREATE TABLE VInt32ParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL,
	sequenceNb 	SMALLINT UNSIGNED NOT NULL,
	value      	BIGINT            NOT NULL,
	UNIQUE(paramId,sequenceNb),
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'UInt32ParamValues'
CREATE TABLE UInt32ParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL UNIQUE,
	value      	BIGINT UNSIGNED   NOT NULL,
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'VUInt32ParamValues'
CREATE TABLE VUInt32ParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL,
	sequenceNb 	SMALLINT UNSIGNED NOT NULL,
	value      	BIGINT UNSIGNED   NOT NULL,
	UNIQUE(paramId,sequenceNb),
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'BoolParamValues'
CREATE TABLE BoolParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL UNIQUE,
	value      	BOOLEAN           NOT NULL,
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'DoubleParamValues'
CREATE TABLE DoubleParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL UNIQUE,
	value      	REAL              NOT NULL,
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'VDoubleParamValues'
CREATE TABLE VDoubleParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL,
	sequenceNb 	SMALLINT UNSIGNED NOT NULL,
	value      	REAL              NOT NULL,
	UNIQUE(paramId,sequenceNb),
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'StringParamValues'
CREATE TABLE StringParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL UNIQUE,
	value      	VARCHAR(256)      NOT NULL,
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'VStringParamValues'
CREATE TABLE VStringParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL,
	sequenceNb 	SMALLINT UNSIGNED NOT NULL,
	value      	VARCHAR(256)      NOT NULL,
	UNIQUE(paramId,sequenceNb),
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'InputTagParamValues'
CREATE TABLE InputTagParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL,
	value      	VARCHAR(64)       NOT NULL,
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'VInputTagParamValues'
CREATE TABLE VInputTagParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL,
	sequenceNb 	SMALLINT UNSIGNED NOT NULL,
	value      	VARCHAR(64)       NOT NULL,
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'EventIDParamValues'
CREATE TABLE EventIDParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL,
	value      	VARCHAR(32)       NOT NULL,
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;

-- TABLE 'VEventIDParamValues'
CREATE TABLE VEventIDParamValues
(
	paramId    	BIGINT UNSIGNED   NOT NULL,
	sequenceNb 	SMALLINT UNSIGNED NOT NULL,
	value      	VARCHAR(32)       NOT NULL,
	FOREIGN KEY(paramId) REFERENCES Parameters(paramId)
) ENGINE=INNODB;


--
-- INSERTs
-- 

-- INSERT root directory
INSERT INTO Directories (parentDirId,dirName,created) VALUES(null,"/",NOW());

-- INSERT valid module types
INSERT INTO ModuleTypes (type) VALUES("EDProducer");
INSERT INTO ModuleTypes (type) VALUES("EDFilter");
INSERT INTO ModuleTypes (type) VALUES("EDAnalyzer");
INSERT INTO ModuleTypes (type) VALUES("HLTProducer");
INSERT INTO ModuleTypes (type) VALUES("HLTFilter");
INSERT INTO ModuleTypes (type) VALUES("ESProducer");
INSERT INTO ModuleTypes (type) VALUES("OutputModule");


-- INSERT valid parameter types
INSERT INTO ParameterTypes (paramType) VALUES("bool");
INSERT INTO ParameterTypes (paramType) VALUES("int32");
INSERT INTO ParameterTypes (paramType) VALUES("vint32");
INSERT INTO ParameterTypes (paramType) VALUES("uint32");
INSERT INTO ParameterTypes (paramType) VALUES("vuint32");
INSERT INTO ParameterTypes (paramType) VALUES("double");
INSERT INTO ParameterTypes (paramType) VALUES("vdouble");
INSERT INTO ParameterTypes (paramType) VALUES("string");
INSERT INTO ParameterTypes (paramType) VALUES("vstring");
INSERT INTO ParameterTypes (paramType) VALUES("InputTag");
INSERT INTO ParameterTypes (paramType) VALUES("VInputTag");
INSERT INTO ParameterTypes (paramType) VALUES("EventID");
INSERT INTO ParameterTypes (paramType) VALUES("VEventID");
