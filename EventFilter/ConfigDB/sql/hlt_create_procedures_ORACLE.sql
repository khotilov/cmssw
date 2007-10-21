--
-- Oracle Stored Procedure
--


--
-- create global temporary tables
--
CREATE GLOBAL TEMPORARY TABLE tmp_template_table
(
  template_id	   NUMBER,
  template_type    VARCHAR2(64),
  template_name    VARCHAR2(128),
  template_cvstag  VARCHAR2(32)
) ON COMMIT PRESERVE ROWS;

CREATE GLOBAL TEMPORARY TABLE tmp_instance_table
(
  instance_id       NUMBER,
  template_id       NUMBER,
  instance_type     VARCHAR2(64),
  instance_name     VARCHAR2(128),
  flag              NUMBER(1)
) ON COMMIT PRESERVE ROWS;

CREATE GLOBAL TEMPORARY TABLE tmp_parameter_table
(
  parameter_id      NUMBER,
  parameter_type    VARCHAR2(64),
  parameter_name    VARCHAR2(128),
  parameter_trkd    NUMBER(1),
  parameter_seqnb   NUMBER,
  parent_id         NUMBER
) ON COMMIT PRESERVE ROWS;

CREATE GLOBAL TEMPORARY TABLE tmp_boolean_table
(
  parameter_id		NUMBER,
  parameter_value   	NUMBER(1)
) ON COMMIT PRESERVE ROWS;

CREATE GLOBAL TEMPORARY TABLE tmp_int_table
(
  parameter_id	  NUMBER,
  parameter_value NUMBER,
  sequence_nb     NUMBER,
  hex		  NUMBER(1)
) ON COMMIT PRESERVE ROWS;

CREATE GLOBAL TEMPORARY TABLE tmp_real_table
(
  parameter_id	  NUMBER,
  parameter_value FLOAT,
  sequence_nb     NUMBER
) ON COMMIT PRESERVE ROWS;

CREATE GLOBAL TEMPORARY TABLE tmp_string_table
(
  parameter_id	    NUMBER,
  parameter_value   VARCHAR2(512),
  sequence_nb       NUMBER
) ON COMMIT PRESERVE ROWS;

CREATE GLOBAL TEMPORARY TABLE tmp_path_entries
(
  path_id           NUMBER,
  entry_id          NUMBER,
  sequence_nb       NUMBER,
  entry_type        VARCHAR2(64)
) ON COMMIT PRESERVE ROWS;

CREATE GLOBAL TEMPORARY TABLE tmp_sequence_entries
(
  sequence_id       NUMBER,
  entry_id          NUMBER,
  sequence_nb       NUMBER,
  entry_type        VARCHAR2(64)
) ON COMMIT PRESERVE ROWS;

CREATE GLOBAL TEMPORARY TABLE tmp_stream_entries
(
  stream_id         NUMBER,
  path_id           NUMBER
) ON COMMIT PRESERVE ROWS;



--
-- TYPE declaration(s)
--
CREATE PACKAGE types
AS
  TYPE ref_cursor IS REF CURSOR;
END;
/


--
-- PROCDEDURE load_parameter_value
--
CREATE PROCEDURE load_parameter_value(parameter_id   NUMBER,
                                      parameter_type CHAR)
AS
  v_bool_value   NUMBER(1);
  v_int32_value  PLS_INTEGER;
  v_int32_hex    NUMBER(1);
  v_double_value FLOAT;
  v_string_value VARCHAR2(512);
  v_sequence_nb  PLS_INTEGER;

  /* declare cursors */
  CURSOR cur_bool IS
    SELECT value FROM BoolParamValues
    WHERE paramId=parameter_id;
  CURSOR cur_int32 IS
    SELECT value,hex FROM Int32ParamValues
    WHERE paramId=parameter_id;
  CURSOR cur_vint32 IS
    SELECT value,sequenceNb,hex FROM VInt32ParamValues
    WHERE paramId=parameter_id
    ORDER BY sequenceNb ASC;
  CURSOR cur_uint32 IS
    SELECT value,hex FROM UInt32ParamValues
    WHERE paramId=parameter_id;
  CURSOR cur_vuint32 IS
    SELECT value,sequenceNb,hex FROM VUInt32ParamValues
    WHERE paramId=parameter_id
    ORDER BY sequenceNb ASC;
  CURSOR cur_double IS
    SELECT value FROM DoubleParamValues
    WHERE paramId=parameter_id;
  CURSOR cur_vdouble IS
    SELECT value,sequenceNb FROM VDoubleParamValues
    WHERE paramId=parameter_id
    ORDER BY sequenceNb ASC;
  CURSOR cur_string IS
    SELECT value FROM StringParamValues
    WHERE paramId=parameter_id;
  CURSOR cur_vstring IS
    SELECT value,sequenceNb FROM VStringParamValues
    WHERE paramId=parameter_id
    ORDER BY sequenceNb ASC;
  CURSOR cur_inputtag IS
    SELECT value FROM InputTagParamValues
    WHERE paramId=parameter_id;
  CURSOR cur_vinputtag IS
    SELECT value,sequenceNb FROM VInputTagParamValues
    WHERE paramId=parameter_id
    ORDER BY sequenceNb ASC;
  CURSOR cur_eventid IS
    SELECT value FROM EventIDParamValues
    WHERE paramId=parameter_id;
  CURSOR cur_veventid IS
    SELECT value,sequenceNb FROM VEventIDParamValues
    WHERE paramId=parameter_id
    ORDER BY sequenceNb ASC;
  CURSOR cur_fileinpath IS
    SELECT value FROM FileInPathParamValues
    WHERE paramId=parameter_id;

BEGIN

  /** load bool values */
  IF parameter_type='bool'
  THEN
    OPEN cur_bool;
    FETCH cur_bool INTO v_bool_value;
    IF cur_bool%FOUND THEN
      INSERT INTO tmp_boolean_table VALUES(parameter_id,v_bool_value);
    END IF;
    CLOSE cur_bool;
  /** load int32 values */
  ELSIF parameter_type='int32'
  THEN
    OPEN cur_int32;
    FETCH cur_int32 INTO v_int32_value,v_int32_hex;
    IF cur_int32%FOUND THEN
      INSERT INTO tmp_int_table
        VALUES(parameter_id,v_int32_value,NULL,v_int32_hex);
    END IF;
    CLOSE cur_int32;
  /** load vint32 values */
  ELSIF parameter_type='vint32'
  THEN
    OPEN cur_vint32;
    LOOP 
      FETCH cur_vint32 INTO v_int32_value,v_sequence_nb,v_int32_hex;
      EXIT WHEN cur_vint32%NOTFOUND;
      INSERT INTO tmp_int_table
        VALUES(parameter_id,v_int32_value,v_sequence_nb,v_int32_hex);
    END LOOP;
    CLOSE cur_vint32;
  /** load uint32 values */
  ELSIF parameter_type='uint32'
  THEN
    OPEN cur_uint32;
    FETCH cur_uint32 INTO v_int32_value,v_int32_hex;
    IF cur_uint32%FOUND THEN
      INSERT INTO tmp_int_table
        VALUES(parameter_id,v_int32_value,NULL,v_int32_hex);
    END IF;
    CLOSE cur_uint32;
  /** load vuint32 values */
  ELSIF parameter_type='vuint32'
  THEN
    OPEN cur_vuint32;
    LOOP
      FETCH cur_vuint32 INTO v_int32_value,v_sequence_nb,v_int32_hex;
      EXIT WHEN cur_vuint32%NOTFOUND;
      INSERT INTO tmp_int_table
        VALUES(parameter_id,v_int32_value,v_sequence_nb,v_int32_hex);
    END LOOP;
    CLOSE cur_vuint32;
  /** load double values */
  ELSIF parameter_type='double'
  THEN
    OPEN cur_double;
    FETCH cur_double INTO v_double_value;
    IF cur_double%FOUND THEN
      INSERT INTO tmp_real_table VALUES(parameter_id,v_double_value,NULL);
    END IF;
    CLOSE cur_double;
  /** load vdouble values */
  ELSIF parameter_type='vdouble'
  THEN
    OPEN cur_vdouble;
    LOOP 
      FETCH cur_vdouble INTO v_double_value,v_sequence_nb;
      EXIT WHEN cur_vdouble%NOTFOUND;
      INSERT INTO tmp_real_table
        VALUES(parameter_id,v_double_value,v_sequence_nb);
    END LOOP;
    CLOSE cur_vdouble;
  /** load string values */
  ELSIF parameter_type='string'
  THEN
    OPEN cur_string;
    FETCH cur_string INTO v_string_value;
    IF cur_string%FOUND THEN
      INSERT INTO tmp_string_table VALUES(parameter_id,v_string_value,NULL);
    END IF;
    CLOSE cur_string;
  /** load vstring values */
  ELSIF parameter_type='vstring'
  THEN
    OPEN cur_vstring;
    LOOP
      FETCH cur_vstring INTO v_string_value,v_sequence_nb;
      EXIT WHEN cur_vstring%NOTFOUND;
      INSERT INTO tmp_string_table
        VALUES(parameter_id,v_string_value,v_sequence_nb);
    END LOOP;
    CLOSE cur_vstring;
  /** load inputtag values */
  ELSIF parameter_type='InputTag'
  THEN
    OPEN cur_inputtag;
    FETCH cur_inputtag INTO v_string_value;
    IF cur_inputtag%FOUND THEN
      INSERT INTO tmp_string_table VALUES(parameter_id,v_string_value,NULL);
    END IF;
    CLOSE cur_inputtag;
  /** load vinputtag values */
  ELSIF parameter_type='VInputTag'
  THEN
    OPEN cur_vinputtag;
    LOOP
      FETCH cur_vinputtag INTO v_string_value,v_sequence_nb;
      EXIT WHEN cur_vinputtag%NOTFOUND;
      INSERT INTO tmp_string_table
        VALUES(parameter_id,v_string_value,v_sequence_nb);
    END LOOP;
    CLOSE cur_vinputtag;
  /** load eventid values */
  ELSIF parameter_type='EventID'
  THEN
    OPEN cur_eventid;
    FETCH cur_eventid INTO v_string_value;
    IF cur_eventid%FOUND THEN
      INSERT INTO tmp_string_table VALUES(parameter_id,v_string_value,NULL);
    END IF;
    CLOSE cur_eventid;
  /** load veventid values */
  ELSIF parameter_type='VEventID'
  THEN
    OPEN cur_veventid;
    LOOP
      FETCH cur_veventid INTO v_string_value,v_sequence_nb;
      EXIT WHEN cur_veventid%NOTFOUND;
      INSERT INTO tmp_string_table
        VALUES(parameter_id,v_string_value,v_sequence_nb);
    END LOOP;
    CLOSE cur_veventid;
  /** load fileinpath values */
  ELSIF parameter_type='FileInPath'
  THEN
    OPEN cur_fileinpath;
    FETCH cur_fileinpath INTO v_string_value;
    IF cur_fileinpath%FOUND THEN
      INSERT INTO tmp_string_table VALUES(parameter_id,v_string_value,NULL);
    END IF;
    CLOSE cur_fileinpath;
  END IF;
END;
/


--
-- PROCEDURE load_parameters
--
CREATE PROCEDURE load_parameters(parent_id IN NUMBER)
AS
  v_parameter_id    PLS_INTEGER;
  v_parameter_type  VARCHAR2(64);
  v_parameter_name  VARCHAR2(128);
  v_parameter_trkd  NUMBER(1);
  v_parameter_seqnb PLS_INTEGER;

  /* cursor for parameters */
  CURSOR cur_parameters IS
    SELECT Parameters.paramId,
           ParameterTypes.paramType,
           Parameters.name,
           Parameters.tracked,
           SuperIdParameterAssoc.sequenceNb
    FROM Parameters
    JOIN ParameterTypes
    ON ParameterTypes.paramTypeId = Parameters.paramTypeId
    JOIN SuperIdParameterAssoc
    ON SuperIdParameterAssoc.paramId = Parameters.paramId
    WHERE SuperIdParameterAssoc.superId = parent_id
    ORDER BY SuperIdParameterAssoc.sequenceNb ASC;

  /* cursor for psets */
  CURSOR cur_psets IS
    SELECT ParameterSets.superId,
           ParameterSets.name,
           ParameterSets.tracked,
           SuperIdParamSetAssoc.sequenceNb
    FROM ParameterSets
    JOIN SuperIdParamSetAssoc
    ON SuperIdParamSetAssoc.psetId = ParameterSets.superId
    WHERE SuperIdParamSetAssoc.superId = parent_id
    ORDER BY SuperIdParamSetAssoc.sequenceNb ASC;

  /* cursor for vpsets */
  CURSOR cur_vpsets IS
    SELECT VecParameterSets.superId,
           VecParameterSets.name,
           VecParameterSets.tracked,
           SuperIdVecParamSetAssoc.sequenceNb
    FROM VecParameterSets
    JOIN SuperIdVecParamSetAssoc
    ON SuperIdVecParamSetAssoc.vpsetId = VecParameterSets.superId
    WHERE SuperIdVecParamSetAssoc.superId = parent_id
    ORDER BY SuperIdVecParamSetAssoc.sequenceNb ASC;

BEGIN

  /* load the parameters and fill them into temporary param table */
  OPEN cur_parameters;
  LOOP
    FETCH cur_parameters
      INTO v_parameter_id,v_parameter_type,
           v_parameter_name,v_parameter_trkd,v_parameter_seqnb;
    EXIT WHEN cur_parameters%NOTFOUND;
    INSERT INTO tmp_parameter_table
      VALUES(v_parameter_id,v_parameter_type,
             v_parameter_name,v_parameter_trkd,v_parameter_seqnb,parent_id);
    load_parameter_value(v_parameter_id,v_parameter_type);
  END LOOP;
  CLOSE cur_parameters;

  /* load psets and fill them into temporary param table */
  OPEN cur_psets;
  LOOP
    FETCH cur_psets
      INTO v_parameter_id,v_parameter_name,v_parameter_trkd,v_parameter_seqnb;
    EXIT WHEN cur_psets%NOTFOUND;
    INSERT INTO tmp_parameter_table
      VALUES(v_parameter_id,'PSet',
             v_parameter_name,v_parameter_trkd,v_parameter_seqnb,parent_id);
    load_parameters(v_parameter_id);
  END LOOP;
  CLOSE cur_psets;

  /* load vpsets and fill them into temporary param table */
  OPEN cur_vpsets;
  LOOP
    FETCH cur_vpsets
      INTO v_parameter_id,
           v_parameter_name,v_parameter_trkd,v_parameter_seqnb;
    EXIT WHEN cur_vpsets%NOTFOUND;
    INSERT INTO tmp_parameter_table
      VALUES(v_parameter_id,'VPSet',
             v_parameter_name,v_parameter_trkd,v_parameter_seqnb,parent_id);
    load_parameters(v_parameter_id);
  END LOOP;
  CLOSE cur_vpsets;
END;
/


--
-- FUNCTION load_template
--
CREATE FUNCTION load_template(release_id IN NUMBER,
                              template_name IN CHAR)
  RETURN types.ref_cursor
AS
  template_cursor   types.ref_cursor;
  v_template_id     PLS_INTEGER;
  v_template_type   VARCHAR2(64);
  v_template_name   VARCHAR2(128);
  v_template_cvstag VARCHAR2(32);
  template_found    BOOLEAN := FALSE;

  /* cursor for edsource templates */
  CURSOR cur_edsource_templates IS
    SELECT EDSourceTemplates.superId,
           EDSourceTemplates.name,
           EDSourceTemplates.cvstag
    FROM EDSourceTemplates
    JOIN SuperIdReleaseAssoc
    ON EDSourceTemplates.superId = SuperIdReleaseAssoc.superId
    JOIN SoftwareReleases
    ON SoftwareReleases.releaseId = SuperIdReleaseAssoc.releaseId
    WHERE SoftwareReleases.releaseId = release_id
    AND   EDSourceTemplates.name = template_name;

  /* cursor for essource templates */
  CURSOR cur_essource_templates IS
    SELECT ESSourceTemplates.superId,
           ESSourceTemplates.name,
           ESSourceTemplates.cvstag
    FROM ESSourceTemplates
    JOIN SuperIdReleaseAssoc
    ON ESSourceTemplates.superId = SuperIdReleaseAssoc.superId
    JOIN SoftwareReleases
    ON SoftwareReleases.releaseId = SuperIdReleaseAssoc.releaseId
    WHERE SoftwareReleases.releaseId = release_id
    AND   ESSourceTemplates.name = template_name;

  /* cursor for esmodule templates */
  CURSOR cur_esmodule_templates IS
    SELECT ESModuleTemplates.superId,
           ESModuleTemplates.name,
           ESModuleTemplates.cvstag
    FROM ESModuleTemplates
    JOIN SuperIdReleaseAssoc
    ON ESModuleTemplates.superId = SuperIdReleaseAssoc.superId
    JOIN SoftwareReleases
    ON SoftwareReleases.releaseId = SuperIdReleaseAssoc.releaseId
    WHERE SuperIdReleaseAssoc.releaseId = release_id
    AND   ESModuleTemplates.name = template_name;

  /* cursor for service templates */
  CURSOR cur_service_templates IS
    SELECT ServiceTemplates.superId,
           ServiceTemplates.name,
           ServiceTemplates.cvstag
    FROM ServiceTemplates
    JOIN SuperIdReleaseAssoc
    ON ServiceTemplates.superId = SuperIdReleaseAssoc.superId
    JOIN SoftwareReleases
    ON SoftwareReleases.releaseId = SuperIdReleaseAssoc.releaseId
    WHERE SuperIdReleaseAssoc.releaseId = release_id
    AND   ServiceTemplates.name = template_name;

  /* cursor for module templates */
  CURSOR cur_module_templates IS
    SELECT ModuleTemplates.superId,
           ModuleTemplates.name,
           ModuleTemplates.cvstag,
           ModuleTypes.type
    FROM ModuleTemplates
    JOIN ModuleTypes
    ON   ModuleTemplates.typeId = ModuleTypes.typeId
    JOIN SuperIdReleaseAssoc
    ON ModuleTemplates.superId = SuperIdReleaseAssoc.superId
    JOIN SoftwareReleases
    ON SoftwareReleases.releaseId = SuperIdReleaseAssoc.releaseId
    WHERE SuperIdReleaseAssoc.releaseId = release_id
    AND   ModuleTemplates.name = template_name;

BEGIN

  /* prepare temporary tables */
  execute immediate 'DELETE FROM tmp_template_table';
  execute immediate 'DELETE FROM tmp_parameter_table';
  execute immediate 'DELETE FROM tmp_boolean_table';
  execute immediate 'DELETE FROM tmp_int_table';
  execute immediate 'DELETE FROM tmp_real_table';
  execute immediate 'DELETE FROM tmp_string_table';

  /* load edsource templates */
  OPEN cur_edsource_templates;
  LOOP
    FETCH cur_edsource_templates
      INTO v_template_id,v_template_name,v_template_cvstag;
    EXIT WHEN cur_edsource_templates%NOTFOUND;
    INSERT INTO tmp_template_table
      VALUES(v_template_id,'EDSource',v_template_name,v_template_cvstag);
    load_parameters(v_template_id);
    template_found := TRUE;
  END LOOP;
  CLOSE cur_edsource_templates;


  /* load essource templates */
  IF template_found=FALSE THEN
    OPEN cur_essource_templates;
    LOOP
      FETCH cur_essource_templates
        INTO v_template_id,v_template_name,v_template_cvstag;
      EXIT WHEN cur_essource_templates%NOTFOUND;
      INSERT INTO tmp_template_table
        VALUES(v_template_id,'ESSource',v_template_name,v_template_cvstag);
      load_parameters(v_template_id);
      template_found := TRUE;
    END LOOP;
    CLOSE cur_essource_templates;
  END IF;

  /* load esmodule templates */
  IF template_found=FALSE THEN
    OPEN cur_esmodule_templates;
    LOOP
      FETCH cur_esmodule_templates
        INTO v_template_id,v_template_name,v_template_cvstag;
      EXIT WHEN cur_esmodule_templates%NOTFOUND;
      INSERT INTO tmp_template_table
         VALUES(v_template_id,'ESModule',v_template_name,v_template_cvstag);
      load_parameters(v_template_id);
      template_found := TRUE;
    END LOOP;
    CLOSE cur_esmodule_templates;
  END IF;

  /* load service templates */
  IF template_found=FALSE THEN
    OPEN cur_service_templates;
    LOOP
      FETCH cur_service_templates
        INTO v_template_id,v_template_name,v_template_cvstag;
      EXIT WHEN cur_service_templates%NOTFOUND;
      INSERT INTO tmp_template_table
        VALUES(v_template_id,'Service',v_template_name,v_template_cvstag);
      load_parameters(v_template_id);
      template_found := TRUE;
    END LOOP;
    CLOSE cur_service_templates;
  END IF;

  /* load module templates */
  IF template_found=FALSE THEN
    OPEN cur_module_templates;
    LOOP
      FETCH cur_module_templates
        INTO v_template_id,v_template_name,v_template_cvstag,v_template_type;
      EXIT WHEN cur_module_templates%NOTFOUND;
      INSERT INTO tmp_template_table
        VALUES(v_template_id,
               v_template_type,v_template_name,v_template_cvstag);
      load_parameters(v_template_id);
      template_found := TRUE;
    END LOOP;
    CLOSE cur_module_templates;
  END IF;
   
  /* generate the final result set by selecting the temporary table */
  OPEN template_cursor FOR
    SELECT template_id,template_type,template_name,template_cvstag
    FROM tmp_template_table;
  RETURN template_cursor;
END;  
/



--
-- FUNCTION load_templates
--
--CREATE FUNCTION load_templates(release_id IN NUMBER)
--  RETURN types.ref_cursor
CREATE PROCEDURE load_templates(release_id IN NUMBER)
AS
--  template_cursor   types.ref_cursor;
  v_template_id     PLS_INTEGER;
  v_template_type   VARCHAR2(64);
  v_template_name   VARCHAR2(128);
  v_template_cvstag VARCHAR2(32);

  /* cursor for edsource templates */
  CURSOR cur_edsource_templates IS
    SELECT EDSourceTemplates.superId,
           EDSourceTemplates.name,
           EDSourceTemplates.cvstag
    FROM EDSourceTemplates
    JOIN SuperIdReleaseAssoc
    ON EDSourceTemplates.superId = SuperIdReleaseAssoc.superId
    WHERE SuperIdReleaseAssoc.releaseId = release_id;

  /* cursor for essource templates */
  CURSOR cur_essource_templates IS
    SELECT ESSourceTemplates.superId,
           ESSourceTemplates.name,
           ESSourceTemplates.cvstag
    FROM ESSourceTemplates
    JOIN SuperIdReleaseAssoc
    ON ESSourceTemplates.superId = SuperIdReleaseAssoc.superId
    WHERE SuperIdReleaseAssoc.releaseId = release_id;

  /* cursor for esmodule templates */
  CURSOR cur_esmodule_templates IS
    SELECT ESModuleTemplates.superId,
           ESModuleTemplates.name,
           ESModuleTemplates.cvstag
    FROM ESModuleTemplates
    JOIN SuperIdReleaseAssoc
    ON ESModuleTemplates.superId = SuperIdReleaseAssoc.superId
    WHERE SuperIdReleaseAssoc.releaseId = release_id;

  /* cursor for service templates */
  CURSOR cur_service_templates IS
    SELECT ServiceTemplates.superId,
           ServiceTemplates.name,
           ServiceTemplates.cvstag
    FROM ServiceTemplates
    JOIN SuperIdReleaseAssoc
    ON ServiceTemplates.superId = SuperIdReleaseAssoc.superId
    WHERE SuperIdReleaseAssoc.releaseId = release_id;

  /* cursor for module templates */
  CURSOR cur_module_templates IS
    SELECT ModuleTemplates.superId,
           ModuleTemplates.name,
           ModuleTemplates.cvstag,
           ModuleTypes.type
    FROM ModuleTemplates
    JOIN ModuleTypes
    ON   ModuleTemplates.typeId = ModuleTypes.typeId
    JOIN SuperIdReleaseAssoc
    ON ModuleTemplates.superId = SuperIdReleaseAssoc.superId
    WHERE SuperIdReleaseAssoc.releaseId = release_id;

BEGIN
  /* prepare temporary tables */
  execute immediate 'DELETE FROM tmp_template_table';
  execute immediate 'DELETE FROM tmp_parameter_table';
  execute immediate 'DELETE FROM tmp_boolean_table';
  execute immediate 'DELETE FROM tmp_int_table';
  execute immediate 'DELETE FROM tmp_real_table';
  execute immediate 'DELETE FROM tmp_string_table';
  
  /* load edsource templates */
  OPEN cur_edsource_templates;
  LOOP
    FETCH cur_edsource_templates
      INTO v_template_id,v_template_name,v_template_cvstag;
    EXIT WHEN cur_edsource_templates%NOTFOUND;
    INSERT INTO tmp_template_table
      VALUES(v_template_id,'EDSource',v_template_name,v_template_cvstag);
    load_parameters(v_template_id);   
  END LOOP;
  CLOSE cur_edsource_templates;

  /* load essource templates */
  OPEN cur_essource_templates;
  LOOP
    FETCH cur_essource_templates
      INTO v_template_id,v_template_name,v_template_cvstag;
    EXIT WHEN cur_essource_templates%NOTFOUND;
    INSERT INTO tmp_template_table
      VALUES(v_template_id,'ESSource',v_template_name,v_template_cvstag);
    load_parameters(v_template_id);
  END LOOP;
  CLOSE cur_essource_templates;

  /* load esmodule templates */
  OPEN cur_esmodule_templates;
  LOOP
    FETCH cur_esmodule_templates
      INTO v_template_id,v_template_name,v_template_cvstag;
    EXIT WHEN cur_esmodule_templates%NOTFOUND;
    INSERT INTO tmp_template_table
      VALUES(v_template_id,'ESModule',v_template_name,v_template_cvstag);
    load_parameters(v_template_id);
  END LOOP;
  CLOSE cur_esmodule_templates;

  /* load service templates */
  OPEN cur_service_templates;
  LOOP 
    FETCH cur_service_templates
      INTO v_template_id,v_template_name,v_template_cvstag;
    EXIT WHEN cur_service_templates%NOTFOUND;
    INSERT INTO tmp_template_table
      VALUES(v_template_id,'Service',v_template_name,v_template_cvstag);
    load_parameters(v_template_id);
  END LOOP;
  CLOSE cur_service_templates;

  /* load module templates */
  OPEN cur_module_templates;
  LOOP
    FETCH cur_module_templates
      INTO v_template_id,v_template_name,
           v_template_cvstag,v_template_type;
    EXIT WHEN cur_module_templates%NOTFOUND;
    INSERT INTO tmp_template_table
      VALUES(v_template_id,v_template_type,
             v_template_name,v_template_cvstag);
    load_parameters(v_template_id);
  END LOOP;
  CLOSE cur_module_templates;

  /* generate the final result set by selecting the temporary table */
--  OPEN template_cursor FOR
--    SELECT template_id,template_type,template_name,template_cvstag
--    FROM tmp_template_table;
--  RETURN template_cursor;
END;  
/


--
-- FUNCTION load_templates_for_config
--
CREATE FUNCTION load_templates_for_config(config_id IN NUMBER)
  RETURN types.ref_cursor
AS
  template_cursor   types.ref_cursor;
  v_template_id     PLS_INTEGER;
  v_template_type   VARCHAR2(64);
  v_template_name   VARCHAR2(128);
  v_template_cvstag VARCHAR2(32);

  /* cursor for edsource templates */
  CURSOR cur_edsource_templates IS
    SELECT DISTINCT
      EDSourceTemplates.superId,
      EDSourceTemplates.name,
      EDSourceTemplates.cvstag
    FROM EDSourceTemplates
    JOIN EDSources 
    ON EDSources.templateId = EDSourceTemplates.superId
    JOIN ConfigurationEDSourceAssoc
    ON EDSources.superId = ConfigurationEDSourceAssoc.edsourceId
    WHERE ConfigurationEDSourceAssoc.configId = config_id;

  /* cursor for essource templates */
  CURSOR cur_essource_templates IS
    SELECT DISTINCT
      ESSourceTemplates.superId,
      ESSourceTemplates.name,
      ESSourceTemplates.cvstag
    FROM ESSourceTemplates
    JOIN ESSources
    ON ESSources.templateId = ESSourceTemplates.superId
    JOIN ConfigurationESSourceAssoc
    ON ESSources.superId = ConfigurationESSourceAssoc.essourceId
    WHERE ConfigurationESSourceAssoc.configId = config_id;

  /* cursor for esmodule templates */
  CURSOR cur_esmodule_templates IS
    SELECT DISTINCT
      ESModuleTemplates.superId,
      ESModuleTemplates.name,
      ESModuleTemplates.cvstag
    FROM ESModuleTemplates
    JOIN ESModules
    ON ESModules.templateId = ESModuleTemplates.superId
    JOIN ConfigurationESModuleAssoc
    ON ESModules.superId = ConfigurationESModuleAssoc.esmoduleId
    WHERE ConfigurationESModuleAssoc.configId = config_id;

  /* cursor for service templates */
  CURSOR cur_service_templates IS
    SELECT DISTINCT
      ServiceTemplates.superId,
      ServiceTemplates.name,
      ServiceTemplates.cvstag
    FROM ServiceTemplates
    JOIN Services
    ON   Services.templateId = ServiceTemplates.superId
    JOIN ConfigurationServiceAssoc
    ON   Services.superId = ConfigurationServiceAssoc.serviceId
    WHERE ConfigurationServiceAssoc.configId = config_id;

  /* cursor for module templates from configuration *paths* */
  CURSOR cur_module_templates_paths IS
    SELECT DISTINCT
      ModuleTemplates.superId,
      ModuleTemplates.name,
      ModuleTemplates.cvstag,
      ModuleTypes.type
    FROM ModuleTemplates
    JOIN ModuleTypes
    ON   ModuleTemplates.typeId = ModuleTypes.typeId
    JOIN Modules
    ON   Modules.templateId = ModuleTemplates.superId
    JOIN PathModuleAssoc
    ON   PathModuleAssoc.moduleId = Modules.superId
    JOIN ConfigurationPathAssoc
    ON   PathModuleAssoc.pathId = ConfigurationPathAssoc.pathId
    WHERE ConfigurationPathAssoc.configId = config_id;

  /* cursor for module templates from configuration *sequences* */
  CURSOR cur_module_templates_sequences IS
    SELECT DISTINCT
      ModuleTemplates.superId,
      ModuleTemplates.name,
      ModuleTemplates.cvstag,
      ModuleTypes.type
    FROM ModuleTemplates
    JOIN ModuleTypes
    ON ModuleTemplates.typeId = ModuleTypes.typeId
    JOIN Modules
    ON Modules.templateId = ModuleTemplates.superId
    JOIN SequenceModuleAssoc
    ON SequenceModuleAssoc.moduleId = Modules.superId
    JOIN ConfigurationSequenceAssoc
    ON SequenceModuleAssoc.sequenceId=ConfigurationSequenceAssoc.sequenceId
    WHERE ConfigurationSequenceAssoc.configId = config_id;

BEGIN

  /* prepare temporary tables */
  execute immediate 'DELETE FROM tmp_template_table';
  execute immediate 'DELETE FROM tmp_parameter_table';
  execute immediate 'DELETE FROM tmp_boolean_table';
  execute immediate 'DELETE FROM tmp_int_table';
  execute immediate 'DELETE FROM tmp_real_table';
  execute immediate 'DELETE FROM tmp_string_table';

  /* load edsource templates */
  OPEN cur_edsource_templates;
  LOOP
    FETCH cur_edsource_templates
      INTO v_template_id,v_template_name,v_template_cvstag;
    EXIT WHEN cur_edsource_templates%NOTFOUND;
    INSERT INTO tmp_template_table
      VALUES(v_template_id,'EDSource',v_template_name,v_template_cvstag);
    load_parameters(v_template_id);
  END LOOP;
  CLOSE cur_edsource_templates;

  /* load essource templates */
  OPEN cur_essource_templates;
  LOOP
    FETCH cur_essource_templates
      INTO v_template_id,v_template_name,v_template_cvstag;
    EXIT WHEN cur_essource_templates%NOTFOUND;
    INSERT INTO tmp_template_table
      VALUES(v_template_id,'ESSource',v_template_name,v_template_cvstag);
    load_parameters(v_template_id);
  END LOOP;
  CLOSE cur_essource_templates;

  /* load esmodule templates */
  OPEN cur_esmodule_templates;
  LOOP
    FETCH cur_esmodule_templates
      INTO v_template_id,v_template_name,v_template_cvstag;
    EXIT WHEN cur_esmodule_templates%NOTFOUND;
    INSERT INTO tmp_template_table
      VALUES(v_template_id,'ESModule',v_template_name,v_template_cvstag);
    load_parameters(v_template_id);
  END LOOP;
  CLOSE cur_esmodule_templates;

  /* load service templates */
  OPEN cur_service_templates;
  LOOP
    FETCH cur_service_templates
      INTO v_template_id,v_template_name,v_template_cvstag;
    EXIT WHEN cur_service_templates%NOTFOUND;
    INSERT INTO tmp_template_table
      VALUES(v_template_id,'Service',v_template_name,v_template_cvstag);
    load_parameters(v_template_id);
  END LOOP;
  CLOSE cur_service_templates;

  /* load module templates from *paths* */
  OPEN cur_module_templates_paths;
  LOOP
    FETCH cur_module_templates_paths
      INTO v_template_id,v_template_name,
           v_template_cvstag,v_template_type;
    EXIT WHEN cur_module_templates_paths%NOTFOUND;
    INSERT INTO tmp_template_table
      VALUES(v_template_id,v_template_type,
             v_template_name,v_template_cvstag);
    load_parameters(v_template_id);
  END LOOP;
  CLOSE cur_module_templates_paths;

  /* load module templates from *sequences* */
  OPEN cur_module_templates_sequences;
  LOOP
    FETCH cur_module_templates_sequences
      INTO v_template_id,v_template_name,
           v_template_cvstag,v_template_type;
    EXIT WHEN cur_module_templates_sequences%NOTFOUND;
    INSERT INTO tmp_template_table
      VALUES(v_template_id,v_template_type,
             v_template_name,v_template_cvstag);
    load_parameters(v_template_id);
  END LOOP;
  CLOSE cur_module_templates_sequences;

  /* generate the final result set by selecting the temporary table */
  OPEN template_cursor FOR
    SELECT template_id,template_type,template_name,template_cvstag
    FROM tmp_template_table;
  RETURN template_cursor;
END;  
/


--
-- FUNCTION load_configuration
--
CREATE FUNCTION load_configuration(config_id IN NUMBER)
  RETURN types.ref_cursor
AS
  instance_cursor   types.ref_cursor;
  v_instance_id     PLS_INTEGER;
  v_template_id     PLS_INTEGER;
  v_instance_type   VARCHAR2(64);
  v_instance_name   VARCHAR2(128);
  v_pset_is_trkd    NUMBER(1);
  v_endpath         NUMBER(1);
  v_prefer          NUMBER(1);
  v_parent_id       PLS_INTEGER;
  v_child_id        PLS_INTEGER;
  v_sequence_nb     PLS_INTEGER;

  /* cursor for global psets */
  CURSOR cur_global_psets IS
    SELECT
      ParameterSets.superId,
      ParameterSets.name,
      ParameterSets.tracked
    FROM ParameterSets
    JOIN ConfigurationParamSetAssoc
    ON ParameterSets.superId = ConfigurationParamSetAssoc.psetId
    WHERE ConfigurationParamSetAssoc.configId = config_id
    ORDER BY ConfigurationParamSetAssoc.sequenceNb;

  /* cursor for edsources */
  CURSOR cur_edsources IS
    SELECT
      EDSources.superId,
      EDSources.templateId
    FROM EDSources
    JOIN ConfigurationEDSourceAssoc
    ON EDSources.superId = ConfigurationEDSourceAssoc.edsourceId
    WHERE ConfigurationEDSourceAssoc.configId = config_id
    ORDER BY ConfigurationEDSourceAssoc.sequenceNb ASC;

  /* cursor for essources */
  CURSOR cur_essources IS
    SELECT
      ESSources.superId,
      ESSources.templateId,
      ESSources.name,
      ConfigurationESSourceAssoc.prefer
    FROM ESSources
    JOIN ConfigurationESSourceAssoc
    ON ESSources.superId = ConfigurationESSourceAssoc.essourceId
    WHERE ConfigurationESSourceAssoc.configId = config_id
    ORDER BY ConfigurationESSourceAssoc.sequenceNb;

  /* cursor for esmodules */
  CURSOR cur_esmodules IS
    SELECT
      ESModules.superId,
      ESModules.templateId,
      ESModules.name,
      ConfigurationESModuleAssoc.prefer
    FROM ESModules
    JOIN ConfigurationESModuleAssoc
    ON ESModules.superId = ConfigurationESModuleAssoc.esmoduleId
    WHERE ConfigurationESModuleAssoc.configId = config_id
    ORDER BY ConfigurationESModuleAssoc.sequenceNb;

  /* cursor for services */
  CURSOR cur_services IS
    SELECT
      Services.superId,
      Services.templateId
    FROM Services
    JOIN ConfigurationServiceAssoc
    ON   Services.superId = ConfigurationServiceAssoc.serviceId
    WHERE ConfigurationServiceAssoc.configId = config_id
    ORDER BY ConfigurationServiceAssoc.sequenceNb;

  /* cursor for modules from configuration *paths* */
  CURSOR cur_modules_from_paths IS
    SELECT
      Modules.superId,
      Modules.templateId,
      Modules.name
    FROM Modules
    JOIN PathModuleAssoc
    ON   PathModuleAssoc.moduleId = Modules.superId
    JOIN ConfigurationPathAssoc
    ON   PathModuleAssoc.pathId = ConfigurationPathAssoc.pathId
    WHERE ConfigurationPathAssoc.configId = config_id;

  /* cursor for modules from configuration *sequences* */
  CURSOR cur_modules_from_sequences IS
    SELECT
      Modules.superId,
      Modules.templateId,
      Modules.name
    FROM Modules
    JOIN SequenceModuleAssoc
    ON SequenceModuleAssoc.moduleId = Modules.superId
    JOIN ConfigurationSequenceAssoc
    ON SequenceModuleAssoc.sequenceId=ConfigurationSequenceAssoc.sequenceId
    WHERE ConfigurationSequenceAssoc.configId = config_id;

  /* cursor for paths */
  CURSOR cur_paths IS
    SELECT
      Paths.pathId,
      Paths.name,
      Paths.isEndPath
    FROM Paths
    JOIN ConfigurationPathAssoc
    ON Paths.pathId = ConfigurationPathAssoc.pathId
    WHERE ConfigurationPathAssoc.configId = config_id
    ORDER BY ConfigurationPathAssoc.sequenceNb ASC;

  /* cursor for sequences */
  CURSOR cur_sequences IS
    SELECT
      Sequences.sequenceId,
      Sequences.name
    FROM Sequences
    JOIN ConfigurationSequenceAssoc
    ON Sequences.sequenceId = ConfigurationSequenceAssoc.sequenceId
    WHERE ConfigurationSequenceAssoc.configId = config_id
    ORDER BY ConfigurationSequenceAssoc.sequenceNb ASC;

  /* cursor for streams */
  CURSOR cur_streams IS
    SELECT
      Streams.streamId,
      Streams.streamLabel
    FROM Streams
    WHERE Streams.configId = config_id;

  /* cursor for path-path associations */
  CURSOR cur_path_path IS
    SELECT
      PathInPathAssoc.parentPathId,
      PathInPathAssoc.childPathId,
      PathInPathAssoc.sequenceNb
    FROM PathInPathAssoc
    JOIN ConfigurationPathAssoc
    ON PathInPathAssoc.parentPathId = ConfigurationPathAssoc.pathId
    WHERE ConfigurationPathAssoc.configId = config_id;

  /* cursor for path-sequence associations */
  CURSOR cur_path_sequence IS
    SELECT
      PathSequenceAssoc.pathId,
      PathSequenceAssoc.sequenceId,
      PathSequenceAssoc.sequenceNb
    FROM PathSequenceAssoc
    JOIN ConfigurationPathAssoc
    ON PathSequenceAssoc.pathId = ConfigurationPathAssoc.pathId
    WHERE ConfigurationPathAssoc.configId = config_id;

  /* cursor for path-module associations */
  CURSOR cur_path_module IS
    SELECT
      PathModuleAssoc.pathId,
      PathModuleAssoc.moduleId,
      PathModuleAssoc.sequenceNb
    FROM PathModuleAssoc
    JOIN ConfigurationPathAssoc
    ON PathModuleAssoc.pathId = ConfigurationPathAssoc.pathId
    WHERE ConfigurationPathAssoc.configId = config_id;

  /* cursor for sequence-sequence associations */
  CURSOR cur_sequence_sequence IS
    SELECT
      SequenceInSequenceAssoc.parentSequenceId,
      SequenceInSequenceAssoc.childSequenceId,
      SequenceInSequenceAssoc.sequenceNb
    FROM SequenceInSequenceAssoc
    JOIN ConfigurationSequenceAssoc
    ON SequenceInSequenceAssoc.parentSequenceId =
       ConfigurationSequenceAssoc.sequenceId
    WHERE ConfigurationSequenceAssoc.configId = config_id;

  /* cursor for sequence-module associations */
  CURSOR cur_sequence_module IS
    SELECT
      SequenceModuleAssoc.sequenceId,
      SequenceModuleAssoc.moduleId,
      SequenceModuleAssoc.sequenceNb
    FROM SequenceModuleAssoc
    JOIN ConfigurationSequenceAssoc
    ON SequenceModuleAssoc.sequenceId =
       ConfigurationSequenceAssoc.sequenceId
    WHERE ConfigurationSequenceAssoc.configId = config_id;

  /* cursor for stream-path associations */
  CURSOR cur_stream_path IS
    SELECT
      StreamPathAssoc.streamId,
      StreamPathAssoc.pathId
    FROM StreamPathAssoc
    JOIN Streams
    ON StreamPathAssoc.streamId = Streams.streamId
    JOIN Paths
    ON StreamPathAssoc.pathId = Paths.pathId
    WHERE Streams.configId = config_id;

BEGIN

  execute immediate 'DELETE FROM tmp_instance_table';
  execute immediate 'DELETE FROM tmp_parameter_table';
  execute immediate 'DELETE FROM tmp_boolean_table';
  execute immediate 'DELETE FROM tmp_int_table';
  execute immediate 'DELETE FROM tmp_real_table';
  execute immediate 'DELETE FROM tmp_string_table';
  execute immediate 'DELETE FROM tmp_path_entries';
  execute immediate 'DELETE FROM tmp_sequence_entries';
  execute immediate 'DELETE FROM tmp_stream_entries';

  /* load global psets */
  OPEN cur_global_psets;
  LOOP
    FETCH cur_global_psets
      INTO v_instance_id,v_instance_name,v_pset_is_trkd;
    EXIT WHEN cur_global_psets%NOTFOUND;
    INSERT INTO tmp_instance_table
      VALUES(v_instance_id,NULL,'PSet',v_instance_name,v_pset_is_trkd);    
    load_parameters(v_instance_id);
  END LOOP;
  CLOSE cur_global_psets;
 
  /* load edsources */
  OPEN cur_edsources;
  LOOP
    FETCH cur_edsources INTO v_instance_id,v_template_id;
    EXIT WHEN cur_edsources%NOTFOUND;
    INSERT INTO tmp_instance_table
      VALUES(v_instance_id,v_template_id,'EDSource',NULL,NULL);
    load_parameters(v_instance_id);
  END LOOP;
  CLOSE cur_edsources;

  /* load essources */
  OPEN cur_essources;
  LOOP
    FETCH cur_essources
      INTO v_instance_id,v_template_id,v_instance_name,v_prefer;
    EXIT WHEN cur_essources%NOTFOUND;
    INSERT INTO tmp_instance_table
    VALUES(v_instance_id,v_template_id,'ESSource',v_instance_name,v_prefer);
    load_parameters(v_instance_id);
  END LOOP;
  CLOSE cur_essources;

  /* load esmodules */
  OPEN cur_esmodules;
  LOOP
    FETCH cur_esmodules
      INTO v_instance_id,v_template_id,v_instance_name,v_prefer;
    EXIT WHEN cur_esmodules%NOTFOUND;
    INSERT INTO tmp_instance_table
    VALUES(v_instance_id,v_template_id,'ESModule',v_instance_name,v_prefer);
    load_parameters(v_instance_id);
  END LOOP;
  CLOSE cur_esmodules;

  /* load services */
  OPEN cur_services;
  LOOP
    FETCH cur_services INTO v_instance_id,v_template_id;
    EXIT WHEN cur_services%NOTFOUND;
    INSERT INTO tmp_instance_table
      VALUES(v_instance_id,v_template_id,'Service',NULL,NULL);
    load_parameters(v_instance_id);
  END LOOP;
  CLOSE cur_services;

  /* load modules from *paths* */
  OPEN cur_modules_from_paths;
  LOOP
    FETCH cur_modules_from_paths
      INTO v_instance_id,v_template_id,v_instance_name;
    EXIT WHEN cur_modules_from_paths%NOTFOUND;
    INSERT INTO tmp_instance_table
      VALUES(v_instance_id,v_template_id,'Module',v_instance_name,NULL);
    load_parameters(v_instance_id);
  END LOOP;
  CLOSE cur_modules_from_paths;

  /* load modules from *sequences* */
  OPEN cur_modules_from_sequences;
  LOOP
    FETCH cur_modules_from_sequences
      INTO v_instance_id,v_template_id,v_instance_name;
    EXIT WHEN cur_modules_from_sequences%NOTFOUND;
    INSERT INTO tmp_instance_table
      VALUES(v_instance_id,v_template_id,'Module',v_instance_name,NULL);
    load_parameters(v_instance_id);
  END LOOP;
  CLOSE cur_modules_from_sequences;

  /* load paths */
  OPEN cur_paths;
  LOOP
    FETCH cur_paths INTO v_instance_id,v_instance_name,v_endpath;
    EXIT WHEN cur_paths%NOTFOUND;
    INSERT INTO tmp_instance_table
      VALUES(v_instance_id,NULL,'Path',v_instance_name,v_endpath);
  END LOOP;
  CLOSE cur_paths;

  /* load sequences */
  OPEN cur_sequences;
  LOOP
    FETCH cur_sequences INTO v_instance_id,v_instance_name;
    EXIT WHEN cur_sequences%NOTFOUND;
    INSERT INTO tmp_instance_table
      VALUES(v_instance_id,NULL,'Sequence',v_instance_name,NULL);
  END LOOP;
  CLOSE cur_sequences;

  /* load streams */
  OPEN cur_streams;
  LOOP
    FETCH cur_streams INTO v_instance_id,v_instance_name;
    EXIT WHEN cur_streams%NOTFOUND;
    INSERT INTO tmp_instance_table
      VALUES(v_instance_id,NULL,'Stream',v_instance_name,NULL);
  END LOOP;
  CLOSE cur_streams;
  
  /* load path-path associations */
  OPEN cur_path_path;
  LOOP
    FETCH cur_path_path INTO v_parent_id,v_child_id,v_sequence_nb;
    EXIT WHEN cur_path_path%NOTFOUND;
    INSERT INTO tmp_path_entries
      VALUES(v_parent_id,v_child_id,v_sequence_nb,'Path');
  END LOOP;
  CLOSE cur_path_path;

  /* load path-sequence associations */
  OPEN cur_path_sequence;
  LOOP
    FETCH cur_path_sequence INTO v_parent_id,v_child_id,v_sequence_nb;
    EXIT WHEN cur_path_sequence%NOTFOUND;
    INSERT INTO tmp_path_entries
      VALUES(v_parent_id,v_child_id,v_sequence_nb,'Sequence');
  END LOOP;
  CLOSE cur_path_sequence;

  /* load path-module associations */
  OPEN cur_path_module;
  LOOP
    FETCH cur_path_module INTO v_parent_id,v_child_id,v_sequence_nb;
    EXIT WHEN cur_path_module%NOTFOUND;
    INSERT INTO tmp_path_entries
      VALUES(v_parent_id,v_child_id,v_sequence_nb,'Module');
  END LOOP;
  CLOSE cur_path_module;

  /* load sequence-sequence associations */
  OPEN cur_sequence_sequence;
  LOOP
    FETCH cur_sequence_sequence INTO v_parent_id,v_child_id,v_sequence_nb;
    EXIT WHEN cur_sequence_sequence%NOTFOUND;
    INSERT INTO tmp_sequence_entries
      VALUES(v_parent_id,v_child_id,v_sequence_nb,'Sequence');
  END LOOP;
  CLOSE cur_sequence_sequence;

  /* load sequence-module associations */
  OPEN cur_sequence_module;
  LOOP
    FETCH cur_sequence_module INTO v_parent_id,v_child_id,v_sequence_nb;
    EXIT WHEN cur_sequence_module%NOTFOUND;
    INSERT INTO tmp_sequence_entries
      VALUES(v_parent_id,v_child_id,v_sequence_nb,'Module');
  END LOOP;
  CLOSE cur_sequence_module;

  /* load stream-path associations*/
  OPEN cur_stream_path;
  LOOP
    FETCH cur_stream_path INTO v_parent_id,v_child_id;
    EXIT WHEN cur_stream_path%NOTFOUND;
    INSERT INTO tmp_stream_entries VALUES(v_parent_id,v_child_id);
  END LOOP;
  CLOSE cur_stream_path;

  /* generate the final result set by selecting the temporary table */
  OPEN instance_cursor FOR
    SELECT DISTINCT instance_id,template_id,instance_type,instance_name,flag
    FROM tmp_instance_table;
  RETURN instance_cursor;  
END;  
/



--
-- FUNCTION get_parameters
--
CREATE FUNCTION get_parameters
  RETURN types.ref_cursor
AS
  parameter_cursor types.ref_cursor;
BEGIN
  OPEN parameter_cursor FOR
    SELECT parameter_id,parameter_type,
           parameter_name,parameter_trkd,parameter_seqnb,parent_id
    FROM tmp_parameter_table;
  RETURN parameter_cursor;
END;
/


--
-- FUNCTION get_boolean_values
--
CREATE FUNCTION get_boolean_values
  RETURN types.ref_cursor
AS
  value_cursor types.ref_cursor;
BEGIN
  OPEN value_cursor FOR
    SELECT DISTINCT parameter_id,parameter_value FROM tmp_boolean_table;
  RETURN value_cursor;
END;
/


--
-- FUNCTION get_int_values
--
CREATE FUNCTION get_int_values
  RETURN types.ref_cursor
AS
  value_cursor types.ref_cursor;
BEGIN
  OPEN value_cursor FOR
    SELECT DISTINCT parameter_id,parameter_value,sequence_nb,hex
      FROM tmp_int_table;
  RETURN value_cursor;
END;
/


--
-- FUNCTION get_real_values
--
CREATE FUNCTION get_real_values
  RETURN types.ref_cursor
AS
  value_cursor types.ref_cursor;
BEGIN
  OPEN value_cursor FOR
    SELECT DISTINCT parameter_id,parameter_value,sequence_nb FROM tmp_real_table;
  RETURN value_cursor;
END;
/


--
-- FUNCTION get_string_values
--
CREATE FUNCTION get_string_values
  RETURN types.ref_cursor
AS
  value_cursor types.ref_cursor;
BEGIN
  OPEN value_cursor FOR
    SELECT DISTINCT parameter_id,parameter_value,sequence_nb FROM tmp_string_table;
  RETURN value_cursor;
END;
/


--
-- FUNCTION get_path_entries
--
CREATE FUNCTION get_path_entries
  RETURN types.ref_cursor
AS
  entry_cursor types.ref_cursor;
BEGIN
  OPEN entry_cursor FOR
    SELECT path_id, entry_id, sequence_nb, entry_type FROM tmp_path_entries
    ORDER BY path_id ASC, sequence_nb ASC;
  RETURN entry_cursor;
END;
/


--
-- FUNCTION get_sequence_entries
--
CREATE FUNCTION get_sequence_entries
  RETURN types.ref_cursor
AS
  entry_cursor types.ref_cursor;
BEGIN
  OPEN entry_cursor FOR
    SELECT sequence_id, entry_id, sequence_nb, entry_type FROM tmp_sequence_entries
    ORDER BY sequence_id ASC, sequence_nb ASC;
  RETURN entry_cursor;
END;
/


--
-- FUNCTION get_stream_entries
--
CREATE FUNCTION get_stream_entries
  RETURN types.ref_cursor
AS
  entry_cursor types.ref_cursor;
BEGIN
  OPEN entry_cursor FOR
    SELECT stream_id, path_id FROM tmp_stream_entries;
  RETURN entry_cursor;
END;
/


COMMIT;
