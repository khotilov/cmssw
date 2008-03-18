package confdb.converter;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import confdb.converter.ConverterBase;
import confdb.converter.ConverterEngine;
import confdb.converter.ConverterException;
import confdb.converter.ConverterFactory;
import confdb.converter.DbProperties;
import confdb.converter.OfflineConverter;
import confdb.converter.IConfigurationWriter.WriteProcess;
import confdb.data.ConfigurationModifier;
import confdb.data.IConfiguration;
import confdb.data.ModifierInstructions;
import confdb.db.ConfDB;
import confdb.db.ConfDBSetups;

public class BrowserConverter extends OfflineConverter
{
    static private HashMap<Integer,BrowserConverter> map = new HashMap<Integer,BrowserConverter>();
    
    private BrowserConverter(String dbType,String dbUrl,
			     String dbUser,String dbPwrd) throws ConverterException
    {
    	super( "HTML", dbType, dbUrl, dbUser, dbPwrd );	
    }
    
    public String getConfigString(int configId,
				  String format,
				  ModifierInstructions modifications,
				  boolean asFragment)
	throws ConverterException,
	       ClassNotFoundException,
	       InstantiationException,
	       IllegalAccessException
    {
    	IConfiguration config = getConfiguration(configId);
    	ConfigurationModifier modifier = new ConfigurationModifier(config);
    	modifier.modify(modifications);
	
    	ConverterEngine engine = ConverterFactory.getConverterEngine( format );
    	if (asFragment)
	    return engine.getConfigurationWriter().toString(modifier,
							    WriteProcess.NO);
    	else
	    return engine.getConfigurationWriter().toString(modifier,
							    WriteProcess.YES);
    }

	
	protected void finalize() throws Throwable
	{
		super.finalize();
		ConfDB db = getDatabase();
		if ( db != null )
			db.disconnect();
	}
	
	static public BrowserConverter getConverter( int dbIndex ) throws ConverterException 
	{
	    ConfDBSetups dbs = new ConfDBSetups();
	    DbProperties dbProperties = new DbProperties( dbs, dbIndex, "convertme!" );
	    if ( dbProperties.getDbUser().endsWith( "_w" ) )
		   	dbProperties.setDbUser( "cms_hlt_r" );
	    else
	    	dbProperties.setDbUser( "cms_hlt_reader" );
		BrowserConverter converter = map.get( new Integer( dbIndex ) );
		if ( converter == null )
		{
			converter = new BrowserConverter( dbs.type( dbIndex ), dbProperties.getDbURL(), dbProperties.getDbUser(), "convertme!" );		
			map.put( new Integer( dbIndex ), converter );
		}
		return converter;
	}

	static public void deleteConverter( ConverterBase converter )
	{
		Set<Map.Entry<Integer,BrowserConverter>> entries = map.entrySet();
		for ( Map.Entry<Integer,BrowserConverter> entry : entries )
		{
			if ( entry.getValue() == converter )
			{
				map.remove( entry.getKey() );
				return;
			}
		}
	}
	
	static public void clearCache()
	{
		map = new HashMap<Integer,BrowserConverter>();
	}
	
    static public String[] listDBs()
    {
    	ConfDBSetups dbs = new ConfDBSetups();
    	ArrayList<String> list = new ArrayList<String>();
     	for ( int i = 0; i < dbs.setupCount(); i++ )
     	{
    		String name = dbs.name(i);
    		if ( name != null && name.length() > 0  )
    		{
    			String host = dbs.host(i);
    			if (     host != null 
       				 && !host.equalsIgnoreCase("localhost") 
   				     && !host.endsWith( ".cms") )
   				     {
   				    	 list.add( dbs.labelsAsArray()[i] );
   				     }
    		}
    	}
    	return list.toArray( new String[ list.size() ] );
    }
    
    static public String[] getAnchors( int dbIndex, int configKey )
    {
		ArrayList<String> list = new ArrayList<String>();
		ConverterBase converter = null;
    	try {
			converter = BrowserConverter.getConverter( dbIndex );
			IConfiguration conf = converter.getConfiguration( configKey );
			if ( conf == null )
				list.add( "??" );
			else
			{
				if ( conf.pathCount() > 0 )
					list.add( "paths" );
				if ( conf.sequenceCount() > 0 )
					list.add( "sequences" );
				if ( conf.moduleCount() > 0 )
					list.add( "modules" );
				if ( conf.edsourceCount() > 0 )
					list.add( "ed_sources" );
				if ( conf.essourceCount() > 0 )
					list.add( "es_sources" );
				if ( conf.esmoduleCount() > 0 )
					list.add( "es_modules" );
				if ( conf.serviceCount() > 0 )
					list.add( "services" );
			}
		} catch (ConverterException e) {
			list.add( e.toString() );
			if ( converter != null )
				BrowserConverter.deleteConverter( converter );
		}
		return list.toArray( new String[ list.size() ] );
    }
    

	
}
