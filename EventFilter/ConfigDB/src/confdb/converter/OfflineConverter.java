package confdb.converter;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import java.sql.Connection;

import confdb.data.*;

import confdb.converter.IConfigurationWriter.WriteProcess;


/**
 * OfflineConverter
 * ----------------
 * @author Philipp Schieferdecker
 *
 * Handle conversion of configuraiton stored in the database for
 * offline use, e.g. the GUI, webserver, or command line tool.
 */
public class OfflineConverter extends ConverterBase
{
    //
    // construction
    //
    
    /** constructor based on explicit connection information */
    public OfflineConverter(String format,
			    String dbType,String dbUrl,String dbUser,String dbPwrd)
	throws ConverterException
    {
	super(format,dbType,dbUrl,dbUser,dbPwrd);
    }

    
    //
    // member functions
    //
    
    /** retrieve the configuration string for the given configId */
    public String getConfigString(int configId,
				  ModifierInstructions modifications,
				  boolean asFragment)
	throws ConverterException
    {
	IConfiguration config = getConfiguration(configId);
	return getConfigString(config,modifications,asFragment);
    }
    
    /** retrieve the configuration string for an IConfiguration object */
    public String getConfigString(IConfiguration config,
				  ModifierInstructions modifications,
				  boolean asFragment)
	throws ConverterException
    {
	ConfigurationModifier modifier = new ConfigurationModifier(config);
	modifier.modify(modifications);
	
	if (asFragment)
	    return getConverterEngine()
		.getConfigurationWriter().toString(modifier,WriteProcess.NO);
	else
	    return getConverterEngine()
		.getConfigurationWriter().toString(modifier,WriteProcess.YES);
    }
    

    //
    // main method, for testing
    //
    public static void main(String[] args)
    {
	String  configId    =          "";
	String  format      =     "Ascii";
	boolean asFragment  =       false;

	String  dbType      =     "mysql";
	String  dbHost      = "localhost";
	String  dbPort      =      "3306";
	String  dbName      =     "hltdb";
	String  dbUser      =          "";
	String  dbPwrd      =          "";

	HashMap<String,String> cnvArgs = new HashMap<String,String>();
	
	for (int iarg=0;iarg<args.length;iarg++) {
	    String arg = args[iarg];
	    if (arg.equals("-c")||arg.equals("--config")) {
		iarg++; configId = args[iarg];
	    }
	    else if (arg.equals("-f")||arg.equals("--format")) {
		iarg++; format = args[iarg];
	    }
	    else if (arg.equals("--cff")) {
		asFragment = true;
		cnvArgs.put("noedsources","");
		cnvArgs.put("nooutput","");
	    }
	    else if (arg.equals("-t")||arg.equals("--dbtype")) {
		iarg++; dbType = args[iarg];
	    }
	    else if (arg.equals("-h")||arg.equals("--dbhost")) {
		iarg++; dbHost = args[iarg];
	    }
	    else if (arg.equals("-p")||arg.equals("--dbport")) {
		iarg++; dbPort = args[iarg];
	    }
	    else if (arg.equals("-d")||arg.equals("--dbname")) {
		iarg++; dbName = args[iarg];
	    }
	    else if (arg.equals("-u")||arg.equals("--dbuser")) {
		iarg++; dbUser = args[iarg];
	    }
	    else if (arg.equals("-s")||arg.equals("--dbpwrd")) {
		iarg++; dbPwrd = args[iarg];
	    }
	    else if (arg.startsWith("--no")) {
		String key = arg.substring(2);
		String val = "";
		cnvArgs.put(key,val);
	    }
	    else if (arg.startsWith("--")) {
		String key = arg.substring(2);
		String val = args[++iarg];
		cnvArgs.put(key,val);
	    }
	    else {
		System.err.println("ERROR: invalid option '" + arg + "'!");
		System.exit(0);
	    }
	}
	
	if (!format.equals("Ascii")&&
	    !format.equals("Python")&&
	    !format.equals("Html")) {
	    System.err.println("ERROR: Invalid format '"+format+"'");
	}
    
	String dbUrl = "";
	if (dbType.equalsIgnoreCase("mysql")) {
	    dbUrl  = "jdbc:mysql://"+dbHost+":"+dbPort+"/"+dbName;
	}
	else if (dbType.equalsIgnoreCase("oracle")) {
	    dbUrl = "jdbc:oracle:thin:@//"+dbHost+":"+dbPort+"/"+dbName;
	}
	else {
	    System.err.println("ERROR: Unknwown db type '"+dbType+"'");
	    System.exit(0);
	}
	
	System.out.println("dbURl  = " + dbUrl);
	System.out.println("dbUser = " + dbUser);
	System.out.println("dbPwrd = " + dbPwrd);

	try {
	    ModifierInstructions modifications = new ModifierInstructions();
	    modifications.interpretArgs(cnvArgs);
	    OfflineConverter cnv = 
		new OfflineConverter(format,dbType,dbUrl,dbUser,dbPwrd);
	    System.out.println(cnv.getConfigString(Integer.parseInt(configId),
						   modifications,
						   asFragment));
	}
	catch(DataException e) {
	    System.err.println("ERROR: " + e.getMessage());
	}
	catch(ConverterException e) {
	    System.out.println("ERROR: " + e.getMessage());
	}
	
    }
    
}
