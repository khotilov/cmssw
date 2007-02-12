package confdb.data;

import java.util.ArrayList;


public class TemplateFactory
{
    /** list holding valid module types */
    private static String[] validModuleTypes =
    { 
	"EDProducer","EDFilter","EDAnalyzer",
	"HLTProducer","HLTFilter",
	"ESProducer",
	"OutputModule"
    }; 
    
    /** check if a type is among valid module types */
    private static boolean isValidModuleType(String type)
    {
	for (String s : validModuleTypes) if (type.equals(s)) return true;
	return false;
    }
    
    /** parametrizez factory method to create a specific template type */
    public static Template create(String type,
				  String name,String cvsTag,int dbSuperId,
				  ArrayList<Parameter> parameters)
    {
	if (type.equals("Service"))
	    return new ServiceTemplate(name,cvsTag,dbSuperId,parameters);
	else if (type.equals("EDSource"))
	    return new EDSourceTemplate(name,cvsTag,dbSuperId,parameters);
	else if (type.equals("ESSource"))
	    return new ESSourceTemplate(name,cvsTag,dbSuperId,parameters);
	else if (isValidModuleType(type))
	    return new ModuleTemplate(name,cvsTag,dbSuperId,parameters,type);
	
	System.out.println("TemplateFactory: unknown type '"+type+"'");
	return null;
    }

}
