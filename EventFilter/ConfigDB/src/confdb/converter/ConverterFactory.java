package confdb.converter;


public class ConverterFactory 
{
	static private final String thisPackage = ConverterFactory.class.getPackage().getName();

	private String writerPackage = "";
	private String subPackage = "ascii";
	private String classHeader = "Ascii";
	

	static public ConverterEngine getConverterEngine( String typeOfConverter ) throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		return new ConverterFactory( typeOfConverter ).getConverterEngine();
	}
	
	
	private ConverterFactory( String typeOfConverter )
	{
		if ( typeOfConverter != null )
		{
			subPackage = typeOfConverter.toLowerCase();
			classHeader = subPackage;
			if ( typeOfConverter.indexOf( '.' ) != -1 )
				classHeader = typeOfConverter.substring( typeOfConverter.lastIndexOf( '.' ) + 1 );
			classHeader = classHeader.substring( 0, 1 ).toUpperCase() + classHeader.substring( 1 );
		}
	}
	
	
	public ConverterEngine getConverterEngine() throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		writerPackage = thisPackage + "." + subPackage;
		ConverterEngine converterEngine = new ConverterEngine( subPackage );
		converterEngine.setConfigurationWriter( getConfigurationWriter() );
		converterEngine.setParameterWriter( getParameterWriter() );

		converterEngine.setEDSourceWriter( getEDSourceWriter() );
		converterEngine.setESSourceWriter( getESSourceWriter() );
		converterEngine.setESModuleWriter( getESModuleWriter() );
		converterEngine.setServiceWriter( getServiceWriter() );
		converterEngine.setModuleWriter( getModuleWriter() );
		converterEngine.setPathWriter( getPathWriter() );
		converterEngine.setSequenceWriter( getSequenceWriter() );
		
		return converterEngine;
	}

	public IConfigurationWriter getConfigurationWriter() 
	  throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		return (IConfigurationWriter)getWriter( "ConfigurationWriter" );
	}

	public IPathWriter getPathWriter() 
	  throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		return (IPathWriter)getWriter( "PathWriter" );
	}

	
	public IServiceWriter getServiceWriter() 
	  throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		return (IServiceWriter)getWriter( "ServiceWriter" );
	}
	
	public IParameterWriter getParameterWriter() 
	  throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		return (IParameterWriter)getWriter( "ParameterWriter" );
	}

	public IEDSourceWriter getEDSourceWriter()
	  throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		return (IEDSourceWriter)getWriter( "EDSourceWriter" );
	}

	public IESSourceWriter getESSourceWriter()
	  throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		return (IESSourceWriter)getWriter( "ESSourceWriter" );
	}

	public IESModuleWriter getESModuleWriter()
	  throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		return (IESModuleWriter)getWriter( "ESModuleWriter" );
	}

	public IModuleWriter getModuleWriter()
	  throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		return (IModuleWriter)getWriter( "ModuleWriter" );
	}
	
	public ISequenceWriter getSequenceWriter()
	  throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		return (ISequenceWriter)getWriter( "SequenceWriter" );
	}

	
	private Object getWriter( String type ) throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		String className = writerPackage + "." + classHeader + type;
		Class<?> c = Class.forName( className );
		return c.newInstance();
	}

	
}
