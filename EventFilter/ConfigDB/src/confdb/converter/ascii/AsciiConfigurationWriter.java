package confdb.converter.ascii;

import confdb.converter.ConverterEngine;
import confdb.converter.IConfigurationWriter;
import confdb.converter.IEDSourceWriter;
import confdb.converter.IESSourceWriter;
import confdb.converter.IESModuleWriter;
import confdb.converter.IModuleWriter;
import confdb.converter.IParameterWriter;
import confdb.converter.IPathWriter;
import confdb.converter.ISequenceWriter;
import confdb.converter.IServiceWriter;
import confdb.data.IConfiguration;
import confdb.data.EDSourceInstance;
import confdb.data.ESSourceInstance;
import confdb.data.ESModuleInstance;
import confdb.data.ModuleInstance;
import confdb.data.Parameter;
import confdb.data.Path;
import confdb.data.Sequence;
import confdb.data.ServiceInstance;

public class AsciiConfigurationWriter implements IConfigurationWriter 
{
	protected ConverterEngine converterEngine = null;

	public String toString( IConfiguration conf, WriteProcess writeProcess  )
	{
		String indent = "  ";
		StringBuffer str = new StringBuffer( 100000 );
		str.append( "// " + conf.name() + " V" + conf.version()
  		   + " (" + conf.releaseTag() + ")" + converterEngine.getNewline() + converterEngine.getNewline() );

		if ( writeProcess == WriteProcess.YES )
			str.append( "process " + conf.processName() + " = {" + converterEngine.getNewline() );
		else
			indent = "";

		IPathWriter pathWriter = converterEngine.getPathWriter();
		for ( int i = 0; i < conf.pathCount(); i++ )
		{
			Path path = conf.path(i);
			str.append( pathWriter.toString( path, converterEngine, indent ) );
		}

		ISequenceWriter sequenceWriter = converterEngine.getSequenceWriter();
		for ( int i = 0; i < conf.sequenceCount(); i++ )
		{
			Sequence sequence = conf.sequence(i);
			str.append( sequenceWriter.toString(sequence, converterEngine, indent ) );
		}

		IParameterWriter parameterWriter = converterEngine.getParameterWriter();
		for ( int i = 0; i < conf.psetCount(); i++ )
		{
			Parameter pset = conf.pset(i);
			str.append( parameterWriter.toString( pset, converterEngine, indent ) );
		}


		IEDSourceWriter edsourceWriter = converterEngine.getEDSourceWriter();
		for ( int i = 0; i < conf.edsourceCount(); i++ )
		{
			EDSourceInstance edsource = conf.edsource(i);
			str.append( edsourceWriter.toString(edsource, converterEngine, indent ) );
		}

		IESSourceWriter essourceWriter = converterEngine.getESSourceWriter();
		for ( int i = 0; i < conf.essourceCount(); i++ )
		{
			ESSourceInstance essource = conf.essource(i);
			str.append( essourceWriter.toString(essource, converterEngine, indent ) );
		}


		IESModuleWriter esmoduleWriter = converterEngine.getESModuleWriter();
		for ( int i = 0; i < conf.esmoduleCount(); i++ )
		{
			ESModuleInstance esmodule = conf.esmodule(i);
			str.append( esmoduleWriter.toString( esmodule, converterEngine, indent ) );
		}


		IServiceWriter serviceWriter = converterEngine.getServiceWriter();
		for ( int i = 0; i < conf.serviceCount(); i++ )
		{
			ServiceInstance service = conf.service(i);
			str.append( serviceWriter.toString( service, converterEngine, indent ) );
		}

		IModuleWriter moduleWriter = converterEngine.getModuleWriter();
		for ( int i = 0; i < conf.moduleCount(); i++ )
		{
			ModuleInstance module = conf.module(i);
			str.append( moduleWriter.toString( module ) );
		}

		if ( writeProcess == WriteProcess.YES )
			str.append( converterEngine.getConfigurationTrailer() );
		return str.toString();
	}

	public void setConverterEngine( ConverterEngine converterEngine ) 
	{
		this.converterEngine = converterEngine;
	}
	
}
