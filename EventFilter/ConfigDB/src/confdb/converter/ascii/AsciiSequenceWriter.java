package confdb.converter.ascii;

import confdb.converter.ConverterEngine;
import confdb.converter.ISequenceWriter;
import confdb.data.Sequence;

public class AsciiSequenceWriter implements ISequenceWriter 
{

	public String toString( Sequence sequence, ConverterEngine converterEngine, String indent ) 
	{
		String str = indent + "sequence " + decorateName( sequence.name()) 
			       + " = { ";
		for ( int i = 0; i < sequence.entryCount(); i++  )
		{
			str += decorate( sequence.entry(i).name() );
			if ( i + 1 < sequence.entryCount() )
				str += " & ";
		}
		str += " }" + converterEngine.getNewline();
		return str;
	}
	
	protected String decorateName( String name )
	{
		return name;
	}

	protected String decorate( String name )
	{
		return name;
	}

}
