package confdb.converter.ascii;

import confdb.converter.ConverterEngine;
import confdb.converter.IESSourceWriter;
import confdb.data.ESSourceInstance;

public class AsciiESSourceWriter extends AsciiInstanceWriter implements IESSourceWriter 
{
	public String toString( ESSourceInstance essource, ConverterEngine converterEngine ) 
	{
		return toString( "es_source", essource, converterEngine );
	}

}
