package confdb.data;

import java.util.Iterator;
import java.util.ArrayList;


/**
 * Block
 * -----
 * @author Philipp Schieferdecker
 *
 */
public class Block
{
    //
    // member data
    //
    
    /** name of the block */
    private String name;
    
    /** list of parameters (only accepted content!) */
    private ArrayList<Parameter> parameters = new ArrayList<Parameter>();
    
    
    //
    // construction
    //

    /** standard constructor */
    public Block(Instance instance,String[] paramNames)
    {
	name = "block_"+instance.name();
	for (String paramName : paramNames) {
	    Parameter param = instance.parameter(paramName);
	    if (param!=null) parameters.add(param);
	}
    }
    
    
    //
    // member functions
    //
    
    /** get the name of this block */
    public String name() { return this.name; }
    
    /** get parameter iterator */
    public Iterator<Parameter> parameterIterator()
    {
	return parameters.iterator();
    }
    
}
