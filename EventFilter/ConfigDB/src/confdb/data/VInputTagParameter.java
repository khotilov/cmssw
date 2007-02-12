package confdb.data;

import java.util.ArrayList;


/**
 * VInputTagParameter
 * -----------------
 * @author Philipp Schieferdecker
 *
 * for parameters of type vector<InputTag>.
 */
public class VInputTagParameter extends VectorParameter
{
    //
    // data members
    //

    /** parameter type string */
    private static final String type = "VInputTag";
    
    /** parameter values */
    private ArrayList<InputTag> values = new ArrayList<InputTag>();
    
    
    //
    // construction
    //
    
    /** standard constructor */
    public VInputTagParameter(String name,ArrayList<String> values,
			      boolean isTracked,boolean isDefault)
    {
	super(name,isTracked,isDefault);
	for (String v : values) {
	    try {
		this.values.add(new InputTag(v));
	    }
	    catch (DataException e) {
		System.out.println("VInputTagParameter ctor error: " +
				   e.getMessage());
	    }
	}
	isValueSet = (this.values.size()>0);
    }
    
    /** constructor from a string */
    public VInputTagParameter(String name,String valuesAsString,
			      boolean isTracked,boolean isDefault)
    {
	super(name,isTracked,isDefault);
	setValue(valuesAsString);
    }

    //
    // member functions
    //
    
    /** make a clone of the parameter */
    public Parameter clone(Object parent)
    {
	ArrayList<String> strValues = new ArrayList<String>();
	for (InputTag tag : values) strValues.add(tag.toString());
	VInputTagParameter result = new VInputTagParameter(name,strValues,
							   isTracked,isDefault);
	result.setParent(parent);
	return result;
    }
    
    /** type of the parameter as a string */
    public String type() { return type; }
    
    /** retrieve the values of the parameter as a string */
    public String valueAsString()
    {
	String result = new String();
	if (isValueSet) {
	    for (InputTag tag : values) result += tag.toString() + ", ";
	    result = result.substring(0,result.length()-2);
	}
	return result;
    }

    /** set the parameter values from string */
    public boolean setValue(String valueAsString)
    {
	values.clear();
	if (valueAsString.length()==0) {
	    isValueSet = false;
	}
	else {
	    try {
		String[] strValues = valueAsString.split(",");
		for (int i=0;i<strValues.length;i++) {
		    String s = strValues[i];
		    while (s.startsWith(" ")) s = s.substring(1,s.length());
		    while (s.endsWith(" ")) s = s.substring(0,s.length()-1);
		    values.add(new InputTag(s));
		}
		isValueSet = true;
	    }
	    catch (DataException e) {
		System.out.println("VInputTagParameter.setValue error: " +
				   e.getMessage());
		return false;
	    }
	}
    	return true;
    }
    
    /** number of vector entries */
    public int vectorSize() { return values.size(); }

    /** i-th value of a vector type parameter */
    public Object value(int i) { return values.get(i); }

    /** set i-th value of a vector-type parameter */
    public boolean setValue(int i,String valueAsString)
    {
	try {
	    InputTag tag = new InputTag(valueAsString);
	    values.set(i,tag);
	}
	catch (DataException e) {
	    System.out.println("VInputTagParameter.setValue error: " +
			       e.getMessage());
	    return false;
	}
	return true;
    }

}


/**
 * InputTag
 * -------
 * @author Philipp Schieferdecker
 */
class InputTag
{
    //
    // member data
    //
    
    /** label */
    private String label    = null;

    /** instance */
    private String instance = null;

    /** process */
    private String process  = null;


    //
    // construction
    //

    /** standard constructor */
    public InputTag(String label,String instance,String process)
    {
	this.label    = label;
	this.instance = instance;
	this.process  = process;
    }
    
    /** constructor from string */
    public InputTag(String valueAsString) throws DataException
    {
	label    = "";
	instance = "";
	process  = "";
	String[] strValues = valueAsString.split(":");
	if (strValues.length>0&&strValues.length<4) {
	    label = strValues[0];
	    if (strValues.length>1) instance = strValues[1];
	    if (strValues.length>2) process  = strValues[2];
	}
	else throw new DataException("InputTag format is " +
				     "<label>[:<instance>[:<process>]]");
    }
    

    //
    // member functions
    //

    /** overload toString() */
    public String toString()
    {
	String result = label;
	if (instance.length()>0) {
	    result += ":" + instance;
	    if (process.length()>0) result += ":" + process;
	}
	return result;
    }
    
}
