package confdb.data;

/**
 * UInt32Parameter
 * ---------------
 * @author Philipp Schieferdecker
 *
 * parameter base class for scalar parameters of type uint32.
 */
public class UInt32Parameter extends ScalarParameter
{
    //
    // member data
    //

    /** parameter type string */
    private static final String type = "uint32";
    
    /** parameter values */
    private Long value = null;
         
    /** flag to indicate that this integer is given in hex format */
    private boolean isHex = false;
    
    
    //
    // construction
    //

    /** standard constructor */
    public UInt32Parameter(String name,Long value,
			   boolean isTracked,boolean isDefault)
    {
	super(name,isTracked,isDefault);
	isValueSet = (value!=null);
	if (isValueSet)	this.value = new Long(value.longValue());
    }
    
    /** constructor from string */
    public UInt32Parameter(String name,String valueAsString,
			   boolean isTracked,boolean isDefault)
    {
	super(name,isTracked,isDefault);
	setValue(valueAsString);
    }
    
    //
    // member functions
    //
    
    /** make a clone of the parameter */
    public Parameter clone(Object parent)
    {
	UInt32Parameter result = new UInt32Parameter(name,value,isTracked,isDefault);
	result.setParent(parent);
	return result;
    }
    
    /** type of the parameter as a string */
    public String type() { return type; }
    
    /** retrieve the value of the parameter */
    public Object value() { return value; }

    /** hex format? */
    public boolean isHex() { return isHex; }
    
    /** retrieve the value of the parameter as a string */
    public String valueAsString()
    {
	if (!isValueSet) return new String();
	return (isHex) ? "0x"+Long.toHexString(value) : value.toString();
    }

    /** set the value  the parameter, indicate if default */
    public boolean setValue(String valueAsString)
    {
	if (valueAsString==null||valueAsString.length()==0) {
	    isValueSet = false;
	    value      = null;
	}
	else {
	    if (valueAsString.startsWith("+"))
		valueAsString = valueAsString.substring(1);
	    
	    isHex = false;
	    if (valueAsString.startsWith("0x")) {
		isHex = true;
		valueAsString = valueAsString.substring(2);
	    }
	    
	    try {
		this.value = (isHex) ?
		    new Long(Long.parseLong(valueAsString,16)) :
		    new Long(valueAsString);
		isValueSet = true;
	    }
	    catch (NumberFormatException e) {
		System.err.println("UInt32Parameter.setValue " +
				   "NumberFormatException: "+
				   e.getMessage());
		return false;
	    }
	}
	return true;
    }

}
