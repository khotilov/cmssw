package confdb.data;

/**
 * VectorParameter
 * ---------------
 * @author Philipp Schieferdecker
 *
 * parameter base class for vector-type parameters.
 */
abstract public class VectorParameter extends Parameter
{
    //
    // construction
    //

    /** standard constructor */
    public VectorParameter(String name,boolean isTracked,boolean isDefault)
    {
	super(name,isTracked,isDefault);
    }

    //
    // abstract member functions
    //
    
    /** number of values */
    abstract public int vectorSize();
    
    /** i-th value of a vector type parameter */
    abstract public Object value(int i);
    
    /** set i-th value of a vector-type parameter */
    abstract public boolean setValue(int i,String valueAsString);

    /** remove i-th value from vector type parameter */
    abstract public Object removeValue(int i);
    
    
    //
    // member functions
    //

    /** override Parameter::isValueSet() */
    public boolean isValueSet()
    {
	return (isTracked()) ? true : super.isValueSet();
    }
}
