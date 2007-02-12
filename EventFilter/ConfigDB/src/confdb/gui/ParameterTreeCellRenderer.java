package confdb.gui;

import javax.swing.*;
import javax.swing.tree.*;
import java.awt.*;

import confdb.data.Parameter;
import confdb.data.PSetParameter;
import confdb.data.VPSetParameter;


/**
 * ParameterTreeCellRenderer
 * -------------------------
 * @author Philipp Schieferdecker
 *
 * To display the appropriate icons in a TreeTable showing parameters.
 */
public class ParameterTreeCellRenderer extends DefaultTreeCellRenderer
{
    //
    // member data
    //

    /** node to be displayed */
    private Object    node      = null;

    /** icon for a ParameterSet node */
    private ImageIcon psetIcon  = new ImageIcon("icons/PSetIcon.png");

    /** icon for a vector<ParameterSet> node */
    private ImageIcon vpsetIcon = new ImageIcon("icons/VPSetIcon.png");
    
    /** icon for any other Parameter */
    private ImageIcon parameterIcon = new ImageIcon("icons/ParameterIcon.png");
    
    
    //
    // member functions
    //
    
    /** prepare the appropriate icon */
    public Icon prepareIcon()
    {
	if (node instanceof PSetParameter)  return psetIcon;
	if (node instanceof VPSetParameter) return vpsetIcon;
	return parameterIcon;
    }
    
    /** TreeCellRenderer */
    public Component getTreeCellRendererComponent(JTree   tree,
						  Object  value,
						  boolean sel,
						  boolean expanded,
						  boolean leaf,
						  int     row,
						  boolean hasFocus)
    {
	super.getTreeCellRendererComponent(tree,value,sel,
					   expanded,leaf,row,
					   hasFocus);
	node = value;
	setIcon(prepareIcon());
	return this;
    }
    
}
