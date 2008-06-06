package confdb.gui;

import javax.swing.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;

import java.util.Iterator;
import java.util.ArrayList;
import java.util.HashMap;

import confdb.data.*;


/**
 * PrescaleDialog
 * --------------
 * @author Philipp Schieferdecker
 *
 * Edit the prescale table, which is encoded in the configuration of
 * the PrescaleService.
 */
public class PrescaleDialog extends JDialog
{
    //
    // member data
    //
    
    /** reference to the configuration */
    private Configuration config;
    
    /** GUI components */
    private JTextField jTextFieldHLT    = new javax.swing.JTextField();
    private JTextField jTextFieldLevel1 = new javax.swing.JTextField();
    private JButton    jButtonOK        = new javax.swing.JButton();
    private JButton    jButtonApply     = new javax.swing.JButton();
    private JButton    jButtonCancel    = new javax.swing.JButton();
    private JTable     jTable           = new javax.swing.JTable();
    
    /** model for the prescale table */
    private PrescaleTableModel tableModel;

    /** index of the selected column */
    private int iColumn = -1;
    
    
    //
    // construction
    //

    /** standard constructor */
    public PrescaleDialog(JFrame jFrame,Configuration config)
    {
	super(jFrame,true);
	this.config = config;

	tableModel = new PrescaleTableModel();
	tableModel.initialize(config);
	jTable.setModel(tableModel);
	jTable.setDefaultRenderer(Integer.class,new PrescaleTableCellRenderer());
	jTextFieldHLT.setText(config.toString());
	
	jButtonCancel.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent e) {
		    setVisible(false);
		}
	    });
	jButtonApply.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent e) {
		    updatePrescaleService();
		}
	    });
	jButtonOK.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent e) {
		    updatePrescaleService();
		    setVisible(false);
		}
	    });
	jTable.addMouseListener(new MouseAdapter() {
		public void mousePressed(MouseEvent e) {
		    jTableShowPopup(e);
		}
		public void mouseReleased(MouseEvent e) {
		    jTableShowPopup(e);
		}
	    });
				
	setTitle("Prescale Editor");
	setContentPane(initComponents());
	adjustTableColumnWidths();
    }
    
    //
    // private member functions
    //
    
    /** update the configurations PrescaleService according to table data */
    private void updatePrescaleService()
    {
	tableModel.updatePrescaleService(config);
    }
    
    /** adjust the width of each table column */
    private void adjustTableColumnWidths()
    {
	int tableWidth = jTable.getPreferredSize().width;
	int columnCount = jTable.getColumnModel().getColumnCount();
	for (int i=0;i<columnCount;i++) {
	    int columnWidth = (i==0) ? tableWidth/2 : tableWidth/2/(columnCount-1);
	    jTable.getColumnModel().getColumn(i).setPreferredWidth(columnWidth);
	}
    }
    
    // listener callbacks
    private void jTableShowPopup(MouseEvent e)
    {
	if (!e.isPopupTrigger()) return;
	iColumn = jTable.columnAtPoint(e.getPoint());
	JPopupMenu popup = new JPopupMenu();
	JMenuItem menuItem = new JMenuItem("Add Column");
	menuItem.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae)
		{
		    tableModel.addColumn(iColumn,
					 JOptionPane
					 .showInputDialog("Enter the level1 "+
							  "label for the new "+
							  "column:"));
		    adjustTableColumnWidths();
		}
	    });
	popup.add(menuItem);
	if (iColumn>0) {
	    menuItem = new JMenuItem("Remove Column");
	    menuItem.addActionListener(new ActionListener() {
		    public void actionPerformed(ActionEvent ae)
		    {
			tableModel.removeColumn(iColumn);
			adjustTableColumnWidths();
		    }
		});
	    popup.add(menuItem);
	}
	popup.show(e.getComponent(),e.getX(),e.getY());
    }
    
    /** initialize GUI components */
    private JPanel initComponents()
    {
	JPanel jPanel = new JPanel();
	
        JLabel      jLabel1     = new javax.swing.JLabel();
        JLabel      jLabel2     = new javax.swing.JLabel();
        JScrollPane jScrollPane = new javax.swing.JScrollPane();
	
        jLabel1.setText("HLT:");
        jLabel2.setText("Level1:");
	
        jTextFieldHLT.setEditable(false);
        jTextFieldHLT.setBorder(BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.LOWERED));

        jTextFieldLevel1.setEditable(false);
        jTextFieldLevel1.setBorder(BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.LOWERED));
	

	


        jScrollPane.setViewportView(jTable);

        jButtonOK.setText("OK");
        jButtonApply.setText("Apply");
        jButtonCancel.setText("Cancel");

        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(jPanel);
        jPanel.setLayout(layout);
        layout.setHorizontalGroup(
				  layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
				  .add(layout.createSequentialGroup()
				       .addContainerGap()
				       .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
					    .add(jScrollPane, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 708, Short.MAX_VALUE)
					    .add(layout.createSequentialGroup()
						 .add(jLabel1)
						 .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
						 .add(jTextFieldHLT, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 299, Short.MAX_VALUE)
						 .add(18, 18, 18)
						 .add(jLabel2)
						 .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
						 .add(jTextFieldLevel1, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 298, Short.MAX_VALUE))
					    .add(org.jdesktop.layout.GroupLayout.TRAILING, layout.createSequentialGroup()
						 .add(jButtonCancel)
						 .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
						 .add(jButtonApply)
						 .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
						 .add(jButtonOK)))
				       .addContainerGap())
				  );
	
        layout.linkSize(new java.awt.Component[] {jButtonApply, jButtonCancel, jButtonOK}, org.jdesktop.layout.GroupLayout.HORIZONTAL);
	
        layout.setVerticalGroup(
				layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
				.add(org.jdesktop.layout.GroupLayout.TRAILING, layout.createSequentialGroup()
				     .add(22, 22, 22)
				     .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
					  .add(jLabel1, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
					  .add(jTextFieldHLT, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
					  .add(jLabel2, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
					  .add(jTextFieldLevel1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
				     .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
				     .add(jScrollPane, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 570, Short.MAX_VALUE)
				     .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
				     .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
					  .add(jButtonOK)
					  .add(jButtonApply)
					  .add(jButtonCancel))
				     .addContainerGap())
				);
	return jPanel;
    }

}


//
//PrescaleTableModel
//
class PrescaleTableModel extends AbstractTableModel
{
    /** the presacale table data structure */
    private PrescaleTable prescaleTable;

    /** update the table model according to configuration's PrescaleService */
    public void initialize(IConfiguration config)
    {
	prescaleTable = new PrescaleTable(config);
	fireTableDataChanged();
    }

    /** update the PrescaleService in configuration according to table data */
    public void updatePrescaleService(IConfiguration config)
    {
	ServiceInstance prescaleSvc = config.service("PrescaleService");
	if (prescaleSvc==null) {
	    System.err.println("No PrescaleService found.");
	    return;
	}
	
	StringBuffer labelsAsString = new StringBuffer();
	for (int i=0;i<prescaleTable.prescaleCount();i++) {
	    if (labelsAsString.length()>0) labelsAsString.append(",");
	    labelsAsString.append(prescaleTable.prescaleColumnName(i));
	}
	prescaleSvc.updateParameter("lvl1Labels","vstring",
				    labelsAsString.toString());
	
	
	VPSetParameter vpsetPrescaleTable =
	    (VPSetParameter)prescaleSvc.parameter("prescaleTable","VPSet");
	if (vpsetPrescaleTable==null) {
	    System.err.println("No VPSet prescaleTable found.");
	    return;
	}
	vpsetPrescaleTable.setValue("");
	
	for (int iPath=0;iPath<prescaleTable.pathCount();iPath++) {
	    if (!prescaleTable.isPrescaled(iPath)) continue;
	    String pathName = prescaleTable.pathName(iPath);
	    String prescalesAsString = prescaleTable.prescalesAsString(iPath);
	    ArrayList<Parameter> params    =new ArrayList<Parameter>();
	    StringParameter      sPathName =new StringParameter("pathName",
								 pathName,
								 true,false);
	    VUInt32Parameter     vPrescales=new VUInt32Parameter("prescales",
								 prescalesAsString,
								 true,false);
	    params.add(sPathName);
	    params.add(vPrescales);
	    vpsetPrescaleTable.addParameterSet(new PSetParameter("",params,
								 true,false));
	}
	prescaleSvc.setHasChanged();
    }
    

    /** add an additional column (-> lvl1Label) */
    public void addColumn(int i,String lvl1Label)
    {
	prescaleTable.addPrescaleColumn(i,lvl1Label);
	fireTableStructureChanged();
    }
    
    /** remove a column of prescales */
    public void removeColumn(int i)
    {
	prescaleTable.removePrescaleColumn(i);
	fireTableStructureChanged();
    }

    /** number of columns */
    public int getColumnCount() { return prescaleTable.prescaleCount()+1; }
    
    /** number of rows */
    public int getRowCount() { return prescaleTable.pathCount(); }
    
    /** get column name for colimn 'col' */
    public String getColumnName(int col)
    { return (col==0) ? "Path" : prescaleTable.prescaleColumnName(col-1); }
    
    /** get the value for row,col */
    public Object getValueAt(int row, int col)
    {
	return (col==0) ?
	    prescaleTable.pathName(row) : prescaleTable.prescale(row,col-1);
    }
    
    /** get the class of the column 'c' */
    public Class getColumnClass(int c)
    {
	return (c==0) ? String.class : Integer.class;
    }
    
    /** is a cell editable or not? */
    public boolean isCellEditable(int row, int col) { return col>0; }
    
    /** set the value of a table cell */
    public void setValueAt(Object value,int row, int col)
    {
	prescaleTable.setPrescale(row,col-1,(Integer)value);
    }
    
    
    /** check if a certain path is already in the list of rows */
    private boolean rowsContainPath(String pathName)
    {
	for (int iPath=0;iPath<prescaleTable.pathCount();iPath++)
	    if (pathName.equals(prescaleTable.pathName(iPath))) return true;
	return false;
    }
}


//
// PrescaleTableCellRenderer
//
class PrescaleTableCellRenderer extends DefaultTableCellRenderer
{
    public Component getTableCellRendererComponent(JTable table,
						   Object value,
						   boolean isSelected,
						   boolean hasFocus,
						   int row,int column)
    {
	setText(value.toString());
	if (value instanceof Integer) {
	    setHorizontalAlignment(SwingConstants.CENTER);
	    int valueAsInt = (Integer)value;
	    if (valueAsInt==0) setBackground(Color.RED);
	    else if (valueAsInt==1) setBackground(Color.GREEN);
	    else setBackground(Color.ORANGE);
	}
	return this;
    }
}
