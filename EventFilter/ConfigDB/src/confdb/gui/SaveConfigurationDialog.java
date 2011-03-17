package confdb.gui;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.tree.*;
import java.awt.*;
import java.awt.event.*;

import org.jdesktop.layout.*;

import confdb.data.Directory;
import confdb.data.Configuration;
import confdb.data.ConfigInfo;
    
import confdb.db.ConfDB;
import confdb.db.DatabaseException;


/**
 * SaveConfigurationDialog
 * -----------------------
 * @author Philipp Schieferdecker
 *
 */
public class SaveConfigurationDialog extends JDialog
{
    //
    // member data
    //
    
    /** define symbols which are forbidden in a dir/config name */
    private static final char[] forbiddenSymbols = { ' ','/','*','\'','"',
						     '[',']','{','}','(',')',
						     '&','$','#',':',';','<','>' };
    
    /** configuration to be saved */
    private Configuration config = null;
    
    /** reference to the database */
    private ConfDB database = null;
    
    /** currently selected directory */
    private Directory selectedDir = null;
    
    /** was a valid choice made? */
    private boolean validChoice = false;
    
    /** directory tree model */
    private DirectoryTreeModel treeModel = null;

    /** GUI components */
    private JTree      jTreeDirectories;
    private JTextField jTextFieldConfigName = new JTextField();
    private JTextField jTextFieldComment    = new JTextField();
    private JButton    jButtonOk            = new JButton();
    private JButton    jButtonCancel        = new JButton();
    
    /** action commands */
    private static final String OK            = new String("OK");
    private static final String CANCEL        = new String("Cancel");
    
    
    //
    // construction
    //
    
    /** standard constructor */
    public SaveConfigurationDialog(JFrame frame,ConfDB database,
				   Configuration config,String comment)
    {
	super(frame,true);
	this.config = config;
	this.database = database;

	jTextFieldComment.setText(comment);
	setTitle("Save Configuration");
	
	
	// initialize tree
	try {
	    Directory rootDir = database.loadConfigurationTree();
	    treeModel = new DirectoryTreeModel(rootDir);
	    jTreeDirectories = new JTree(treeModel);
	    jTreeDirectories.setEditable(true);
	    jTreeDirectories
		.setCellRenderer(new DirectoryTreeCellRenderer());
	    jTreeDirectories
		.setCellEditor(new DirectoryTreeCellEditor(jTreeDirectories,
							   new DirectoryTreeCellRenderer()));
	}
	catch (DatabaseException e) {
	    jTreeDirectories = new JTree();
	}
	
	setContentPane(createContentPane());
	
	if (config.version()==0) {
	    jTextFieldConfigName.setText(config.name());
	    jTextFieldConfigName.selectAll();
	}
    
	// register listener callbacks
	jTreeDirectories
	    .addMouseListener(new DirectoryTreeMouseListener(jTreeDirectories,
							     database));

	// register focus listener to save current edit by clicking elsewhere.
	// See: TreeDirectoryFocusListener.java
	jTreeDirectories.addFocusListener(new TreeDirectoryFocusListener(jTreeDirectories, database));
	
	
	jTreeDirectories.addTreeSelectionListener(new TreeSelectionListener() {
		public void valueChanged(TreeSelectionEvent e) {
		    jTreeDirectoriesValueChanged(e);
		}
	    });
	jTextFieldConfigName.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent e) {
		    jTextFieldConfigNameActionPerformed(e);
		}
	    });
	jTextFieldConfigName
	    .getDocument().addDocumentListener(new DocumentListener() {
		    public void insertUpdate(DocumentEvent e) {
			jTextFieldConfigNameInsertUpdate(e);
		    }
		    public void removeUpdate(DocumentEvent e) {
			jTextFieldConfigNameRemoveUpdate(e);
		    }
		    public void changedUpdate(DocumentEvent e) {}
		});
	jButtonOk.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent e) {
		    jButtonOkActionPerformed(e);
		}
	    });
	jButtonCancel.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent e) {
		    jButtonCancelActionPerformed(e);
		}
	    });
    }
    
    
    //
    // member functions
    //
    
    /** indicate if valid choise was made */
    public boolean validChoice() { return validChoice; }
    
    /** retrieve comment */
    public String comment() { return jTextFieldComment.getText(); }
    

    
    //
    // private member functions
    //

    /** check if the string contains a forbidden symbol */
    private boolean containsForbiddenSymbol(String s)
    {
	for (char symbol : forbiddenSymbols)
	    if (s.indexOf(symbol)>=0) return true;
	return false;
    }

    private void jTreeDirectoriesValueChanged(TreeSelectionEvent e)
    {
	Object o = jTreeDirectories.getLastSelectedPathComponent();
	if (o instanceof Directory) {
	    selectedDir = (Directory)o;
	    
	    // Next line was commented since the focus listener is introduced.
	    // raul.jimenez.estupinan@cern.ch
	    //if (selectedDir.name().split(" ").length>1) selectedDir=null;
	    // Instead of returning a null Directory it is returned using underscores. 
	    if(selectedDir.name().split(" ").length>1) selectedDir.setName(selectedDir.name().replace(" ", "_"));
	    
	    
	    updateOkButton();
	    
	}
	else if (o==null||(o instanceof ConfigInfo)) {
	    selectedDir = null;
	    jTreeDirectories.getSelectionModel().clearSelection();
	    jButtonOk.setEnabled(false);
	}
    }
    
    private void jTextFieldConfigNameActionPerformed(ActionEvent e)
    {
	if (jButtonOk.isEnabled()) jButtonOkActionPerformed(e);
    }
    
    private void jTextFieldConfigNameInsertUpdate(DocumentEvent e)
    {
	updateOkButton();
    }
    
    private void jTextFieldConfigNameRemoveUpdate(DocumentEvent e)
    {
	updateOkButton();
    }
    
    private void jButtonOkActionPerformed(ActionEvent e)
    {
	String    configName  = jTextFieldConfigName.getText();
	Directory parentDir   = selectedDir;
	String    releaseTag  = config.releaseTag();
	
	if (configName.length()>0&&parentDir!=null) {
	    ConfigInfo configInfo = new ConfigInfo(configName,parentDir,releaseTag);
	    config.setConfigInfo(configInfo);
	    validChoice = true;
	    setVisible(false);
	}
    }
    private void jButtonCancelActionPerformed(ActionEvent e)
    {
	validChoice = false;
	setVisible(false);
    }
    
    private void updateOkButton()
    {
	if (selectedDir==null) {
	    jButtonOk.setEnabled(false);
	    return;
	}
	String configName = jTextFieldConfigName.getText();
	if (configName.length()==0||containsForbiddenSymbol(configName)) {
	    jButtonOk.setEnabled(false);
	    return;
	}
	
	
	
	
	for (int i=0;i<selectedDir.configInfoCount();i++) {
	    if (selectedDir.configInfo(i).name().equals(configName)) {
		jButtonOk.setEnabled(false);
		return;
	    }
	}
	jButtonOk.setEnabled(true);
    }
    
    /** init GUI components [generated with NetBeans] */
    private JPanel createContentPane()
    {
	JPanel      contentPane  = new JPanel();
        JScrollPane jScrollPane1 = new JScrollPane();
        JLabel      jLabel1      = new JLabel();

        jScrollPane1.setViewportView(jTreeDirectories);

        jLabel1.setText("Configuration Name:");

        jButtonOk.setText("OK");
	jButtonOk.setEnabled(false);
        jButtonCancel.setText("Cancel");
	
        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(contentPane);
        contentPane.setLayout(layout);
        layout
	    .setHorizontalGroup(
				layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
				.add(layout.createSequentialGroup()
				     .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
					  .add(layout.createSequentialGroup()
					       .addContainerGap()
					       .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
						    .add(jScrollPane1,
							 org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 358, Short.MAX_VALUE)
						    .add(layout.createSequentialGroup()
							 .add(jLabel1)
							 .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
							 .add(jTextFieldConfigName,
							      org.jdesktop.layout.GroupLayout.DEFAULT_SIZE,
							      217, Short.MAX_VALUE))))
					  .add(layout.createSequentialGroup()
					       .add(91, 91, 91)
					       .add(jButtonOk, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 88, Short.MAX_VALUE)
					       .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
					       .add(jButtonCancel,
						    org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 86, Short.MAX_VALUE)
					       .add(99, 99, 99)))
				     .addContainerGap())
				);
        layout.setVerticalGroup(
				layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
				.add(layout.createSequentialGroup()
				     .addContainerGap()
				     .add(jScrollPane1,org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 387, Short.MAX_VALUE)
				     .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
				     .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
					  .add(jLabel1)
					  .add(jTextFieldConfigName,
					       org.jdesktop.layout.GroupLayout.PREFERRED_SIZE,
					       org.jdesktop.layout.GroupLayout.DEFAULT_SIZE,
					       org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
				     .add(19, 19, 19)
				     .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
					  .add(jButtonOk)
					  .add(jButtonCancel))
				     .addContainerGap())
				);
	
	return contentPane;
    }
    
}

