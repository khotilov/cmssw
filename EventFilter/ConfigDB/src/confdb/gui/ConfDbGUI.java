package confdb.gui;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.tree.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;

import java.util.ArrayList;

import java.util.concurrent.ExecutionException;

import confdb.data.Configuration;
import confdb.data.ConfigInfo;
import confdb.data.Template;
import confdb.data.Parameter;

import confdb.db.CfgDatabase;
import confdb.db.DatabaseException;

import confdb.gui.treetable.TreeTableTableModel;


/**
 * ConfDbGUI
 * ---------
 * @author Philipp Schieferdecker
 *
 * Graphical User Interface to create and manipulate cmssw
 * configurations stored in the Configuration Database, ConfDB.
 */
public class ConfDbGUI implements TableModelListener
{
    //
    // member data
    //
    
    /** the current configuration */
    private Configuration config = null;
    
    /** handle access to the database */
    private CfgDatabase database = null;
    
    /** main frame of the application */
    private JFrame frame = null; 
    
    /** content pane of the application */
    private JPanel contentPane = null;

    /** database information panel */
    private DatabaseInfoPanel dbInfoPanel = null;
    
    /** panel to display the currently selected instance and its properties */
    private InstancePanel instancePanel = null;

    /** panel to display information about the current configuration */
    private ConfigurationPanel configurationPanel = null;

    /** progress bar for database operations */
    private JProgressBar progressBar = null;
    
    /** tree structure holding the current configuration */
    private JTree tree = null;
    private ConfigurationTreeModel treeModel = null;
    
    /** list of all available Event Data Source templates */
    private ArrayList<Template> edsourceTemplateList = null;
    
    /** list of all available Event Setup Source templates */
    private ArrayList<Template> essourceTemplateList = null;
    
    /** list of all available Service templates */
    private ArrayList<Template> serviceTemplateList = null;
    
    /** list of all available Module templates */
    private ArrayList<Template> moduleTemplateList = null;
    

    //
    // construction
    //
    
    /** standard constructor */
    public ConfDbGUI(JFrame frame)
    {
	this.frame           = frame;
	config               = new Configuration();
	database             = new CfgDatabase();
	serviceTemplateList  = new ArrayList<Template>();
	edsourceTemplateList = new ArrayList<Template>();
	essourceTemplateList = new ArrayList<Template>();
	moduleTemplateList   = new ArrayList<Template>();
    }
    
    
    //
    // member functions
    //

    /** TableModelListener: tableChanged() */
    public void tableChanged(TableModelEvent e)
    {
	Object source = e.getSource();
	if (source instanceof TreeTableTableModel) {
	    TreeTableTableModel tableModel = (TreeTableTableModel)source;
	    Parameter node = (Parameter)tableModel.changedNode();
	    if (node!=null) {
		Object parent = node.parent();
		while (parent instanceof Parameter) {
		    Parameter p = (Parameter)parent;
		    parent = p.parent();
		}
		treeModel.nodeChanged(parent);
		treeModel.nodeChanged(node);
		treeModel.updateLevel1Nodes();
	    }
	}
    }
    
    /** connect to the database */
    public void connectToDatabase()
    {
	// query the database info from the user in a dialog
	DatabaseConnectionDialog dbDialog = new DatabaseConnectionDialog(frame);
	dbDialog.pack();
	dbDialog.setLocationRelativeTo(frame);
	dbDialog.setVisible(true);
	
	// retrieve the database parameters from the dialog
	if (!dbDialog.validChoice()) return;
	String dbType = dbDialog.getDbType();
	String dbHost = dbDialog.getDbHost();
	String dbPort = dbDialog.getDbPort();
	String dbName = dbDialog.getDbName();
	String dbUrl  = dbDialog.getDbUrl();
	String dbUser = dbDialog.getDbUser();
	String dbPwrd = dbDialog.getDbPassword();
	
	try {
	    database.connect(dbType,dbUrl,dbUser,dbPwrd);
	    dbInfoPanel.connectedToDatabase(dbType,dbHost,dbPort,dbName,dbUser);
	}
	catch (DatabaseException e) {
	    String msg = "Failed to connect to DB: " + e.getMessage();
	    JOptionPane.showMessageDialog(frame,msg,"",JOptionPane.ERROR_MESSAGE);
	}
    }
    
    // connect to the database
    public void disconnectFromDatabase()
    {
	try {
	    database.disconnect();
	    dbInfoPanel.disconnectedFromDatabase();
	    serviceTemplateList.clear();
	    edsourceTemplateList.clear();
	    essourceTemplateList.clear();
	    moduleTemplateList.clear();
	}
	catch (DatabaseException e) {
	    String msg = "Failed to disconnect from DB: " + e.getMessage();
	    JOptionPane.showMessageDialog(frame,msg,"",JOptionPane.ERROR_MESSAGE);
	}
    }
    
    /** new configuration */
    public void newConfiguration()
    {
	closeConfiguration();
	ConfigurationNameDialog dialog = new ConfigurationNameDialog(frame,database);
	dialog.pack();
	dialog.setLocationRelativeTo(frame);
	dialog.setVisible(true);
	
	if (dialog.validChoice()) {
	    String cfgName    = dialog.name();
	    String releaseTag = dialog.releaseTag();
	    
	    NewConfigurationThread worker =
		new NewConfigurationThread(cfgName,releaseTag);
	    worker.start();
	    progressBar.setIndeterminate(true);
	    progressBar.setVisible(true);
	    progressBar.setString("Loading Templates for Release " +
				  dialog.releaseTag() + " ... ");
	}
    }
    
    /** open configuration */
    public void openConfiguration()
    {
	closeConfiguration();
	
	OpenConfigurationDialog dialog =
	    new OpenConfigurationDialog(frame,database);
	dialog.pack();
	dialog.setLocationRelativeTo(frame);
	dialog.setVisible(true);
	
	if (dialog.validChoice()) {
	    OpenConfigurationThread worker =
		new OpenConfigurationThread(dialog.configInfo());
	    worker.start();
	    progressBar.setIndeterminate(true);
	    progressBar.setVisible(true);
	    progressBar.setString("Loading Configuration ...");
	}
    }

    /** close configuration */
    public void closeConfiguration()
    {
	config.reset();
	treeModel.setConfiguration(config);
	configurationPanel.update(config);
	instancePanel.clear();
    }
    
    /** save configuration */
    public void saveConfiguration()
    {
	if (!config.hasChanged()) return;
	if (config.version()==0) { saveAsConfiguration(); return; }
	if (!checkConfiguration()) return;	
	
	SaveConfigurationThread worker =
	    new SaveConfigurationThread();
	worker.start();
	progressBar.setIndeterminate(true);
	progressBar.setString("Save Configuration ...");
	progressBar.setVisible(true);
    }
    
    /** saveAs configuration */
    public void saveAsConfiguration()
    {
	if (!checkConfiguration()) return;

	SaveConfigurationDialog dialog =
	    new SaveConfigurationDialog(frame,database,config);
	dialog.pack();
	dialog.setLocationRelativeTo(frame);
	dialog.setVisible(true);
	
	if (dialog.validChoice()) {
	    SaveConfigurationThread worker =
		new SaveConfigurationThread();
	    worker.start();
	    progressBar.setIndeterminate(true);
	    progressBar.setString("Save Configuration ...");
	    progressBar.setVisible(true);
	}
    }
    
    /** check if configuration is in a storable state */
    public boolean checkConfiguration()
    {
	int unsetParamCount = config.unsetTrackedParameterCount();
	if (unsetParamCount>0) {
	    String msg =
		"The configuration contains " + unsetParamCount +
		" unset tracked parameters. They must be set before saving!";
	    JOptionPane.showMessageDialog(frame,msg,"",
					  JOptionPane.ERROR_MESSAGE);
	    return false;
	}
	return true;
    }

    /** show 'create templates' dialog */
    public void createTemplates()
    {
	CreateTemplateDialog dialog = new CreateTemplateDialog(frame,database);
	dialog.pack();
	dialog.setLocationRelativeTo(frame);
	dialog.setVisible(true);
	String releaseTag = config.releaseTag();
	
	UpdateTemplatesThread worker =
	    new UpdateTemplatesThread(releaseTag);
	worker.start();
	progressBar.setIndeterminate(true);
	progressBar.setVisible(true);
	progressBar.setString("Updating Templates for Release "+releaseTag+" ... ");
    }

    /** create the content pane */
    private JPanel createContentPane()
    {
	// the contentPane of the main frame
	JPanel contentPane = new JPanel(new GridBagLayout());
	contentPane.setOpaque(true);
	
	GridBagConstraints c = new GridBagConstraints();
	c.fill = GridBagConstraints.BOTH;
	c.weightx = 0.5;
	
	c.gridx=0;c.gridy=0; c.weighty=0.01;
	contentPane.add(createDbInfoView(),c);
	
	// create the instance view
	JPanel instanceView = createInstanceView(new Dimension(500,800));
	
	// create the tree view
	JPanel treeView = createTreeView(new Dimension(500,800));
	
	// create the tree/component panel split pane
	JSplitPane  horizontalSplitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
							 treeView,instanceView);
	horizontalSplitPane.setOneTouchExpandable(true);
	horizontalSplitPane.setResizeWeight(0.01);
	horizontalSplitPane.setDividerLocation(0.5);
	
	// add horizontal split pane to content pane
	c.gridx=0;c.gridy=1; c.weighty=0.98;
	contentPane.add(horizontalSplitPane,c);
	
	// add the status bar at the bottom
	c.gridx=0;c.gridy=2; c.weighty=0.01;
	JPanel statusPanel = new JPanel(new GridLayout());
	progressBar = new JProgressBar(0);
	progressBar.setIndeterminate(true);
	progressBar.setStringPainted(true);
	progressBar.setVisible(false);
	statusPanel.add(progressBar);
	contentPane.add(statusPanel,c);
	
	
	// return content pane, to be set in main frame
	return contentPane;
    }

    /** craete the top panel with the database connection info */
    private JPanel createDbInfoView()
    {
	dbInfoPanel = new DatabaseInfoPanel();
	return dbInfoPanel;
    }
        
    /** create the component panel, to display the currently selected component */
    private JPanel createInstanceView(Dimension dim)
    {
	instancePanel = new InstancePanel(frame,dim);
	instancePanel.addTableModelListener(this);
	return instancePanel;
    }
    
    /** create the configuration tree */
    private JPanel createTreeView(Dimension dim)
    {
	configurationPanel = new ConfigurationPanel();
	treeModel = new ConfigurationTreeModel(config);
	tree      = new JTree(treeModel);
	tree.addTreeSelectionListener(instancePanel);

	ConfigurationTreeMouseListener mouseListener =
	    new ConfigurationTreeMouseListener(tree,frame,
					       serviceTemplateList,
					       edsourceTemplateList,
					       essourceTemplateList,
					       moduleTemplateList);
	tree.addMouseListener(mouseListener);
	treeModel.addTreeModelListener(mouseListener);
	
	tree.setRootVisible(false);
	tree.setEditable(true);
	tree.getSelectionModel()
	    .setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
	
	ConfigurationTreeRenderer renderer = new ConfigurationTreeRenderer();
	tree.setCellRenderer(renderer);
	tree.setCellEditor(new ConfigurationTreeEditor(tree,renderer));
	
	JPanel treeView = new JPanel(new GridBagLayout());
	treeView.setPreferredSize(dim);
	
	JScrollPane spTree = new JScrollPane(tree);
	//spTree.setPreferredSize(??);
	
	GridBagConstraints c = new GridBagConstraints();
	c.fill = GridBagConstraints.BOTH;
	c.weightx = 0.5;
	
	c.gridx=0;c.gridy=0;c.gridwidth=1; c.weighty = 0.01;
	treeView.add(configurationPanel,c);
	
	c.gridx=0;c.gridy=1;c.gridwidth=1; c.weighty=0.99;
	treeView.add(spTree,c);
	
	treeView.setPreferredSize(dim);
	return treeView;
    }
    
    /** create the menu bar */
    private void createMenuBar()
    {
	CfgMenuBar menuBar = new CfgMenuBar(frame,this);
    }
    
    //
    // static member functions
    //
    
    /** create and show the ConfDbGUI */
    private static void createAndShowGUI()
    {
	// create the application's main frame
	JFrame frame = new JFrame("ConfDbGUI");
	
	// create the ConfDbGUI app and set it as the main frame's content pane
	ConfDbGUI app = new ConfDbGUI(frame);
	
	// add the application's content pane to the main frame
	frame.setContentPane(app.createContentPane());
	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	
	// add menu bar to the main frame
	app.createMenuBar();
	
	// display the main frame
	int frameWidth  = (int)(0.75*frame.getToolkit().getScreenSize().getWidth());
	int frameHeight = (int)(0.75*frame.getToolkit().getScreenSize().getHeight());
	int frameX      = (int)(0.125*frame.getToolkit().getScreenSize().getWidth());
	int frameY      = (int)(0.10*frame.getToolkit().getScreenSize().getHeight());
	frame.pack();
	frame.setSize(frameWidth,frameHeight);
	frame.setLocation(frameX,frameY);
	frame.setVisible(true);
		
	// try to etablish a database connection
	app.connectToDatabase();

	// reset configurations
	app.closeConfiguration();
    }
    

    /** main, create and show GUI, thread-safe */
    public static void main(String[] args)
    {
	javax.swing.SwingUtilities
	    .invokeLater(new Runnable() 
		{
		    public void run()
		    {
			createAndShowGUI();
		    }
		});
    }
    

    //
    // threads *not* to be executed on the EDT
    //

    /** load release templates from the database */
    private class NewConfigurationThread extends SwingWorker<String>
    {
	/** name of the new configuration */
	private String cfgName = null;
	
	/** release to be loaded */
	private String releaseTag = null;
	
	/** start time */
	private long startTime;
	
	/** standard constructor */
	public NewConfigurationThread(String cfgName,String releaseTag)
	{
	    this.cfgName    = cfgName;
	    this.releaseTag = releaseTag;
	}
	
	/** SwingWorker: construct() */
	protected String construct() throws InterruptedException
	{
	    startTime = System.currentTimeMillis();
	    database.loadEDSourceTemplates(releaseTag,edsourceTemplateList);
	    database.loadESSourceTemplates(releaseTag,essourceTemplateList);
	    database.loadServiceTemplates(releaseTag,serviceTemplateList);
	    database.loadModuleTemplates(releaseTag,moduleTemplateList);
	    return new String("Done!");
	}
	
	/** SwingWorker: finished */
	protected void finished()
	{
	    try {
		config.initialize(new ConfigInfo(cfgName,null,releaseTag),
				  edsourceTemplateList,
				  essourceTemplateList,
				  serviceTemplateList,
				  moduleTemplateList);
		treeModel.setConfiguration(config);
		configurationPanel.update(config);
		long elapsedTime = System.currentTimeMillis() - startTime;
		progressBar.setString(progressBar.getString() +
				      get() + " (" + elapsedTime + " ms)");
	    }
	    catch (Exception e) {
		System.out.println("EXCEPTION: "+ e.getMessage());
		e.printStackTrace();
		progressBar.setString(progressBar.getString() + "FAILED!");	
	    }
	    progressBar.setIndeterminate(false);
	}
    }
    

    /** load a configuration from the database */
    private class OpenConfigurationThread extends SwingWorker<String>
    {
	/** configuration info */
	private ConfigInfo configInfo = null;
	
	/** start time */
	private long startTime;
	
	/** standard constructor */
	public OpenConfigurationThread(ConfigInfo configInfo)
	{
	    this.configInfo = configInfo;
	}
	
	/** SwingWorker: construct() */
	protected String construct() throws InterruptedException
	{
	    startTime = System.currentTimeMillis();

	    config = database.loadConfiguration(configInfo,
						edsourceTemplateList,
						essourceTemplateList,
						serviceTemplateList,
						moduleTemplateList);
	    return new String("Done!");
	}
	
	/** SwingWorker: finished */
	protected void finished()
	{
	    try {
		treeModel.setConfiguration(config);
		configurationPanel.update(config);
		long elapsedTime = System.currentTimeMillis() - startTime;
		progressBar.setString(progressBar.getString() +
				      get() + " (" + elapsedTime + " ms)");
	    }
	    catch (ExecutionException e) {
		System.out.println("EXECUTION-EXCEPTION: "+ e.getCause());
		e.printStackTrace();
	    }
	    catch (Exception e) {
		System.out.println("EXCEPTION: "+ e.getMessage());
		e.printStackTrace();
		progressBar.setString(progressBar.getString() + "FAILED!");
	    }
	    progressBar.setIndeterminate(false);
	}
    }

    /** save a configuration in the database */
    private class SaveConfigurationThread extends SwingWorker<String>
    {
	/** start time */
	private long startTime;
	
	/** standard constructor */
	public SaveConfigurationThread() {}
	
	/** SwingWorker: construct() */
	protected String construct() throws InterruptedException
	{
	    startTime = System.currentTimeMillis();
	    database.insertConfiguration(config);
	    return new String("Done!");
	}
	
	/** SwingWorker: finished */
	protected void finished()
	{
	    try {
		config.setHasChanged(false);
		configurationPanel.update(config);
		long elapsedTime = System.currentTimeMillis() - startTime;
		progressBar.setString(progressBar.getString() +
				      get() + " (" + elapsedTime + " ms)");
	    }
	    catch (Exception e) {
		System.out.println("EXCEPTION: "+ e.getMessage());
		e.printStackTrace();
		progressBar.setString(progressBar.getString() + "FAILED!");
	    }
	    progressBar.setIndeterminate(false);
	}
    }

    /** load release templates from the database */
    private class UpdateTemplatesThread extends SwingWorker<String>
    {
	/** release to be loaded */
	private String releaseTag = null;
	
	/** start time */
	private long startTime;
	
	/** standard constructor */
	public UpdateTemplatesThread(String releaseTag)
	{
	    this.releaseTag = releaseTag;
	}
	
	/** SwingWorker: construct() */
	protected String construct() throws InterruptedException
	{
	    startTime = System.currentTimeMillis();
	    database.loadEDSourceTemplates(releaseTag,edsourceTemplateList);
	    database.loadESSourceTemplates(releaseTag,essourceTemplateList);
	    database.loadServiceTemplates(releaseTag,serviceTemplateList);
	    database.loadModuleTemplates(releaseTag,moduleTemplateList);
	    return new String("Done!");
	}
	
	/** SwingWorker: finished */
	protected void finished()
	{
	    try {
		config.updateHashMaps(edsourceTemplateList,
				      essourceTemplateList,
				      serviceTemplateList,
				      moduleTemplateList);
		long elapsedTime = System.currentTimeMillis() - startTime;
		progressBar.setString(progressBar.getString() +
				      get() + " (" + elapsedTime + " ms)");
	    }
	    catch (Exception e) {
		System.out.println("EXCEPTION: "+ e.getMessage());
		e.printStackTrace();
		progressBar.setString(progressBar.getString() + "FAILED!");	
	    }
	    progressBar.setIndeterminate(false);
	}
    }
    
}
