package confdb.gui;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.tree.*;
import java.awt.*;
import java.awt.event.*;

import java.util.ArrayList;
import java.util.HashMap;

import confdb.gui.menu.ScrollableMenu;

import confdb.data.*;

/**
 * ConfigurationTreeMouseListener
 * ------------------------------
 * @author Philipp Schieferdecker
 *
 */
public class ConfigurationTreeMouseListener extends    MouseAdapter 
                                            implements TreeModelListener
{
    //
    // member data
    //

    /** reference to the tree to be manipulated */
    private JTree tree = null;
    
    /** reference to the configuration tree-model */
    private ConfigurationTreeModel treeModel = null;
    
    /** list of available service templates */
    private ArrayList<Template> services = null;
    
    /** list of available event data source templates */
    private ArrayList<Template>  edsources= null;
    
    /** list of available event setup source templates */
    private ArrayList<Template> essources = null;
    
    /** list of available module templates */
    private ArrayList<Template> modules = null;

    /** popup menu showing all available service templates */
    private JPopupMenu popupServices = null;

    /** popup menu showing all available event data source templates */
    private JPopupMenu popupEDSources = null;

    /** popup menu showing all available event setup source templates */
    private JPopupMenu popupESSources = null;

    /** popup menu associated with paths */
    private JPopupMenu popupPaths = null;
    
    /** popup menu associated with modules */
    private JPopupMenu popupModules = null;
    
    /** popup menu associated with sequences */
    private JPopupMenu popupSequences = null;
    
    /** action listener for services-menu actions */
    private ServiceMenuListener serviceListener = null;
    
    /** action listener for edsources-menu actions */
    private EDSourceMenuListener edsourceListener = null;

    /** action listener for essources-menu actions */
    private ESSourceMenuListener essourceListener = null;
    
    /** action listener for paths-menu actions */
    private PathMenuListener pathListener = null;

    /** action listener for modules-menu actions */
    private ModuleMenuListener moduleListener = null;

    /** action listener for sequences-menu actions */
    private SequenceMenuListener sequenceListener = null;


    //
    // construction
    //
    
    /** standard constructor */
    public ConfigurationTreeMouseListener(JTree tree, JFrame frame,
					  ArrayList<Template> services,
					  ArrayList<Template> edsources,
					  ArrayList<Template> essources,
					  ArrayList<Template> modules)
    {
	this.tree      = tree;
	this.treeModel = (ConfigurationTreeModel)tree.getModel();
	this.services  = services;
	this.modules   = modules;
	this.edsources = edsources;
	this.essources = essources;
	
	serviceListener  = new ServiceMenuListener(tree);
	edsourceListener = new EDSourceMenuListener(tree);
	essourceListener = new ESSourceMenuListener(tree);
	pathListener     = new PathMenuListener(tree,frame);
	moduleListener   = new ModuleMenuListener(tree,frame);
	sequenceListener = new SequenceMenuListener(tree,frame);
    }
    
    
    //
    // member functions
    //
    
    /** MouseAdapter: mousePressed() */
    public void mousePressed(MouseEvent e)
    {
	maybeShowPopup(e);
    }
    
    /** MouseAdapter: mouseReleased() */
    public void mouseReleased(MouseEvent e)
    {
	maybeShowPopup(e);
    }
    
    /** check if this event should really trigger the menu to be displayed */
    private void maybeShowPopup(MouseEvent e)
    {
	if (!e.isPopupTrigger()) return;
	
	Configuration config = (Configuration)treeModel.getRoot();
	if (config.name().length()==0) return;
	
	tree = (JTree)e.getComponent();
	TreePath treePath = tree.getPathForLocation(e.getX(), e.getY());
	if(treePath==null) return;
	int      depth = treePath.getPathCount();
	if(depth<=1) return;
	
	tree.setSelectionPath(treePath);
	
	Object node = treePath.getLastPathComponent();
	
	// show the 'EDSources' popup?
	if (isInTreePath(treePath,treeModel.edsourcesNode())&&depth<=3) {
	    updateEDSourceMenu();
	    popupEDSources.show(e.getComponent(),e.getX(),e.getY());
	    return;
	}
	
	// show the 'ESSources' popup?
	if (isInTreePath(treePath,treeModel.essourcesNode())&&depth<=3) {
	    updateESSourceMenu();
	    popupESSources.show(e.getComponent(),e.getX(),e.getY());
	    return;
	}
	
	// show the 'Services' popup?
	if (isInTreePath(treePath,treeModel.servicesNode())&&depth<=3) {
	    updateServiceMenu();
	    popupServices.show(e.getComponent(),e.getX(),e.getY());
	    return;
	}
	
	// show the 'Paths' popup?
	if (isInTreePath(treePath,treeModel.pathsNode())&&depth<=4) {
	    updatePathMenu();
	    popupPaths.show(e.getComponent(),e.getX(),e.getY());
	    return;
	}
	
	// show the 'Modules' popup?
	if (node instanceof ModuleInstance) {
	    updateModuleMenu();
	    popupModules.show(e.getComponent(),e.getX(),e.getY());
	    return;
	}
	
	// show the 'Sequences' popup?
	if (isInTreePath(treePath,treeModel.sequencesNode())&&depth<=4) {
	    updateSequenceMenu();
	    popupSequences.show(e.getComponent(),e.getX(),e.getY());
	    return;
	}
	
    }
    
    /** check if a node is in a tree path */
    private boolean isInTreePath(TreePath treePath,Object node)
    {
	for (int i=0;i<treePath.getPathCount();i++)
	    if (treePath.getPathComponent(i).equals(node)) return true;
	return false;
    }
    
    /** update 'EDSource' menu */
    public void updateEDSourceMenu()
    {
	JMenuItem menuItem;
	popupEDSources = new JPopupMenu();
	
	TreePath treePath = tree.getSelectionPath();
	int      depth    = treePath.getPathCount();
	Object   node     = treePath.getLastPathComponent();
	
	JMenu edsourceMenu = new JMenu("Add EDSource");
	if (depth==3||
	    (depth==2&&treeModel.getChildCount(node)>0)) {
	    edsourceMenu.setEnabled(false);
	    popupEDSources.add(edsourceMenu);
	}
	else {
	    for (Template t : edsources) {
		menuItem = new JMenuItem(t.name());
		menuItem.addActionListener(edsourceListener);
		edsourceMenu.add(menuItem);
	    }
	    popupEDSources.add(edsourceMenu);
	}
	if (depth==3) {
	    menuItem = new JMenuItem("Remove EDSource");
	    menuItem.addActionListener(edsourceListener);
	    popupEDSources.add(menuItem);
	}
    }
    
    /** update 'ESSources' menu */
    public void updateESSourceMenu()
    {
	JMenuItem menuItem;
	popupESSources    = new JPopupMenu();
	
	TreePath treePath = tree.getSelectionPath();
	int      depth    = treePath.getPathCount();
	
	// 'ESSources' or specific event setup source
	if (depth>=2&&depth<=3) {
	    JMenu essourceMenu = new JMenu("Add ESSource");
	    for (Template t : essources) {
		if (t.instanceCount()>0) {
		    JMenu instanceMenu = new JMenu(t.name());
		    for (int i=0;i<t.instanceCount();i++) {
			Instance instance = t.instance(i);
			menuItem = new JMenuItem(instance.name());
			menuItem.setEnabled(false);
			instanceMenu.add(menuItem);
		    }
		    instanceMenu.addSeparator();
		    menuItem = new JMenuItem("New Instance");
		    menuItem.addActionListener(essourceListener);
		    menuItem.setActionCommand(t.name());
		    instanceMenu.add(menuItem);
		    essourceMenu.add(instanceMenu);
		}
		else {
		    menuItem = new JMenuItem(t.name());
		    menuItem.addActionListener(essourceListener);
		    essourceMenu.add(menuItem);
		}
	    }
	    popupESSources.add(essourceMenu);
	}
	
	// a specific event setup source is selected
	if (depth==3) {
	    menuItem = new JMenuItem("Remove ESSource");
	    menuItem.addActionListener(essourceListener);
	    popupESSources.add(menuItem);
	}
    }
    
    /** update 'Services' menu */
    public void updateServiceMenu()
    {
	JMenuItem menuItem;
	popupServices = new JPopupMenu();
	int depth = tree.getSelectionPath().getPathCount();
	
	JMenu serviceMenu = new JMenu("Add Service");
	for (Template t : services) {
	    menuItem = new JMenuItem(t.name());
	    if (t.instanceCount()>0) menuItem.setEnabled(false);
	    else menuItem.addActionListener(serviceListener);
	    serviceMenu.add(menuItem);
	}
	popupServices.add(serviceMenu);
	
	if (depth==3) {
	    menuItem = new JMenuItem("Remove Service");
	    menuItem.addActionListener(serviceListener);
	    popupServices.add(menuItem);
	}
	
    }
    
    /** update 'Paths' menu */
    public void updatePathMenu()
    {
	JMenuItem menuItem;
	popupPaths = new JPopupMenu();

	TreePath  treePath = tree.getSelectionPath();
	int       depth    = treePath.getPathCount();
	Object    node     = treePath.getPathComponent(depth-1);
	Object    parent   = treePath.getPathComponent(depth-2);

	Configuration config = (Configuration)treeModel.getRoot();
	
	// 'Paths' selectd
	if (depth==2) {
	    menuItem = new JMenuItem("Add Path");
	    menuItem.addActionListener(pathListener);
	    menuItem.setActionCommand("NEWPATH");
	    popupPaths.add(menuItem);
	    return;
	}
	
	// specific path is selected
	if (depth==3) {
	    Path path = (Path)node;

	    JMenu addModuleMenu = createAddModuleMenu(path,pathListener);
	    popupPaths.add(addModuleMenu);
	    
	    JMenu addPathMenu = createAddPathMenu(path);
	    popupPaths.add(addPathMenu);
	    
	    JMenu addSequenceMenu = createAddSequenceMenu(path);
	    popupPaths.add(addSequenceMenu);
	    
	    popupPaths.addSeparator();
	    
	    menuItem = new JMenuItem("Rename Path");
	    menuItem.addActionListener(pathListener);
	    popupPaths.add(menuItem);
	    
	    menuItem = new JMenuItem("Remove Path");
	    menuItem.addActionListener(pathListener);
	    popupPaths.add(menuItem);

	    return;
	}
	
	// a specific module/path/sequence reference is selected
	if (depth==4) {
	    Path path = (Path)parent;
	    JMenu addModuleMenu = createAddModuleMenu(path,pathListener);
	    popupPaths.add(addModuleMenu);
	    
	    JMenu addPathMenu = createAddPathMenu(path);
	    popupPaths.add(addPathMenu);
		
	    if (config.sequenceCount()>0) {
		JMenu addSequenceMenu = createAddSequenceMenu(path);
		popupPaths.add(addSequenceMenu);
	    }
	    
	    popupPaths.addSeparator();
	    
	    if (node instanceof ModuleReference) {
		menuItem = new JMenuItem("Remove Module");
		menuItem.addActionListener(pathListener);
		popupPaths.add(menuItem);
	    }
	    if (node instanceof PathReference) {
		menuItem = new JMenuItem("Remove Path");
		menuItem.addActionListener(pathListener);
		popupPaths.add(menuItem);
	    }
	    if (node instanceof SequenceReference) {
		menuItem = new JMenuItem("Remove Sequence");
		menuItem.addActionListener(pathListener);
		popupPaths.add(menuItem);
	    }
	}
    }
    
    /** update 'Modules' Menu */
    public void updateModuleMenu()
    {
	popupModules = new JPopupMenu();
	JMenuItem menuItem = new JMenuItem("Rename Module");
	menuItem.addActionListener(moduleListener);
	popupModules.add(menuItem);
    }
    
    /** update 'Sequences' Menu */
    public void updateSequenceMenu()
    {
	JMenuItem menuItem;
	popupSequences = new JPopupMenu();
	 
	TreePath treePath = tree.getSelectionPath();
	int      depth    = treePath.getPathCount();
	Object   node     = treePath.getPathComponent(depth-1);
	Object   parent   = treePath.getPathComponent(depth-2);
	
	if (depth==2) {
	    menuItem = new JMenuItem("Add Sequence");
	    menuItem.addActionListener(sequenceListener);
	    popupSequences.add(menuItem);
	}
	else if (depth==3) {
	    Sequence sequence = (Sequence)node;

	    JMenu addModuleMenu = createAddModuleMenu(sequence,sequenceListener);
	    popupSequences.add(addModuleMenu);

	    menuItem = new JMenuItem("Add Sequence");
	    menuItem.addActionListener(sequenceListener);
	    popupSequences.add(menuItem);
	    
	    popupSequences.addSeparator();
	    
	    menuItem = new JMenuItem("Rename Sequence");
	    menuItem.addActionListener(sequenceListener);
	    popupSequences.add(menuItem);
	    
	    menuItem = new JMenuItem("Remove Sequence");
	    menuItem.addActionListener(sequenceListener);
	    popupSequences.add(menuItem);
	}
	else if (depth==4) {
	    Sequence sequence = (Sequence)parent;

	    JMenu addModuleMenu = createAddModuleMenu(sequence,sequenceListener);
	    popupSequences.add(addModuleMenu);

	    menuItem = new JMenuItem("Remove Module");
	    menuItem.addActionListener(sequenceListener);
	    popupSequences.add(menuItem);
	}
    }
    
    /** create the 'Add Module' submenu */
    private JMenu createAddModuleMenu(ReferenceContainer container,
				      ActionListener     listener)
    {
	JMenu     addModuleMenu = new JMenu("Add Module");
	JMenuItem menuItem;
	
	HashMap<String,JMenu> menuHashMap = new HashMap<String,JMenu>();
	
	for (Template t : modules) {
	    
	    String moduleType = t.type();
	    if (!menuHashMap.containsKey(moduleType)) {
		JMenu menu = new ScrollableMenu(moduleType);
		menuHashMap.put(moduleType,menu);
		addModuleMenu.add(menu);
	    }
	    
	    if (t.instanceCount()>0) {
		JMenu instanceMenu = new JMenu(t.name());
		for (int i=0;i<t.instanceCount();i++) {
		    ModuleInstance instance = (ModuleInstance)t.instance(i);
		    menuItem = new JMenuItem(instance.name());
		    menuItem.addActionListener(listener);
		    menuItem.setActionCommand(t.name());
		    
		    for (int j=0;j<container.entryCount();j++) {
			Reference reference = container.entry(j);
			if (instance.isReferencedBy(reference))
			    menuItem.setEnabled(false);
		    }
		    instanceMenu.add(menuItem);
		}
		instanceMenu.addSeparator();
		menuItem = new JMenuItem("New Instance");
		menuItem.addActionListener(listener);
		menuItem.setActionCommand(t.name());
		instanceMenu.add(menuItem);
		menuHashMap.get(moduleType).add(instanceMenu);
	    }
	    else {
		menuItem = new JMenuItem(t.name());
		menuItem.addActionListener(listener);
		menuHashMap.get(moduleType).add(menuItem);
	    }
	}
	return addModuleMenu;
    }

    
    /**create 'Add Path' submenu */
    private JMenu createAddPathMenu(Object node)
    {
	JMenuItem menuItem;
	JMenu     addPathMenu = new JMenu("Add Path");
	
	ArrayList<Path> forbiddenPaths = new ArrayList<Path>();
	
	if (node instanceof Path) {
	    Path path = (Path)node;
	    forbiddenPaths.add(path);
	    menuItem = new JMenuItem("New Path");
	    menuItem.addActionListener(pathListener);
	    menuItem.setActionCommand("NEWPATH");
	    addPathMenu.add(menuItem);
	    Configuration config = (Configuration)treeModel.getRoot();
	    if (config.pathCount()>1) addPathMenu.addSeparator();
	}
	else if (node instanceof Reference) {
	    Reference reference  = (Reference)node;
	    Path      parentPath = (Path)reference.container();
	    forbiddenPaths.add(parentPath);
	    for (int i=0;i<parentPath.entryCount();i++) {
		Reference r = parentPath.entry(i);
		if (r instanceof PathReference) {
		    PathReference pathreference = (PathReference)r;
		    Path          path          = (Path)reference.parent();
		    forbiddenPaths.add(path);
		}
	    }
	}
	Configuration config = (Configuration)treeModel.getRoot();
	for (int i=0;i<config.pathCount();i++) {
	    Path path = config.path(i);
	    menuItem = new JMenuItem(path.name());
	    menuItem.addActionListener(pathListener);
	    menuItem.setActionCommand("PATHREF");
	    if (forbiddenPaths.contains(path)) menuItem.setEnabled(false);
	    addPathMenu.add(menuItem);
	}
	return addPathMenu;
    }
    
    /** create 'Add Sequence' Menu */
    private JMenu createAddSequenceMenu(Path path)
    {
	JMenu     addSequenceMenu = new JMenu("Add Sequence");
	JMenuItem menuItem;
	
	ArrayList<Sequence> forbiddenSequences = new ArrayList<Sequence>();
	for (int i=0;i<path.entryCount();i++) {
	    Reference reference = path.entry(i);
	    if (reference instanceof SequenceReference) {
		SequenceReference seqreference = (SequenceReference)reference;
		Sequence          sequence     = (Sequence)seqreference.parent();
		forbiddenSequences.add(sequence);
	    }
	}
	
	Configuration config = (Configuration)treeModel.getRoot();
	for (int i=0;i<config.sequenceCount();i++) {
	    Sequence sequence = config.sequence(i);
	    menuItem = new JMenuItem(sequence.name());
	    menuItem.addActionListener(pathListener);
	    menuItem.setActionCommand("SEQREF");
	    if (forbiddenSequences.contains(sequence)) menuItem.setEnabled(false);
	    addSequenceMenu.add(menuItem); 
	}
	return addSequenceMenu;
    }


    /** TreeModelListener: treeNodesChanged() */
    public void treeNodesChanged(TreeModelEvent e)
    {
	TreePath treePath = e.getTreePath(); if (treePath==null) return;
	int      depth    = treePath.getPathCount(); if (depth<2) return;
	int      index    = e.getChildIndices()[0];
	Object   child    = e.getChildren()[0];
	Object   parent   = treePath.getLastPathComponent();
	
	Configuration      config = (Configuration)treeModel.getRoot();
	ReferenceContainer container = null;
	
	if (config==null) return;
	
	if (parent.equals(treeModel.pathsNode())) {
	    Path newPath = config.path(index);
	    for (int i=0;i<config.pathCount();i++) {
		if (i==index) continue;
		Path path = config.path(i);
		if (0==path.name().compareToIgnoreCase(newPath.name()))
		    container = newPath;
	    }
	}
	else if (parent.equals(treeModel.sequencesNode())) {
	    Sequence newSequence = config.sequence(index);
	    for (int i=0;i<config.sequenceCount();i++) {
		if (i==index) continue;
		Sequence sequence = config.sequence(i);
		if (0==sequence.name().compareToIgnoreCase(newSequence.name()))
		    container = newSequence;
	    }
	}
	if (container!=null) {
	    TreePath containerTreePath = treePath.pathByAddingChild(child);
	    container.setName("<ENTER UNIQUE NAME>");
	    tree.startEditingAtPath(containerTreePath);
	}
    }

    /** TreeModelListener: treeNodesInserted() */
    public void treeNodesInserted(TreeModelEvent e) {}
    
    /** TreeModelListener: treeNodesRemoved() */
    public void treeNodesRemoved(TreeModelEvent e) {}

    /** TreeModelListener: treeStructureChanged() */
    public void treeStructureChanged(TreeModelEvent e) {}

    
}


/**
 * listen to actions from the 'EDSource' popup menu
 */
class EDSourceMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;
    
    /** reference to the tree model */
    private ConfigurationTreeModel treeModel = null;
    
    /** standard constructor */
    public EDSourceMenuListener(JTree tree)
    {
	this.tree = tree;
	this.treeModel = (ConfigurationTreeModel)tree.getModel();
    }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	TreePath  treePath = tree.getSelectionPath();
	Object    node     = treePath.getLastPathComponent();
	
	Configuration config = (Configuration)treeModel.getRoot();

	if (cmd.equals("Remove EDSource")) {
	    EDSourceInstance edsource       = (EDSourceInstance)node;
	    TreePath         parentTreePath = treePath.getParentPath();
	    int              index          = config.indexOfEDSource(edsource);
	    config.removeEDSource(edsource);
	    treeModel.nodeRemoved(treeModel.edsourcesNode(),index,edsource);
	    tree.setSelectionPath(parentTreePath);
	}
	else {
	    String templateName = cmd;
	    EDSourceInstance edsource = config.insertEDSource(templateName);
	    treeModel.nodeInserted(treeModel.edsourcesNode(),0);
	    tree.expandPath(treePath);
	}
	treeModel.updateLevel1Nodes();
    }
    
}


/**
 * listen to actions from the 'ESSources' popup menu
 */
class ESSourceMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;

    /** reference to the tree model */
    private ConfigurationTreeModel treeModel = null;
    
    /** standard constructor */
    public ESSourceMenuListener(JTree tree)
    {
	this.tree = tree;
	this.treeModel = (ConfigurationTreeModel)tree.getModel();
    }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	TreePath  treePath = tree.getSelectionPath();
       	int       depth    = treePath.getPathCount();
	Object    node     = treePath.getPathComponent(depth-1);
	Object    parent   = treePath.getPathComponent(depth-2);
	
	Configuration config = (Configuration)treeModel.getRoot();
	
	if (cmd.equals("Remove ESSource")) {
	    ESSourceInstance essource = (ESSourceInstance)node;
	    int              index    = config.indexOfESSource(essource);
	    config.removeESSource(essource);
	    treeModel.nodeRemoved(treeModel.essourcesNode(),index,essource);
	    tree.setSelectionPath(treePath.getParentPath());
 	}
	// add an event setup source TODO
	else if (depth>=2&&depth<=3) {
	    int insertIndex = 0;
	    if (depth==3) {
		ESSourceInstance essource = (ESSourceInstance)node;
		insertIndex = config.indexOfESSource(essource) + 1;
	    }
	    if (cmd.equals("New Instance")) {
		// TODO: dialog
		String templateName = source.getActionCommand();
		config.insertESSource(insertIndex,templateName,templateName);
	    }
	    else {
		String templateName = cmd;
		config.insertESSource(insertIndex,templateName,templateName);
	    }
	    treeModel.nodeInserted(treeModel.essourcesNode(),insertIndex);
	    if (insertIndex==0) tree.expandPath(treePath);
	    TreePath parentTreePath = (depth==2) ? treePath : treePath.getParentPath();
	    TreePath newTreePath=parentTreePath.pathByAddingChild(config.essource(insertIndex));
	    tree.setSelectionPath(newTreePath);
	}
	treeModel.updateLevel1Nodes();
    }
}


/**
 * listen to actions from the 'Services' popup menu
 */
class ServiceMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;
    
    /** reference to the tree model */
    private ConfigurationTreeModel treeModel = null;
    
    /** standard constructor */
    public ServiceMenuListener(JTree tree)
    {
	this.tree = tree;
	this.treeModel = (ConfigurationTreeModel)tree.getModel();
    }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	TreePath  treePath = tree.getSelectionPath();
	int       depth    = treePath.getPathCount();
	Object    node     = treePath.getPathComponent(depth-1);
	Object    parent   = treePath.getPathComponent(depth-2);

	Configuration config = (Configuration)treeModel.getRoot();
	
	if (cmd.equals("Remove Service")) {
	    ServiceInstance service = (ServiceInstance)node;	
	    int             index   = config.indexOfService(service);
	    config.removeService(service);
	    treeModel.nodeRemoved(treeModel.servicesNode(),index,service);
	    tree.setSelectionPath(treePath.getParentPath());
	}
	else {
	    int insertIndex = 0;
	    if (depth==3) {
		ServiceInstance service = (ServiceInstance)node;
		insertIndex = config.indexOfService(service) + 1;
	    }
	    String templateName = cmd;
	    config.insertService(insertIndex,templateName);
	    treeModel.nodeInserted(treeModel.servicesNode(),insertIndex);
	    if (insertIndex==0) tree.expandPath(treePath);
	    TreePath parentTreePath = (depth==2) ? treePath : treePath.getParentPath();
	    TreePath newTreePath = parentTreePath.pathByAddingChild(config.service(insertIndex));
	    tree.setSelectionPath(newTreePath);
	}
	treeModel.updateLevel1Nodes();
    }
    
}


/**
 * listen to actions from the 'Paths' popup menu
 */
class PathMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;

    /** reference to the tree model */
    private ConfigurationTreeModel treeModel = null;

    /** application frame */
    private JFrame frame = null;
    
    /** standard constructor */
    public PathMenuListener(JTree tree,JFrame frame)
    {
	this.tree      = tree;
	this.treeModel = (ConfigurationTreeModel)tree.getModel();
	this.frame     = frame;
    }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	String    action   = source.getActionCommand();
	TreePath  treePath = tree.getSelectionPath();
	int       depth    = treePath.getPathCount();
	Object    node     = treePath.getPathComponent(depth-1);
	Object    parent   = treePath.getPathComponent(depth-2);
	
	Configuration config = (Configuration)treeModel.getRoot();

	if (action.equals("NEWPATH")) {
	    int insertIndex = 0;
	    if (depth>2) {
		if (depth==4) node = parent;
		Path path   = (Path)node;
		insertIndex = config.indexOfPath(path)+1;
	    }
	    config.insertPath(insertIndex,"<ENTER PATH NAME>");
	    treeModel.nodeInserted(treeModel.pathsNode(),insertIndex);
	    if (insertIndex==0) tree.expandPath(treePath);
	    treeModel.updateLevel1Nodes();
	    TreePath parentTreePath = treePath;
	    if (depth>2) parentTreePath = parentTreePath.getParentPath();
	    if (depth>3) parentTreePath = parentTreePath.getParentPath();
	    TreePath newTreePath =
		parentTreePath.pathByAddingChild(config.path(insertIndex));
	    tree.expandPath(newTreePath);
	    editPathName(newTreePath);
	    return;
	}
	else if (action.equals("PATHREF")) {
	    TreePath parentTreePath = (depth==3) ? treePath:treePath.getParentPath();
	    Path     parentPath  = (depth==3) ? (Path)node : (Path)parent;
	    int      insertIndex = (depth==3) ? 0 : parentPath.indexOfEntry((Reference)node)+1;
	    String   pathName    = cmd;
	    for (int i=0;i<config.pathCount();i++) {
		Path path = config.path(i);
		if (path.name().equals(pathName)) {
		    config.insertPathReference(parentPath,insertIndex,path);
		    treeModel.nodeInserted(parentPath,insertIndex);
		    treeModel.updateLevel1Nodes();
		    return;
		}
	    }
	}
	else if (action.equals("SEQREF")) {
	    TreePath parentTreePath = (depth==3) ? treePath : treePath.getParentPath();
	    Path     parentPath  = (depth==3) ? (Path)node : (Path)parent;
	    int      insertIndex = (depth==3) ? 0 : parentPath.indexOfEntry((Reference)node)+1;
	    String   sequenceName= cmd;
	    for (int i=0;i<config.sequenceCount();i++) {
		Sequence sequence = config.sequence(i);
		if (sequence.name().equals(sequenceName)) {
		    config.insertSequenceReference(parentPath,insertIndex,sequence);
		    treeModel.nodeInserted(parentPath,insertIndex);
		    treeModel.updateLevel1Nodes();
		    return;
		}
	    }
	}
	else if (cmd.equals("Rename Path")) {
	    editPathName(treePath);
	    return;
	}
	else if (cmd.equals("Remove Path")) {
	    if (depth==3) {
		Path path = (Path)node;
		int  index = config.indexOfPath(path);
		config.removePath(path);
		treeModel.nodeRemoved(treeModel.pathsNode(),index,path);
	    }
	    else {
		PathReference reference  = (PathReference)node;
		Path          parentPath = (Path)parent;
		int           index = parentPath.indexOfEntry(reference);
		parentPath.removeEntry(reference);
		treeModel.nodeRemoved(parentPath,index,reference);
	    }
	    tree.setSelectionPath(treePath.getParentPath());
	}
	else if (cmd.equals("Rename Module")) {
	    ModuleReference reference = (ModuleReference)node;
	    ModuleInstance  instance  = (ModuleInstance)reference.parent();
	    ModuleTemplate  template  = (ModuleTemplate)instance.template();
	    ValidatedNameDialog dialog = new ValidatedNameDialog(frame,template);
	    dialog.pack();
	    dialog.setLocationRelativeTo(frame);
	    dialog.setVisible(true);
	    if (dialog.success()) {
		instance.setName(dialog.instanceName());
		config.setHasChanged(true);
	    }
	}
	else if (cmd.equals("Remove Module")) {
	    ModuleReference reference  = (ModuleReference)node;
	    Path            parentPath = (Path)reference.container();
	    int             index      = parentPath.indexOfEntry(reference);
	    config.removeModuleReference(reference);
	    treeModel.nodeRemoved(parentPath,index,reference);
	    tree.setSelectionPath(treePath.getParentPath());
	}
	else if (cmd.equals("Remove Sequence")) {
	    SequenceReference reference  = (SequenceReference)node;
	    Path              parentPath = (Path)parent;
	    int               index      = parentPath.indexOfEntry(reference);
	    parentPath.removeEntry(reference);
	    config.setHasChanged(true);
	    treeModel.nodeRemoved(parentPath,index,reference);
	    tree.setSelectionPath(treePath.getParentPath());
	}
 	// add a module to the currently selected path
	else if (depth>=3&&depth<=4) {
	    TreePath parentTreePath = (depth==3) ? treePath : treePath.getParentPath();
	    Path     parentPath   = (depth==3) ? (Path)node : (Path)parent;
	    int      insertIndex  = (depth==3) ? 0 : parentPath.indexOfEntry((Reference)node)+1;
	    String   sequenceName = cmd;
	    String   templateName = source.getActionCommand();
	    String   instanceName = cmd;
	    if (templateName.length()==0) templateName = cmd;
	    if (instanceName.equals("New Instance")) {
		Template template = config.moduleTemplate(templateName);
		ValidatedNameDialog dialog = new ValidatedNameDialog(frame,template);
		dialog.pack();
		dialog.setLocationRelativeTo(frame);
		dialog.setVisible(true);
		if (dialog.success()) {
		    instanceName = dialog.instanceName();
		    ArrayList<Parameter> parameters = dialog.instanceParameters();
		    ModuleReference reference =
			config.insertModuleReference(parentPath,insertIndex,
						     templateName,instanceName);
		    ModuleInstance instance = (ModuleInstance)reference.parent();
		    instance.setParameters(parameters);
		}
	    }
	    else {
		config.insertModuleReference(parentPath,insertIndex,templateName,instanceName);
	    }
	    treeModel.nodeInserted(parentPath,insertIndex);
	}
	treeModel.updateLevel1Nodes();
    }
    
    /** edit a path name in place */
    private void editPathName(TreePath treePath)
    {
	tree.expandPath(treePath);
	tree.scrollPathToVisible(treePath);
	tree.setSelectionPath(treePath);
	tree.startEditingAtPath(treePath);
    }
    
}

/**
 * listen to actions from the 'Modules' popup menu
 */
class ModuleMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;
    
    /** tree model */
    private ConfigurationTreeModel treeModel = null;
    
    /** the application frame */
    private JFrame frame = null;
    
    /** standard constructor */
    public ModuleMenuListener(JTree tree,JFrame frame)
    {
	this.tree  = tree;
	this.treeModel = (ConfigurationTreeModel)tree.getModel();
	this.frame = frame;
    }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	TreePath  treePath = tree.getSelectionPath();
	Object    node     = treePath.getLastPathComponent();
	
	Configuration config = (Configuration)treeModel.getRoot();
	
	if (cmd.equals("Rename Module")) {
	    ModuleInstance  instance   = (ModuleInstance)node;
	    ModuleTemplate  template   = (ModuleTemplate)instance.template();
	    ValidatedNameDialog dialog = new ValidatedNameDialog(frame,template);
	    dialog.pack();
	    dialog.setLocationRelativeTo(frame);
	    dialog.setVisible(true);
	    if (dialog.success()) {
		instance.setName(dialog.instanceName());
		treeModel.nodeChanged(instance);
		for (int i=0;i<instance.referenceCount();i++) {
		    Reference r= instance.reference(i);
		    treeModel.nodeChanged(r);
		}
		config.setHasChanged(true);
	    }
 	}
    }
}


/**
 * listen to actions from the 'Sequences' popup menu
 */
class SequenceMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;

    /** tree model */
    private ConfigurationTreeModel treeModel = null;
    
    /** application frame */
    private JFrame frame = null;

    /** standard constructor */
    public SequenceMenuListener(JTree tree,JFrame frame)
    {
	this.tree = tree;
	this.treeModel = (ConfigurationTreeModel)tree.getModel();
	this.frame = frame;
    }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	TreePath  treePath = tree.getSelectionPath();
	int       depth    = treePath.getPathCount();
	Object    node     = treePath.getPathComponent(depth-1);
	Object    parent   = treePath.getPathComponent(depth-2);
	
	Configuration config = (Configuration)treeModel.getRoot();
	
	if (cmd.equals("Add Sequence")) {
	    int insertIndex = 0;
	    if (depth==3) {
		Sequence seq = (Sequence)node;
		insertIndex = config.indexOfSequence(seq)+1;
	    }
	    config.insertSequence(insertIndex,"<ENTER SEQUENCE NAME>");
	    treeModel.nodeInserted(treeModel.sequencesNode(),insertIndex);
	    if (insertIndex==0) tree.expandPath(treePath);
	    treeModel.updateLevel1Nodes();
	    TreePath parentTreePath = treePath;
	    if (depth>2) parentTreePath = parentTreePath.getParentPath();
	    TreePath newTreePath = parentTreePath.pathByAddingChild(config.sequence(insertIndex));
	    editSequenceName(newTreePath);
	    return;
	}
	else if (cmd.equals("Rename Sequence")) {
	    editSequenceName(treePath);
	    return;
	}
	else if (cmd.equals("Remove Sequence")) {
	    Sequence sequence = (Sequence)node;
	    int      index    = config.indexOfSequence(sequence);
	    config.removeSequence(sequence);
	    treeModel.nodeRemoved(parent,index,sequence);
	    tree.setSelectionPath(treePath.getParentPath());
 	}
	else if (cmd.equals("Remove Module")) {
	    ModuleReference reference = (ModuleReference)node;
	    Sequence        sequence  = (Sequence)reference.container();
	    int             index     = sequence.indexOfEntry(reference);
	    config.removeModuleReference(reference);
	    treeModel.nodeRemoved(sequence,index,reference);
	    tree.setSelectionPath(treePath.getParentPath());
	}
	else {
	    TreePath parentTreePath = (depth==3) ? treePath : treePath.getParentPath();
	    Sequence sequence     = (depth==3) ? (Sequence)node : (Sequence)parent;
	    int      insertIndex  = (depth==3) ? 0 : sequence.indexOfEntry((Reference)node)+1;
	    String   templateName = source.getActionCommand();
	    String   instanceName = cmd;
	    if (templateName.length()==0) templateName = cmd;
	    if (instanceName.equals("New Instance")) {
		Template template = config.moduleTemplate(templateName);
		ValidatedNameDialog dialog = new ValidatedNameDialog(frame,template);
		dialog.pack();
		dialog.setLocationRelativeTo(frame);
		dialog.setVisible(true);
		if (dialog.success()) {
		    instanceName = dialog.instanceName();
		    ArrayList<Parameter> parameters = dialog.instanceParameters();
		    ModuleReference reference =
			config.insertModuleReference(sequence,insertIndex,
						     templateName,instanceName);
		    ModuleInstance instance = (ModuleInstance)reference.parent();
		    instance.setParameters(parameters);
		}
	    }
	    else {
		config.insertModuleReference(sequence,insertIndex,templateName,instanceName);
	    }
	    treeModel.nodeInserted(sequence,insertIndex);
	}
	treeModel.updateLevel1Nodes();
    }
    
    /** edit a sequence name in place */
    private void editSequenceName(TreePath treePath)
    {
	tree.expandPath(treePath);
	tree.scrollPathToVisible(treePath);
	tree.setSelectionPath(treePath);
	tree.startEditingAtPath(treePath);
    }
    
}
