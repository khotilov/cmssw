package confdb.gui;

import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.*;
import javax.swing.tree.*;
import java.awt.event.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import confdb.gui.menu.ScrollableMenu;

import confdb.data.*;
import java.util.Collections;


/**
 * ConfigurationTreeMouseListener
 * ------------------------------
 * @author Philipp Schieferdecker
 *
 */
public class ConfigurationTreeMouseListener extends MouseAdapter 
{
    //
    // member data
    //

    /** reference to the tree to be manipulated */
    private JTree tree = null;
    
    /** reference to the configuration tree-model */
    private ConfigurationTreeModel treeModel = null;
    
    /** popup mneu associated with global psets */
    private JPopupMenu popupPSets = null;

    /** popup menu showing all available service templates */
    private JPopupMenu popupServices = null;

    /** popup menu showing all available event data source templates */
    private JPopupMenu popupEDSources = null;

    /** popup menu showing all available event setup source templates */
    private JPopupMenu popupESSources = null;

    /** popup menu showing all available event setup module templates */
    private JPopupMenu popupESModules = null;

    /** popup menu associated with paths */
    private JPopupMenu popupPaths = null;
    
    /** popup menu associated with sequences */
    private JPopupMenu popupSequences = null;
    
    /** popup menu associated with modules */
    private JPopupMenu popupModules = null;
    
    /** popup menu associated with event contents */
    private JPopupMenu popupContents = null;
    
    /** popup menu associated with streams */
    private JPopupMenu popupStreams = null;
    
    /** popup menu associated with datasets */
    private JPopupMenu popupDatasets = null;
    
    /** action listener for psets-menu actions */
    private PSetMenuListener psetListener = null;
    
    /** action listener for services-menu actions */
    private ServiceMenuListener serviceListener = null;
    
    /** action listener for edsources-menu actions */
    private EDSourceMenuListener edsourceListener = null;

    /** action listener for essources-menu actions */
    private ESSourceMenuListener essourceListener = null;
    
    /** action listener for esmodules-menu actions */
    private ESModuleMenuListener esmoduleListener = null;
    
    /** action listener for paths-menu actions */
    private PathMenuListener pathListener = null;

    /** action listener for sequences-menu actions */
    private SequenceMenuListener sequenceListener = null;

    /** action listener for modules-menu actions */
    private ModuleMenuListener moduleListener = null;

    /** action listern for contents-menu actions */
    private ContentMenuListener contentListener = null;

    /** action listern for streams-menu actions */
    private StreamMenuListener streamListener = null;

    /** action listern for datasets-menu actions */
    private DatasetMenuListener datasetListener = null;

    /** enable the ability to sort components */
    private boolean enableSort = true;

    
    //
    // construction
    //
    
    /** standard constructor */
    public ConfigurationTreeMouseListener(JTree tree, JFrame frame)
    {
	this.tree      = tree;
	this.treeModel = (ConfigurationTreeModel)tree.getModel();
	
	psetListener     = new PSetMenuListener(tree,frame);
	edsourceListener = new EDSourceMenuListener(tree);
	essourceListener = new ESSourceMenuListener(tree);
	esmoduleListener = new ESModuleMenuListener(tree);
	serviceListener  = new ServiceMenuListener(tree);
	pathListener     = new PathMenuListener(tree);
	sequenceListener = new SequenceMenuListener(tree);
	moduleListener   = new ModuleMenuListener(tree);
	contentListener  = new ContentMenuListener(tree);
	streamListener   = new StreamMenuListener(tree,frame);
	datasetListener  = new DatasetMenuListener(tree,frame);
    }
    
    
    //
    // member functions
    //
    
    /** MouseAdapter: mousePressed() */
    public void mousePressed(MouseEvent e) { maybeShowPopup(e); }
    
    /** MouseAdapter: mouseReleased() */
    public void mouseReleased(MouseEvent e) { maybeShowPopup(e); }
    
    /** check if this event should really trigger the menu to be displayed */
    private void maybeShowPopup(MouseEvent e)
    {
	if (!e.isPopupTrigger()) return;
	
	if (!tree.isEditable()) return;

	if (treeModel.getRoot() instanceof ConfigurationModifier) return;

	IConfiguration config = (IConfiguration)treeModel.getRoot();
	if (config.name().length()==0) return;
	
	tree = (JTree)e.getComponent();
	TreePath treePath = tree.getPathForLocation(e.getX(), e.getY());
	if(treePath==null) return;
	int      depth = treePath.getPathCount();
	if(depth<=1) return;
	
	tree.setSelectionPath(treePath);
	
	Object node = treePath.getLastPathComponent();
	
	// show the 'PSets' popup?
	if (isInTreePath(treePath,treeModel.psetsNode())&&depth<=3) {
	    updatePSetMenu();
	    popupPSets.show(e.getComponent(),e.getX(),e.getY());
	    return;
	}
	
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
	
	// show the 'ESModules' popup?
	if (isInTreePath(treePath,treeModel.esmodulesNode())&&depth<=3) {
	    updateESModuleMenu();
	    popupESModules.show(e.getComponent(),e.getX(),e.getY());
	    return;
	}
	
	// show the 'Services' popup?
	if (isInTreePath(treePath,treeModel.servicesNode())&&depth<=3) {
	    updateServiceMenu();
	    popupServices.show(e.getComponent(),e.getX(),e.getY());
	    return;
	}
	
	// show the 'Paths' popup?
	if (isInTreePath(treePath,treeModel.pathsNode())/*&&depth<=4*/) {
	    updatePathMenu();
	    popupPaths.show(e.getComponent(),e.getX(),e.getY());
	    return;
	}
	
	// show the 'Sequences' popup?
	if (isInTreePath(treePath,treeModel.sequencesNode())/*&&depth<=4*/) {
	    updateSequenceMenu();
	    popupSequences.show(e.getComponent(),e.getX(),e.getY());
	    return;
	}
	
	// show the 'Modules' popup?
	if (isInTreePath(treePath,treeModel.modulesNode())&&depth<=3) {
	    updateModuleMenu();
	    popupModules.show(e.getComponent(),e.getX(),e.getY());
	    return;
	}

	// show 'Contents' popup?
	if (isInTreePath(treePath,treeModel.contentsNode())) {
	    updateContentMenu();
	    popupContents.show(e.getComponent(),e.getX(),e.getY());
	}

	// show 'Streams' popup?
	if (isInTreePath(treePath,treeModel.streamsNode())) {
	    updateStreamMenu();
	    popupStreams.show(e.getComponent(),e.getX(),e.getY());
	}

	// show 'Datasets' popup?
	if (isInTreePath(treePath,treeModel.datasetsNode())) {
	    updateDatasetMenu();
	    popupDatasets.show(e.getComponent(),e.getX(),e.getY());
	}


    }
    
    /** check if a node is in a tree path */
    private boolean isInTreePath(TreePath treePath,Object node)
    {
	for (int i=0;i<treePath.getPathCount();i++)
	    if (treePath.getPathComponent(i).equals(node)) return true;
	return false;
    }
    
    /** update 'PSets' menu */
    public void updatePSetMenu()
    {
	JMenuItem menuItem;
	popupPSets = new JPopupMenu();
	
	TreePath treePath = tree.getSelectionPath();
	int      depth    = treePath.getPathCount();
	Object   node     = treePath.getLastPathComponent();
	
	menuItem = new JMenuItem("Add PSet");
	menuItem.addActionListener(psetListener);
	menuItem.setActionCommand("NEWPSET");
	popupPSets.add(menuItem);
    
	if (depth==3) {
	    menuItem = new JMenuItem("Remove PSet");
	    menuItem.addActionListener(psetListener);
	    popupPSets.add(menuItem);
	}
	
	//if (depth==2&&enableSort) {
	//popupPSets.addSeparator();
	//menuItem = new JMenuItem("Sort");
	//menuItem.addActionListener(psetListener);
	//popupPSets.add(menuItem);
	//}
    }
    
    /** update 'EDSource' menu */
    public void updateEDSourceMenu()
    {
	JMenuItem menuItem;
	popupEDSources = new JPopupMenu();
	
	TreePath treePath = tree.getSelectionPath();
	int      depth    = treePath.getPathCount();
	Object   node     = treePath.getLastPathComponent();
	
	JMenu edsourceMenu = new ScrollableMenu("Add EDSource");
	if (depth==3||
	    (depth==2&&treeModel.getChildCount(node)>0)) {
	    edsourceMenu.setEnabled(false);
	    popupEDSources.add(edsourceMenu);
	}
	else {
	    TreeModel       model   = tree.getModel();
	    Configuration   config  = (Configuration)model.getRoot();
	    SoftwareRelease release = config.release();
	    Iterator<EDSourceTemplate> it = release.edsourceTemplateIterator();
	    while (it.hasNext()) {
		EDSourceTemplate t = it.next();
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

	if (depth==2&&enableSort) {
	    popupEDSources.addSeparator();
	    menuItem = new JMenuItem("Sort");
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
	    JMenu essourceMenu = new ScrollableMenu("Add ESSource");

	    TreeModel       model   = tree.getModel();
	    Configuration   config  = (Configuration)model.getRoot();
	    SoftwareRelease release = config.release();
	    Iterator<ESSourceTemplate> it = release.essourceTemplateIterator();
	    
	    while (it.hasNext()) {
		ESSourceTemplate t = it.next();
		if (t.instanceCount()>0) {
		    JMenu instanceMenu = new ScrollableMenu(t.name());
		    menuItem = new JMenuItem("New Instance");
		    menuItem.addActionListener(essourceListener);
		    menuItem.setActionCommand(t.name());
		    instanceMenu.add(menuItem);
		    
		    JMenu copyMenu = new ScrollableMenu("Copy");
		    instanceMenu.add(copyMenu);

		    for (int i=0;i<t.instanceCount();i++) {
			Instance instance = t.instance(i);
			menuItem = new JMenuItem(instance.name());
			menuItem.addActionListener(essourceListener);
			menuItem.setActionCommand("copy:"+
						  t.name()+":"+instance.name());
			copyMenu.add(menuItem);
		    }
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
	    menuItem = new JMenuItem("Rename ESSource");
	    menuItem.addActionListener(essourceListener);
	    popupESSources.add(menuItem);
	    popupESSources.addSeparator();

	    ESPreferable essource = (ESPreferable)treePath.getLastPathComponent();
	    JCheckBoxMenuItem cbMenuItem = new JCheckBoxMenuItem("preferred");
	    cbMenuItem.setState(essource.isPreferred());
	    cbMenuItem.addItemListener(new ESPreferableItemListener(tree));
	    popupESSources.add(cbMenuItem);
	}

	if (depth==2&&enableSort) {
	    popupESSources.addSeparator();
	    menuItem = new JMenuItem("Sort");
	    menuItem.addActionListener(essourceListener);
	    popupESSources.add(menuItem);
	}
    }
    
    /** update 'ESModules' menu */
    public void updateESModuleMenu()
    {
	JMenuItem menuItem;
	popupESModules    = new JPopupMenu();
	
	TreePath treePath = tree.getSelectionPath();
	int      depth    = treePath.getPathCount();
	
	// 'ESModules' or specific event setup module
	if (depth>=2&&depth<=3) {
	    JMenu esmoduleMenu = new ScrollableMenu("Add ESModule");
	    
	    TreeModel       model   = tree.getModel();
	    Configuration   config  = (Configuration)model.getRoot();
	    SoftwareRelease release = config.release();
	    Iterator<ESModuleTemplate>  it = release.esmoduleTemplateIterator();

	    while (it.hasNext()) {
		ESModuleTemplate t = it.next();
		if (t.instanceCount()>0) {
		    JMenu instanceMenu = new ScrollableMenu(t.name());
		    menuItem = new JMenuItem("New Instance");
		    menuItem.addActionListener(esmoduleListener);
		    menuItem.setActionCommand(t.name());
		    instanceMenu.add(menuItem);

		    JMenu copyMenu = new ScrollableMenu("Copy");
		    instanceMenu.add(copyMenu);
		    
		    for (int i=0;i<t.instanceCount();i++) {
			Instance instance = t.instance(i);
			menuItem = new JMenuItem(instance.name());
			menuItem.addActionListener(esmoduleListener);
			menuItem.setActionCommand("copy:"+
						  t.name()+":"+instance.name());
			copyMenu.add(menuItem);
		    }
		    esmoduleMenu.add(instanceMenu);
		}
		else {
		    menuItem = new JMenuItem(t.name());
		    menuItem.addActionListener(esmoduleListener);
		    esmoduleMenu.add(menuItem);
		}
	    }
	    popupESModules.add(esmoduleMenu);
	}
	
	// a specific event setup source is selected
	if (depth==3) {
	    menuItem = new JMenuItem("Remove ESModule");
	    menuItem.addActionListener(esmoduleListener);
	    popupESModules.add(menuItem);
	    menuItem = new JMenuItem("Rename ESModule");
	    menuItem.addActionListener(esmoduleListener);
	    popupESModules.add(menuItem);

	    ESPreferable esmodule = (ESPreferable)treePath.getLastPathComponent();
	    JCheckBoxMenuItem cbMenuItem = new JCheckBoxMenuItem("preferred");
	    cbMenuItem.setState(esmodule.isPreferred());
	    cbMenuItem.addItemListener(new ESPreferableItemListener(tree));
	    popupESModules.add(cbMenuItem);
	}

	if (depth==2&&enableSort) {
	    popupESModules.addSeparator();
	    menuItem = new JMenuItem("Sort");
	    menuItem.addActionListener(esmoduleListener);
	    popupESModules.add(menuItem);
	}
    }
    
    /** update 'Services' menu */
    public void updateServiceMenu()
    {
	JMenuItem menuItem;
	popupServices = new JPopupMenu();
	int depth = tree.getSelectionPath().getPathCount();
	
	JMenu serviceMenu = new ScrollableMenu("Add Service");
	
	TreeModel       model   = tree.getModel();
	Configuration   config  = (Configuration)model.getRoot();
	SoftwareRelease release = config.release();
	Iterator<ServiceTemplate> it = release.serviceTemplateIterator();
	
	while (it.hasNext()) {
	    ServiceTemplate t = it.next();
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
	
	if (depth==2&&enableSort) {
	    popupServices.addSeparator();
	    menuItem = new JMenuItem("Sort");
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

	IConfiguration config = (IConfiguration)treeModel.getRoot();
	
	// 'Paths' selectd
	if (depth==2) {
	    menuItem = new JMenuItem("Add Path");
	    menuItem.addActionListener(pathListener);
	    menuItem.setActionCommand("NEWPATH");
	    popupPaths.add(menuItem);
	}
	
	// specific path is selected
	if (depth==3) {
	    Path path = (Path)node;

	    JMenu addModuleMenu = createAddModuleMenu(path,pathListener);
	    popupPaths.add(addModuleMenu);
	    
	    JMenu addPathMenu = createAddPathMenu(path);
	    popupPaths.add(addPathMenu);
	    
	    JMenu addSequenceMenu = createAddSequenceMenu(path,pathListener,false);
	    popupPaths.add(addSequenceMenu);
	    
	    popupPaths.addSeparator();
	    
	    menuItem = new JMenuItem("Rename Path");
	    menuItem.addActionListener(pathListener);
	    popupPaths.add(menuItem);
	    
	    menuItem = new JMenuItem("Remove Path");
	    menuItem.addActionListener(pathListener);
	    popupPaths.add(menuItem);

	    
	    popupPaths.addSeparator();
	    JCheckBoxMenuItem cbMenuItem = new JCheckBoxMenuItem("endpath");
	    cbMenuItem.setState(path.isEndPath());
	    if (path.hasOutputModule()) cbMenuItem.setEnabled(false);
	    cbMenuItem.addItemListener(new PathItemListener(tree));
	    popupPaths.add(cbMenuItem);
	    
	    return;
	}
	
	// a specific module/path/sequence reference is selected
	if (depth==4) {
	    Path path = (Path)parent;
	    JMenu addModuleMenu = createAddModuleMenu(path,pathListener);
	    popupPaths.add(addModuleMenu);

	    JMenu addPathMenu = createAddPathMenu(path);
	    popupPaths.add(addPathMenu);
	    
	    JMenu addSequenceMenu = createAddSequenceMenu(path,pathListener,false);
	    popupPaths.add(addSequenceMenu);
	    
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

	if (depth==2&&enableSort) {
	    popupPaths.addSeparator();
	    menuItem = new JMenuItem("Sort");
	    menuItem.addActionListener(pathListener);
	    popupPaths.add(menuItem);
	}

	if (node instanceof Reference) {
	    if (depth==4) popupPaths.addSeparator();
	    Reference reference = (Reference)node;
	    menuItem = new JMenuItem("GOTO "+reference.name());
	    menuItem.setActionCommand("GOTO");
	    menuItem.addActionListener(pathListener);
	    popupPaths.add(menuItem);
	}
	
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
	
	IConfiguration config = (IConfiguration)treeModel.getRoot();
	
	if (depth==2) {
	    menuItem = new JMenuItem("Add Sequence");
	    menuItem.addActionListener(sequenceListener);
	    popupSequences.add(menuItem);
	    popupSequences.addSeparator();
	    menuItem = new JMenuItem("Remove Unreferenced Sequences");
	    menuItem.addActionListener(sequenceListener);
	    menuItem.setActionCommand("UNREF");
	    popupSequences.add(menuItem);
	    menuItem = new JMenuItem("Resolve Unnecessary Sequences");
	    menuItem.addActionListener(sequenceListener);
	    menuItem.setActionCommand("RESOLVE");
	    popupSequences.add(menuItem);
	}
	else if (depth==3) {
	    Sequence sequence = (Sequence)node;

	    JMenu addModuleMenu = createAddModuleMenu(sequence,
						      sequenceListener);
	    popupSequences.add(addModuleMenu);
	    
	    JMenu addSequenceMenu = createAddSequenceMenu(sequence,
							  sequenceListener,
							  false);
	    popupSequences.add(addSequenceMenu);
	    
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

	    JMenu addModuleMenu = createAddModuleMenu(sequence,
						      sequenceListener);
	    popupSequences.add(addModuleMenu);

	    JMenu addSequenceMenu = createAddSequenceMenu(sequence,
							  sequenceListener,
							  true);
	    popupSequences.add(addSequenceMenu);
	    
	    if (node instanceof ModuleReference) {
		menuItem = new JMenuItem("Remove Module");
		menuItem.addActionListener(sequenceListener);
		popupSequences.add(menuItem);
	    }
	    if (node instanceof PathReference) {
		menuItem = new JMenuItem("Remove Path");
		menuItem.addActionListener(sequenceListener);
		popupSequences.add(menuItem);
	    }
	    if (node instanceof SequenceReference) {
		menuItem = new JMenuItem("Remove Sequence");
		menuItem.addActionListener(sequenceListener);
		popupSequences.add(menuItem);
	    }
	}

	if (depth==2&&enableSort) {
	    popupSequences.addSeparator();
	    menuItem = new JMenuItem("Sort");
	    menuItem.addActionListener(sequenceListener);
	    popupSequences.add(menuItem);
	}

	if (node instanceof Reference) {
	    if (depth==4) popupSequences.addSeparator();
	    Reference reference = (Reference)node;
	    menuItem = new JMenuItem("GOTO "+reference.name());
	    menuItem.setActionCommand("GOTO");
	    menuItem.addActionListener(sequenceListener);
	    popupSequences.add(menuItem);
	}
	
    }
    
    /** update 'Modules' Menu */
    public void updateModuleMenu()
    {

	TreePath treePath = tree.getSelectionPath();
	int      depth    = treePath.getPathCount();

	JMenuItem menuItem;
	popupModules = new JPopupMenu();	

	if (depth==3) {
	    menuItem = new JMenuItem("Rename Module");
	    menuItem.addActionListener(moduleListener);
	    popupModules.add(menuItem);
	}
	
	if (depth==2&&enableSort) {
	    popupModules.addSeparator();
	    menuItem = new JMenuItem("Sort");
	    menuItem.addActionListener(moduleListener);
	    popupModules.add(menuItem);
	}
    }

    /** update 'EventContents' menu */
    public void updateContentMenu()
    {
	JMenuItem menuItem;
	popupContents = new JPopupMenu();
	
	TreePath treePath = tree.getSelectionPath();
	int      depth    = treePath.getPathCount();
	Object   node     = treePath.getLastPathComponent();
	
	if (depth>3) return;
	
	menuItem = new JMenuItem("Add EventContent");
	menuItem.addActionListener(contentListener);
	menuItem.setActionCommand("ADD");
	popupContents.add(menuItem);
	
	if (depth==3) {
	    EventContent content = (EventContent)node;

	    //menuItem = new JMenuItem("Add Stream");
	    //menuItem.addActionListener(contentListener);
	    //popupContents.add(menuItem);
	    
	    //JMenu addPathMenu = new ScrollableMenu("Add Path");
	    //popupContents.add(addPathMenu);

	    menuItem = new JMenuItem("<html>Rename <i>" + content.label() +
				     "</i></html>");
	    menuItem.addActionListener(contentListener);
	    menuItem.setActionCommand("RENAME");
	    popupContents.add(menuItem);
	    
	    menuItem = new JMenuItem("<html>Remove <i>" + content.label()
				     + "</i></html>");
	    menuItem.addActionListener(contentListener);
	    menuItem.setActionCommand("REMOVE");
	    popupContents.add(menuItem);
	    
	    popupContents.addSeparator();
	    
	    //menuItem = new JMenuItem("Edit " + content.label());
	    //menuItem.addActionListener(contentListener);
	    //popupContents.add(menuItem);
	    
	    //TreeModel      model  = tree.getModel();
	    //IConfiguration config = (IConfiguration)model.getRoot();
	    //ArrayList<Path> paths = new ArrayList<Path>();
	    //Iterator<Path> itP = config.pathIterator();
	    //while (itP.hasNext()) {
	    //Path path = itP.next();
	    //if (content.indexOfPath(path)<0) paths.add(path);
	    //}
	    //Collections.sort(paths);
	    //itP = paths.iterator();
	    //while (itP.hasNext()) {
	    //menuItem = new JMenuItem(itP.next().name());
	    //menuItem.addActionListener(contentListener);
	    //addPathMenu.add(menuItem);
	    //}
	}
    }
    
    /** update 'Streams' menu */
    public void updateStreamMenu()
    {
	JMenuItem menuItem;
	popupStreams = new JPopupMenu();
	
	TreePath treePath = tree.getSelectionPath();
	int      depth    = treePath.getPathCount();
	Object   node     = treePath.getLastPathComponent();
	
	if (depth>3) return;
	
	if (depth==2) {
	    menuItem = new JMenuItem("Add Stream");
	    menuItem.addActionListener(streamListener);
	    menuItem.setActionCommand("ADD");
	    popupStreams.add(menuItem);
	}
	else if (depth==3) {
	    Stream stream = (Stream)node;

	    //menuItem = new JMenuItem("Add Dataset");
	    //menuItem.addActionListener(streamListener);
	    //popupStreams.add(menuItem);
	    
	    //JMenu addPathMenu = new ScrollableMenu("Add Path");
	    //popupStreams.add(addPathMenu);

	    menuItem = new JMenuItem("<html>Rename <i>" + stream.label() +
				     "</i></html>");
	    menuItem.addActionListener(streamListener);
	    menuItem.setActionCommand("RENAME");
	    popupStreams.add(menuItem);
	    
	    menuItem = new JMenuItem("<html>Remove <i>" + stream.label() +
				     "</i></html>");
	    menuItem.addActionListener(streamListener);
	    menuItem.setActionCommand("REMOVE");
	    popupStreams.add(menuItem);
	    
	    //popupStreams.addSeparator();
	    
	    //TreeModel      model  = tree.getModel();
	    //IConfiguration config = (IConfiguration)model.getRoot();
	    //ArrayList<Path> paths = new ArrayList<Path>();
	    //Iterator<Path> itP = config.pathIterator();
	    //while (itP.hasNext()) {
	    //Path path = itP.next();
	    //if (stream.indexOfPath(path)<0) paths.add(path);
	    //}
	    //Collections.sort(paths);
	    //itP = paths.iterator();
	    //while (itP.hasNext()) {
	    //menuItem = new JMenuItem(itP.next().name());
	    //	menuItem.addActionListener(streamListener);
	    //	addPathMenu.add(menuItem);
	    //}
	}
    }
    
    /** update 'PrimaryDatasets' menu */
    public void updateDatasetMenu()
    {
	JMenuItem menuItem;
	popupDatasets = new JPopupMenu();
	
	TreePath treePath = tree.getSelectionPath();
	int      depth    = treePath.getPathCount();
	Object   node     = treePath.getLastPathComponent();
	
	TreeModel      model  = tree.getModel();
	IConfiguration config = (IConfiguration)model.getRoot();

	if (depth>3) return;
	
	if (depth==2) {
	    menuItem = new JMenuItem("Add Primary Dataset");
	    menuItem.addActionListener(datasetListener);
	    menuItem.setActionCommand("ADD");
	    if (config.streamCount()==0) menuItem.setEnabled(false);
	    popupDatasets.add(menuItem);
	}
	else if (depth==3) {
	    PrimaryDataset dataset = (PrimaryDataset)node;
	    
	    //menuItem = new JMenuItem("Add Path");
	    //menuItem.addActionListener(datasetListener);
	    //popupDatasets.add(menuItem);
	    
	    //JMenu addPathMenu = new ScrollableMenu("Add Path");
	    //popupDatasets.add(addPathMenu);
	    
	    menuItem = new JMenuItem("<html>Rename <i>" + dataset.label() +
				     "</i></html>");
	    menuItem.addActionListener(datasetListener);
	    menuItem.setActionCommand("RENAME");
	    popupDatasets.add(menuItem);
	    
	    menuItem = new JMenuItem("<html>Remove <i>" + dataset.label() +
				     "</i></html>");
	    menuItem.addActionListener(datasetListener);
	    menuItem.setActionCommand("REMOVE");
	    popupDatasets.add(menuItem);
	    
	    //popupDatasets.addSeparator();
	    
	    //TreeModel      model  = tree.getModel();
	    //IConfiguration config = (IConfiguration)model.getRoot();
	    //ArrayList<Path> paths =
	    //dataset.parentStream().listOfUnassignedPaths();
	    //Collections.sort(paths);
	    //Iterator<Path> itP = paths.iterator();
	    //while (itP.hasNext()) {
	    //menuItem = new JMenuItem(itP.next().name());
	    //menuItem.addActionListener(datasetListener);
	    //addPathMenu.add(menuItem);
	    //}
	}
    }
    
    
    //
    // private member functions
    //
    
    /** create the 'Add Module' submenu */
    private JMenu createAddModuleMenu(ReferenceContainer container,
				      ActionListener     listener)
    {
	JMenu     addModuleMenu = new JMenu("Add Module");
	JMenuItem menuItemAll;
	JMenuItem menuItem;
	JMenuItem copyItemAll;
	JMenuItem copyItem;

	TreeModel       model   = tree.getModel();
	Configuration   config  = (Configuration)model.getRoot();
	SoftwareRelease release = config.release();

	// explicitely add OutputModule menu
	JMenu outputMenu = new ScrollableMenu("OutputModule");
	addModuleMenu.add(outputMenu);
	Iterator<OutputModule> itOM = config.outputIterator();
	while (itOM.hasNext()) {
	    OutputModule om = itOM.next();
	    menuItem = new JMenuItem(om.name());
	    menuItem.addActionListener(listener);
	    menuItem.setActionCommand("OutputModule");
	    outputMenu.add(menuItem);
	}
	
	// dynamically fill menus to add remaining module types
	HashMap<String,JMenu> menuHashMap = new HashMap<String,JMenu>();
	Iterator<ModuleTemplate> it = release.moduleTemplateIterator();
	while (it.hasNext()) {
	    ModuleTemplate t = it.next();
	    
	    JMenu moduleTypeMenu;
	    JMenu moduleTypeAllMenu;
	    JMenu moduleTypeAndLetterMenu;
	    
	    String moduleType = t.type();
	    if (moduleType.equals("OutputModule")) continue;
	    
	    String moduleTypeAll = moduleType + "All";
	    if (!menuHashMap.containsKey(moduleType)) {
		moduleTypeMenu = new ScrollableMenu(moduleType);
		moduleTypeAllMenu = new ScrollableMenu("All");
		menuHashMap.put(moduleType,moduleTypeMenu);
		menuHashMap.put(moduleTypeAll,moduleTypeAllMenu);
		addModuleMenu.add(moduleTypeMenu);
		moduleTypeMenu.add(moduleTypeAllMenu);
	    }
	    else {
		moduleTypeMenu = menuHashMap.get(moduleType);
		moduleTypeAllMenu = menuHashMap.get(moduleTypeAll);
	    }
	    
	    String moduleLetter = t.name().substring(0,1);
	    String moduleTypeAndLetter = t.type() + moduleLetter;
	    if (!menuHashMap.containsKey(moduleTypeAndLetter)) {
		moduleTypeAndLetterMenu = new ScrollableMenu(moduleLetter);
		menuHashMap.put(moduleTypeAndLetter,moduleTypeAndLetterMenu);
		moduleTypeMenu.add(moduleTypeAndLetterMenu);
	    }
	    else {
		moduleTypeAndLetterMenu = menuHashMap.get(moduleTypeAndLetter);
	    }
	    
	    if (t.instanceCount()>0) {
		JMenu instanceMenuAll = new ScrollableMenu(t.name());
		JMenu instanceMenu = new ScrollableMenu(t.name());
		menuItemAll = new JMenuItem("New Instance");
		menuItem = new JMenuItem("New Instance");
		menuItemAll.addActionListener(listener);
		menuItem.addActionListener(listener);
		menuItemAll.setActionCommand(t.name());
		menuItem.setActionCommand(t.name());
		instanceMenuAll.add(menuItemAll);
		instanceMenu.add(menuItem);

		JMenu copyMenuAll = new ScrollableMenu("Copy");
		JMenu copyMenu = new ScrollableMenu("Copy");
		instanceMenuAll.add(copyMenuAll);
		instanceMenu.add(copyMenu);
		
		instanceMenuAll.addSeparator();
		instanceMenu.addSeparator();

		ArrayList<String> sortedInstanceNames = new ArrayList<String>();
		for (int i=0;i<t.instanceCount();i++) {
			ModuleInstance instance = (ModuleInstance)t.instance(i);
			sortedInstanceNames.add(instance.name());
		}
		Collections.sort(sortedInstanceNames);

		Iterator<String> itI = sortedInstanceNames.iterator();
		while (itI.hasNext()) {
		    String instanceName = itI.next();
			ModuleInstance instance = null;
			try {
				instance = (ModuleInstance)t.instance(instanceName);
			}
			catch (DataException ex) {}
			
			menuItemAll = new JMenuItem(instance.name());
		    menuItem    = new JMenuItem(instance.name());
		    copyItemAll = new JMenuItem(instance.name());
		    copyItem    = new JMenuItem(instance.name());
		    menuItemAll.addActionListener(listener);
		    menuItem.addActionListener(listener);
		    copyItemAll.addActionListener(listener);
		    copyItem.addActionListener(listener);
		    menuItemAll.setActionCommand(t.name()+":"+instance.name());
		    menuItem.setActionCommand(t.name()+":"+instance.name());
		    copyItemAll.setActionCommand("copy:"+t.name()+":"+instance.name());
		    copyItem.setActionCommand("copy:"+t.name()+":"+instance.name());
		    
		    for (int j=0;j<container.entryCount();j++) {
				Reference reference = container.entry(j);
				if (instance.isReferencedBy(reference)) {
					menuItemAll.setEnabled(false);
					menuItem.setEnabled(false);
					copyItemAll.setEnabled(false);
					copyItem.setEnabled(false);
					break;
				}
		    }
		    instanceMenuAll.add(menuItemAll);
		    instanceMenu.add(menuItem);
		    copyMenuAll.add(copyItemAll);
		    copyMenu.add(copyItem);
		}
		moduleTypeAllMenu.add(instanceMenuAll);
		moduleTypeAndLetterMenu.add(instanceMenu);
		}
	    else {
		menuItemAll = new JMenuItem(t.name());
		menuItem = new JMenuItem(t.name());
		menuItemAll.setActionCommand(t.name());
		menuItem.setActionCommand(t.name());
		menuItemAll.addActionListener(listener);
		menuItem.addActionListener(listener);
		moduleTypeAllMenu.add(menuItemAll);
		moduleTypeAndLetterMenu.add(menuItem);
	    }
	}
	return addModuleMenu;
    }

    
    /**create 'Add Path' submenu */
    private JMenu createAddPathMenu(Object node)
    {
	JMenuItem menuItem;
	JMenu     addPathMenu = new ScrollableMenu("Add Path");
	
	ArrayList<Path> forbiddenPaths = new ArrayList<Path>();
	
	if (node instanceof Path) {
	    Path path = (Path)node;
	    forbiddenPaths.add(path);
	    menuItem = new JMenuItem("New Path");
	    menuItem.addActionListener(pathListener);
	    menuItem.setActionCommand("NEWPATH");
	    addPathMenu.add(menuItem);
	    addPathMenu.addSeparator();
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
	IConfiguration config = (IConfiguration)treeModel.getRoot();
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
    private JMenu createAddSequenceMenu(ReferenceContainer pathOrSequence,
					ActionListener     listener,
					boolean            isSeqRef)
    {
	JMenu     addSequenceMenu = new ScrollableMenu("Add Sequence");
	JMenuItem menuItem;
	ArrayList<Sequence> forbiddenSequences = new ArrayList<Sequence>();
	
	if (pathOrSequence instanceof Sequence) {
	    Sequence sequence = (Sequence)pathOrSequence;
	    forbiddenSequences.add(sequence);
	    if (!isSeqRef) {
		menuItem = new JMenuItem("New Sequence");
		menuItem.addActionListener(listener);
		menuItem.setActionCommand("NEWSEQ");
		addSequenceMenu.add(menuItem);
		addSequenceMenu.addSeparator();
	    }
	}
	
	for (int i=0;i<pathOrSequence.entryCount();i++) {
	    Reference reference = pathOrSequence.entry(i);
	    if (reference instanceof SequenceReference) {
		SequenceReference seqreference = (SequenceReference)reference;
		Sequence          sequence     = (Sequence)seqreference.parent();
		forbiddenSequences.add(sequence);
	    }
	}
	
	IConfiguration config = (IConfiguration)treeModel.getRoot();
	for (int i=0;i<config.sequenceCount();i++) {
	    Sequence sequence = config.sequence(i);
	    menuItem = new JMenuItem(sequence.name());
	    menuItem.addActionListener(listener);
	    menuItem.setActionCommand("SEQREF");
	    if (forbiddenSequences.contains(sequence)) menuItem.setEnabled(false);
	    addSequenceMenu.add(menuItem); 
	}
	return addSequenceMenu;
    }

}


/**
 * listen to actions from the 'PSets' popup menu
 */
class PSetMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;
    
    /** reference to the parent frame */
    private JFrame frame = null;
    
    /** standard constructor */
    public PSetMenuListener(JTree tree,JFrame frame)
    {
	this.tree = tree;
	this.frame = frame;
    }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source = (JMenuItem)(e.getSource());
	String    cmd    = source.getText();
	
	/*if (cmd.equals("Sort")) {
	    ConfigurationTreeActions.sortPSets(tree);
	}
	else */
	if (cmd.equals("Remove PSet")) {
	    PSetParameter pset =
		(PSetParameter)tree.getSelectionPath().getLastPathComponent();
	    ConfigurationTreeActions.removePSet(tree,pset);
	}
	else if (cmd.equals("Add PSet")) {
	    AddParameterDialog dlg = new AddParameterDialog(frame,true);
	    dlg.addParameterSet();
	    dlg.pack();
	    dlg.setLocationRelativeTo(frame);
	    dlg.setVisible(true);
	    if (dlg.validChoice()) {
		PSetParameter pset =
		    (PSetParameter)ParameterFactory.create(dlg.type(),
							   dlg.name(),
							   dlg.valueAsString(),
							   dlg.isTracked());
		ConfigurationTreeActions.insertPSet(tree,pset);
	    }
	}
    }
    
}


/**
 * listen to actions from the 'EDSource' popup menu
 */
class EDSourceMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;
    
    /** standard constructor */
    public EDSourceMenuListener(JTree tree) { this.tree = tree; }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	TreePath  treePath = tree.getSelectionPath();
	Object    node     = treePath.getLastPathComponent();
	
	if (cmd.equals("Sort")) {
	    ConfigurationTreeActions.sortEDSources(tree);
	}
	else if (cmd.equals("Remove EDSource")) {
	    EDSourceInstance edsource = (EDSourceInstance)node;
	    ConfigurationTreeActions.removeEDSource(tree,edsource);
	}
	else { // cmd = template name of new EDSource !
	    ConfigurationTreeActions.insertEDSource(tree,cmd);
	}
    }
    
}


/**
 * listen to actions from the 'ESSources' popup menu
 */
class ESSourceMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;

    /** standard constructor */
    public ESSourceMenuListener(JTree tree) { this.tree = tree; }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	TreePath  treePath = tree.getSelectionPath();
	Object    node     = treePath.getLastPathComponent();
	
	if (cmd.equals("Sort")) {
	    ConfigurationTreeActions.sortESSources(tree);
	}
	else if (cmd.equals("Remove ESSource")) {
	    ESSourceInstance essource = (ESSourceInstance)node;
	    ConfigurationTreeActions.removeESSource(tree,essource);
 	}
	else if (cmd.equals("Rename ESSource")) {
	    ConfigurationTreeActions.editNodeName(tree);
 	}
	else {
	    String templateName = source.getActionCommand();
	    ConfigurationTreeActions.insertESSource(tree,templateName);
	}
    }
}


/**
 * listen to actions from the 'ESModules' popup menu
 */
class ESModuleMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;

    /** standard constructor */
    public ESModuleMenuListener(JTree tree) { this.tree = tree; }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	TreePath  treePath = tree.getSelectionPath();
	Object    node     = treePath.getLastPathComponent();
	
	if (cmd.equals("Sort")) {
	    ConfigurationTreeActions.sortESModules(tree);
	}
	else if (cmd.equals("Remove ESModule")) {
	    ESModuleInstance esmodule = (ESModuleInstance)node;
	    ConfigurationTreeActions.removeESModule(tree,esmodule);
 	}
	else if (cmd.equals("Rename ESModule")) {
	    ConfigurationTreeActions.editNodeName(tree);
 	}
	else {
	    String templateName = source.getActionCommand();
	    ConfigurationTreeActions.insertESModule(tree,templateName);
	}
    }
}


/**
 * listen to actions from the 'Services' popup menu
 */
class ServiceMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;
    
    /** standard constructor */
    public ServiceMenuListener(JTree tree) { this.tree = tree; }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	TreePath  treePath = tree.getSelectionPath();
	Object    node     = treePath.getLastPathComponent();
	
	if (cmd.equals("Sort")) {
	    ConfigurationTreeActions.sortServices(tree);
	}
	else if (cmd.equals("Remove Service")) {
	    ServiceInstance service = (ServiceInstance)node;	
	    ConfigurationTreeActions.removeService(tree,service);
	}
	else {
	    String templateName = cmd;
	    ConfigurationTreeActions.insertService(tree,templateName);
	}
    }
    
}


/**
 * listen to actions from the 'Paths' popup menu
 */
class PathMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;

    /** standard constructor */
    public PathMenuListener(JTree tree) { this.tree = tree; }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	String    action   = source.getActionCommand();
	TreePath  treePath = tree.getSelectionPath();
	Object    node     = treePath.getLastPathComponent();
	
	if (action.equals("NEWPATH")) {
	    ConfigurationTreeActions.insertPath(tree);
	}
	else if (action.equals("PATHREF")) {
	    ConfigurationTreeActions.insertReference(tree,"Path",cmd);
	}
	else if (action.equals("SEQREF")) {
	    ConfigurationTreeActions.insertReference(tree,"Sequence",cmd);
	}
	else if (action.equals("GOTO")) {
	    ConfigurationTreeActions.scrollToInstance(tree);
	}
	else if (cmd.equals("Sort")) {
	    ConfigurationTreeActions.sortPaths(tree);
	}
	else if (cmd.equals("Rename Path")) {
	    ConfigurationTreeActions.editNodeName(tree);
	}
	else if (cmd.equals("Remove Path")) {
	    if (node instanceof Path)
		ConfigurationTreeActions.removeReferenceContainer(tree);
	    else if (node instanceof PathReference)
		ConfigurationTreeActions.removeReference(tree);
	}
	else if (cmd.equals("Remove Module")) {
	    ConfigurationTreeActions.removeReference(tree);
	}
	else if (cmd.equals("Remove Sequence")) {
	    ConfigurationTreeActions.removeReference(tree);
	}
	else if (action.equals("OutputModule")) {
	    ConfigurationTreeActions.insertReference(tree,"OutputModule",cmd);
	}
 	// add a module(-reference) to the currently selected path
	else {
	    ConfigurationTreeActions.insertReference(tree,"Module",action);
	}
    }

}


/**
 * listen to actions from the 'Sequences' popup menu
x */
class SequenceMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;

    /** standard constructor */
    public SequenceMenuListener(JTree tree) { this.tree = tree; }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	String    action   = source.getActionCommand();
	TreePath  treePath = tree.getSelectionPath();
	Object    node     = treePath.getLastPathComponent();
	
	if (action.equals("RMUNREF")) {
	    ConfigurationTreeActions.removeUnreferencedSequences(tree);
	}
	else if (action.equals("RESOLVE")) {
	    ConfigurationTreeActions.resolveUnnecessarySequences(tree);
	}
	else if (action.equals("NEWSEQ")) {
	    ConfigurationTreeActions.insertSequence(tree);
	}
	else if (action.equals("SEQREF")) {
	    ConfigurationTreeActions.insertReference(tree,"Sequence",cmd);
	}
	else if (action.equals("GOTO")) {
	    ConfigurationTreeActions.scrollToInstance(tree);
	}
	else if (cmd.equals("Sort")) {
	    ConfigurationTreeActions.sortSequences(tree);
	}
	else if (cmd.equals("Add Sequence")) {
	    ConfigurationTreeActions.insertSequence(tree);
	}
	else if (cmd.equals("Rename Sequence")) {
	    ConfigurationTreeActions.editNodeName(tree);
	}
	else if (cmd.equals("Remove Sequence")) {
	    if (node instanceof Sequence) {
		Sequence sequence = (Sequence)node;
		if (sequence.referenceCount()>0) {
		    StringBuffer warning = new StringBuffer();
		    warning
			.append("Do you really want to remove '")
			.append(sequence.name())
			.append("', which is referenced ")
			.append(sequence.referenceCount()).append(" times?");
		    if (JOptionPane.CANCEL_OPTION==
			JOptionPane.showConfirmDialog(null,warning.toString(),"",
						      JOptionPane.OK_CANCEL_OPTION))
			return;
		}
		ConfigurationTreeActions.removeReferenceContainer(tree);
	    }
	    else if (node instanceof SequenceReference)
		ConfigurationTreeActions.removeReference(tree);
 	}
	else if (cmd.equals("Remove Module")) {
	    ConfigurationTreeActions.removeReference(tree);
	}
	else if (action.equals("OutputModule")) {
	    ConfigurationTreeActions.insertReference(tree,"OutputModule",cmd);
	}
	// add a module to the selected sequence
	else {
	    ConfigurationTreeActions.insertReference(tree,"Module",action);
	}
    }

}


/**
 * listen to actions from the 'Modules' popup menu
 */
class ModuleMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;
    
    /** standard constructor */
    public ModuleMenuListener(JTree tree) { this.tree = tree; }
    
    /** ActionListener interface */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	
	if (cmd.equals("Sort")) {
	    ConfigurationTreeActions.sortModules(tree);
	}
	else if (cmd.equals("Rename Module")) {
	    ConfigurationTreeActions.editNodeName(tree);
 	}
    }
}

/**
 * listen to actions from the 'Contents' popup menu
 */
class ContentMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;

    /** standard constructor */
    public ContentMenuListener(JTree tree) { this.tree = tree; }

    /** ActionListener: actionPerformed() */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	String    action   = source.getActionCommand();
	TreePath  treePath = tree.getSelectionPath();
	Object    node     = treePath.getLastPathComponent();

	if (action.equals("ADD")) {
	    ConfigurationTreeActions.insertContent(tree);
	}
	else if (action.equals("REMOVE")) {
	    ConfigurationTreeActions.removeContent(tree);
	}
	else if (action.equals("RENAME")) {
	    ConfigurationTreeActions.editNodeName(tree);
	}
    }
}

/**
 * listen to actions from the 'Streams' popup menu
 */
class StreamMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;

    /** frame */
    private JFrame frame = null;
    
    /** standard constructor */
    public StreamMenuListener(JTree tree,JFrame frame)
    {
	this.tree = tree;
	this.frame = frame;
    }

    /** ActionListener: actionPerformed() */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	String    action   = source.getActionCommand();
	TreePath  treePath = tree.getSelectionPath();
	Object    node     = treePath.getLastPathComponent();

    	if (action.equals("ADD")) {
	    ConfigurationTreeModel model =
		(ConfigurationTreeModel)tree.getModel();
	    Configuration config = (Configuration)model.getRoot();
	    CreateStreamDialog dlg = new CreateStreamDialog(frame,config);
	    dlg.pack(); dlg.setLocationRelativeTo(frame);
	    dlg.setVisible(true);
	    if (dlg.isSuccess())
		ConfigurationTreeActions.insertStream(tree,dlg.stream());
	}
	else if (action.equals("REMOVE")) {
	    ConfigurationTreeActions.removeStream(tree);
	}
	else if (action.equals("RENAME")) {
	    ConfigurationTreeActions.editNodeName(tree);
	}
    }
}

/**
 * listen to actions from the 'PrimaryDatasets' popup menu
 */
class DatasetMenuListener implements ActionListener
{
    /** reference to the tree to be manipulated */
    private JTree tree = null;

    /** frame */
    private JFrame frame = null;

    /** standard constructor */
    public DatasetMenuListener(JTree tree,JFrame frame)
    {
	this.tree = tree;
	this.frame = frame;
    }

    /** ActionListener: actionPerformed() */
    public void actionPerformed(ActionEvent e)
    {
	JMenuItem source   = (JMenuItem)(e.getSource());
	String    cmd      = source.getText();
	String    action   = source.getActionCommand();
	TreePath  treePath = tree.getSelectionPath();
	Object    node     = treePath.getLastPathComponent();
	
    	if (action.equals("ADD")) {
	    ConfigurationTreeModel model =
		(ConfigurationTreeModel)tree.getModel();
	    Configuration config = (Configuration)model.getRoot();
	    CreateDatasetDialog dlg = new CreateDatasetDialog(frame,config);
	    dlg.pack(); dlg.setLocationRelativeTo(frame);
	    dlg.setVisible(true);
	    if (dlg.isSuccess())
		ConfigurationTreeActions.insertPrimaryDataset(tree,
							      dlg.dataset());
	}
	else if (action.equals("REMOVE")) {
	    ConfigurationTreeActions.removePrimaryDataset(tree);
	}
	else if (action.equals("RENAME")) {
	    ConfigurationTreeActions.editNodeName(tree);
	}
    }
}


/**
 * listen to item events from the ESSource / ESModule menus
 */
class ESPreferableItemListener implements ItemListener
{
    /** reference to  tree */
    private JTree tree = null;

    /** constructor */
    public ESPreferableItemListener(JTree tree) { this.tree = tree; }
    
    /** ItemListener.itemStateChanged() */
    public void itemStateChanged(ItemEvent e)
    {
	JCheckBoxMenuItem item = (JCheckBoxMenuItem)e.getItemSelectable();
	ConfigurationTreeActions.setPreferred(tree,item.isSelected());
    }
}


/**
 * listen to item events of the Paths menu
 */
class PathItemListener implements ItemListener
{
    /** reference to tree */
    private JTree tree = null;

    /** constructor */
    public PathItemListener(JTree tree) { this.tree = tree; }
    
    /** ItemListener.itemStateChanged() */
    public void itemStateChanged(ItemEvent e)
    {
	JCheckBoxMenuItem item = (JCheckBoxMenuItem)e.getItemSelectable();
	ConfigurationTreeActions.setPathAsEndpath(tree,item.isSelected());
    }

}
