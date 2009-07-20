package confdb.gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;


/**
 * MenuBar
 * -------
 * @author Philipp Schieferdecker
 */
public class MenuBar
{
    //
    // member data
    //
    
    /** the menu bar */
    private JMenuBar jMenuBar = null;

    /** menu bar item names: confdbMenu */
    private static final String confdbMenuAbout = "About";
    private static final String confdbMenuQuit  = "Quit";
    
    /** menu bar item names: configMenu */
    private static final String configMenuNew         = "New";
    private static final String configMenuParse       = "Parse";
    private static final String configMenuOpen        = "Open";
    private static final String configMenuClose       = "Close";
    private static final String configMenuSave        = "Save";
    private static final String configMenuCommentSave = "Comment&Save";
    private static final String configMenuSaveAs      = "Save As";
    private static final String configMenuImport      = "Import";
    private static final String configMenuMigrate     = "Migrate";
    private static final String configMenuConvert     = "Convert";
    
    /** menu bar item names: toolMenu */
    private static final String toolMenuReplace       = "Search&Replace";
    private static final String toolMenuDiff          = "Compare (Diff)";
    private static final String toolMenuOMEditor      = "Edit OutputModules";
    private static final String toolMenuPSEditor      = "Edit Prescales";
    private static final String toolMenuMLEditor      = "Edit MessageLogger";

    /** menu bar item names: options */
    private static final String optionsMenuTrack      = "Track InputTags";

    /** menu bar item names: dbMenu */
    private static final String dbMenuConnectToDB      = "Connect to DB";
    private static final String dbMenuDisconnectFromDB = "Disconnect from DB";
    private static final String dbMenuExportConfig     = "Export Configuration";
    
    /** meny bar items */
    private JMenuItem confdbMenuAboutItem       = null;
    private JMenuItem confdbMenuQuitItem        = null;
    
    private JMenuItem configMenuNewItem         = null;
    private JMenuItem configMenuParseItem       = null;
    private JMenuItem configMenuOpenItem        = null;
    private JMenuItem configMenuCloseItem       = null;
    private JMenuItem configMenuCommentSaveItem = null;
    private JMenuItem configMenuSaveItem        = null;
    private JMenuItem configMenuSaveAsItem      = null;

    private JMenuItem configMenuImportItem      = null;
    private JMenuItem configMenuMigrateItem     = null;
    private JMenuItem configMenuConvertItem     = null;

    private JMenuItem toolMenuDiffItem          = null;
    private JMenuItem toolMenuReplaceItem       = null;
    private JMenuItem toolMenuOMEditorItem      = null;
    private JMenuItem toolMenuPSEditorItem      = null;
    private JMenuItem toolMenuMLEditorItem      = null;
    
    private JCheckBoxMenuItem optionsMenuTrackItem = null;

    private JMenuItem dbMenuConnectItem         = null;
    private JMenuItem dbMenuDisconnectItem      = null;
    private JMenuItem dbMenuExportItem          = null;
    
    

    //
    // construction
    //
    public MenuBar(JMenuBar jMenuBar,ConfDbGUI app)
    {
	this.jMenuBar = jMenuBar;
	populateMenuBar(new CommandActionListener(app));
    }
    
    
    //
    // memeber functions
    //
    
    /** a configuration was opened */
    public void configurationIsOpen()
    {
	dbConnectionIsEstablished();
	configMenuNewItem.setEnabled(true);
	if (configMenuParseItem!=null) configMenuParseItem.setEnabled(true);
	configMenuOpenItem.setEnabled(true);
	configMenuCloseItem.setEnabled(true);
	configMenuSaveItem.setEnabled(true);
	configMenuCommentSaveItem.setEnabled(true);
	configMenuSaveAsItem.setEnabled(true);
	configMenuImportItem.setEnabled(true);
	configMenuMigrateItem.setEnabled(true);
	configMenuConvertItem.setEnabled(true);
	toolMenuReplaceItem.setEnabled(true);
	toolMenuDiffItem.setEnabled(true);
	toolMenuOMEditorItem.setEnabled(true);
	toolMenuPSEditorItem.setEnabled(true);
	toolMenuMLEditorItem.setEnabled(true);
	optionsMenuTrackItem.setEnabled(true);
	dbMenuExportItem.setEnabled(true);
    }

    /** no configuration is open */
    public void configurationIsNotOpen()
    {
	configMenuCloseItem.setEnabled(false);
	configMenuSaveItem.setEnabled(false);
	configMenuCommentSaveItem.setEnabled(false);
	configMenuSaveAsItem.setEnabled(false);
	configMenuImportItem.setEnabled(false);
	configMenuMigrateItem.setEnabled(false);
	configMenuConvertItem.setEnabled(false);
	//toolMenuDiffItem.setEnabled(false);
	toolMenuReplaceItem.setEnabled(false);
	toolMenuOMEditorItem.setEnabled(false);
	toolMenuPSEditorItem.setEnabled(false);
	toolMenuMLEditorItem.setEnabled(false);
	optionsMenuTrackItem.setEnabled(false);
	dbMenuExportItem.setEnabled(false);
    }
    
    /** database connection is established */
    public void dbConnectionIsEstablished()
    {
	configMenuNewItem.setEnabled(true);
	if (configMenuParseItem!=null) configMenuParseItem.setEnabled(true);
	configMenuOpenItem.setEnabled(true);
	toolMenuDiffItem.setEnabled(true);
	dbMenuDisconnectItem.setEnabled(true);
    }
    
    /** no database connection is established */
    public void dbConnectionIsNotEstablished()
    {
	configurationIsNotOpen();
	configMenuNewItem.setEnabled(false);
	if (configMenuParseItem!=null) configMenuParseItem.setEnabled(false);
	configMenuOpenItem.setEnabled(false);
	toolMenuDiffItem.setEnabled(false);
	dbMenuDisconnectItem.setEnabled(false);
    }

    /** populate the menu bar with all menus and their items */
    private void populateMenuBar(CommandActionListener listener)
    {
	JMenuItem menuItem;
	
	JMenu confdbMenu = new JMenu("ConfDbGUI");
	confdbMenu.setMnemonic(KeyEvent.VK_E);
	jMenuBar.add(confdbMenu);
	confdbMenuAboutItem = new JMenuItem(confdbMenuAbout,KeyEvent.VK_A);
	confdbMenuAboutItem.setActionCommand(confdbMenuAbout);
	confdbMenuAboutItem.addActionListener(listener);
	confdbMenu.add(confdbMenuAboutItem);
	confdbMenuQuitItem = new JMenuItem(confdbMenuQuit,KeyEvent.VK_Q);
	confdbMenuQuitItem.setActionCommand(confdbMenuQuit);
	confdbMenuQuitItem.addActionListener(listener);
	confdbMenu.add(confdbMenuQuitItem);
	
	JMenu configMenu = new JMenu("Configurations");
	configMenu.setMnemonic(KeyEvent.VK_C);
	jMenuBar.add(configMenu);
	configMenuNewItem = new JMenuItem(configMenuNew,KeyEvent.VK_N);
	configMenuNewItem.setActionCommand(configMenuNew);
	configMenuNewItem.addActionListener(listener);
	configMenu.add(configMenuNewItem);
	configMenuParseItem = new JMenuItem(configMenuParse,KeyEvent.VK_P);
	configMenuParseItem.setActionCommand(configMenuParse);
	configMenuParseItem.addActionListener(listener);
	configMenu.add(configMenuParseItem);
	configMenuOpenItem = new JMenuItem(configMenuOpen,KeyEvent.VK_O);
	configMenuOpenItem.setActionCommand(configMenuOpen);
	configMenuOpenItem.addActionListener(listener);
	configMenu.add(configMenuOpenItem);
	configMenuCloseItem = new JMenuItem(configMenuClose,KeyEvent.VK_C);
	configMenuCloseItem.setActionCommand(configMenuClose);
	configMenuCloseItem.addActionListener(listener);
	configMenu.add(configMenuCloseItem);
	configMenu.addSeparator();
	configMenuSaveItem = new JMenuItem(configMenuSave,KeyEvent.VK_S);
	configMenuSaveItem.addActionListener(listener);
	configMenu.add(configMenuSaveItem);
	configMenuCommentSaveItem=new JMenuItem(configMenuCommentSave,KeyEvent.VK_R);
	configMenuCommentSaveItem.setActionCommand(configMenuCommentSave);
	configMenuCommentSaveItem.addActionListener(listener);
	configMenu.add(configMenuCommentSaveItem);
	configMenuSaveAsItem = new JMenuItem(configMenuSaveAs,KeyEvent.VK_A);
	configMenuSaveAsItem.setActionCommand(configMenuSaveAs);
	configMenuSaveAsItem.addActionListener(listener);
	configMenu.add(configMenuSaveAsItem);
	configMenu.addSeparator();
	configMenuImportItem = new JMenuItem(configMenuImport,KeyEvent.VK_I);
	configMenuImportItem.setActionCommand(configMenuImport);
	configMenuImportItem.addActionListener(listener);
	configMenu.add(configMenuImportItem);
	configMenuMigrateItem = new JMenuItem(configMenuMigrate,KeyEvent.VK_M);
	configMenuMigrateItem.setActionCommand(configMenuMigrate);
	configMenuMigrateItem.addActionListener(listener);
	configMenu.add(configMenuMigrateItem);
	configMenuConvertItem = new JMenuItem(configMenuConvert,KeyEvent.VK_N);
	configMenuConvertItem.setActionCommand(configMenuConvert);
	configMenuConvertItem.addActionListener(listener);
	configMenu.add(configMenuConvertItem);
	configMenu.addSeparator();
	
	JMenu toolMenu = new JMenu("Tools");
	toolMenu.setMnemonic(KeyEvent.VK_T);
	jMenuBar.add(toolMenu);
	toolMenuDiffItem = new JMenuItem(toolMenuDiff,KeyEvent.VK_D);
	toolMenuDiffItem.setActionCommand(toolMenuDiff);
	toolMenuDiffItem.addActionListener(listener);
	toolMenu.add(toolMenuDiffItem);
	toolMenuReplaceItem = new JMenuItem(toolMenuReplace,KeyEvent.VK_R);
	toolMenuReplaceItem.setActionCommand(toolMenuReplace);
	toolMenuReplaceItem.addActionListener(listener);
	toolMenu.add(toolMenuReplaceItem);
	toolMenuOMEditorItem = new JMenuItem(toolMenuOMEditor,KeyEvent.VK_O);
	toolMenuOMEditorItem.setActionCommand(toolMenuOMEditor);
	toolMenuOMEditorItem.addActionListener(listener);
	toolMenu.add(toolMenuOMEditorItem);
	toolMenuPSEditorItem = new JMenuItem(toolMenuPSEditor,KeyEvent.VK_P);
	toolMenuPSEditorItem.setActionCommand(toolMenuPSEditor);
	toolMenuPSEditorItem.addActionListener(listener);
	toolMenu.add(toolMenuPSEditorItem);
	toolMenuMLEditorItem = new JMenuItem(toolMenuMLEditor,KeyEvent.VK_Q);
	toolMenuMLEditorItem.setActionCommand(toolMenuMLEditor);
	toolMenuMLEditorItem.addActionListener(listener);
	toolMenu.add(toolMenuMLEditorItem);
	


	JMenu optionsMenu = new JMenu("Options");
	optionsMenu.setMnemonic(KeyEvent.VK_O);
	jMenuBar.add(optionsMenu);
	optionsMenuTrackItem = new JCheckBoxMenuItem(optionsMenuTrack);
	optionsMenuTrackItem.addActionListener(listener);
	optionsMenu.add(optionsMenuTrackItem);

	JMenu dbMenu = new JMenu("Database");
	dbMenu.setMnemonic(KeyEvent.VK_D);
	jMenuBar.add(dbMenu);
	dbMenuConnectItem = new JMenuItem(dbMenuConnectToDB,KeyEvent.VK_C);
	dbMenuConnectItem.setActionCommand(dbMenuConnectToDB);
	dbMenuConnectItem.addActionListener(listener);
	dbMenu.add(dbMenuConnectItem);
	dbMenuDisconnectItem = new JMenuItem(dbMenuDisconnectFromDB,KeyEvent.VK_D);
	dbMenuDisconnectItem.setActionCommand(dbMenuDisconnectFromDB);
	dbMenuDisconnectItem.addActionListener(listener);
	dbMenu.add(dbMenuDisconnectItem);
	dbMenuExportItem = new JMenuItem(dbMenuExportConfig,KeyEvent.VK_E);
	dbMenuExportItem.setActionCommand(dbMenuExportConfig);
	dbMenuExportItem.addActionListener(listener);
	dbMenu.add(dbMenuExportItem);
    }
    

}
