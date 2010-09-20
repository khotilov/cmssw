package confdb.gui;

import java.util.Iterator;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.event.*;

import confdb.data.EventContent;
import confdb.data.Stream;
import confdb.data.Configuration;

/**
 * CreateStreamDialog
 * ------------------
 * @author Philipp Schieferdecker
 *
 * Let the user create a new stream and either link it an existing
 * content or create a new one.
 */
public class CreateStreamDialog extends JDialog
{
    //
    // member data
    //
    
    /** reference to the configuration */
    private Configuration config = null;

    /** created stream */
    private Stream stream = null;

    /** GUI components */
    private JButton    jButtonCancel;
    private JButton    jButtonOK;
    private JComboBox  jComboBoxEventContent;
    private JLabel     jLabelEventContent;
    private JLabel     jLabelStreamName;
    private JTextField jTextFieldStreamName;

    //
    // construction
    //
    
    /** Creates new form CreateStreamDialog */
    public CreateStreamDialog(JFrame frame,Configuration config)
    {
	super(frame,true);
        this.config = config;
        setContentPane(initComponents());
        updateEventContentList();
	setTitle("Create New Stream");
    }

    //
    // public member functions
    //

    /** indicate if a stream was successfully created */
    public boolean isSuccess() { return stream!=null; }

    /** retrieve the created stream */
    public Stream stream() { return stream; }


    //
    // private member functions
    //

    /** update the list of available event context, put into combo box */
    private void updateEventContentList()
    {
        DefaultComboBoxModel cbm =
           (DefaultComboBoxModel)jComboBoxEventContent.getModel();
        cbm.removeAllElements();
        cbm.addElement(new String("<NEW>"));
        Iterator<EventContent> itC = config.contentIterator();
        while (itC.hasNext()) {
            EventContent content = itC.next();
            cbm.addElement(content.name());
        }
    }

    //
    // ACTIONLISTENER CALLBACKS
    //
    private void jButtonOKActionPerformed(ActionEvent evt)
    {
        String streamName = jTextFieldStreamName.getText();
        String contentName = (String)jComboBoxEventContent.getSelectedItem();
        EventContent content = config.content(contentName);
        if (content==null) {
            content = config.insertContent("hltEventContent" + streamName);
        }
        stream = content.insertStream(streamName);
        setVisible(false);
    }
    private void jButtonCancelActionPerformed(ActionEvent evt)
    {
        setVisible(false);
    }

    //
    // DOCUMENTLISTENER CALLBACKS
    //
    private void jTextFieldStreamNameInsertUpdate(DocumentEvent e)
    {
	String streamName = jTextFieldStreamName.getText();
	if (config.stream(streamName)==null) jButtonOK.setEnabled(true);
	else jButtonOK.setEnabled(false);
    }
    public void jTextFieldStreamNameRemoveUpdate(DocumentEvent e)
    {
	String streamName = jTextFieldStreamName.getText();
	if (config.stream(streamName)==null) jButtonOK.setEnabled(true);
	else jButtonOK.setEnabled(false);
    }
    

    /** generate the graphical components */
    private JPanel initComponents() {

	JPanel jPanel = new JPanel();
	
        jLabelStreamName = new javax.swing.JLabel();
        jTextFieldStreamName = new javax.swing.JTextField();
        jComboBoxEventContent = new javax.swing.JComboBox();
        jLabelEventContent = new javax.swing.JLabel();
        jButtonOK = new javax.swing.JButton();
        jButtonCancel = new javax.swing.JButton();

        jLabelStreamName.setText("Stream Name:");

	jTextFieldStreamName.getDocument()
	    .addDocumentListener(new DocumentListener() {
		    public void insertUpdate(DocumentEvent e) {
			jTextFieldStreamNameInsertUpdate(e);
		    }
		    public void removeUpdate(DocumentEvent e) {
			jTextFieldStreamNameRemoveUpdate(e);
		    }
		    public void changedUpdate(DocumentEvent e) {}
		});

        jComboBoxEventContent.setModel(new javax.swing.DefaultComboBoxModel(new String[] { }));

        jLabelEventContent.setText("Event Content:");

        jButtonOK.setText("OK");
        jButtonCancel.setText("Cancel");
        jButtonOK.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButtonOKActionPerformed(evt);
            }
        });
        jButtonCancel.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButtonCancelActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(jPanel);
        jPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .addContainerGap()
                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(layout.createSequentialGroup()
                        .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(jTextFieldStreamName, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 338, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                            .add(jLabelStreamName))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(jLabelEventContent)
                            .add(jComboBoxEventContent, 0, 324, Short.MAX_VALUE)))
                    .add(org.jdesktop.layout.GroupLayout.TRAILING, layout.createSequentialGroup()
                        .add(jButtonCancel)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(jButtonOK)))
                .addContainerGap())
        );

        layout.linkSize(new java.awt.Component[] {jButtonCancel, jButtonOK}, org.jdesktop.layout.GroupLayout.HORIZONTAL);

        layout.setVerticalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .addContainerGap()
                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabelStreamName)
                    .add(jLabelEventContent))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jTextFieldStreamName, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jComboBoxEventContent, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jButtonOK)
                    .add(jButtonCancel))
                .addContainerGap())
        );
	
	return jPanel;
    }

}
