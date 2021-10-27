package edu.iastate.metnet.metaomgraph.ui;

import java.awt.EventQueue;
import java.awt.FlowLayout;
import java.awt.Font;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JLabel;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.JTextPane;
import javax.swing.JOptionPane;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;

import org.apache.commons.collections4.CollectionUtils;

import edu.iastate.metnet.metaomgraph.AnimatedSwingWorker;
import edu.iastate.metnet.metaomgraph.MetaOmAnalyzer;
import edu.iastate.metnet.metaomgraph.MetaOmGraph;
import edu.iastate.metnet.metaomgraph.MetaOmProject;
import edu.iastate.metnet.metaomgraph.MetadataHybrid;
import edu.iastate.metnet.metaomgraph.CalculateLogFC;
import edu.iastate.metnet.metaomgraph.DifferentialExpResults;
import edu.iastate.metnet.metaomgraph.FrameModel;
import edu.iastate.metnet.metaomgraph.Metadata.MetadataQuery;
import edu.iastate.metnet.metaomgraph.logging.ActionProperties;

import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.ActionEvent;
import javax.swing.JCheckBox;


/**
 *
 * UI class to choose which features to use in the Limma Analysis
 *
 */

public class LimmaLowExpGeneFrame extends TaskbarInternalFrame {

    private JComboBox comboBox;
    private JComboBox comboBox_1;

    JLabel lblN2;
    JLabel lblN1;

    private JTextField txtGroup1;
    private JTextField txtGroup2;

    private JScrollPane jscp1;
    private JTable tableGrp1;

    private JScrollPane jscp2;
    private JTable tableGrp2;

    private MetadataHybrid mdob;
    private MetaOmProject myProject;

    private boolean[] excludedCopy;

    private JCheckBox chckbxSaveResultsWith;

    /**
     * Default Properties
     */

    private Color SELECTIONBCKGRND = MetaOmGraph.getTableSelectionColor();
    private Color BCKGRNDCOLOR1 = MetaOmGraph.getTableColor1();
    private Color BCKGRNDCOLOR2 = MetaOmGraph.getTableColor2();
    private Color HIGHLIGHTCOLOR = MetaOmGraph.getTableHighlightColor();
    private Color HYPERLINKCOLOR = MetaOmGraph.getTableHyperlinkColor();

    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                try {
                    LimmaLowExpGeneFrame frame = new LimmaLowExpGeneFrame();
                    frame.setVisible(true);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }

    /**
     * Create the frame.
     */
    public LimmaLowExpGeneFrame() {
        setBounds(100, 100, 450, 300);
        getContentPane().setLayout(new BorderLayout(0, 0));
        setTitle("Limma Analysis");

        // init objects
        myProject = MetaOmGraph.getActiveProject();
        mdob = myProject.getMetadataHybrid();
        if (mdob == null) {
            JOptionPane.showMessageDialog(null, "Error. No metadata found", "Error", JOptionPane.ERROR_MESSAGE);
            dispose();
        }

        // get excluded
        boolean[] excluded = MetaOmAnalyzer.getExclude();
        if (excluded != null) {
            excludedCopy = new boolean[excluded.length];
            System.arraycopy(excluded, 0, excludedCopy, 0, excluded.length);
        }

        JPanel panel_1 = new JPanel();
        getContentPane().add(panel_1, BorderLayout.SOUTH);

        JButton btnOk = new JButton("Ok");
        btnOk.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {

                LimmaLowExpGeneFrame lframe = new LimmaLowExpGeneFrame();
                lframe.setSize(lframe.getWidth(), MetaOmGraph.getMainWindow().getHeight() / 2);
                MetaOmGraph.getDesktop().add(lframe);
                lframe.setVisible(true);

                new AnimatedSwingWorker("Working...", true) {
                    @Override
                    public Object construct() {
                        EventQueue.invokeLater(new Runnable() {
                            @Override
                            public void run() {
                                try {

                                    logFCResultsFrame frame = null;
                                    frame.setSize(MetaOmGraph.getMainWindow().getWidth() / 2, MetaOmGraph.getMainWindow().getHeight() / 2);

                                    frame.setEnabled(true);


                                } catch (Exception e) {
                                    StringWriter sw = new StringWriter();
                                    PrintWriter pw = new PrintWriter(sw);
                                    e.printStackTrace(pw);
                                    String sStackTrace = sw.toString();

                                    JDialog jd = new JDialog();
                                    JTextPane jt = new JTextPane();
                                    jt.setText(sStackTrace);
                                    jt.setBounds(10, 10, 300, 100);
                                    jd.getContentPane().add(jt);
                                    jd.setBounds(100, 100, 500, 200);
                                    jd.setVisible(true);
                                }
                            }
                        });
                        return null;
                    }
                }.start();



                //Harsha - reproducibility log

                HashMap<String,Object> actionMap = new HashMap<String,Object>();
                HashMap<String,Object> dataMap = new HashMap<String,Object>();
                HashMap<String,Object> result = new HashMap<String,Object>();

            }
        });
        panel_1.add(btnOk);

        JPanel panel_2 = new JPanel();
        getContentPane().add(panel_2, BorderLayout.CENTER);
        panel_2.setLayout(new BorderLayout(0, 0));

        JSplitPane splitPane = new JSplitPane();
        panel_2.add(splitPane, BorderLayout.CENTER);
        splitPane.setResizeWeight(0.5);

        JPanel panel_3 = new JPanel();
        splitPane.setRightComponent(panel_3);
        panel_3.setLayout(new BorderLayout(0, 0));

        txtGroup2 = new JTextField();
        txtGroup2.setText("Group2");
        txtGroup2.setColumns(10);
        JPanel topbtnPnl2 = new JPanel();
        JButton sendLeft = new JButton("<<");
        sendLeft.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
            }
        });
        topbtnPnl2.add(sendLeft);

        JLabel lblGroupName_1 = new JLabel("Group name:");
        topbtnPnl2.add(lblGroupName_1);
        topbtnPnl2.add(txtGroup2);
        panel_3.add(topbtnPnl2, BorderLayout.NORTH);
        //topbtnPnl2.setLayout(new BoxLayout(topbtnPnl2, BoxLayout.LINE_AXIS));
        topbtnPnl2.setLayout(new FlowLayout());
        lblN2 = new JLabel("n=0");
        lblN2.setFont(new Font("Tahoma", Font.BOLD, 13));
        topbtnPnl2.add(lblN2);

        // add table2
        jscp2 = new JScrollPane();
        jscp2.setViewportView(tableGrp2);
        panel_3.add(jscp2, BorderLayout.CENTER);

        JButton btnAdd2 = new JButton("Add");
        btnAdd2.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
            }
        });
        JPanel btnPnl2 = new JPanel(new FlowLayout());
        btnPnl2.add(btnAdd2);
        JButton btnRem2 = new JButton("Remove");
        btnRem2.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
            }
        });
        btnPnl2.add(btnRem2);
        JButton btnSearch2 = new JButton("Search");
        btnSearch2.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // search in list by metadatata
            }
        });
        btnPnl2.add(btnSearch2);
        panel_3.add(btnPnl2, BorderLayout.SOUTH);

        //btnPnl2.setLayout(new BoxLayout(btnPnl2, BoxLayout.LINE_AXIS));
        btnPnl2.setLayout(new FlowLayout());
        JPanel panel_4 = new JPanel();
        splitPane.setLeftComponent(panel_4);
        panel_4.setLayout(new BorderLayout(0, 0));

        txtGroup1 = new JTextField();
        txtGroup1.setText("Group1");
        txtGroup1.setColumns(10);
        JPanel topbtnPnl1 = new JPanel();
        JButton sendRight = new JButton(">>");
        sendRight.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent arg0) {
            }
        });

        lblN1 = new JLabel("n=0");
        lblN1.setFont(new Font("Tahoma", Font.BOLD, 13));
        topbtnPnl1.add(lblN1);

        JLabel lblGroupName = new JLabel("Group name:");
        topbtnPnl1.add(lblGroupName);
        topbtnPnl1.add(txtGroup1);
        topbtnPnl1.add(sendRight);

        panel_4.add(topbtnPnl1, BorderLayout.NORTH);

        //topbtnPnl1.setLayout(new BoxLayout(topbtnPnl1, BoxLayout.LINE_AXIS));
        topbtnPnl1.setLayout(new FlowLayout());

        // add table1
        jscp1 = new JScrollPane();
        jscp1.setViewportView(tableGrp1);
        panel_4.add(jscp1, BorderLayout.CENTER);

        JButton btnAdd1 = new JButton("Add");
        btnAdd1.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {

            }
        });
        JPanel btnPnl1 = new JPanel();
        btnPnl1.add(btnAdd1);
        JButton btnRem1 = new JButton("Remove");
        btnRem1.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent arg0) {
            }
        });
        btnPnl1.add(btnRem1);
        JButton btnSearch1 = new JButton("Search");
        btnSearch1.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // search in list by metadatata
            }
        });
        btnPnl1.add(btnSearch1);
        panel_4.add(btnPnl1, BorderLayout.SOUTH);

        //btnPnl1.setLayout(new BoxLayout(btnPnl1, BoxLayout.LINE_AXIS));
        btnPnl1.setLayout(new FlowLayout());

        // frame properties
        this.setClosable(true);
        pack();

        int defaultWidth = this.getWidth();
        int defaultHeight = MetaOmGraph.getMainWindow().getHeight() / 2;
        LimmaLowExpGeneFrame thisFrame = this;

        this.addComponentListener(new ComponentListener() {

            @Override
            public void componentShown(ComponentEvent e) {
                // TODO Auto-generated method stub

            }

            @Override
            public void componentResized(ComponentEvent e) {
                // TODO Auto-generated method stub

                if(thisFrame.getWidth() < defaultWidth || thisFrame.getHeight() < defaultHeight) {
                    thisFrame.setSize(defaultWidth, defaultHeight);
                }
            }

            @Override
            public void componentMoved(ComponentEvent e) {
                // TODO Auto-generated method stub

            }

            @Override
            public void componentHidden(ComponentEvent e) {
                // TODO Auto-generated method stub

            }
        });

        putClientProperty("JInternalFrame.frameType", "normal");
        setResizable(true);
        setMaximizable(true);
        setIconifiable(true);
        setClosable(true);

        FrameModel diffExpFrameModel = new FrameModel("DEA","Differential Expression Analysis",12);
        setModel(diffExpFrameModel);
    }
}
