package edu.iastate.metnet.metaomgraph.ui;

import edu.iastate.metnet.metaomgraph.*;
import org.apache.commons.collections4.CollectionUtils;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;
import java.util.List;

public class LimmaGroupPanel extends JPanel  {
    private Color SELECTIONBCKGRND = MetaOmGraph.getTableSelectionColor();
    private Color BCKGRNDCOLOR1 = MetaOmGraph.getTableColor1();
    private Color BCKGRNDCOLOR2 = MetaOmGraph.getTableColor2();
    private Color HIGHLIGHTCOLOR = MetaOmGraph.getTableHighlightColor();
    private Color HYPERLINKCOLOR = MetaOmGraph.getTableHyperlinkColor();
    MetaOmProject myProject = MetaOmGraph.getActiveProject();
    MetadataHybrid mdob = myProject.getMetadataHybrid();

    public LimmaGroupPanel(int id) {
        JPanel panel = new JPanel();
        getRootPane().getContentPane().add(panel, BorderLayout.CENTER);
        panel.setLayout(new BorderLayout(0, 0));

        JSplitPane splitPane = new JSplitPane();
        panel.add(splitPane, BorderLayout.CENTER);
        splitPane.setResizeWeight(0.5);

        JPanel panel_3 = new JPanel();
        splitPane.setRightComponent(panel_3);
        panel_3.setLayout(new BorderLayout(0, 0));

        JTextField txtGroup = new JTextField();
        txtGroup.setText("Group" + id);
        txtGroup.setColumns(10);
        JPanel topbtnPnl2 = new JPanel();
        JButton sendLeft = new JButton("<<");
        sendLeft.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                moveSelectedtoLeft();
            }
        });
        topbtnPnl2.add(sendLeft);

        JLabel lblGroupName_1 = new JLabel("Group name:");
        topbtnPnl2.add(lblGroupName_1);
        topbtnPnl2.add(txtGroup);
        panel_3.add(topbtnPnl2, BorderLayout.NORTH);
        //topbtnPnl2.setLayout(new BoxLayout(topbtnPnl2, BoxLayout.LINE_AXIS));
        topbtnPnl2.setLayout(new FlowLayout());
        JLabel lblN = new JLabel("n=0");
        lblN.setFont(new Font("Tahoma", Font.BOLD, 13));
        topbtnPnl2.add(lblN);

        // add table2
        JScrollPane jscp = new JScrollPane();
        JTable tableGrp = initTableModel();
        // updateTableData(tableGrp, mdob.getMetadataCollection().getAllDataCols());
        updateTableData(tableGrp, null);
        jscp.setViewportView(tableGrp);
        panel_3.add(jscp, BorderLayout.CENTER);

        JButton btnAdd2 = new JButton("Add");
        btnAdd2.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                java.util.List<String> queryRes = showSearchMetadataPanel();
                if (queryRes == null || queryRes.size() < 1) {
                    return;
                }
                // JOptionPane.showConfirmDialog(null, "match:" + queryRes.toString());
                addRows(tableGrp, queryRes);
            }
        });
        JPanel btnPnl2 = new JPanel(new FlowLayout());
        btnPnl2.add(btnAdd2);
        JButton btnRem2 = new JButton("Remove");
        btnRem2.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                removeSelectedRows(tableGrp);
            }
        });
        btnPnl2.add(btnRem2);
        JButton btnSearch2 = new JButton("Search");
        btnSearch2.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // search in list by metadatata
                java.util.List<String> queryRes = showSearchMetadataPanel();
                if (queryRes == null || queryRes.size() < 1) {
                    return;
                }
                // get intersection
                java.util.List<String> allRows = getAllRows(tableGrp);
                java.util.List<String> res = (List<String>) CollectionUtils.intersection(queryRes, allRows);
                // set selected and bring to top
                setSelectedRows(res, tableGrp);
            }
        });
        btnPnl2.add(btnSearch2);
        panel_3.add(btnPnl2, BorderLayout.SOUTH);
    }

    private JTable initTableModel() {
        JTable table = new JTable() {
            @Override
            public boolean getScrollableTracksViewportWidth() {
                return getPreferredSize().width < getParent().getWidth();
            }

            @Override
            public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
                Component c = super.prepareRenderer(renderer, row, column);

                if (!isRowSelected(row)) {
                    c.setBackground(getBackground());
                    int modelRow = convertRowIndexToModel(row);

                    if (row % 2 == 0) {
                        c.setBackground(BCKGRNDCOLOR1);
                    } else {
                        c.setBackground(BCKGRNDCOLOR2);
                    }

                } else {
                    c.setBackground(SELECTIONBCKGRND);
                }

                return c;
            }

        };
        // disable colum drag
        table.getTableHeader().setReorderingAllowed(false);
        DefaultTableModel model = new DefaultTableModel() {
            private static final long serialVersionUID = 1L;

            @Override
            public Class<?> getColumnClass(int column) {
                switch (column) {
                    case 0:
                        return String.class;
                    default:
                        return String.class;
                }
            }

            @Override
            public boolean isCellEditable(int row, int column) {
                return false;
            }
        };
        table.setModel(model);
        // set properties
        table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        table.setAutoCreateRowSorter(true);
        table.setPreferredScrollableViewportSize(table.getPreferredSize());
        table.setFillsViewportHeight(true);
        table.getTableHeader().setFont(new Font("Garamond", Font.BOLD, 14));
        return table;
    }

    private void updateTableData(JTable table, List<String> rows) {
        DefaultTableModel tablemodel = (DefaultTableModel) table.getModel();
        tablemodel.setRowCount(0);
        tablemodel.setColumnCount(0);
        // add data
        String dcName = mdob.getDataColName();
        tablemodel.addColumn(dcName);

        if (rows == null) {
            return;
        }

        // convert list to set remove duplicates
        Set<String> tempSet = new TreeSet<String>(rows);
        rows = new ArrayList<>();
        rows.addAll(tempSet);
        Vector temp = null;
        for (String s : rows) {
            temp = new Vector<>();
            temp.add(s);
            tablemodel.addRow(temp);
        }
    }

    private void initComboBoxes() {
        JComboBox comboBox = new JComboBox(MetaOmGraph.getActiveProject().getGeneListNames());
        String[] methods = new String[] { "M-W U test", "Student's t-test", "Welch's t-test", "Permutation test",
                "Paired t-test", "Wilcoxon Signed Rank Test", "Permutation test (paired samples)" };
        JComboBox comboBox_1 = new JComboBox(methods);
    }

    /**
     * move selected rows from table 1 to table 2
     */
    private void moveSelectedtoRight() {
//        List<String> selected1 = getSelectedRows(tableGrp1);
//        addRows(tableGrp2, selected1);
//        removeSelectedRows(tableGrp1);
//        updateLabelN();
    }

    private void moveSelectedtoLeft() {
//        List<String> selected2 = getSelectedRows(tableGrp2);
//        addRows(tableGrp1, selected2);
//        removeSelectedRows(tableGrp2);
//        updateLabelN();

    }

    /**
     * get selected rows from a table
     *
     * @param table
     * @return
     */
    private List<String> getSelectedRows(JTable table) {
        return getSelectedRows(table, false);
    }

    private List<String> getSelectedRows(JTable table, boolean invert) {
        int selected[] = table.getSelectedRows();
        List<String> res = new ArrayList<>();
        for (int i = 0; i < selected.length; i++) {
            String thisRow = "";
            thisRow = (String) table.getValueAt(selected[i], 0);
            res.add(thisRow);
        }
        // JOptionPane.showMessageDialog(null, "sel:" + res);

        if (invert) {
            List<String> temp = new ArrayList<>();
            for (int i = 0; i < table.getRowCount(); i++) {
                String thisRow = "";
                thisRow = (String) table.getValueAt(i, 0);
                if (!res.contains(thisRow)) {
                    temp.add(thisRow);
                }
            }
            res = temp;
        }

        return res;
    }

    private void removeSelectedRows(JTable table) {
        List<String> toKeep = getSelectedRows(table, true);
        // JOptionPane.showMessageDialog(null, "tokeep:" + toKeep.toString());
        updateTableData(table, toKeep);
        updateLabelN();

    }

    private void addRows(JTable table, List<String> toAdd) {
        toAdd.addAll(getAllRows(table));
        // JOptionPane.showMessageDialog(null, "toAdd:" + toAdd.toString());
        updateTableData(table, toAdd);
        updateLabelN();
    }

    private List<String> getAllRows(JTable table) {
        // get existing rows
        List<String> temp = new ArrayList<>();
        for (int i = 0; i < table.getRowCount(); i++) {
            String thisRow = "";
            thisRow = (String) table.getValueAt(i, 0);
            temp.add(thisRow);
        }
        return temp;
    }

    private void updateLabelN() {

//        int n1 = getAllRows(tableGrp1).size();
//        lblN1.setText("n=" + String.valueOf(n1));
//        int n2 = getAllRows(tableGrp2).size();
//        lblN2.setText("n=" + String.valueOf(n2));

    }

    /**
     * @author urmi bring the matched items to top and set them as selected
     * @param res
     * @param tab
     */
    private void setSelectedRows(List<String> res, JTable tab) {
        DefaultTableModel model = (DefaultTableModel) tab.getModel();
        List<String> newVals = new ArrayList<>();
        // bring matched values at top
        String temp;
        int total_matches = 0;
        for (int c = 0; c < model.getRowCount(); c++) {
            temp = model.getValueAt(c, 0).toString();
            if (res.contains(temp)) {
                newVals.add("\t:::" + temp);
                total_matches++;
            } else {
                newVals.add("~:::" + temp);
            }
        }
        java.util.Collections.sort(newVals);
        // JOptionPane.showMessageDialog(null, newVals.toString());
        model.setRowCount(0);
        for (int i = 0; i < newVals.size(); i++) {
            Vector<String> v = new Vector<>();
            v.add(newVals.get(i).split(":::")[1]);
            model.addRow(v);

        }
        if (total_matches > 0) {
            tab.setRowSelectionInterval(0, total_matches - 1);
        }
    }

    private List<String> showSearchMetadataPanel() {
        // search datacolumns by metadata and add results
        // display query panel
        try {
            final TreeSearchQueryConstructionPanel tsp = new TreeSearchQueryConstructionPanel(myProject, false);
            final Metadata.MetadataQuery[] queries;
            queries = tsp.showSearchDialog();
            if (tsp.getQueryCount() <= 0) {
                // System.out.println("Search dialog cancelled");
                // User didn't enter any queries
                return null;
            }
            // final int[] result = new int[myProject.getDataColumnCount()];
            Collection<Integer> result = new ArrayList<>();
            new AnimatedSwingWorker("Searching...", true) {
                @Override
                public Object construct() {
                    try {
                        ArrayList<Integer> toAdd = new ArrayList<Integer>();
                        for (int i = 0; i < myProject.getDataColumnCount(); i++) {
                            toAdd.add(i);
                        }
                        Integer[] hits = myProject.getMetadataHybrid().search(queries, tsp.matchAll());

                        // remove excluded cols from list
                        // urmi
//                        boolean[] excluded = excludedCopy;
//                        if (excluded != null) {
//                            List<Integer> temp = new ArrayList<>();
//                            for (Integer i : hits) {
//                                if (!excluded[i]) {
//                                    temp.add(i);
//                                }
//                            }
//                            hits = new Integer[temp.size()];
//                            hits = temp.toArray(hits);
//                        }
                        int index;
                        for (index = 0; index < hits.length; index++) {
                            result.add(hits[index]);
                            toAdd.remove(hits[index]);
                        }

                        return null;
                    }
                    catch(Exception e) {
                        return null;
                    }
                }
            }.start();

            // category
            if (result.size() < 1) {
                JOptionPane.showMessageDialog(null, "No hits found", "No hits", JOptionPane.INFORMATION_MESSAGE);
                return null;
            } else {
                List<String> hitsColumns = new ArrayList<>();
                // get datacolumn names
                for (int i : result) {
                    hitsColumns.add(myProject.getDataColumnHeader(i));
                }
                return hitsColumns;
            }

        }
        catch(Exception e) {
            JOptionPane.showMessageDialog(null, "No hits found", "No hits", JOptionPane.INFORMATION_MESSAGE);
            return null;
        }
    }
}
