package edu.iastate.metnet.metaomgraph.ui;

import edu.iastate.metnet.metaomgraph.MetaOmGraph;
import edu.iastate.metnet.metaomgraph.MetaOmProject;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class ExpressionFilterConstructionPanel extends JPanel
        implements ActionListener {

    private JFrame window = new JFrame("Expression Filter");
    private JPanel queryViewport;
    private JTextField max;
    private JLabel maxDesc;
    private JTextField min;
    private JLabel minDesc;
    private JPanel maxPanel;
    private JPanel minPanel;
    private JButton filter;
    private JButton cancel;
    private JPanel buttonPanel;
    private boolean minPopulated = false;
    private boolean maxPopulated = false;

    private MetaOmProject project;

    public ExpressionFilterConstructionPanel(MetaOmProject project) {
        this.project = project;
        Container mypane = window.getContentPane();
        window.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();
        c.gridx = 0;
        c.gridy = 0;
        c.anchor = GridBagConstraints.CENTER;

        // Align Labels
        JPanel entryPanel = new JPanel(new GridBagLayout());
        GridBagConstraints c2 = new GridBagConstraints();
        c2.insets = new Insets(0,0,10,0);
        setMyConstraints(c2,0,0,GridBagConstraints.EAST);
        maxDesc = new JLabel("Maximum Value: ");
        maxDesc.setPreferredSize(new Dimension(100, 30));
        entryPanel.add(maxDesc,c2);

        c2.insets = new Insets(0,0,10,0);
        setMyConstraints(c2,1,0,GridBagConstraints.WEST);
        max = new JTextField();
        max.setPreferredSize(new Dimension(200, 30));
        entryPanel.add(max,c2);

        setMyConstraints(c2,0,1,GridBagConstraints.EAST);
        minDesc = new JLabel("Minimum Value: ");
        minDesc.setPreferredSize(new Dimension(100, 30));
        entryPanel.add(minDesc, c2);

        setMyConstraints(c2,1,1,GridBagConstraints.WEST);
        min = new JTextField();
        min.setPreferredSize(new Dimension(200, 30));
        entryPanel.add(min,c2);

        mypane.add(entryPanel, c);

        // Align Buttons
        buttonPanel = new JPanel(new GridBagLayout());
        GridBagConstraints c3 = new GridBagConstraints();

        c3.insets = new Insets(0,0,0,5);
        setMyConstraints(c3,0,0,GridBagConstraints.EAST);
        filter = new JButton("Filter");
        filter.setEnabled(false);
        filter.setActionCommand("filter");
        filter.addActionListener(this);
        filter.setPreferredSize(new Dimension(90,50));
        buttonPanel.add(filter,c3);

        setMyConstraints(c3,1,0,GridBagConstraints.WEST);
        cancel = new JButton("Cancel");
        cancel.setActionCommand("cancel");
        cancel.setPreferredSize(new Dimension(90,50));
        cancel.addActionListener(this);
        buttonPanel.add(cancel,c3);

        c.gridx = 0;
        c.gridy = 1;
        c.insets = new Insets(10, 0,0,0);
        c.anchor = GridBagConstraints.CENTER;
        mypane.add(buttonPanel, c);

        max.getDocument().addDocumentListener(new DocumentListener() {
            public void changedUpdate(DocumentEvent e) {
                changed();
            }
            public void removeUpdate(DocumentEvent e) {
                changed();
            }
            public void insertUpdate(DocumentEvent e) {
                changed();
            }

            public void changed() {
                if (!max.getText().equals("")) {
                    maxPopulated = true;
                    if (minPopulated) {
                        filter.setEnabled(true);
                    }
                }
                else {
                    maxPopulated = false;
                    filter.setEnabled(false);
                }
            }
        });

        min.getDocument().addDocumentListener(new DocumentListener() {
            public void changedUpdate(DocumentEvent e) {
                changed();
            }
            public void removeUpdate(DocumentEvent e) {
                changed();
            }
            public void insertUpdate(DocumentEvent e) {
                changed();
            }

            public void changed() {
                if (!min.getText().equals("")) {
                    minPopulated = true;
                    if (maxPopulated) {
                        filter.setEnabled(true);
                    }
                }
                else {
                    minPopulated = false;
                    filter.setEnabled(false);
                }

            }
        });

        // Align Frame
        window.setSize(450, 225);
        int width = MetaOmGraph.getMainWindow().getWidth();
        int height = MetaOmGraph.getMainWindow().getHeight();
        // align window to the middle of the screen
        window.setLocation((width - window.getWidth()) / 2, (height - window.getHeight()) / 2);
        window.setVisible(true);
    }
    private static void setMyConstraints(GridBagConstraints c, int gridx, int gridy, int anchor) {
        c.gridx = gridx;
        c.gridy = gridy;
        c.anchor = anchor;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        boolean maxError = false;
        boolean minError = false;
        if ("filter".equals(e.getActionCommand())) {
            try {
                Double maxValue = Double.parseDouble(max.getText());
            } catch (NumberFormatException ne) {
                maxError = true;

            }
            try {
                Double minValue = Double.parseDouble(min.getText());
            } catch (NumberFormatException ne) {
                minError = true;

            }
            String errorMsg = "";
            if (maxError) {
                errorMsg += "Maximum field must be a number\n";
            }
            if (minError) {
                errorMsg += "Minimum field must be a number";
            }
            if (maxError || minError) {
                JOptionPane.showMessageDialog(window, errorMsg);
            }
            //filterWithExpressionBounds(maxValue, minValue);
            return;
        }
        if ("cancel".equals(e.getActionCommand())) {
            window.dispose();
            return;
        }
    }
}
