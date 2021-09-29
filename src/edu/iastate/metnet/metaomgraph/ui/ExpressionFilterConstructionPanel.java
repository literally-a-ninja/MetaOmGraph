package edu.iastate.metnet.metaomgraph.ui;

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

        max = new JTextField();
        max.setBounds(150, 25, 200, 30);
        maxDesc = new JLabel("Maximum Value:", SwingConstants.RIGHT);
        maxDesc.setBounds(25, 25, 100, 30);
        min = new JTextField();
        min.setBounds(150, 75, 200, 30);
        minDesc = new JLabel("Minimum Value:", SwingConstants.RIGHT);
        minDesc.setBounds(25, 75, 100, 30);
        filter = new JButton("Filter");
        filter.setBounds(125, 130, 75, 50);
        filter.setEnabled(false);
        filter.setActionCommand("filter");
        filter.addActionListener(this);
        cancel = new JButton("Cancel");
        cancel.setBounds(225, 130, 75, 50);
        cancel.setActionCommand("cancel");
        cancel.addActionListener(this);

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

        window.add(max);
        window.add(maxDesc);
        window.add(min);
        window.add(minDesc);
        window.add(filter);
        window.add(cancel);
        window.getContentPane().setLayout(null);
        window.setSize(450, 225);
        window.setVisible(true);
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
