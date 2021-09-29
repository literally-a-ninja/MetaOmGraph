package edu.iastate.metnet.metaomgraph.ui;

import edu.iastate.metnet.metaomgraph.MetaOmProject;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class ExpressionFilterConstructionPanel extends JPanel
        implements ActionListener {

    private JFrame window = new JFrame("Expression Filter");
    private JTextField max;
    private JLabel maxDesc;
    private JTextField min;
    private JLabel minDesc;
    private JPanel queryPanel;
    private JButton filter;
    private JButton cancel;

    private MetaOmProject project;

    public ExpressionFilterConstructionPanel(MetaOmProject project) {
        this.project = project;
        setLayout(new BorderLayout());
        queryPanel = new JPanel();
        queryPanel.setLayout(new BoxLayout(queryPanel, 1));
        JPanel queryButtonPanel = new JPanel();

        max = new JTextField();
        maxDesc = new JLabel("Maximum Value:", SwingConstants.RIGHT);
        min = new JTextField();
        minDesc = new JLabel("Minimum Value:", SwingConstants.RIGHT);
        filter = new JButton("Filter");
        filter.setEnabled(false);
        filter.setActionCommand("filter");
        filter.addActionListener(this);
        cancel = new JButton("Cancel");
        cancel.setActionCommand("cancel");
        cancel.addActionListener(this);

        window.add(max);
        window.add(maxDesc);
        window.add(min);
        window.add(minDesc);
        window.add(filter);
        window.add(cancel);
        window.getContentPane().setLayout(null);
        window.setSize(450, 250);
        window.setVisible(true);
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        if ("filter".equals(e.getActionCommand())) {
            try {
                Double maxValue = Double.parseDouble(max.getText());
            } catch (NumberFormatException ne) {
                JOptionPane.showMessageDialog(window, "Maximum field must be a number");
            }
            try {
                Double minValue = Double.parseDouble(min.getText());
            } catch (NumberFormatException ne) {
                JOptionPane.showMessageDialog(window, "Minimum field must be a number");
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
