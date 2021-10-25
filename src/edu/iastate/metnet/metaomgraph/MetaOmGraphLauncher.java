package edu.iastate.metnet.metaomgraph;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class MetaOmGraphLauncher implements ActionListener {

    // launcher variables
    private final String jarPath = "metaomgraph4-1.8.2beta.jar"; // replace with jar for current version of MOG

    //JFrame
    private JFrame frame;
    private JCheckBox extraMemoryCheckBox;
    private JScrollPane messagePane;
    private JTextArea messageLog;
    private JButton runButton;

    //commands
    private final String HEAP_MEMORY = "switch heap memory";
    private final String RUN = "run";

    //launch variables
    private static boolean extraMemory;

    private MetaOmGraphLauncher () {
        // build window
        frame = new JFrame("MetaOmGraph Launcher");
        frame.setLayout(new GridBagLayout());

        // gridbag formatting
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0; // x pos
        gbc.gridy = 0; // y pos
        gbc.weightx = 0.5; // how much extra space in the x dir should be given
        gbc.weighty = 0; // how much extra space in the y dir should be given
        gbc.gridwidth = 2; // number of grids to fill horizontally
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.NORTHWEST;


        // extra heap memory checkbox
        extraMemoryCheckBox = new JCheckBox("Run with higher heap memory");
        extraMemoryCheckBox.setActionCommand(HEAP_MEMORY);
        extraMemoryCheckBox.addActionListener(this);
        gbc.weighty = 1;
        frame.add(extraMemoryCheckBox, gbc);

        // message log
        messageLog = new JTextArea("=== Message Log ===");
        messageLog.setRows(5);
        messagePane = new JScrollPane(messageLog);
        gbc.weighty = 0;
        gbc.gridy++;
        gbc.anchor = GridBagConstraints.SOUTHWEST;
        frame.add(messagePane, gbc);

        // run program button
        runButton = new JButton("Run");
        runButton.setActionCommand(RUN);
        runButton.addActionListener(this);
        gbc.weighty = 0;
        gbc.gridy++;
        gbc.anchor = GridBagConstraints.SOUTHWEST;
        frame.add(runButton, gbc);

        // init launch variables
        initVariables();

        // open launcher
        init();
    }

    private void initVariables() {
        checkExtraHeapMemory();
    }

    // start the program
    private void init() {
        frame.setSize(800, 600);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        String com = e.getActionCommand();

        switch (com) {
            case HEAP_MEMORY:
            {
                checkExtraHeapMemory();
                break;
            }
            case RUN:
            {
                runProgram();
                break;
            }
            default:
            {
                break;
            }
        }
    }


    private void checkExtraHeapMemory() {
        if (extraMemoryCheckBox.isSelected()) {
            extraMemory = true;
        } else {
            extraMemory = false;
        }
        log("extra heap memory: " + Boolean.toString(extraMemory));
    }

    // method which launches MOG
    private void runProgram() {
        log("Starting MetaOmGraph...");

        File file = new File("target"); // jar directory
        String cmd = ""; // cmd for launching MOG

        String os = System.getProperty("os.name"); // get operating system
        if (os.toLowerCase().startsWith("windows")) { // check if windows
            cmd += "cmd /c ";
        } else {
            cmd += "sh -c ";
        }
        cmd += "java ";

        ArrayList<String> flags = new ArrayList<String>();
        // add logic for adding flags here
        if (extraMemory) {
            flags.add("-Xmx4g");
        }

        for (String f: flags) {
            cmd += f + " ";
        }
        cmd += "-jar " + jarPath; // add jar path to command
        try {
            // log for debugging
            // log(cmd);
            Process pr = Runtime.getRuntime().exec(cmd, null, file); // run command on command line in target dir
            // printResults(pr); // print run results for debugging

        } catch (IOException e) {
            e.printStackTrace();
            log(e.getMessage());
        }
    }

    // log messages in launcher
    private void log(String message) {
        this.messageLog.append("\n" + message);

        // scroll to the bottom of the pane
        messageLog.setCaretPosition(messageLog.getDocument().getLength());
    }

    // main method
    public static void main(String[] args) {
        new MetaOmGraphLauncher();
    }

    // check if process if running correctly
    public void printResults(Process process) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        String line = "";
        while ((line = reader.readLine()) != null) {
            System.out.println(line);
        }
    }
}
