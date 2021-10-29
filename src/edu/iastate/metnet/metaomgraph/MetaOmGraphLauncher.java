package edu.iastate.metnet.metaomgraph;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
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
    private JLabel mogLogo;
    private JRadioButton defaultMemoryRadio;
    private JRadioButton extraMemoryRadio;
    private ButtonGroup memoryGroup;
    private JTextArea minMemoryBox;
    private JTextArea maxMemoryBox;
    private JScrollPane messagePane;
    private JTextArea messageLog;
    private JButton runButton;

    //commands
    private final String EXTRA_MEMORY = "extra memory";
    private final String DEFAULT_MEMORY = "default memory";
    private final String RUN = "run";

    //launch variables
    private boolean extraMemory;
    private int minMemory;
    private int maxMemory;

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

        // logo
        try {
            BufferedImage image = ImageIO.read(new File("/resource/MetaOmicon.png"));
            mogLogo = new JLabel(new ImageIcon(image));
            frame.add(mogLogo, gbc);
        } catch (IOException e) {
            e.printStackTrace();
            log(e.getMessage());
        }

        // default memory radio button
        defaultMemoryRadio = new JRadioButton("Use Default Memory");
        defaultMemoryRadio.setActionCommand(DEFAULT_MEMORY);
        defaultMemoryRadio.addActionListener(this);
        gbc.gridy++;
        gbc.gridx = 0;
        gbc.anchor = GridBagConstraints.WEST;
        frame.add(defaultMemoryRadio, gbc);

        // extra memory radio button
        extraMemoryRadio = new JRadioButton("Use Extra Memory");
        extraMemoryRadio.setActionCommand(EXTRA_MEMORY);
        extraMemoryRadio.addActionListener(this);
        gbc.gridy++;
        gbc.gridx = 0;
        gbc.anchor = GridBagConstraints.WEST;
        frame.add(extraMemoryRadio, gbc);

        // group buttons
        memoryGroup = new ButtonGroup();
        memoryGroup.add(defaultMemoryRadio);
        memoryGroup.add(extraMemoryRadio);

        // min memory
        minMemoryBox = new JTextArea("Min", 1, 5);
        gbc.gridy++;
        gbc.anchor = GridBagConstraints.NORTH;
        frame.add(minMemoryBox, gbc);

        // max memory
        maxMemoryBox = new JTextArea("Max", 1, 5);
        gbc.gridy++;
        gbc.anchor = GridBagConstraints.NORTH;
        frame.add(maxMemoryBox, gbc);

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
        gbc.gridx = 1;
        gbc.gridwidth = 1;
        gbc.gridy++;
        gbc.anchor = GridBagConstraints.SOUTHWEST;
        frame.add(runButton, gbc);

        // init launch variables
        initVariables();

        // open launcher
        init();
    }

    private void initVariables() {
        defaultMemory();
        minMemory = 0;
        maxMemory = 0;
    }

    // start the program
    private void init() {
        frame.setSize(400, 300);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        String com = e.getActionCommand();

        switch (com) {
            case DEFAULT_MEMORY:
            {
                defaultMemory();
                break;
            }
            case EXTRA_MEMORY:
            {
                extraMemory();
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

    private void defaultMemory() {
        extraMemory = false;
    }

    private void extraMemory() {
        extraMemory = true;
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
