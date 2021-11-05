package edu.iastate.metnet.metaomgraph;

import com.formdev.flatlaf.FlatDarkLaf;
import com.formdev.flatlaf.FlatLightLaf;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.plaf.ColorUIResource;
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
    private final String logoPath = "/resource/MetaOmicon.png";

    //JFrame
    private JFrame frame;
    private JLabel mogLogo;
    private JRadioButton defaultMemoryRadio;
    private JRadioButton extraMemoryRadio;
    private ButtonGroup memoryGroup;
    private JTextField minMemoryBox;
    private JTextField maxMemoryBox;
    private JScrollPane messagePane;
    private JTextArea messageLog;
    private JButton runButton;
    private JLabel text;

    //commands
    private final String EXTRA_MEMORY = "extra memory";
    private final String DEFAULT_MEMORY = "default memory";
    private final String RUN = "run";
    private final String SET_MIN_MEMORY = "set min memory";
    private final String SET_MAX_MEMORY = "set max memory";

    //launch variables
    private boolean extraMemory;
    private int minMemory;
    private int maxMemory;

    //debugger log on/off
    private boolean loggerOn = false;

    private MetaOmGraphLauncher () {
        // build window
        frame = new JFrame("MetaOmGraph Launcher");
        frame.setLayout(new GridBagLayout());
        try {
            UIManager.setLookAndFeel(new FlatLightLaf());
        } catch (UnsupportedLookAndFeelException e) {
            e.printStackTrace();
        }

        // gridbag formatting
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.weighty = 1;
        gbc.weightx = 1;

        //logo constraints
        gbc.gridx = 3; // x pos
        gbc.gridy = 0; // y pos
        gbc.gridwidth = 1; // number of grids to fill horizontally
        gbc.anchor = GridBagConstraints.CENTER;

        // logo
        try {
            BufferedImage image = ImageIO.read(this.getClass().getResourceAsStream(logoPath));
            mogLogo = new JLabel(new ImageIcon(image));
            frame.add(mogLogo, gbc);
            frame.setIconImage(image);
        } catch (IOException e) {
            e.printStackTrace();
            log(e.getMessage());
        }
        // Logo Text
        text = new JLabel("MetaOmGraph");
        gbc.gridy++;
        gbc.gridx = 3;
        gbc.gridwidth = 2;
        frame.add(text, gbc);

        // default memory radio button
        defaultMemoryRadio = new JRadioButton("Use Default Memory");
        defaultMemoryRadio.setActionCommand(DEFAULT_MEMORY);
        defaultMemoryRadio.addActionListener(this);
        gbc.gridy++;
        gbc.gridx = 0;
        gbc.gridwidth = 7;
        gbc.anchor = GridBagConstraints.WEST;
        frame.add(defaultMemoryRadio, gbc);

        // extra memory radio button
        extraMemoryRadio = new JRadioButton("Use Extra Memory");
        extraMemoryRadio.setActionCommand(EXTRA_MEMORY);
        extraMemoryRadio.addActionListener(this);
        gbc.gridy++;
        gbc.gridx = 0;
        gbc.gridwidth = 7;
        gbc.anchor = GridBagConstraints.WEST;
        frame.add(extraMemoryRadio, gbc);

        // group buttons
        memoryGroup = new ButtonGroup();
        memoryGroup.add(defaultMemoryRadio);
        memoryGroup.add(extraMemoryRadio);

        // min label
        text = new JLabel("Min:");
        gbc.gridy++;
        gbc.gridx = 0;
        gbc.gridwidth = 1;
        gbc.anchor = GridBagConstraints.EAST;
        frame.add(text, gbc);
        // min memory
        minMemoryBox = new JTextField(5);
        minMemoryBox.setHorizontalAlignment(JTextField.CENTER);
        minMemoryBox.setActionCommand(SET_MIN_MEMORY);
        minMemoryBox.addActionListener(this);
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.CENTER;
        frame.add(minMemoryBox, gbc);
        // min label GB
        text = new JLabel("GB");
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.WEST;
        frame.add(text, gbc);

        // max label
        text = new JLabel("Max:");
        gbc.gridx = 4;
        gbc.anchor = GridBagConstraints.EAST;
        frame.add(text, gbc);
        // max memory
        maxMemoryBox = new JTextField(5);
        maxMemoryBox.setHorizontalAlignment(JTextField.CENTER);
        maxMemoryBox.setActionCommand(SET_MAX_MEMORY);
        maxMemoryBox.addActionListener(this);
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.CENTER;
        frame.add(maxMemoryBox, gbc);
        // mac label GB
        text = new JLabel("GB");
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.EAST;
        frame.add(text, gbc);

        // message log
        if (loggerOn) {
            messageLog = new JTextArea("=== Message Log ===");
            messageLog.setRows(5);
            messagePane = new JScrollPane(messageLog);
            gbc.gridx = 0;
            gbc.gridy++;
            gbc.gridwidth = 6;
            gbc.anchor = GridBagConstraints.CENTER;
            frame.add(messagePane, gbc);
        }

        // run program button
        runButton = new JButton("Run");
        runButton.setActionCommand(RUN);
        runButton.addActionListener(this);
        gbc.gridy++;
        gbc.gridx = 3;
        gbc.gridwidth = 1;
        gbc.anchor = GridBagConstraints.CENTER;
        frame.add(runButton, gbc);

        // init launch variables
        initVariables();

        // open launcher
        init();
    }

    private void initVariables() {
        extraMemory = false;
        minMemory = 0;
        maxMemory = 0;
        defaultMemoryRadio.setSelected(true);
        defaultMemory();
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
            case EXTRA_MEMORY: {
                extraMemory();
                break;
            }
            case RUN:
            {
                runProgram();
                break;
            }
            case SET_MIN_MEMORY: {
                setMinMemory();
                break;
            }
            case SET_MAX_MEMORY: {
                setMaxMemory();
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
        minMemoryBox.setEnabled(false);
        maxMemoryBox.setEnabled(false);
        minMemory = 2;
        maxMemory = 6;
        minMemoryBox.setText(Integer.toString(minMemory));
        maxMemoryBox.setText(Integer.toString(maxMemory));
        log("extra heap memory: " + Boolean.toString(extraMemory));
    }

    private void extraMemory() {
        extraMemory = true;
        minMemoryBox.setEnabled(true);
        maxMemoryBox.setEnabled(true);
        minMemoryBox.setText(Integer.toString(minMemory));
        maxMemoryBox.setText(Integer.toString(maxMemory));
        log("extra heap memory: " + Boolean.toString(extraMemory));
    }

    // sets min memory
    // needs error checking - aorgler
    private void setMinMemory() {
        String input = minMemoryBox.getText();
        for(char c : input.toCharArray()) {
            if (!Character.isDigit(c)) {
                return;
            }
        }
        minMemory = Integer.parseInt(input);
        log("Set min: " +  minMemory + "gb");
    }

    //
    private void setMaxMemory() {
        String input = maxMemoryBox.getText();
        for(char c : input.toCharArray()) {
            if (!Character.isDigit(c)) {
                return;
            }
        }
        maxMemory = Integer.parseInt(input);
        log("Set min: " +  maxMemory + "gb");
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
            setMinMemory();
            String minHeap = "-Xms" + Integer.toString(minMemory) + "g";
            flags.add(minHeap);

            setMaxMemory();
            String maxHeap = "-Xmx" + Integer.toString(maxMemory) + "g";
            flags.add(maxHeap);
        }

        for (String f: flags) {
            cmd += f + " ";
        }
        cmd += "-jar " + jarPath; // add jar path to command
        try {
            // log for debugging
            log(cmd);
            Process pr = Runtime.getRuntime().exec(cmd, null, file); // run command on command line in target dir
            // printResults(pr); // print run results for debugging

        } catch (IOException e) {
            e.printStackTrace();
            log(e.getMessage());
        }
    }

    // log messages in launcher
    private void log(String message) {
        if (loggerOn) {
            this.messageLog.append("\n" + message);

            // scroll to the bottom of the pane
            messageLog.setCaretPosition(messageLog.getDocument().getLength());
        }
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
