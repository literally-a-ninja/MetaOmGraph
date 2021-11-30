package edu.iastate.metnet.metaomgraph.ui.VersionFrame;

import edu.iastate.metnet.metaomgraph.ui.VersionFrame.layouts.MainLayout.MainLayout;
import edu.iastate.metnet.simpleui.LayoutFactory;
import org.apache.logging.log4j.core.Layout;

import javax.swing.*;
import java.awt.*;

public class Controller {
    final protected LayoutFactory factory;

    Controller() {
        factory = new LayoutFactory(new MainLayout());
    }

    public void dialog() {
        factory.make("Java version mismatch!");
    }

    public static void main(String[] args) {
        Controller controller = new Controller();
        controller.dialog();
    }

}
