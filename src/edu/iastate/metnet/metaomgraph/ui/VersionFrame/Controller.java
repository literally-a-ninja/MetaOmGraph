package edu.iastate.metnet.metaomgraph.ui.VersionFrame;

import edu.iastate.metnet.metaomgraph.ui.VersionFrame.layouts.MainLayout.Manifest;
import edu.iastate.metnet.simpleui.LayoutFactory;

public class Controller {
    final protected LayoutFactory factory;

    Controller() {
        factory = new LayoutFactory(new Manifest());
    }

    public void dialog() {
        factory.make("Java version mismatch!");
    }

    public static void main(String[] args) {
        Controller controller = new Controller();
        controller.dialog();
    }

}
