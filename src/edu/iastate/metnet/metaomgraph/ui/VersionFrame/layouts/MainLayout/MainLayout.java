package edu.iastate.metnet.metaomgraph.ui.VersionFrame.layouts.MainLayout;

import edu.iastate.metnet.simpleui.*;
import edu.iastate.metnet.simpleui.components.*;

public class MainLayout extends AbstractLayout {
    public ISimpleConstraint constraint() {
        return new StackConstraint(0, 10);
    }

    @Override
    public AbstractComponent[] components() {
        return new AbstractComponent[]{
                new Header("UNSUPPORTED JAVA VERSION"),
                new Label("You using JDK 1, but JDK 17+ is the oldest Java framework we support."),
                new Label("We can download and setup JDK 17 for you automatically or you can choose to do \n" +
                        "so manually. If you choose to allow us to manually setup JDK 17, the binaries will be\n" +
                        "stored relative to MetaOmGraph.")
        };
    }
}
