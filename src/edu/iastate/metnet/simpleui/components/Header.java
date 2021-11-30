package edu.iastate.metnet.simpleui.components;

import edu.iastate.metnet.simpleui.AbstractComponent;

import javax.swing.*;
import java.awt.*;

public class Header extends Label {

    public Header(String title) {
        super(title);
    }

    @Override
    public Container create() {
        return new edu.iastate.metnet.metaomgraph.ui.Header(this.m_textContent);
    }
}
