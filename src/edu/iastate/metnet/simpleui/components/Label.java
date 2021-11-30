package edu.iastate.metnet.simpleui.components;

import edu.iastate.metnet.simpleui.AbstractComponent;

import javax.swing.*;
import java.awt.*;

public class Label extends AbstractComponent {
    protected String m_textContent;

    public Label(String title) {
        this.m_textContent = title;
    }

    @Override
    public Container create() {
        return new JLabel(this.m_textContent);
    }
}
