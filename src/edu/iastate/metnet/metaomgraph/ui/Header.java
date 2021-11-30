package edu.iastate.metnet.metaomgraph.ui;

import javax.swing.*;
import javax.swing.border.Border;
import java.awt.*;

public class Header extends JLabel {
    public Header(String text) {
        this(text, 18.0f);
    }

    public Header(String text, float size) {
        super(text);
        this.setFont(this.getFont().deriveFont(size));
    }

    public Header(String text, Font font) {
        super(text);
        this.setFont(font);
    }

    @Override
    public Dimension getPreferredSize() {
        FontMetrics metrics        = this.getFontMetrics(this.getFont());
        int         typeWidth      = metrics.stringWidth(this.getText());
        int         componentWidth = this.getWidth();

        Insets margin = this.getBorder().getBorderInsets(this);

        int height = Math.min(this.getFont().getSize() * (componentWidth / typeWidth), this.getHeight())
                     + metrics.getMaxAscent()
                     + margin.bottom + margin.top;

        return new Dimension(componentWidth, height);
    }
}
