package edu.iastate.metnet.metaomgraph.ui;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Rectangle2D;

public class Header extends JLabel {
    String text;

    public Header(String text) {
        this(text, 18.0f);
    }

    public Header(String text, float size) {
        this.text = text;
        this.setFont(this.getFont().deriveFont(size));
    }

    @Override
    public Dimension getPreferredSize() {
        Graphics g = this.getGraphics();

        if (null == g)
            return new Dimension(0, 0);

        FontMetrics metrics = this.getComponentGraphics(g).getFontMetrics(this.getFont());

        Dimension bbox = metrics.getMaxCharBounds(g).getBounds().getSize();

        bbox.setSize(bbox.width, bbox.height + 10);
        return bbox;
    }

    @Override
    protected void paintComponent(Graphics g) {
        Dimension bbox = this.getPreferredSize();
//        GradientPaint gradient = new GradientPaint(0.0F, 0.0F, Color.WHITE, 0.0F,
//                bbox.width, new Color(255, 255, 255, 0));

        Graphics2D g2 = (Graphics2D) g.create();
        g2.setPaint(Color.BLACK);
        g2.drawString(text, 0, this.getHeight());
    }
}
