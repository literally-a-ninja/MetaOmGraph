package edu.iastate.metnet.simpleui;

import org.jetbrains.annotations.NotNull;

import java.awt.*;

public class StackConstraint extends GridBagConstraints implements ISimpleConstraint {

    public StackConstraint(int gutterWidth, int gutterHeight) {
        super(0, 0,
                0, 0,
                1.5f, 1.5f,
                GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL,
                new Insets(0, 0, 0, 0),
                gutterWidth, gutterHeight);
    }

    @Override
    public void beforeInsert(final Container container, final Container preInsertItem) {
    }

    @Override
    public void afterInsert(final Container container, @NotNull final Container postInsertedItem) {

        Dimension size = postInsertedItem.getPreferredSize();

        if (this.ipadx > 0) this.gridx++;
        if (this.ipady > 0) this.gridy++;

        int height = (0 < this.ipady ? size.height : 0);
        int width = (0 < this.ipadx ? size.width : 0);

        Insets insets = this.insets;
        this.insets = new Insets(
                height + insets.top + this.ipady,
                width + insets.left + this.ipadx,
                height + insets.bottom + this.ipady,
                width + insets.right + this.ipadx
        );
    }
}
