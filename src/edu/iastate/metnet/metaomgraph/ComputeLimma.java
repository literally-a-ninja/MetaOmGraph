package edu.iastate.metnet.metaomgraph;

import javax.script.ScriptEngine;
import javax.script.ScriptException;
import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.*;

import org.renjin.primitives.matrix.DoubleMatrixBuilder;
import org.renjin.primitives.matrix.StringMatrixBuilder;
import org.renjin.script.*;

public class ComputeLimma {

    private final String selectedList;
    private final HashMap<String, String> groups;
    private final MetaOmProject myProject;
    private final ArrayList<Integer> limmaCountsInd;

    public ComputeLimma(String selectedList, HashMap<String, String> groups, MetaOmProject myProject,  ArrayList<Integer> limmaCountsInd) {
        this.selectedList = selectedList;
        this.groups = groups;
        this.myProject = myProject;
        this.limmaCountsInd = limmaCountsInd;
    }

    public void calc() throws ScriptException, IOException, URISyntaxException {
        int[] selected = myProject.getGeneListRowNumbers(this.selectedList);
        double[][] thisData;
        String[] rowHeaders;
        String[] columnHeaders;

        StringMatrixBuilder rGroups = new StringMatrixBuilder(2, limmaCountsInd.size());

        // create a script engine manager:
        RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
        // create a Renjin engine:
        ScriptEngine engine = factory.getScriptEngine();

        int[] countsInd = limmaCountsInd.stream().mapToInt(i->i).toArray();

        // Start of Limma analysis
        thisData = myProject.getSelectedListOnlyRowData(countsInd, selectedList);

        rowHeaders = myProject.getDefaultRowNames(selected);
        columnHeaders = myProject.getDataColumnHeaders(countsInd);

        DoubleMatrixBuilder rData = new DoubleMatrixBuilder(selected.length, limmaCountsInd.size());
        for (int row = 0; row < selected.length; row++) {
            for (int col = 0; col < limmaCountsInd.size(); col++) {
                rData.setValue(row, col, thisData[col][row]);
            }
        }
        rData.setColNames(Arrays.asList(columnHeaders));
        rData.setRowNames(Arrays.asList(rowHeaders));

        for (int i = 0; i < limmaCountsInd.size(); i++) {
            rGroups.setValue(0, i, columnHeaders[i]);
            rGroups.setValue(1, i, groups.get(columnHeaders[i]));
        }
        rGroups.setRowNames(Arrays.asList("Samples", "Groups"));

        // Preprocessing
        engine.eval("library(edgeR)");

        engine.put("counts", rData.build());

        // CONVERTING THE COUNT MATRIX TO A DGELIST OBJECT AND CALCULATING THE NORMALIZATION FACTORS

        engine.eval("d <- DGEList(counts)");
        engine.eval("d <- calcNormFactors(d)");

        // READING GROUPS AND TODO PLOTTING MDS CHART

        String homePath = System.getProperty("user.home");
        engine.eval("snames <- colnames(counts)");
        engine.put("groups", rGroups.build());
        engine.eval("group2 <- interaction(groups[2,])");
        String mdsPath = "png(\"" + homePath + "/metaomgraph/mds.png\")";
//            engine.eval(mdsPath);
//            engine.eval("plotMDS(d, col = as.numeric(group))");
//            engine.eval("dev.off()");

        // PERFORMING THE LIMMA VOOM OPERATION
        String voomPath = "png(\"" + homePath + "/metaomgraph/voom.png\")";
        engine.eval("mm <- model.matrix(~0 + group2)");
        engine.eval(voomPath);
        engine.eval("y <- voom(d, mm, plot = T)");
        engine.eval("dev.off()");

        // Finding differentially expressed genes
        engine.eval("fit <- lmFit(y, mm)");
        engine.eval("fitbit <- colnames(coef(fit))");

        engine.eval("contr <- makeContrasts(group21 - group22, levels = colnames(coef(fit)))");

        engine.eval("tmp <- contrasts.fit(fit, contr)");
        engine.eval("tmp <- eBayes(tmp)");

        engine.eval("top.table <- topTable(tmp, coef=1, sort.by = \"P\", n = Inf)");
        engine.eval("top.table$Gene <- rownames(top.table)");
        engine.eval("top.table <- top.table[,c(\"Gene\", names(top.table)[1:6])]");

        String writeDETable = "write.table(top.table, file = \"" + homePath + "/metaomgraph/differentialexpression.tsv" + "\", row.names = F, sep = \"\t\", quote = F)";
        engine.eval(writeDETable);

//            StringVector genes = (StringVector)engine.eval("top.table$Gene");
//            DoubleVector logFC = (DoubleVector)engine.eval("top.table$logFC");
//            DoubleVector aveExpr = (DoubleVector)engine.eval("top.table$AveExpr");
//            DoubleVector t = (DoubleVector)engine.eval("top.table$t");
//            DoubleVector pVal = (DoubleVector)engine.eval("top.table$P.Value");
//            DoubleVector adjPVal = (DoubleVector)engine.eval("top.table$adj.P.Val");
    }
}
