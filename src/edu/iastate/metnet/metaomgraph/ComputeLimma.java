package edu.iastate.metnet.metaomgraph;

import com.github.rcaller.datatypes.DataFrame;
import com.github.rcaller.rstuff.RCaller;
import com.github.rcaller.rstuff.RCode;
import javax.script.ScriptEngine;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URISyntaxException;
import java.util.*;
import java.util.stream.Collectors;

import org.renjin.primitives.matrix.DoubleMatrixBuilder;
import org.renjin.primitives.matrix.Matrix;
import org.renjin.primitives.matrix.StringMatrixBuilder;
import org.renjin.script.*;
import org.renjin.sexp.DoubleArrayVector;
import org.renjin.sexp.DoubleVector;
import org.renjin.sexp.ListVector;
import org.renjin.sexp.StringVector;


public class ComputeLimma {

    private String selectedList;
    private HashMap<String, String> groups;
    private MetaOmProject myProject;
    private ArrayList<Integer> limmaCountsInd;

    public ComputeLimma(String selectedList, HashMap<String, String> groups, MetaOmProject myProject,  ArrayList<Integer> limmaCountsInd) {
        this.selectedList = selectedList;
        this.groups = groups;
        this.myProject = myProject;
        this.limmaCountsInd = limmaCountsInd;
    }

    public void calc() {
        int[] selected = myProject.getGeneListRowNumbers(this.selectedList);
        double[][] thisData;
        ArrayList<ArrayList<Double>> data = new ArrayList<>();
        String[] rowHeaders = new String[selected.length];
        String[] columnHeaders = new String[limmaCountsInd.size()];

        DoubleMatrixBuilder rData = new DoubleMatrixBuilder(limmaCountsInd.size(), selected.length);
        StringMatrixBuilder rGroups = new StringMatrixBuilder(2, limmaCountsInd.size());

        // create a script engine manager:
        RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
        // create a Renjin engine:
        ScriptEngine engine = factory.getScriptEngine();

        StringWriter outputWriter = new StringWriter();
        engine.getContext().setWriter(outputWriter);

        int[] example = limmaCountsInd.stream().mapToInt(i->i).toArray();

        // Start of Limma analysis
        try {
            thisData = myProject.getSelectedListOnlyRowData(example, selectedList);
            rowHeaders = myProject.getDefaultRowNames(selected);
            columnHeaders = myProject.getDataColumnHeaders(example);

            for (int row = 0; row < limmaCountsInd.size(); row++) {
                for (int col = 0; col < selected.length; col++) {
                    rData.setValue(row, col, thisData[row][col]);
                }
            }
            rData.setColNames(Arrays.asList(rowHeaders));
            rData.setRowNames(Arrays.asList(columnHeaders));

            for (int i = 0; i < limmaCountsInd.size(); i++) {
                rGroups.setValue(0, i, columnHeaders[i]);
                rGroups.setValue(1, i, groups.get(columnHeaders[i]));
            }
            rGroups.setRowNames(Arrays.asList("Samples", "Groups"));

            // Preprocessing
            engine.eval("library(edgeR)");

            engine.put("counts", rData.build());

            engine.eval("d <- DGEList(counts)");
            engine.eval("d <- calcNormFactors(d)");

            engine.eval("snames <- colnames(counts)");

            engine.put("groups", rGroups.build());
            engine.eval("group2 <- interaction(groups[2,])");

            /////////////
            engine.eval("countsx <- read.table(\"/home/fahmi/3032_3_sd6_mog/rscripts/limma_counts.csv\", header = TRUE, sep = \",\")");
            engine.eval("counts2x <- countsx[,-1]");
            engine.eval("rownames(counts2x) <- countsx[,1]");
            engine.eval("countsx <- data.matrix(counts2x)");
            engine.eval("d0x <- DGEList(countsx)");
            engine.eval("d0x <- calcNormFactors(d0x)");
            engine.eval("snamesx <- colnames(countsx) # Sample names");
            engine.eval("groupsx <- read.table(\"/home/fahmi/3032_3_sd6_mog/rscripts/limma_groups.csv\", header = TRUE, sep = \",\")");
            engine.eval("groupx <- interaction(groupsx['Groups'])");


            /////////////

            engine.eval("mm <- model.matrix(~0 + group2)");
            engine.eval("y <- voom(d, mm, plot = F)");

            engine.eval("fit <- lmFit(y, mm)");
            engine.eval("contr <- makeContrasts(groupgroup1 - groupgroup2, levels = colnames(coef(fit)))");

            engine.eval("tmp <- contrasts.fit(fit, contr)");
            engine.eval("tmp <- eBayes(tmp)");

            engine.eval("top.table <- topTable(tmp, coeff=1, sort.by = \"P\", n = Inf)");
            engine.eval("top.table$Gene <- rownames(top.table)");
            engine.eval("top.table <- top.table[,c(\"Gene\", names(top.table)[1:6])]");

            String currentPath = System.getProperty("user.dir");
            String dePath = currentPath + "/rscripts/de.tsv";
            String writeDETable = "write.table(top.table, file = \"" + dePath + "\", row.names = F, sep = \"\t\", quote = F)";
            engine.eval(writeDETable);

            StringVector genes = (StringVector)engine.eval("top.table$Gene");
            DoubleVector logFC = (DoubleVector)engine.eval("top.table$logFC");
            DoubleVector aveExpr = (DoubleVector)engine.eval("top.table$AveExpr");
            DoubleVector t = (DoubleVector)engine.eval("top.table$t");
            DoubleVector pVal = (DoubleVector)engine.eval("top.table$P.Value");
            DoubleVector adjPVal = (DoubleVector)engine.eval("top.table$adj.P.Val");

            System.out.println("Done!");

        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }
}
