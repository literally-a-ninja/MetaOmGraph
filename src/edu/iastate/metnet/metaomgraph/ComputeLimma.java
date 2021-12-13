package edu.iastate.metnet.metaomgraph;

import com.github.rcaller.datatypes.DataFrame;
import com.github.rcaller.rstuff.RCaller;
import com.github.rcaller.rstuff.RCode;
import javax.script.ScriptEngine;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.List;

import org.renjin.script.*;
import org.renjin.sexp.DoubleVector;
import org.renjin.sexp.StringVector;


public class ComputeLimma {

    private HashMap<String, String> groups;
    private List<String> names;
    private MetaOmProject myProject;
    private Double lowly;

    public ComputeLimma(String selectedList, HashMap<String, String> groups, MetaOmProject myProject) {
        this.groups = groups;
    }

    public static void calc() {
        // create a script engine manager:
        RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
        // create a Renjin engine:
        ScriptEngine engine = factory.getScriptEngine();
        // Start of Limma analysis
        try {
            // Preprocessing
            String currentPath = System.getProperty("user.dir");
            String countsPath = currentPath + "/rscripts/limma_counts.csv";
            String groupsPath = currentPath + "/rscripts/limma_groups.csv";
            String mdsPath = currentPath + "/rscripts/mds-plot.png";
            String voomPath = currentPath + "/rscripts/voom-plot.png";
            String dePath = currentPath + "/rscripts/de.tsv";

            String readFileCounts = "counts <- read.table(\"" + countsPath + "\", header = TRUE, sep = \",\")";
            String readFileGroups = "groups <- read.table(\"" + groupsPath + "\", header = TRUE, sep = \",\")";
            String writeDETable = "write.table(top.table, file = \"" + dePath + "\", row.names = F, sep = \"\t\", quote = F)";

            engine.eval("library(edgeR)");
            engine.eval("library(limma)");
            engine.eval(readFileCounts);
            engine.eval("counts2 <- counts[,-1]");
            engine.eval("rownames(counts2) <- counts[,1]");
            engine.eval("counts <- data.matrix(counts2)");

            engine.eval("d0 <- DGEList(counts)");
            engine.eval("d0 <- calcNormFactors(d0)");
            engine.eval("d <- d0 ");

            engine.eval("snames <- colnames(counts)");
            engine.eval(readFileGroups);
            engine.eval("group <- interaction(groups['Groups'])");

            engine.eval("mm <- model.matrix(~0 + group)");
            engine.eval("y <- voom(d, mm, plot = F)");

            engine.eval("fit <- lmFit(y, mm)");
            engine.eval("contr <- makeContrasts(groupgroup1 - groupgroup2, levels = colnames(coef(fit)))");

            engine.eval("tmp <- contrasts.fit(fit, contr)");
            engine.eval("tmp <- eBayes(tmp)");

            engine.eval("top.table <- topTable(tmp, coeff=1, sort.by = \"P\", n = Inf)");
            engine.eval("top.table$Gene <- rownames(top.table)");
            engine.eval("top.table <- top.table[,c(\"Gene\", names(top.table)[1:6])]");

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

    public void doCalc() {
        RCaller caller = RCaller.create();
        RCode code = RCode.create();
        // Start of Limma analysis
        try {
            // Preprocessing
            String currentPath = System.getProperty("user.dir");
            String countsPath = currentPath + "/rscripts/limma_counts.csv";
            String groupsPath = currentPath + "/rscripts/limma_groups.csv";
            String mdsPath = currentPath + "/rscripts/mds-plot.png";
            String voomPath = currentPath + "/rscripts/voom-plot.png";
            String dePath = currentPath + "/rscripts/de.tsv";

            String readFileCounts = "counts <- read.table(\"" + countsPath + "\", header = TRUE, sep = \",\")";
            String readFileGroups = "groups <- read.table(\"" + groupsPath + "\", header = TRUE, sep = \",\")";
            String writeDETable = "write.table(top.table, file = \"" + dePath + "\", row.names = F, sep = \"\t\", quote = F)";

            code.R_require("edgeR");
            code.addRCode(readFileCounts);
            code.addRCode("counts2 <- counts[,-1]");
            code.addRCode("rownames(counts2) <- counts[,1]");
            code.addRCode("counts <- data.matrix(counts2)");

            code.addRCode("d0 <- DGEList(counts)");
            code.addRCode("d0 <- calcNormFactors(d0)");

            code.addRCode("snames <- colnames(counts)");
            code.addRCode(readFileGroups);

            code.addRCode("group <- interaction(groups['Groups'])");

            code.addRCode("png" + "(" + mdsPath + ")");
            code.addRCode("plotMDS(d, col = as.numeric(group))");
            code.addRCode("dev.off()");

            code.addRCode("mm <- model.matrix(~0 + group)");
            code.addRCode("png" + "(" + voomPath + ")");
            code.addRCode("y <- voom(d, mm, plot = T)");
            code.addRCode("dev.off()");

            code.addRCode("fit <- lmFit(y, mm)");
            code.addRCode("contr <- makeContrasts(groupgroup1 - groupgroup2, levels = colnames(coef(fit)))");

            code.addRCode("tmp <- contrasts.fit(fit, contr)");
            code.addRCode("tmp <- eBayes(tmp)");

            code.addRCode("top.table <- topTable(tmp, sort.by = \"P\", n = Inf)");

            code.addRCode("length(which(top.table$adj.P.Val < 0.05))");
            code.addRCode("top.table$Gene <- rownames(top.table)");
            code.addRCode("top.table <- top.table[,c(\"Gene\", names(top.table)[1:6])]");
            code.addRCode(writeDETable);
            //        caller.runAndReturnResult("result");
            System.out.println("Done!");

        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }

    public static void main(String[] args) {
        calc();
    }
}
