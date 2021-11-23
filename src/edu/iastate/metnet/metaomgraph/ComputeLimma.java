package edu.iastate.metnet.metaomgraph;

import org.openxmlformats.schemas.spreadsheetml.x2006.main.STIconSetType;
import org.renjin.script.RenjinScriptEngineFactory;
import org.renjin.script.*;
//import org.renjin.bioconductor.limma.limma;
import org.renjin.bioconductor.edgeR.edgeR;

import javax.script.ScriptEngine;
import java.io.IOException;
import java.net.URISyntaxException;


public class ComputeLimma {

    public static void doCalc() throws IOException {
        RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
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

            engine.eval("library(grid)");
            engine.eval("library(lattice)");
            engine.eval("library(locfit)");
            engine.eval("library(limma)");
            engine.eval("library(edgeR)");
            engine.eval(readFileCounts);
            engine.eval("counts2 <- counts[,-1]");
            engine.eval("rownames(counts2) <- counts[,1]");
            engine.eval("counts <- data.matrix(counts2)");

            engine.eval("d0 <- DGEList(counts)");
            engine.eval("d0 <- calcNormFactors(d0)");
            engine.eval("dim(d0)");

            engine.eval("cutoff <- 2");
            engine.eval("drop <- which(apply(cpm(d0), 1, max) < cutoff)");
            engine.eval("print(d0)");
            engine.eval("d <- d0[-drop,]");
            engine.eval("dim(d)");

            engine.eval("snames <- colnames(counts)");
            engine.eval(readFileGroups);

            engine.eval("group <- interaction(groups['Groups'])");

//            engine.eval("png" + "(" + mdsPath + ")");
//            engine.eval("plotMDS(d, col = as.numeric(group))");
//            engine.eval("dev.off()");

//            engine.eval("mm <- model.matrix(~0 + group)");
//            engine.eval("png" + "(" + voomPath + ")");
//            engine.eval("y <- voom(d, mm, plot = T)");
//            engine.eval("dev.off()");

            engine.eval("fit <- lmFit(y, mm)");
            engine.eval("contr <- makeContrasts(groupgroup1 - groupgroup2, levels = colnames(coef(fit)))");

            engine.eval("tmp <- contrasts.fit(fit, contr)");
            engine.eval("tmp <- eBayes(tmp)");

            engine.eval("top.table <- topTable(tmp, sort.by = \"P\", n = Inf)");

            engine.eval("length(which(top.table$adj.P.Val < 0.05))");
            engine.eval("top.table$Gene <- rownames(top.table)");
            engine.eval("top.table <- top.table[,c(\"Gene\", names(top.table)[1:6])]");
            engine.eval(writeDETable);
            System.out.println("Done!");

        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }

    //    public static void doCalcRCaller() throws IOException, URISyntaxException {
//        String fileContent = RUtils.getMeanScriptContent();
//        RCode code = RCode.create();
//        code.addRCode(fileContent);
//        code.addIntArray("input", values);
//        code.addRCode("result <- customMean(input)");
//        RCaller caller = RCaller.create(code, RCallerOptions.create());
//        caller.runAndReturnResult("result");
//    }
    public static void main(String[] args) throws IOException {
        doCalc();
    }
}
