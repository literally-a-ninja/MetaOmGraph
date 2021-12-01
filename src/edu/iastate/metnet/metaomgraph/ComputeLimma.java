package edu.iastate.metnet.metaomgraph;

import com.github.rcaller.datatypes.DataFrame;
import com.github.rcaller.rstuff.RCaller;
import com.github.rcaller.rstuff.RCode;
import javax.script.ScriptEngine;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.List;


public class ComputeLimma {

    private HashMap<String, String> groups;
    private List<String> names;
    private MetaOmProject myProject;
    private double lowly;

    public ComputeLimma(HashMap<String, String> groups, List<String> names, MetaOmProject myProject, double lowly) {
        this.groups = groups;
        this.names = names;
        this.lowly = lowly;
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

            code.addRCode("cutoff <- 2");
            code.addRCode("drop <- which(apply(cpm(d0), 1, max) < cutoff)");
            code.addRCode("print(d0)");
            code.addRCode("d <- d0[-drop,]");

            code.addRCode("snames <- colnames(counts)");
            code.addRCode(readFileGroups);

            code.addRCode("group <- interaction(groups['Groups'])");

//            code.addRCode("png" + "(" + mdsPath + ")");
//            code.addRCode("plotMDS(d, col = as.numeric(group))");
//            code.addRCode("dev.off()");

//            code.addRCode("mm <- model.matrix(~0 + group)");
//            code.addRCode("png" + "(" + voomPath + ")");
//            code.addRCode("y <- voom(d, mm, plot = T)");
//            code.addRCode("dev.off()");

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

    public static void main(String[] args) throws IOException {
//        doCalc();
    }
}
