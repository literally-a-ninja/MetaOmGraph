package edu.iastate.metnet.metaomgraph;

import org.objenesis.ObjenesisException;
import org.renjin.script.RenjinScriptEngineFactory;

import javax.script.ScriptEngine;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * This class uses Renjin to perform Limma analysis
 */
public class ComputeLimma {

    private String selectedList;
    private String grpID;
    private MetaOmProject myProject;
    private Map<String, Collection<Integer>> splitIndex;
    private boolean[] excluded;
    private int testMethod;
    private int filterLowExpGene;
    private RenjinScriptEngineFactory factory;
    private ScriptEngine engine;

    private List<String> featureNames;
    Collection<Collection<Integer>> grpInds;
    Collection<String> grpName;

    public List<String> getFeatureNames() {
        return this.featureNames;
    }

    public ComputeLimma(String selectedList, String grpID, MetaOmProject myProject, boolean tflag, int filterLowExpGene) {
        this.selectedList = selectedList;
        this.grpID = grpID;
        this.myProject = myProject;
        excluded = MetaOmAnalyzer.getExclude();
        this.testMethod = 0;
        this.filterLowExpGene = filterLowExpGene;

        factory = new RenjinScriptEngineFactory();
        engine = factory.getScriptEngine();
    }

    public void doCalc() throws IOException {

        int[] selected = myProject.getGeneListRowNumbers(this.selectedList);

        for (int i : selected) {
            double[] thisDataRaw = null; // untransformed data
            double[] thisData = null;

            thisDataRaw = myProject.getAllData(i, true);
            featureNames.add(myProject.getDefaultRowNames(i));
        }


            // Start of Limma analysis
        try {
            // Preprocessing
            engine.eval("d0 <- DGEList(counts)");
            engine.eval("d0 <- calcNormFactors(d0)");

            // Filter low expression gene (Optional, if filterLowExpGene > 0)
            if (filterLowExpGene > 0) {
                String cutoff = "cutoff <- " + filterLowExpGene;
                String drop = "drop <- which(apply(cpm(d0), " + filterLowExpGene + ", max) < cutoff)";
                engine.eval(cutoff);
                engine.eval(drop);
                engine.eval("d <- d0[-drop,");
            }

            // Voom
            engine.eval("mm <- model.matrix(~0 + group)");

        } catch (Exception ignored) {

        }
    }
}
