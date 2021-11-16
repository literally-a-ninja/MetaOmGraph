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

    }

    public void doCalc() throws IOException {

        int[] selected = myProject.getGeneListRowNumbers(this.selectedList);

        for (int i : selected) {
            double[] thisDataRaw = null; // untransformed data
            double[] thisData = null;

            thisDataRaw = myProject.getAllData(i, true);
            featureNames.add(myProject.getDefaultRowNames(i));
        }
        } catch (Exception ignored) {

        }
    }
}
