package edu.iastate.metnet.metaomgraph;

import org.renjin.script.RenjinScriptEngine;

import java.util.Collection;
import java.util.Map;

public class ComputeLimma {

    private String selectedList;
    private String grpID;
    private MetaOmProject myProject;
    private Map<String, Collection<Integer>> splitIndex;
    private boolean[] excluded;
    private int testMethod;
    private RenjinScriptEngine engine;

    public ComputeLimma(String selectedList, String grpID, MetaOmProject myProject, boolean tflag) {
        this.selectedList = selectedList;
        this.grpID = grpID;
        this.myProject = myProject;
        excluded = MetaOmAnalyzer.getExclude();
        this.testMethod = 0;

        engine = new RenjinScriptEngine();
    }
}
