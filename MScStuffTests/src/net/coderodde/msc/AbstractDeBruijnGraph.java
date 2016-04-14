package net.coderodde.msc;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public abstract class AbstractDeBruijnGraph {
    
    protected final int k;
    protected final Map<String, Set<String>> childrenMap = new HashMap<>();
    protected final Map<String, Set<String>> parentsMap  = new HashMap<>();
    
    AbstractDeBruijnGraph(int k) {
        checkKMerSize(k);
        this.k = k;
    }
    
    public abstract Set<String> getAllNodes();
    public abstract Set<String> getChildrenOf(String kmer);
    public abstract Set<String> getParentsOf(String kmer);
    
    protected void checkKMerSize(int k) {
        if (k < 2) {
            throw new IllegalArgumentException(
                    "The k is too small: " + k + ". Must be at least 2.");
        }
    }
    
    protected void checkStringIsKmer(String string) {
        if (string.length() != k) {
            throw new IllegalArgumentException(
                    "The input string \"" + string + "\" is not a " +
                    k + "-mer.");
        }
    }
}
