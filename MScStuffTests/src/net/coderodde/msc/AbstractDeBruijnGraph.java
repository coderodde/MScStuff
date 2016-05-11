package net.coderodde.msc;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public abstract class AbstractDeBruijnGraph {
    
    protected final int k;
    protected final Map<Kmer, Set<Kmer>> childrenMap = new HashMap<>();
    protected final Map<Kmer, Set<Kmer>> parentsMap  = new HashMap<>();
    protected int arcs;
    
    AbstractDeBruijnGraph(int k) {
        checkKMerSize(k);
        this.k = k;
    }
    
    public abstract Set<Kmer> getAllNodes();
    public abstract Set<Kmer> getChildrenOf(Kmer kmer);
    public abstract Set<Kmer> getParentsOf(Kmer kmer);
    
    public int getNumberOfArcs() {
        return arcs;
    }
    
    protected void checkKMerSize(int k) {
        if (k < 2) {
            throw new IllegalArgumentException(
                    "The k is too small: " + k + ". Must be at least 2.");
        }
    }
    
    protected void checkStringIsKmer(Kmer kmer) {
        if (kmer.length() != k) {
            throw new IllegalArgumentException(
                    "The input string \"" + kmer + "\" is not a " +
                    k + "-mer.");
        }
    }
}
