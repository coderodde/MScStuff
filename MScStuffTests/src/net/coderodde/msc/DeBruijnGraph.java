package net.coderodde.msc;

import java.util.Set;

public interface DeBruijnGraph {
    public Set<String> getAllNodes();
    public Set<String> getChildrenOf(String kmer);
    public Set<String> getParentsOf(String kmer);
}
