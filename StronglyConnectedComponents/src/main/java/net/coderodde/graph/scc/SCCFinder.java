package net.coderodde.graph.scc;

import java.util.List;
import net.coderodde.graph.DirectedGraph;

/**
 * This interface defines the API for algorithms finding strongly connected 
 * components.
 * 
 * @author Rodion "rodde" Efremov
 * @version 1.6 (May 3, 2016)
 */
public interface SCCFinder {
   
    /**
     * Returns a list of strongly connected components in the input graph
     * {@code digraph}.
     * 
     * @param digraph
     * @return 
     */
    public List<List<Integer>> 
        findStronglyConnectedCmponents(final DirectedGraph digraph);
}
