package net.coderodde.graph.scc;


import java.util.Random;
import net.coderodde.graph.DirectedGraph;
import net.coderodde.graph.scc.support.IterativeKosarajuSCCFinder;
import net.coderodde.graph.scc.support.KosarajuSCCFinder;
import net.coderodde.graph.scc.support.TarjanSCCFinder;
import net.coderodde.graph.scc.support.IterativeTarjanSCCFinder;

/**
 * This class contains various utility methods.
 *
 * @author Rodion "rodde" Efremov
 * @version 1.6 (May 3, 2016)
 */
public class Util {

    public static DirectedGraph createRandomDigraph(final Random random,
                                                    final int numberOfNodes,
                                                    final int numberOfArcs) {
        final DirectedGraph digraph = new DirectedGraph();

        for (int i = 0; i < numberOfNodes; ++i) {
            digraph.addNode(i);
        }

        for (int i = 0; i < numberOfArcs; ++i) {
            final int tail = random.nextInt(numberOfNodes);
            final int head = random.nextInt(numberOfNodes);
            digraph.addEdge(tail, head);
        }

        return digraph;
    }
    
    public static void warmup(final DirectedGraph digraph, 
                              final int iterations) {
        final SCCFinder scc1 = new IterativeKosarajuSCCFinder();
        final SCCFinder scc2 = new KosarajuSCCFinder();
        final SCCFinder scc3 = new IterativeTarjanSCCFinder();
        final SCCFinder scc4 = new TarjanSCCFinder();
        
        for (int i = 0; i < iterations; ++i) {
            scc1.findStronglyConnectedCmponents(digraph);
            scc2.findStronglyConnectedCmponents(digraph);
            scc3.findStronglyConnectedCmponents(digraph);
            scc4.findStronglyConnectedCmponents(digraph);
        }
    }
}
