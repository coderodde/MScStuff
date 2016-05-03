import java.util.Random;
import net.coderodde.graph.DirectedGraph;
import net.coderodde.graph.scc.SCCFinder;
import net.coderodde.graph.scc.Util;
import net.coderodde.graph.scc.support.KosarajuSCCFinder;
import net.coderodde.graph.scc.support.RecursiveKosarajuSCCFinder;
import net.coderodde.graph.scc.support.RecursiveTarjanSCCFinder;
import net.coderodde.graph.scc.support.TarjanSCCFinder;

public class RecursiveExperiment {

    private static final int NUMBER_OF_NODES = 100000;
    private static final int NUMBER_OF_ARCS  = 150000;
    
    private static final int WARMUP_NODES = 1000;
    private static final int WARMUP_ARCS = 1500;
    private static final int WARMUP_ITERATIONS = 300;
    
    public static void main(String[] args) {
        final long seed = System.nanoTime();
        final Random random = new Random(seed);
        final DirectedGraph digraph = Util.createRandomDigraph(random,
                                                               NUMBER_OF_NODES,
                                                               NUMBER_OF_ARCS);
        final DirectedGraph warmupDigraph = 
                Util.createRandomDigraph(random, WARMUP_NODES, WARMUP_ARCS);
        
        Util.warmup(warmupDigraph, WARMUP_ITERATIONS);
        
        System.out.println("Seed = " + seed);
        
        final SCCFinder kosaraju          = new KosarajuSCCFinder();
        final SCCFinder recursiveKosaraju = new RecursiveKosarajuSCCFinder();
//        final SCCFinder tarjan            = new TarjanSCCFinder();
        final SCCFinder recursiveTarjan   = new RecursiveTarjanSCCFinder();
        
        long startTime = System.nanoTime();
//        tarjan.findStronglyConnectedCmponents(digraph);
        long endTime = System.nanoTime();
        
        System.out.printf("Tarjan's algorithm in %.2f milliseconds.\n",
                          (endTime - startTime) / 1e6);
        
        startTime = System.nanoTime();
        recursiveTarjan.findStronglyConnectedCmponents(digraph);
        endTime = System.nanoTime();
        
        System.out.printf(
                "Recursive Tarjan's algorithm in %.2f milliseconds.\n",
                (endTime - startTime) / 1e6);
        
        startTime = System.nanoTime();
        kosaraju.findStronglyConnectedCmponents(digraph);
        endTime = System.nanoTime();
        
        System.out.printf("Kosaraju's algorithm in %.2f milliseconds.\n",
                          (endTime - startTime) / 1e6);
        
        startTime = System.nanoTime();
        recursiveKosaraju.findStronglyConnectedCmponents(digraph);
        endTime = System.nanoTime();
        
        System.out.printf(
                "Recursive Kosaraju's algorithm in %.2f milliseconds.\n",
                (endTime - startTime) / 1e6);
    }
}
