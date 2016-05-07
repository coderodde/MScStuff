
import java.util.Random;
import net.coderodde.graph.DirectedGraph;
import net.coderodde.graph.scc.SCCFinder;
import net.coderodde.graph.scc.Util;
import net.coderodde.graph.scc.support.KosarajuSCCFinder;
import net.coderodde.graph.scc.support.PathBasedSCCFinder;
import net.coderodde.graph.scc.support.TarjanSCCFinder;

public class ExperimentDataGenerator {
    
    private static final int WARMUP_ITERATIONS = 100;
    private static final int WARMUP_NODES = 10_000;
    private static final int WARMUP_ARCS = 12_000;
    private static final String STATUS = "[STATUS] ";
    
    private static final SCCFinder KOSARAJU = new KosarajuSCCFinder();
    private static final SCCFinder PATHBASED = new PathBasedSCCFinder();
    private static final SCCFinder TARJAN = new TarjanSCCFinder();
    private static final int ITERATIONS = 10;
    
    /**
     * The number of arcs will approximately equal to the number of nodes times
     * this constant.
     */
    private static final double ARC_FACTOR = 1.4;
    
    public static void main(String[] args) {
        final long seed = System.nanoTime();
        final Random random = new Random(seed);
        
        System.out.println("Seed = " + seed);
        System.out.println(STATUS + "Warming up...");
        
        final DirectedGraph warmupDigraph = 
                Util.createRandomDigraph(random, WARMUP_NODES, WARMUP_ARCS);
        
        Util.warmup(warmupDigraph, WARMUP_ITERATIONS);
        System.out.println(STATUS + "Warming up done!");
        
        for (int size = 10_000; size <= 200_000; size += 10_000) {
            System.out.println(generateDataRow(size, random));
        }
    }
    
    private static String generateDataRow(final int size, final Random random) {
        final DirectedGraph graph = 
                Util.createRandomDigraph(random, 
                                         size, 
                                         (int)(ARC_FACTOR * size));
        
        final StringBuilder sb = new StringBuilder();
        sb.append(String.format("%3dk ", size / 1000));
        
        long startTime = System.nanoTime();
        
        for (int i = 0; i < ITERATIONS; ++i) {
            KOSARAJU.findStronglyConnectedCmponents(graph);
        }
            
        long endTime = System.nanoTime();
        long duration = (endTime - startTime) / ITERATIONS;
        
        sb.append(String.format("%4.0f ", duration / 1e6));
        
        startTime = System.nanoTime();
        
        for (int i = 0; i < ITERATIONS; ++i) {
            PATHBASED.findStronglyConnectedCmponents(graph);
        }
            
        endTime = System.nanoTime();
        duration = (endTime - startTime) / ITERATIONS;
        
        sb.append(String.format("%4.0f ", duration / 1e6));
        
        startTime = System.nanoTime();
        
        for (int i = 0; i < ITERATIONS; ++i) {
            TARJAN.findStronglyConnectedCmponents(graph);
        }
            
        endTime = System.nanoTime();
        duration = (endTime - startTime) / ITERATIONS;
        
        sb.append(String.format("%4.0f", duration / 1e6));
        return sb.toString();
    }
}
