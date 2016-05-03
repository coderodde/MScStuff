package net.coderodde.graph.scc.support;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import net.coderodde.graph.scc.Util;
import net.coderodde.graph.DirectedGraph;

import org.junit.Test;
import static org.junit.Assert.*;

public class SCCFinderTest {

    private static final int ITERATIONS = 100;
    private static final int MAXIMUM_NUMBER_OF_NODES = 500;
    private static final float MINIMUM_ARC_LOAD_FACTOR = 0.001f;
    
    @Test
    public void testCorrectness() {
        final long seed = System.nanoTime();
        final Random random = new Random(seed);
        System.out.println("Seed = " + seed);
        
        for (int iteration = 0; iteration < ITERATIONS; ++iteration) {
            testCorrectnessOnce(random);
        }
    }
    
    private void testCorrectnessOnce(final Random random) {
        final int numberOfNodes = random.nextInt(MAXIMUM_NUMBER_OF_NODES + 1);
        final float arcLoadFactor = Math.max(MINIMUM_ARC_LOAD_FACTOR,
                                             random.nextFloat());
        final int numberOfArcs = (int)(arcLoadFactor * 
                                       Math.pow(numberOfNodes, 2.0));
        final DirectedGraph digraph = Util.createRandomDigraph(random,
                                                          numberOfNodes,
                                                          numberOfArcs);
        
        final List<List<Integer>> scc1 = 
                new RecursiveKosarajuSCCFinder()
                        .findStronglyConnectedCmponents(digraph);
        
        final List<List<Integer>> scc2 = 
                new KosarajuSCCFinder().findStronglyConnectedCmponents(digraph);
        
        final List<List<Integer>> scc3 = 
                new RecursiveTarjanSCCFinder()
                .findStronglyConnectedCmponents(digraph);
        
        for (final List<Integer> component : scc1) {
            Collections.sort(component);
        }
        
        for (final List<Integer> component : scc2) {
            Collections.sort(component);
        }
        
        for (final List<Integer> component : scc3) {
            Collections.sort(component);
        }
        
        final Set<List<Integer>> scc1set = new HashSet<>();
        final Set<List<Integer>> scc2set = new HashSet<>();
        final Set<List<Integer>> scc3set = new HashSet<>();
        
        for (final List<Integer> scc : scc1) {
            scc1set.add(scc);
        }
        
        for (final List<Integer> scc : scc2) {
            scc2set.add(scc);
        }
        
        for (final List<Integer> scc : scc3) {
            scc3set.add(scc);
        }
        
        assertEquals(scc1set, scc2set);
        assertEquals(scc2set, scc3set);
    }
    
    
}
