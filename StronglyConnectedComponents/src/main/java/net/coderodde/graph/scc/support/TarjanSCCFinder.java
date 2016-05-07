package net.coderodde.graph.scc.support;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import net.coderodde.graph.DirectedGraph;
import net.coderodde.graph.scc.SCCFinder;

/**
 * This class implements 
 * <a href="https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm">Tarjan's strongly connected components algorithm</a> using 
 * recursive depth-first search.
 * 
 * @author Rodion "rodde" Efremov
 * @version 1.6 (May 3, 2016)
 */
public final class TarjanSCCFinder implements SCCFinder {

    private DirectedGraph digraph;
    private int index;
    private Deque<Integer> stack;
    private Set<Integer> onStackSet;
    private Map<Integer, Integer> indexMap;
    private Map<Integer, Integer> lowLinkMap;
    private List<List<Integer>> solution;
    
    public TarjanSCCFinder() {}
    
    private TarjanSCCFinder(final DirectedGraph digraph) {
        Objects.requireNonNull(digraph, "The input digraph is null.");
        this.digraph = digraph;
        this.stack = new ArrayDeque<>();
        this.indexMap = new HashMap<>();
        this.lowLinkMap = new HashMap<>();
        this.onStackSet = new HashSet<>();
        this.solution = new ArrayList<>();
    }
    
    @Override
    public List<List<Integer>> 
        findStronglyConnectedCmponents(final DirectedGraph digraph) {
        return new TarjanSCCFinder(digraph).compute();
    }
    
    private List<List<Integer>> compute() {
        Objects.requireNonNull(digraph, "The input directed graph is null.");
        
        for (final Integer node : digraph.getAllNodes()) {
            if (!indexMap.containsKey(node)) {
//                System.out.println("TOP: !indexMap.containsKey(" + node + ")");
                strongConnect(node);
            }
        }
        
        return this.solution;
    }
    
    private void strongConnect(final Integer node) {
//        System.out.println("strongConnect(" + node + ")");
        indexMap.put(node, index);
        lowLinkMap.put(node, index);
        ++index;
        stack.push(node);
        onStackSet.add(node);
        
        for (final Integer child : digraph.getChildrenOf(node)) {
//            System.out.println("for child " + child);
            
            if (!indexMap.containsKey(child)) {
//                System.out.println("!indexMap.containsKey(" + child + ")");
                strongConnect(child);
//                System.out.println("after recursive call for child " + child);
                lowLinkMap.put(node, Math.min(lowLinkMap.get(node), 
                                              lowLinkMap.get(child)));
            } else if (onStackSet.contains(child)) {
//                System.out.println("onStackSet.contains(" + child + ")");
                lowLinkMap.put(node, Math.min(lowLinkMap.get(node), 
                                              indexMap.get(child)));
            }
        }
        
        if (lowLinkMap.get(node).equals(indexMap.get(node))) {
//            System.out.println("lowLinkMap.get(" + node + ").equals(indexMap.get(" + node + "))");
            final List<Integer> newStronglyConnectedComponent = 
                    new ArrayList<>();
            
            Integer top;
            
            do {
                top = stack.pop();
                onStackSet.remove(top);
                newStronglyConnectedComponent.add(top);
            } while (!top.equals(node));
            
            this.solution.add(newStronglyConnectedComponent);
        }
    }
    
    private static void main(String[] args) {
        final int a = 0;
        final int b = 1;
        final int c = 2;
        final int d = 3;
        final int e = 4;
        final int f = 5;
        final int g = 6;
        final int h = 7;
        
        final DirectedGraph digraph = new DirectedGraph();
        
        digraph.addNode(a);
        digraph.addNode(b);
        digraph.addNode(c);
        digraph.addNode(d);
        digraph.addNode(e);
        digraph.addNode(f);
        digraph.addNode(g);
        digraph.addNode(h);
        
        digraph.addEdge(a, b);
        
        digraph.addEdge(b, e);
        digraph.addEdge(b, f);
        digraph.addEdge(b, c);
    
        digraph.addEdge(c, d);
        digraph.addEdge(c, g);
        
        digraph.addEdge(d, c);
        digraph.addEdge(d, h);
        
        digraph.addEdge(e, a);
        digraph.addEdge(e, f);
        
        digraph.addEdge(f, g);
        
        digraph.addEdge(g, f);
        digraph.addEdge(g, h);
        
        digraph.addEdge(h, h);  
        
        final SCCFinder finder = new TarjanSCCFinder();
//        System.out.println(finder.findStronglyConnectedCmponents(digraph));
        
        digraph.clear();
        
        digraph.addNode(0);
        digraph.addNode(1);
        digraph.addNode(2);
        digraph.addNode(3);
        
        digraph.addEdge(0, 2);
        digraph.addEdge(2, 0);
        digraph.addEdge(0, 1);
        digraph.addEdge(2, 3);
        
        System.out.println(finder.findStronglyConnectedCmponents(digraph));
    }
}
