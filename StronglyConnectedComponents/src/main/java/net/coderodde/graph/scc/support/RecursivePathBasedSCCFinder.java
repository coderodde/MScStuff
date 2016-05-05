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
 * This class implements a <a href="">path-based strong component algorithm</a>.
 * 
 * @author Rodion "rodde" Efremov
 * @version 1.6
 */
public class RecursivePathBasedSCCFinder implements SCCFinder {

    private DirectedGraph digraph;
    private int counter;
    private Set<Integer> assignedNodeSet;
    private Map<Integer, Integer> preorderNumberMap;
    private Deque<Integer> stackP;
    private Deque<Integer> stackS;
    private List<List<Integer>> solution;
    
    public RecursivePathBasedSCCFinder() {}
    
    private RecursivePathBasedSCCFinder(final DirectedGraph digraph) {
        Objects.requireNonNull(digraph, "The input directed graph is null.");
        this.digraph = digraph;
        this.assignedNodeSet = new HashSet<>();
        this.preorderNumberMap = new HashMap<>();
        this.stackP = new ArrayDeque<>();
        this.stackS = new ArrayDeque<>();
        this.solution = new ArrayList<>();
    }
    
    @Override
    public List<List<Integer>> 
    findStronglyConnectedCmponents(DirectedGraph digraph) {
        return new RecursivePathBasedSCCFinder(digraph).compute();
    }
    
    private List<List<Integer>> compute() {
        for (final Integer node : digraph.getAllNodes()) {
            if (!preorderNumberMap.containsKey(node)) {
                visit(node);
            }
        }
        
        return this.solution;
    }
    
    private void visit(final Integer node) {
        preorderNumberMap.put(node, counter++);
        stackP.addLast(node);
        stackS.addLast(node);
        
        for (final Integer child : digraph.getChildrenOf(node)) {
            if (!preorderNumberMap.containsKey(child)) {
                visit(child);
            } else if (!assignedNodeSet.contains(child)) {
                while (preorderNumberMap.get(stackP.getLast()) > 
                       preorderNumberMap.get(child)) {
                    stackP.removeLast();
                }
            }
        }
        
        if (node.equals(stackP.getLast())) {
            stackP.removeLast();
            Integer topOfStackS;
            final List<Integer> component = new ArrayList<>();
            
            do {
                topOfStackS = stackS.removeLast();
                component.add(topOfStackS);
                assignedNodeSet.add(topOfStackS);
            } while (!topOfStackS.equals(node));
            
            this.solution.add(component);
        }
    }
    
    private static void main(String[] args) {
        System.out.println("Path-based");
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

        digraph.addEdge(b, c);
        digraph.addEdge(b, e);
        digraph.addEdge(b, f);

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

        final SCCFinder finder = new RecursivePathBasedSCCFinder();
        System.out.println(finder.findStronglyConnectedCmponents(digraph));
    }
}
