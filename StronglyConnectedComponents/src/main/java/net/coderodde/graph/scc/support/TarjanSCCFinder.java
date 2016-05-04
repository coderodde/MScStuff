package net.coderodde.graph.scc.support;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
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
                strongConnect(node);
            }
        }
        
        return this.solution;
    }
    
    private void strongConnect(final Integer node) {
        final Deque<Integer> nodeStack = new ArrayDeque<>();
        final Deque<Integer> childStack = new ArrayDeque<>();
        final Deque<Iterator<Integer>> nodeChildIteratorStack = 
                new ArrayDeque<>();
        
        nodeStack.addLast(node);
        nodeChildIteratorStack.add(digraph.getChildrenOf(node).iterator());
        
        outer:
        while (!nodeStack.isEmpty()) {
            indexMap.put(node, index);
            lowLinkMap.put(node, index);
            ++index;
            stack.push(node);
            onStackSet.add(node);
            
            final Iterator<Integer> currentNodeChildIterator = 
                    nodeChildIteratorStack.getLast();
            final Integer current = nodeStack.getLast();
            
            while (currentNodeChildIterator.hasNext()) {
                final Integer child = currentNodeChildIterator.next();
                childStack.addLast(child);
                
                if (!indexMap.containsKey(child)) {
                    nodeStack.addLast(child);
                    nodeChildIteratorStack.addLast(
                            digraph.getChildrenOf(child).iterator());
                    
                    continue outer;
                } else if (onStackSet.contains(child)) {
                    lowLinkMap.put(current, 
                                   Math.min(lowLinkMap.get(current),
                                              indexMap.get(child)));
                }
                
                childStack.removeLast();
            }
            
            if (lowLinkMap.get(current).equals(indexMap.get(current))) {
                final List<Integer> stronglyConnectedComponent = 
                        new ArrayList<>();

                Integer w;
                
                do {
                    w = stack.pop();
                    onStackSet.remove(w);
                    stronglyConnectedComponent.add(w);
                } while (!w.equals(current));
                
                this.solution.add(stronglyConnectedComponent);
            }
            
            while (!nodeChildIteratorStack.isEmpty()
                    && !nodeChildIteratorStack.getLast().hasNext()) {
                final Integer tmp = nodeStack.removeLast();
                nodeChildIteratorStack.removeLast();
                final Integer child = childStack.removeLast();
                lowLinkMap.put(tmp, Math.min(lowLinkMap.get(tmp),
                                             lowLinkMap.get(child)));
            }
        }
    }
    
    public static void main(String[] args) {
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
        
        final SCCFinder finder = new KosarajuSCCFinder();
        System.out.println(finder.findStronglyConnectedCmponents(digraph));
    }
}
