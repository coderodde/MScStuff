package net.coderodde.graph.scc.support;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import net.coderodde.graph.DirectedGraph;
import net.coderodde.graph.scc.SCCFinder;

/**
 * This class implements the recursive
 * <a href="https://en.wikipedia.org/wiki/Kosaraju%27s_algorithm">Kosarajus's
 * algorithm</a>
 * for finding strongly connected components in an input directed graph.
 *
 * @author Rodion "rodde" Efremov
 * @version 1.6 (May 3, 2016)
 */
public final class KosarajuSCCFinder implements SCCFinder {

    private DirectedGraph digraph;
    private List<Integer> nodeList;
    private Set<Integer> visitedSet;
    private Map<Integer, Integer> assignmentMap;

    public KosarajuSCCFinder() {}

    private KosarajuSCCFinder(final DirectedGraph digraph) {
        Objects.requireNonNull(digraph, "The input directed graph is null.");
        this.digraph = digraph;
        this.nodeList = new ArrayList<>(digraph.size());
        this.visitedSet = new HashSet<>(digraph.size());
        this.assignmentMap = new HashMap<>(digraph.size());
    }

    @Override
    public List<List<Integer>>
            findStronglyConnectedCmponents(DirectedGraph digraph) {
        return new KosarajuSCCFinder(digraph).compute();
    }

    private List<List<Integer>> compute() {
        for (final Integer node : digraph.getAllNodes()) {
            visit(node);
        }
        
        Collections.<Integer>reverse(nodeList);

        for (final Integer node : nodeList) {
            assign(node, node);
        }

        final Map<Integer, List<Integer>> map = new HashMap<>();

        for (final Map.Entry<Integer, Integer> entry
                : assignmentMap.entrySet()) {
            final Integer component = entry.getValue();

            if (!map.containsKey(component)) {
                map.put(component, new ArrayList<>());
            }

            map.get(component).add(entry.getKey());
        }

        return new ArrayList<>(map.values());
    }

    private void visit(final Integer node) {
        if (visitedSet.contains(node)) {
            return;
        }
        
        visitedSet.add(node);
        
        for (final Integer child : digraph.getChildrenOf(node)) {
            visit(child);
        }
        
        nodeList.add(node);
    }

    private void assign(final Integer node, final Integer root) {
        if (!assignmentMap.containsKey(node)) {
            assignmentMap.put(node, root);

            for (final Integer parent : digraph.getParentsOf(node)) {
                assign(parent, root);
            }
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

        final SCCFinder finder = new KosarajuSCCFinder();
        System.out.println(finder.findStronglyConnectedCmponents(digraph));
    }
}
