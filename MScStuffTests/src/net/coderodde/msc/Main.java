package net.coderodde.msc;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class Main {

    private static final int K_MER_SIZE = 31;
    
    public static void main(final String... args) {
        readFilesDemo(args);
//        smallDemo();
    }
    
    private static void readFilesDemo(final String... args) {
        final List<String> genomeList = new ArrayList<>();
        final long startTimeTotal = System.currentTimeMillis();
        
        System.out.println("Args: " + args[0]);
        
        final File[] files = new File(args[0]).listFiles();
        
        for (final File file : files) {
            System.out.println("---");
            System.out.println("Reading \"" + file.getAbsolutePath() + "\"...");
            
            long startTime = System.currentTimeMillis();
            final String genome = GenomeReader.readFile(file);
            genomeList.add(genome);
            long endTime = System.currentTimeMillis();
            
            final long duration = endTime - startTime;
            
            System.out.println("Read in " + duration + " millisecons.");
            System.out.println("Length: " + genome.length());
        }
        
        final long endTimeTotal = System.currentTimeMillis();
        System.out.println("[RESULT] Done in " + (endTimeTotal - startTimeTotal)
                                               + " milliseconds.");
        
        final long startTime = System.currentTimeMillis();
        final NodeCentricDeBruijnGraph graph = 
                new NodeCentricDeBruijnGraph(genomeList, K_MER_SIZE);
        final long endTime = System.currentTimeMillis();
        
        System.out.println("Graph built in " + (endTime - startTime) 
                                             + " milliseconds.");
        
        System.out.println("Nodes: " + graph.getAllNodes().size());
        System.out.println("Arcs:  " + graph.getNumberOfArcs());
    }
    
    private static void smallDemo() {
        List<String> readList = Arrays.asList("ACCGCTA", 
                                              "TTACGG", 
                                              "GTTA", 
                                              "AATAG");
        System.out.println("=== Node-centric de Bruijn graph ===");
        
        NodeCentricDeBruijnGraph graph = 
                new NodeCentricDeBruijnGraph(readList, 3);
        
        System.out.println(graphToString(graph));
        
        System.out.println();
        System.out.println("=== Edge-centric de Bruijn graph ===");
        
        EdgeCentricDeBruijnGraph graph2 = 
                new EdgeCentricDeBruijnGraph(readList, 3);
        
        System.out.println(graphToString(graph2));
    }
    
    private static String graphToString(AbstractDeBruijnGraph graph) {
        List<Kmer> nodeList = new ArrayList<>(graph.getAllNodes());
        String tmp = Integer.toString(nodeList.size());
        int fieldLength = tmp.length();
        int lineNumber = 1;
        Collections.<Kmer>sort(nodeList);
        StringBuilder sb = new StringBuilder();
        String lineNumberFormatToken = "%" + fieldLength + "d: ";
        
        for (Kmer node : nodeList) {
            sb.append(String.format(lineNumberFormatToken, lineNumber++));
            sb.append(node);
            sb.append(", children: [");
            
            if (!graph.getChildrenOf(node).isEmpty()) {
                for (Kmer child : graph.getChildrenOf(node)) {
                    sb.append(child);
                    sb.append(" ");
                }

                sb.deleteCharAt(sb.length() - 1);
            }
                
            sb.append("], parents[");
            
            if (!graph.getParentsOf(node).isEmpty()) {
                for (Kmer parent : graph.getParentsOf(node)) {
                    sb.append(parent);
                    sb.append(" ");
                }

                sb.deleteCharAt(sb.length() - 1);
            }
               
            sb.append("]\n");
        }
        
        sb.deleteCharAt(sb.length() - 1);
        return sb.toString();
    }
}
