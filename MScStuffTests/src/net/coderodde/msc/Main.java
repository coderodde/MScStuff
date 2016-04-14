package net.coderodde.msc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class Main {

    public static void main(final String... args) {
        List<String> readList = Arrays.asList("ACCGCTA", 
                                              "TTACGG", 
                                              "GTTA", 
                                              "AATAG");
        
        NodeCentricDeBruijnGraph graph = 
                new NodeCentricDeBruijnGraph(readList, 3);
        System.out.println(graphToString(graph));
    }
    
    private static String graphToString(AbstractDeBruijnGraph graph) {
        List<String> nodeList = new ArrayList<>(graph.getAllNodes());
        String tmp = Integer.toString(nodeList.size());
        int fieldLength = tmp.length();
        int lineNumber = 1;
        Collections.<String>sort(nodeList);
        StringBuilder sb = new StringBuilder();
        String lineNumberFormatToken = "%" + fieldLength + "d: ";
        
        for (String node : nodeList) {
            sb.append(String.format(lineNumberFormatToken, lineNumber++));
            sb.append(node);
            sb.append(", children: [");
            
            if (!graph.getChildrenOf(node).isEmpty()) {
                for (String child : graph.getChildrenOf(node)) {
                    sb.append(child);
                    sb.append(" ");
                }

                sb.deleteCharAt(sb.length() - 1);
            }
                
            sb.append("], parents[");
            
            if (!graph.getParentsOf(node).isEmpty()) {
                for (String parent : graph.getParentsOf(node)) {
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
