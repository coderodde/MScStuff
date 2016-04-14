package net.coderodde.msc;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

public class EdgeCentricDeBruijnGraph {
    
    private final Map<String, Set<String>> childrenMap = new HashMap<>();
    private final Map<String, Set<String>> parentMap   = new HashMap<>();
    private final int k;
    
    public EdgeCentricDeBruijnGraph(Collection<String> reads, int k) {
        Objects.requireNonNull(reads, "The collection of reads is null.");
        checkKMerSize(k);
        this.k = k;
        
        buildGraph(reads);
    }
    
    private void buildGraph(Collection<String> reads) {
        // Maps each kmer suffix to all the kmers that contain it.
        final Map<String, Set<String>> mapSuffixToKmers = new HashMap<>();
        // Maps each kmer prefix to all the kmers that contain it.
        final Map<String, Set<String>> mapPrefixToKmers = new HashMap<>();
        
        for (String string : reads) {
            if (string.length() < k) {
                throw new IllegalArgumentException(
                        "The length of the string \"" + string + "\" is too " +
                        "small (" + string.length() + "). Must be at least " +
                        k);
            }
            
            final int kmers = string.length() - k + 1;
            
            // Create all kmers of string 'string'.
            for (int i = 0; i < kmers; ++i) {
                String kmer = string.substring(i, i + k);
                
                if (childrenMap.containsKey(kmer)) {
                    continue;
                }
                
                if (!childrenMap.containsKey(kmer)) {
                    childrenMap.put(kmer, new HashSet<>());
                }
                
                if (!parentMap.containsKey(kmer)) {
                    parentMap.put(kmer, new HashSet<>());
                }
                
                String kmerPrefix = kmer.substring(0, k - 1);
                String kmerSuffix = kmer.substring(1);
                                
                if (!mapPrefixToKmers.containsKey(kmerPrefix)) {
                    mapPrefixToKmers.put(kmerPrefix, new HashSet<>());
                }
                
                mapPrefixToKmers.get(kmerPrefix).add(kmer);
                
                if (!mapSuffixToKmers.containsKey(kmerSuffix)) {
                    mapSuffixToKmers.put(kmerSuffix, new HashSet<>());
                }
                
                mapSuffixToKmers.get(kmerSuffix).add(kmer);
            }
        }
        
        // Create edges.
        childrenMap.keySet().stream().forEach((kmer) -> {
            String kmerPrefix = kmer.substring(0, k - 1);
            String kmerSuffix = kmer.substring(1);
            
            Set<String> parentKmers = mapSuffixToKmers.get(kmerSuffix);
            Set<String> childKmers  = mapPrefixToKmers.get(kmerPrefix);
            
            parentKmers.stream().forEach((parentKmer) -> {
                parentMap.get(kmer).add(parentKmer);
            });
            
            childKmers.stream().forEach((childKmer) -> {
                childrenMap.get(kmer).add(childKmer);
            });
        });
    }
    
    private void checkKMerSize(int k) {
        if (k < 2) {
            throw new IllegalArgumentException(
                    "The k is too small: " + k + ". Must be at least 2.");
        }
    }
}
