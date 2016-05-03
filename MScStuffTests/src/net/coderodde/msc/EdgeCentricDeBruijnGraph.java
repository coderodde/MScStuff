package net.coderodde.msc;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

public class EdgeCentricDeBruijnGraph extends AbstractDeBruijnGraph {
    
    public EdgeCentricDeBruijnGraph(Collection<String> reads, int k) {
        super(k);
        Objects.requireNonNull(reads, "The collection of reads is null.");
        buildGraph(reads);
    }
    
    @Override
    public Set<String> getAllNodes() {
        return Collections.<String>unmodifiableSet(childrenMap.keySet());
    }
    
    @Override
    public Set<String> getChildrenOf(String kmer) {
        Objects.requireNonNull(kmer, "The input kmer is null.");
        checkStringIsKmer(kmer);
        return Collections.<String>unmodifiableSet(childrenMap.get(kmer));
    }
    
    @Override
    public Set<String> getParentsOf(String kmer) {
        Objects.requireNonNull(kmer, "The input kmer is null.");
        checkStringIsKmer(kmer);
        return Collections.<String>unmodifiableSet(parentsMap.get(kmer));
    }
    
    private void buildGraph(Collection<String> reads) {
        // Maps each kmer suffix to all the kmers that contain it.
        final Map<String, Set<String>> mapSuffixToKmers = new HashMap<>();
        // Maps each kmer prefix to all the kmers that contain it.
        final Map<String, Set<String>> mapPrefixToKmers = new HashMap<>();
        final Set<String> allKmers = new HashSet<>();
        
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
                
                allKmers.add(kmer);
                
                if (childrenMap.containsKey(kmer)) {
                    continue;
                }
                
                if (!childrenMap.containsKey(kmer)) {
                    childrenMap.put(kmer, new HashSet<>());
                }
                
                if (!parentsMap.containsKey(kmer)) {
                    parentsMap.put(kmer, new HashSet<>());
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
        for (String kmer : childrenMap.keySet()) {
            String kmerPrefix = kmer.substring(0, k - 1);
            String kmerSuffix = kmer.substring(1);
            
            Set<String> parentKmers = mapSuffixToKmers.get(kmerPrefix);
            Set<String> childKmers  = mapPrefixToKmers.get(kmerSuffix);
            StringBuilder sb = new StringBuilder(k);
            
            if (parentKmers != null) {
                for (String parentKmer : parentKmers) {
                    sb.delete(0, sb.length());
                    sb.append(parentKmer);
                    sb.deleteCharAt(sb.length() - 1);
                    sb.append(kmer.charAt(kmer.length() - 1));
                    String tmp = sb.toString();
                    
                    if (allKmers.contains(tmp)) {
                        parentsMap.get(kmer).add(parentKmer);
                    }
                }
            }
            
            if (childKmers != null) {
                for (String childKmer : childKmers) {
                    sb.delete(0, sb.length());
                    sb.append(kmer);
                    sb.deleteCharAt(sb.length() - 1);
                    sb.append(childKmer.charAt(childKmer.length() - 1));
                    String tmp = sb.toString();
                    
                    if (allKmers.contains(tmp)) {
                        childrenMap.get(kmer).add(childKmer);
                    }
                }
            }
        }
    }
}
