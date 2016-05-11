package net.coderodde.msc;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

public class NodeCentricDeBruijnGraph extends AbstractDeBruijnGraph {

    public NodeCentricDeBruijnGraph(Collection<String> reads, int k) {
        super(k);
        Objects.requireNonNull(reads, "The collection of reads is null.");
        checkReadsAgainstKmerSize(reads, k);
        buildGraph(reads);
    }
    
    private void checkReadsAgainstKmerSize(final Collection<String> reads,
                                           final int k) {
        for (final String read : reads) {
            if (read.length() - 1 < k) {
                throw new IllegalArgumentException(
                        "A read is too short, length = " + read.length() + 
                        ". Must be at least " + (k + 1));
            }
        }
    }
    
    @Override
    public Set<Kmer> getAllNodes() {
        return Collections.<Kmer>unmodifiableSet(childrenMap.keySet());
    }
    
    @Override
    public Set<Kmer> getChildrenOf(Kmer kmer) {
        Objects.requireNonNull(kmer, "The input kmer is null.");
        checkStringIsKmer(kmer);
        return Collections.<Kmer>unmodifiableSet(childrenMap.get(kmer));
    }
    
    @Override
    public Set<Kmer> getParentsOf(Kmer kmer) {
        Objects.requireNonNull(kmer, "The input kmer is null.");
        checkStringIsKmer(kmer);
        return Collections.<Kmer>unmodifiableSet(parentsMap.get(kmer));
    }
    
    private void buildGraph(Collection<String> reads) {
        // Maps each kmer suffix to all the kmers that contain it.
        final Map<Kmer, Set<Kmer>> mapSuffixToKmers = new HashMap<>();
        // Maps each kmer prefix to all the kmers that contain it.
        final Map<Kmer, Set<Kmer>> mapPrefixToKmers = new HashMap<>();
        
        int index = 0;
        for (String string : reads) {
            if (string.length() < k) {
                throw new IllegalArgumentException(
                        "The length of the string \"" + string + "\" is too " +
                        "small (" + string.length() + "). Must be at least " +
                        k);
            }
            
            System.out.println("Processing read number " + index++);
            final int kmers = string.length() - k + 1;
            int previousKmers = 0;
            
            // Create all kmers of string 'string'.
            for (int i = 0; i < kmers; ++i) {
                Kmer kmer = new Kmer(string, i, k);
                
                if (childrenMap.containsKey(kmer)) {
                    // We already created a node for the 'kmer'.
                    continue;
                } else {
                    childrenMap.put(kmer, new HashSet<Kmer>());
                    parentsMap.put(kmer, new HashSet<Kmer>());
                }
                
                Kmer kmerPrefix = kmer.substring(0, k - 1);
                Kmer kmerSuffix = kmer.substring(1);
                                
                if (!mapPrefixToKmers.containsKey(kmerPrefix)) {
                    mapPrefixToKmers.put(kmerPrefix, new HashSet<Kmer>());
                }
                
                if (!mapSuffixToKmers.containsKey(kmerSuffix)) {
                    mapSuffixToKmers.put(kmerSuffix, new HashSet<Kmer>());
                }
                
                mapPrefixToKmers.get(kmerPrefix).add(kmer);
                mapSuffixToKmers.get(kmerSuffix).add(kmer);
                
                int ii = i / 1_000;
                
                if (previousKmers != ii) {
                    System.out.println(ii);
                    previousKmers = ii;
                }
            }
        }
        
        System.out.println("--- Building arcs ---");
        int prevArcs = 0;
        
        // Create arcs.
        for (Kmer kmer : childrenMap.keySet()) {
            Kmer kmerPrefix = kmer.substring(0, k - 1);
            Kmer kmerSuffix = kmer.substring(1);
            
            Set<Kmer> parentKmers = mapSuffixToKmers.get(kmerPrefix);
            Set<Kmer> childKmers  = mapPrefixToKmers.get(kmerSuffix);
            
            if (parentKmers != null) {
                for (Kmer parentKmer : parentKmers) {
                    parentsMap.get(kmer).add(parentKmer);
                    arcs++;
                }
            }
            
            if (childKmers != null) {
                for (Kmer childKmer : childKmers) {
                    childrenMap.get(kmer).add(childKmer);
                    arcs++;
                }
            }
            
            int a = arcs / 10_000;
            
            if (prevArcs != a) {
                prevArcs = a;
                System.out.println("Arcs: " + a);
            }
        }
    }
}
