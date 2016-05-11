package net.coderodde.msc;

/**
 * This class implements a <tt>k</tt>-mer, a substring of length <tt>k</tt>.
 * 
 * @author Rodion "rodde" Efremov
 * @version 1.6 (May 11, 2016)
 */
public final class Kmer implements Comparable<Kmer> {
    
    private final String string;
    private final int startIndex;
    private final int length;
    private final int hashCode;
    
    public Kmer(final String string, final int startIndex, final int length) {
        this.string     = string;
        this.startIndex = startIndex;
        this.length     = length;
        this.hashCode   = string.substring(startIndex, 
                                           startIndex + length).hashCode();
    }
    
    @Override
    public int hashCode() {
        return hashCode;
    }
    
    @Override
    public boolean equals(Object o) {
        final Kmer other = (Kmer) o;
        
        if (hashCode != other.hashCode) {
            return false;
        }
        
        return string.substring(startIndex, 
                                startIndex + length)
                     .equals(other.string
                                  .substring(other.startIndex,
                                             other.startIndex + other.length));
    }
    
    @Override
    public String toString() {
        return string.substring(startIndex, startIndex + length);
    }
    
    public int length() {
        return length;
    }
    
    public Kmer substring(final int startIndex, final int length) {
        final int actualStartIndex = this.startIndex + startIndex;
        return new Kmer(string, actualStartIndex, length);
    }
    
    public Kmer substring(final int startIndex) {
        return substring(startIndex, length - startIndex);
    }
    
    public char charAt(final int index) {
        return string.charAt(this.startIndex + index);
    }

    @Override
    public int compareTo(Kmer o) {
        return toString().compareTo(o.toString());
    }
}
