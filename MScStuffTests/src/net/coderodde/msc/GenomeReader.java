package net.coderodde.msc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * This class is responsible for reading a genome from a file and storing it in
 * a string over character <tt>A, C, G, T, N</tt>.
 * 
 * @author Rodion "rodde" Efremov
 * @version 1.6 (May 10, 2016)
 */
public class GenomeReader {
   
    private static final String FNA_FILE_EXTENSION = "fna";
    private static final String TOKEN_COMPLETE_GENOME = "complete genome";
    
    private GenomeReader() {}
    
    public static String readFile(final File file) {
        if (!file.getName().endsWith("." + FNA_FILE_EXTENSION)) {
            Logger.getLogger(GenomeReader.class.getName())
                  .log(Level.SEVERE, "The input file does not have the " +
                                     "extension \"" + FNA_FILE_EXTENSION + 
                                     "\".");
            return null;
        }
        
        Scanner scanner;
        
        try {
            scanner = new Scanner(new FileReader(file));
            
            final String header = scanner.nextLine();
            
//            if (!header.contains(TOKEN_COMPLETE_GENOME)) {
//                Logger.getLogger(GenomeReader.class.getName())
//                      .log(Level.SEVERE, "The reader expected a header, " + 
//                                         "but none was found.");
//                return null;
//            }
            
            final StringBuilder all = new StringBuilder((int) file.length());
            final StringBuilder work = new StringBuilder();
            
            while (scanner.hasNextLine()) {
                final String line = scanner.nextLine();
                work.delete(0, work.length());
                work.append(line);
                
                final int workBuilderLength = work.length();
                
                for (int i = 0; i < workBuilderLength; ++i) {
                    final char c = Character.toUpperCase(work.charAt(i));
                    work.setCharAt(i, fix(c));
                }
                
                all.append(work.toString());
            }
            
            return all.toString();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(GenomeReader.class.getName())
                  .log(Level.SEVERE, "", ex);
            return null;
        }
    }
    
    private static char fix(final char c) {
        switch (c) {
            case 'a':
                System.out.println("a");
            case 'A':
                return 'A';
                
            case 'c':
                System.out.println("c");
            case 'C':
                return 'C';
                
            case 'g':
                System.out.println("g");
            case 'G': 
                return 'G';
                
            case 't':
                System.out.println("t");
            case 'T': 
                return 'T';
                
            default:
                return 'N';
        }
    }
}
