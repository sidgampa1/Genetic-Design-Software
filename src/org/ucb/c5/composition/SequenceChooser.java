package org.ucb.c5.composition;

import com.sun.org.apache.xpath.internal.SourceTree;
import com.sun.xml.internal.bind.v2.TODO;

import java.util.HashMap;
import java.util.Map;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;



/**
 * This reverse translates and chooses the best codons to use for a given protein
 * using the GeneOptimizer algorithm (enumerating and ranking a sliding window)
 *
 * @author Siddharth Gampa
 *
 */
public class SequenceChooser {
    private Map<Character, String[]> aminoAcidToCodon;
    private List<String> aa_windows;
    private SequenceChecker seqCheck;
    

    public void initiate() throws Exception {
        aminoAcidToCodon = new HashMap<>();
        aa_windows = new ArrayList<>();
        seqCheck = new SequenceChecker();
        seqCheck.initiate();
        
        aminoAcidToCodon.put('A', new String[]{"GCG", "GCA", "GCC", "GCT"});
        aminoAcidToCodon.put('C', new String[]{"TGC", "TGT"});
        aminoAcidToCodon.put('D', new String[]{"GAT", "GAC"});
        aminoAcidToCodon.put('E', new String[]{"GAA", "GAG"});
        aminoAcidToCodon.put('F', new String[]{"TTC", "TTT"});
        aminoAcidToCodon.put('G', new String[]{"GGT", "GGC", "GGA", "GGG"});
        aminoAcidToCodon.put('H', new String[]{"CAC", "CAT"});
        aminoAcidToCodon.put('I', new String[]{"ATC", "ATT", "ATA"});
        aminoAcidToCodon.put('K', new String[]{"AAA", "AAG"});
        aminoAcidToCodon.put('L', new String[]{"CTG", "CTA", "CTC", "CTT", "TTA", "TTG"});
        aminoAcidToCodon.put('M', new String[]{"ATG"});
        aminoAcidToCodon.put('N', new String[]{"AAC", "AAT"});
        aminoAcidToCodon.put('P', new String[]{"CCG", "CCA", "CCC", "CCT"});
        aminoAcidToCodon.put('Q', new String[]{"CAG", "CAA"});
        aminoAcidToCodon.put('R', new String[]{"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"});
        aminoAcidToCodon.put('S', new String[]{"TCT", "TCC", "TCA", "TCG", "AGC", "AGT"});
        aminoAcidToCodon.put('T', new String[]{"ACC", "ACT", "ACA", "ACG"});
        aminoAcidToCodon.put('V', new String[]{"GTT", "GTC", "GTA", "GTG"});
        aminoAcidToCodon.put('W', new String[]{"TGG"});
        aminoAcidToCodon.put('Y', new String[]{"TAC", "TAT"});
    }

    public String[] run(String peptide) {
        String[] codons = new String[peptide.length()];

        // add sliding windows in peptide seq
        // TODO this wont work for peptide length that's not a multiple of 3
        for (int i = 0; i < peptide.length(); i += 3) {
            aa_windows.add(peptide.substring(i, i+3));
        }

        // for each window of 3 aa, grab best codons (no fb), retain middle one
        int codon_count = 0;
        for (String aas: aa_windows) {
            int codon_index = 0;
            String dna_window = "";
            Random rand = new Random();
            boolean done = false;
            while(!done) {
                dna_window = "";
                for (int j = 0; j < 3; j++) {
                    char aa = aas.charAt(j);
                    int index = rand.nextInt(aminoAcidToCodon.get(aa).length);
                    dna_window += (aminoAcidToCodon.get(aa)[index]);
                }
                System.out.println(dna_window);
                done = seqCheck.run(dna_window);
            }
             // if seq w/ no fb found, add middle codon to codon list
            /** error here if peptide length not multiple of 3
             * the last dna window will be < 3  because of 3 aa windows
             * TODO fix potentional error
            **/

            if (dna_window.length() == 9){
                System.out.println("Success, new codon match found!");
                codon_count += 1;
                System.out.println("codon count: " + codon_count);
                codons[codon_index] = dna_window.substring(3,6);
            }
        }

        return codons;
    }
}
