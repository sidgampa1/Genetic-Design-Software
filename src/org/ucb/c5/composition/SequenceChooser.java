package org.ucb.c5.composition;

import com.sun.org.apache.xpath.internal.SourceTree;
import com.sun.xml.internal.bind.v2.TODO;

import java.util.*;


/**
 * This reverse translates and chooses the best codons to use for a given protein
 * using the GeneOptimizer algorithm (enumerating and ranking a sliding window)
 *
 * @author Siddharth Gampa
 *
 */
public class SequenceChooser {
    private List<String> aa_windows;
    private SequenceChecker seqCheck;
    private AminoAcidToCodon translate;
    

    public void initiate() throws Exception {
        aa_windows = new ArrayList<>();
        seqCheck = new SequenceChecker();
        translate = new AminoAcidToCodon();
        seqCheck.initiate();
        translate.initiate();
    }

    private double GCCheck(String seq) {
        double GC_content = 0;
        for (int i = 0; i < seq.length(); i++) {
            char nuc = seq.charAt(i);
            if (nuc == 'G' || nuc == 'C') {
                GC_content += 1;
            }
        }
        return GC_content / seq.length();
    }

    private boolean isValidGC(String seq) {
        double GC = GCCheck(seq);
//        System.out.println("GC content is " + GC);
        if (GC > 0.40 && GC < 0.60) {
//            System.out.println("GC content is " + GC);
            return true;
        }
        return false;
    }

    /**
     * helper method to create permutations of aa seq
     * @param aa_window
     * @return perms = 100 random, valid (no forbidden seqs and good GC content) permutations
     * of dna seqs for aa window
     * TODO I dont think random enumeration is best, maybe try something else
     */
    private DNAPermutation[] getValidPerms(String preamble, String aa_window) {
        DNAPermutation[] perms = new DNAPermutation[100];
        Random rand = new Random(100); //seeded to ensure consistency
        int counter = 0;
        while(counter < 100) {
            String perm = "";
            for (int j = 0; j < aa_window.length(); j++) {
                char aa = aa_window.charAt(j);
                if (!translate.table.containsKey(aa)) {
                    throw new IllegalArgumentException("amino acid is not valid");
                }
                String[] possible_codons = translate.table.get(aa);
                int codon_ind = rand.nextInt(possible_codons.length);
                String codon = possible_codons[codon_ind];
                perm += codon;
            }

            double perm_GC = GCCheck(perm);
            String total_perm = preamble + perm;
            if (seqCheck.run(total_perm) && isValidGC(total_perm)) {
                DNAPermutation good_perm = new DNAPermutation(perm, perm_GC);
                perms[counter] = good_perm;
                counter++;
            }
        }
        return perms;
    }


    public String[] run(String peptide) {
        String[] codons = new String[peptide.length()];
        /**
         * for each window of 3 aa:
         *  grab preamble + 3 aa of interest + 6 aa downstream
         *  construct 100 permutations of nucs with aa seq (preamble already chosen)
         *  throw out forbidden seq permutations, choose best from rest
         **/
        String preamble = "";
        for (int i = 0; i < peptide.length(); i += 3) {
            String target_aas;
            int downstream_start = i + 3;
            int downstream_end = i+ 9;

            // account for edge cases at end of protein windows
            if ((i + 3) >= peptide.length()) {
                target_aas = peptide.substring(i, peptide.length());
                downstream_end = 0;
                downstream_start = 0;
            }
            else {
                target_aas = peptide.substring(i, i + 3);
            }
            if (i + 9 >= peptide.length() && !((i+3) >= peptide.length())) {
                downstream_end = peptide.length();
            }

            String aa_sub_window = target_aas + peptide.substring(downstream_start, downstream_end);
            DNAPermutation[] dna_perms = getValidPerms(preamble, aa_sub_window);
            Arrays.sort(dna_perms);

            // TODO create weighting system to score perms (not just choosing first perm)
            String best_perm = dna_perms[dna_perms.length - 1].getSeq();
            preamble += best_perm.substring(0,target_aas.length() * 3);
//            System.out.println("Found a 3 codon seq for index" + i);
        }

        System.out.println("GC content of chosen DNA seq: " + GCCheck(preamble));
        // turn complete dna seq into codon list to return
        int counter = 0;
        for (int j = 0; j < preamble.length(); j += 3) {
            String codon = preamble.substring(j, j + 3);
            codons[counter] = codon;
            counter++;
        }

        return codons;
    }
}
