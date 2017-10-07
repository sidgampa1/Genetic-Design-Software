package org.ucb.c5.composition;

import org.ucb.c5.sequtils.HairpinCounter;

import java.util.*;


/**
 * This reverse translates and chooses the best codons to use for a given protein
 * using the GeneOptimizer algorithm (enumerating and ranking a sliding window)
 *
 * @author Siddharth Gampa
 *
 */
public class SequenceChooser {
    private SequenceChecker seqCheck;
    private AminoAcidToCodon translator;
    private HairpinCounter hairpin;
    

    public void initiate() throws Exception {
        // initiate tools and tables
        seqCheck = new SequenceChecker();
        translator = new AminoAcidToCodon();
        hairpin = new HairpinCounter();

        seqCheck.initiate();
        translator.initiate();
        hairpin.initiate();
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
        return (GC > 0.40 && GC < 0.60);
    }

    /**
     * helper method to create permutations of aa seq
     * @param aa_window - sliding window of aa's to optimize codons for
     * @return perms = 100 random, valid (no forbidden seqs) permutations
     *                  of dna seqs for aa window
     */
    private Set<DNAPermutation> getValidPerms(String preamble, String aa_window) throws Exception{
        //count number of possible codon permutations for window
        int num_poss = 1;
        for(int i = 0; i < aa_window.length() && num_poss < 100; i++) {
            char aa = aa_window.charAt(i);
            if (!translator.table.containsKey(aa)) {
                throw new IllegalArgumentException("amino acid is not valid: " + aa);
            }
            num_poss *= translator.table.get(aa).length;
        }

        // limit permutations to 100 to speed up process
        if (num_poss > 100) {
            num_poss = 100;
        }

        // generate DNA permutations from RNG, keeping ones that have no forbidden seqs and good GC content and hairpins
        // if possible
        Set<DNAPermutation> perms = new HashSet<>(); //Hashset to ensure no duplicate perms
        Random rand = new Random(100); //seeded to ensure consistency
        int num_good_perms = 0;
        int total_tried_perms = 0;

        // iterate through permutations until we find either all possible perms or 100 perms if large perm possiblities
        while(num_good_perms < num_poss) {
            StringBuilder perm = new StringBuilder();

            //grab one random codon per amino acid to create a permutation
            for (int j = 0; j < aa_window.length(); j++) {
                char aa = aa_window.charAt(j);
                if (!translator.table.containsKey(aa)) {
                    throw new IllegalArgumentException("amino acid is not valid");
                }
                String[] possible_codons = translator.table.get(aa);
                int codon_ind = rand.nextInt(possible_codons.length);
                String codon = possible_codons[codon_ind];
                perm.append(codon);
            }

            total_tried_perms++;
            String perm_seq = perm.toString();
            double perm_GC = GCCheck(perm_seq);
            double perm_hp = hairpin.run(perm_seq);
            String total_perm = preamble + perm_seq;
            DNAPermutation good_perm = new DNAPermutation(perm_seq, perm_GC, perm_hp);

            boolean isFBfree = seqCheck.run(total_perm);
            boolean isValidGC = isValidGC(total_perm);

            // only add to perms list if it meets min reqs for structure, or if there are few possible perms or if
            // we have already iterated through many perms without finding an optimal one
            if (num_poss == 100 && isFBfree && isValidGC){
                perms.add(good_perm);
                num_good_perms++; }
            else if ((num_poss < 100 || total_tried_perms > 1000) && isFBfree) {
                perms.add(good_perm);
                num_good_perms++;
            }
        }
        return perms;
    }


    public String[] run(String peptide) throws Exception{
        String[] codons = new String[peptide.length()];
        /**
         * for each window of 3 aa (GeneOptimizer algorithm):
         *  grab preamble + 3 aa of interest + 6 aa downstream
         *  construct 100 permutations of nucs with aa seq (preamble already chosen)
         *  throw out forbidden seq permutations, choose best from rest
         **/
        StringBuilder preamble = new StringBuilder();
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
            // convert to list to allow sorting of permutations by attributes
            List<DNAPermutation> dna_perms = new ArrayList(getValidPerms(preamble.toString(), aa_sub_window));

            //sort dna permutations for this window by checking hairpin count first, then good GC
            Collections.sort(dna_perms, new Comparator<DNAPermutation>() {
                @Override
                public int compare(DNAPermutation t1, DNAPermutation t2) {
                    double GC1 = t1.getGC_content();
                    double GC2 = t2.getGC_content();
                    double hp1 = t1.getHairpin();
                    double hp2 = t2.getHairpin();

                    if (hp1 == hp2 && hp1 == 0) {
                        if (GC1 > GC2)
                            return 1;
                        else if (GC2 > GC1) {
                            return -1;
                        }
                    }
                    else if (hp1 > hp2) {
                        return -1;
                    }
                    else if (hp2 > hp1) {
                        return 1;
                    }
                    return 0;
                }
            });

            // keep best permutation and append it to the preamble which will be considered when optimizing next window
            String best_perm = dna_perms.get(dna_perms.size() - 1).getSeq();
            preamble.append(best_perm.substring(0,target_aas.length() * 3));
        }

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
