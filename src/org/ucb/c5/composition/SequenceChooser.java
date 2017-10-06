package org.ucb.c5.composition;

import com.sun.org.apache.xpath.internal.SourceTree;
import com.sun.xml.internal.bind.v2.TODO;
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
    private List<String> aa_windows;
    private SequenceChecker seqCheck;
    private AminoAcidToCodon translate;
    private HairpinCounter hairpin;
    

    public void initiate() throws Exception {
        aa_windows = new ArrayList<>();
        seqCheck = new SequenceChecker();
        translate = new AminoAcidToCodon();
        hairpin = new HairpinCounter();
        seqCheck.initiate();
        translate.initiate();
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
     */
    private List<DNAPermutation> getValidPerms(String preamble, String aa_window) throws Exception{
        int num_poss = 1;
        for(int i = 0; i < aa_window.length() && num_poss < 100; i++) {
            num_poss *= translate.table.get(aa_window.charAt(i)).length;
            if (num_poss > 100) {
                num_poss = 100;
            }
        }
        System.out.println("amino acid window is " + aa_window);
        System.out.println("num possiblities is " + num_poss);
        List<DNAPermutation> perms = new ArrayList<>();
        Random rand = new Random(100); //seeded to ensure consistency
        int counter = 0;
        int tot_tried_perms = 0;
        while(counter < num_poss) {
            StringBuilder perm = new StringBuilder();
            for (int j = 0; j < aa_window.length(); j++) {
                char aa = aa_window.charAt(j);
                if (!translate.table.containsKey(aa)) {
                    throw new IllegalArgumentException("amino acid is not valid");
                }
                String[] possible_codons = translate.table.get(aa);
                int codon_ind = rand.nextInt(possible_codons.length);
                String codon = possible_codons[codon_ind];
                perm.append(codon);
            }

            tot_tried_perms++;
            String perm_seq = perm.toString();
            double perm_GC = GCCheck(perm_seq);
            double perm_hp = hairpin.run(perm_seq);
            String total_perm = preamble + perm_seq;
            DNAPermutation good_perm = new DNAPermutation(perm_seq, perm_GC, perm_hp);

            boolean isUnique = !perms.contains(good_perm);
            boolean isFBfree = seqCheck.run(total_perm);
            boolean isValidGC = isValidGC(total_perm);
            boolean isGoodHairpin = perm_hp < 30;
            if (num_poss == 100 && isUnique && isFBfree && isValidGC && isGoodHairpin){
                perms.add(good_perm);
                counter++;
//                System.out.println("perms length is " + perms.size());
            }
            else if ((num_poss < 100 || tot_tried_perms > 1000) && isFBfree) {
                perms.add(good_perm);
                counter++;
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

            peptide = peptide.toUpperCase();
            String aa_sub_window = target_aas + peptide.substring(downstream_start, downstream_end);
            System.out.println("created target_aas " + i);
            List<DNAPermutation> dna_perms = getValidPerms(preamble.toString(), aa_sub_window);
            System.out.println("I have a collection of good perms for window " + i);
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
//            System.out.println("i have a sorted collection of good perms by GC");

            // TODO create weighting system to score perms (not just choosing first perm)
            String best_perm = dna_perms.get(dna_perms.size() - 1).getSeq();
            preamble.append(best_perm.substring(0,target_aas.length() * 3));
            System.out.println("Found a 3 codon seq for index" + i);
        }

        System.out.println("GC content of chosen DNA seq: " + GCCheck(preamble.toString()));
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
