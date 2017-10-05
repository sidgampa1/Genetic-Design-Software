package org.ucb.c5.composition;

import java.util.HashMap;
import java.util.Map;

/**
 * class to hold an amino acid to codon table
 */
public class AminoAcidToCodon {
    public Map<Character, String[]> table;

    public void initiate() throws Exception {
        table = new HashMap<>();
        table.put('A', new String[]{"GCG", "GCA", "GCC", "GCT"});
        table.put('C', new String[]{"TGC", "TGT"});
        table.put('D', new String[]{"GAT", "GAC"});
        table.put('E', new String[]{"GAA", "GAG"});
        table.put('F', new String[]{"TTC", "TTT"});
        table.put('G', new String[]{"GGT", "GGC", "GGA", "GGG"});
        table.put('H', new String[]{"CAC", "CAT"});
        table.put('I', new String[]{"ATC", "ATT", "ATA"});
        table.put('K', new String[]{"AAA", "AAG"});
        table.put('L', new String[]{"CTG", "CTA", "CTC", "CTT", "TTA", "TTG"});
        table.put('M', new String[]{"ATG"});
        table.put('N', new String[]{"AAC", "AAT"});
        table.put('P', new String[]{"CCG", "CCA", "CCC", "CCT"});
        table.put('Q', new String[]{"CAG", "CAA"});
        table.put('R', new String[]{"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"});
        table.put('S', new String[]{"TCT", "TCC", "TCA", "TCG", "AGC", "AGT"});
        table.put('T', new String[]{"ACC", "ACT", "ACA", "ACG"});
        table.put('V', new String[]{"GTT", "GTC", "GTA", "GTG"});
        table.put('W', new String[]{"TGG"});
        table.put('Y', new String[]{"TAC", "TAT"});
    }

    public void run() throws Exception {
        //does nothing
    }
}
