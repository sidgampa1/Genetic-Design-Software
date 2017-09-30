package org.ucb.c5.composition;

import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.composition.model.Transcript;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.Random;

/**
 * This reverse translates a protein sequence to a DNA and chooses an RBS. It
 * uses the highest CAI codon for each amino acid in the specified Host.
 *
 * @author J. Christopher Anderson
 * @author Siddharth Gampa
 * 
 * TODO:  REPLACE WITH YOUR VERSION
 */
public class TranscriptDesigner {

    private Map<Character, String[]> aminoAcidToCodon;
    private RBSChooser rbsChooser;  //Delete the 2 to use old algorithm
    private SequenceChecker seqCheck;

    public void initiate() throws Exception {
        //Initialize the RBSChooser
        rbsChooser = new RBSChooser();  //Delete the 2 to use old algorithm
        seqCheck = new SequenceChecker();

        rbsChooser.initiate();
        seqCheck.initiate();

        //Construct a map between each amino acid and the highest-CAI codon for E coli
        aminoAcidToCodon = new HashMap<>();

        aminoAcidToCodon.put('A', new String[]{"GCG","GCA","GCC","GCT"});
        aminoAcidToCodon.put('C', new String[]{"TGC","TGT"});
        aminoAcidToCodon.put('D', new String[]{"GAT","GAC"});
        aminoAcidToCodon.put('E', new String[]{"GAA","GAG"});
        aminoAcidToCodon.put('F', new String[]{"TTC","TTT"});
        aminoAcidToCodon.put('G', new String[]{"GGT","GGC","GGA","GGG"});
        aminoAcidToCodon.put('H', new String[]{"CAC","CAT"});
        aminoAcidToCodon.put('I', new String[]{"ATC","ATT","ATA"});
        aminoAcidToCodon.put('K', new String[]{"AAA","AAG"});
        aminoAcidToCodon.put('L', new String[]{"CTG","CTA","CTC","CTT","TTA","TTG"});
        aminoAcidToCodon.put('M', new String[]{"ATG"});
        aminoAcidToCodon.put('N', new String[]{"AAC","AAT"});
        aminoAcidToCodon.put('P', new String[]{"CCG","CCA","CCC","CCT"});
        aminoAcidToCodon.put('Q', new String[]{"CAG","CAA"});
        aminoAcidToCodon.put('R', new String[]{"CGT", "CGC","CGA","CGG","AGA","AGG"});
        aminoAcidToCodon.put('S', new String[]{"TCT","TCC","TCA","TCG","AGC","AGT"});
        aminoAcidToCodon.put('T', new String[]{"ACC","ACT","ACA","ACG"});
        aminoAcidToCodon.put('V', new String[]{"GTT","GTC","GTA","GTG"});
        aminoAcidToCodon.put('W', new String[]{"TGG"});
        aminoAcidToCodon.put('Y', new String[]{"TAC","TAT"});
    }

    public Transcript run(String peptide, Set<RBSOption> ignores) throws Exception {
        //Choose codons for each amino acid
        String[] codons = new String[peptide.length()];
        Random random = new Random();
        for(int i=0; i<peptide.length(); i++) {
            char aa = peptide.charAt(i);
            String[] codon_opt = aminoAcidToCodon.get(aa);
            int rand = random.nextInt(codon_opt.length);
            String codon = aminoAcidToCodon.get(aa)[rand];
            codons[i] = codon;
        }
        
        //Choose an RBS
        StringBuilder cds = new StringBuilder();
        for(String codon : codons) {
            cds.append(codon);
        }
        boolean no_forbidden = seqCheck.run(cds.toString());

        RBSOption selectedRBS = null;
        if(no_forbidden) {
            selectedRBS = rbsChooser.run(cds.toString(), ignores);
        }
        else {
            System.out.println("creating new cds");
            run(peptide,ignores);
        }
        //Construct the Transcript and return it
        Transcript out = new Transcript(selectedRBS, peptide, codons);
        return out;
    }

    public static void main(String[] args) throws Exception {
        TranscriptDesigner td = new TranscriptDesigner();
        td.initiate();
        Set<RBSOption> rbs_exclude = new HashSet<>();
        String peptide = "MLSDTIDTKQQQQQLHVLFIDSYDSFTYNVVRLIEQQTDISPGVNAVHVTTVHSDTFQSMDQLLPLLPLFDAIVVGPGPGNPNNGAQDMGIISELFENANGKLDEVPILGICLGFQAMCLAQGADVSELNTIKHGQVYEMHLNDAARACGLFSGYPDTFKSTRYHSLHVNAEGIDTLLPLCTTEDENGILLMSAQTKNKPWFGVQYHPESCCSELGGLLVSNFLKLSFINNVKTGRWEKKKLNGEFSDILSRLDRTIDRDPIYKVKEKYPKGEDTTYVKQFEVSEDPKLTFEICNIIREEKFVMSSSVISENTGEWSIIALPNSASQVFTHYGAMKKTTVHYWQDSEISYTLLKKCLDGQDSDLPGSLEVIHEDKSQFWITLGKFMENKIIDNHREIPFIGGLVGILGYEIGQYIACGRCNDDENSLVPDAKLVFINNSIVINHKQGKLYCISLDNTFPVALEQSLRDSFVRKKNIKQSLSWPKYLPEEIDFIITMPDKLDYAKAFKKCQDYMHKGDSYEMCLTTQTKVVPSAVIEPWRIFQTLVQRNPAPFSSFFEFKDIIPRQDETPPVLCFLSTSPERFLKWDADTCELRPIKGTVKKGPQMNLAKATRILKTPKEFGENLMILDLIRNDLYELVPDVRVEEFMSVQEYATVYQLVSVVKAHGLTSASKKTRYSGIDVLKHSLPPGSMTGAPKKITVQLLQDKIESKLNKHVNGGARGVYSGVTGYWSVNSNGDWSVNIRCMYSYNGGTSWQLGAGGAITVLSTLDGELEEMYNKLESNLQIFM";
        td.run(peptide, rbs_exclude);
//        int num_codons = 0;
//        for (String[] chars: td.aminoAcidToCodon.values()) {
//            num_codons += chars.length;
//        }
//        System.out.println(num_codons);
    }
}
