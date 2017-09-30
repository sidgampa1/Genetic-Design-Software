package org.ucb.c5.composition;

import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.composition.model.Transcript;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

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

    private SequenceChooser seqChooser;
    private RBSChooser rbsChooser;  //Delete the 2 to use old algorithm
    private SequenceChecker seqCheck;

    public void initiate() throws Exception {
        //Initialize the RBSChooser
        rbsChooser = new RBSChooser();  //Delete the 2 to use old algorithm
        seqChooser = new SequenceChooser();
        seqCheck = new SequenceChecker();

        rbsChooser.initiate();
        seqChooser.initiate();
        seqCheck.initiate();

        //Construct a map between each amino acid and the highest-CAI codon for E coli
    }

    public Transcript run(String peptide, Set<RBSOption> ignores) throws Exception {

//        Random random = new Random();
//        for(int i=0; i<peptide.length(); i++) {
//            char aa = peptide.charAt(i);
//            String[] codon_opt = aminoAcidToCodon.get(aa);
//            int rand = random.nextInt(codon_opt.length);
//            String codon = aminoAcidToCodon.get(aa)[rand];
//            codons[i] = codon;
//        }

        //Choose best codon for each amino acid
        String[] codons = seqChooser.run(peptide);

        
        //Choose an RBS
        StringBuilder cds = new StringBuilder();
        for(String codon : codons) {
            cds.append(codon);
        }

        RBSOption selectedRBS;
        selectedRBS = rbsChooser.run(cds.toString(), ignores);

        System.out.println("made it to the rbs + cds check");
        String total_seq = selectedRBS + cds.toString();
        boolean done = seqCheck.run(total_seq);
        if (done) {
            System.out.println("Woot Woot completed dna seq construction!!!");
        }

        //Construct the Transcript and return it
        Transcript out = new Transcript(selectedRBS, peptide, codons);
        return out;
    }

    public static void main(String[] args) throws Exception {
        TranscriptDesigner td = new TranscriptDesigner();
        td.initiate();
        Set<RBSOption> rbs_exclude = new HashSet<>();
        // removed M from end of peptide to make it multiple of 3
        String peptide = "MLSDTIDTKQQQQQLHVLFIDSYDSFTYNVVRLIEQQTDISPGVNAVHVTTVHSDTFQSMDQLLPLLPLFDAIVVGPGPGNPNNGAQDMGIISELFENANGKLDEVPILGICLGFQAMCLAQGADVSELNTIKHGQVYEMHLNDAARACGLFSGYPDTFKSTRYHSLHVNAEGIDTLLPLCTTEDENGILLMSAQTKNKPWFGVQYHPESCCSELGGLLVSNFLKLSFINNVKTGRWEKKKLNGEFSDILSRLDRTIDRDPIYKVKEKYPKGEDTTYVKQFEVSEDPKLTFEICNIIREEKFVMSSSVISENTGEWSIIALPNSASQVFTHYGAMKKTTVHYWQDSEISYTLLKKCLDGQDSDLPGSLEVIHEDKSQFWITLGKFMENKIIDNHREIPFIGGLVGILGYEIGQYIACGRCNDDENSLVPDAKLVFINNSIVINHKQGKLYCISLDNTFPVALEQSLRDSFVRKKNIKQSLSWPKYLPEEIDFIITMPDKLDYAKAFKKCQDYMHKGDSYEMCLTTQTKVVPSAVIEPWRIFQTLVQRNPAPFSSFFEFKDIIPRQDETPPVLCFLSTSPERFLKWDADTCELRPIKGTVKKGPQMNLAKATRILKTPKEFGENLMILDLIRNDLYELVPDVRVEEFMSVQEYATVYQLVSVVKAHGLTSASKKTRYSGIDVLKHSLPPGSMTGAPKKITVQLLQDKIESKLNKHVNGGARGVYSGVTGYWSVNSNGDWSVNIRCMYSYNGGTSWQLGAGGAITVLSTLDGELEEMYNKLESNLQIF";
        td.run(peptide, rbs_exclude);

    }
}
