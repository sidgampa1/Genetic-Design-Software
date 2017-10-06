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

        if (peptide == "") {
            throw new IllegalArgumentException();
        }
        //Choose best codon for each amino acid
        String[] codons = seqChooser.run(peptide);

        
        //Choose an RBS
        StringBuilder cds = new StringBuilder();
        for(String codon : codons) {
            cds.append(codon);
        }

        RBSOption selectedRBS;
        selectedRBS = rbsChooser.run(cds.toString(), peptide, ignores);

        System.out.println("made it to the rbs + cds check");
        if (seqCheck.run(cds.toString())) {
            System.out.println("cds is forbidden seq free!");
        }
        String total_seq = selectedRBS.getRbs() + 'x' + cds.toString();
        boolean done = seqCheck.run(total_seq);
        if (done) {
            System.out.println("Woot Woot completed dna seq construction!!!");
            System.out.println("peptide: " + peptide);
            System.out.println("rbs + cds: " + total_seq);

        }
        else {
            System.out.println("rbs + cds check failed");
        }

        //Construct the Transcript and return it
        Transcript out = new Transcript(selectedRBS, peptide, codons);
        return out;
    }

    public static void main(String[] args) throws Exception {
        final long startTime = System.currentTimeMillis();
        TranscriptDesigner td = new TranscriptDesigner();
        td.initiate();
        Set<RBSOption> rbs_exclude = new HashSet<>();
        String peptide = "MLSDTIDTKQQQQQLHVLFIDSYDSFTYNVVRLIEQQTDISPGVNAVHVTTVHSDTFQSMDQLLPLLPLFDAIVVGPGPGNPNNGAQDMGIISELFENANGKLDEVPILGICLGFQAMCLAQGADVSELNTIKHGQVYEMHLNDAARACGLFSGYPDTFKSTRYHSLHVNAEGIDTLLPLCTTEDENGILLMSAQTKNKPWFGVQYHPESCCSELGGLLVSNFLKLSFINNVKTGRWEKKKLNGEFSDILSRLDRTIDRDPIYKVKEKYPKGEDTTYVKQFEVSEDPKLTFEICNIIREEKFVMSSSVISENTGEWSIIALPNSASQVFTHYGAMKKTTVHYWQDSEISYTLLKKCLDGQDSDLPGSLEVIHEDKSQFWITLGKFMENKIIDNHREIPFIGGLVGILGYEIGQYIACGRCNDDENSLVPDAKLVFINNSIVINHKQGKLYCISLDNTFPVALEQSLRDSFVRKKNIKQSLSWPKYLPEEIDFIITMPDKLDYAKAFKKCQDYMHKGDSYEMCLTTQTKVVPSAVIEPWRIFQTLVQRNPAPFSSFFEFKDIIPRQDETPPVLCFLSTSPERFLKWDADTCELRPIKGTVKKGPQMNLAKATRILKTPKEFGENLMILDLIRNDLYELVPDVRVEEFMSVQEYATVYQLVSVVKAHGLTSASKKTRYSGIDVLKHSLPPGSMTGAPKKITVQLLQDKIESKLNKHVNGGARGVYSGVTGYWSVNSNGDWSVNIRCMYSYNGGTSWQLGAGGAITVLSTLDGELEEMYNKLESNLQIFM";
        td.run(peptide, rbs_exclude);

        final long endTime = System.currentTimeMillis();
        System.out.println("total runtime: " + (endTime - startTime));

    }
}
