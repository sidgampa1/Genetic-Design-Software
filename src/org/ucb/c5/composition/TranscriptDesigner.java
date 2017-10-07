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
    private RBSChooser rbsChooser;
    private SequenceChecker seqCheck;

    public void initiate() throws Exception {
        //Initialize the RBSChooser
        rbsChooser = new RBSChooser();
        seqChooser = new SequenceChooser();
        seqCheck = new SequenceChecker();

        rbsChooser.initiate();
        seqChooser.initiate();
        seqCheck.initiate();

    }

    public Transcript run(String peptide, Set<RBSOption> ignores) throws Exception {

        if (peptide == "") {
            throw new IllegalArgumentException();
        }
        //Choose best codon for each amino acid
        String[] codons = seqChooser.run(peptide.toUpperCase());

        
        //Choose an RBS
        StringBuilder cds = new StringBuilder();
        for(String codon : codons) {
            cds.append(codon);
        }

        RBSOption selectedRBS;
        selectedRBS = rbsChooser.run(cds.toString().toUpperCase(), peptide.toUpperCase(), ignores);

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


    }
}
