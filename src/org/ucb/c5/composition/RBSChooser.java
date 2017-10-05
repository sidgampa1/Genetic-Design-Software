package org.ucb.c5.composition;

import java.util.HashSet;
import java.util.Set;
import java.util.HashMap;
import java.util.ArrayList;
import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.utils.FileUtils;
import org.ucb.c5.sequtils.CalcEditDistance;
import org.ucb.c5.sequtils.HairpinCounter;
import org.ucb.c5.sequtils.RevComp;
import org.ucb.c5.sequtils.Translate;

/**
 * Second generation RBSChooser algorithm
 *
 * Employs a list of genes and their associated ribosome binding sites for
 * highly-expressed proteins in E. coli.
 *
 * @author J. Christopher Anderson
 */
public class RBSChooser {

    private ArrayList<RBSOption> rbss;
    private CalcEditDistance six_aa_scorer;
    private HairpinCounter second_struct_scorer;
    private Translate translator;
    private RevComp reverse;
    private AminoAcidToCodon translate;

    public void initiate() throws Exception {
        six_aa_scorer = new CalcEditDistance();
        second_struct_scorer = new HairpinCounter();
        translator = new Translate();
        reverse = new RevComp();
        translate = new AminoAcidToCodon();

        translate.initiate();
        six_aa_scorer.initiate();
        second_struct_scorer.initiate();
        translator.initiate();
        reverse.initiate();


        rbss = new ArrayList<>();


        //Read the data file
        String gene_data = FileUtils.readResourceFile("composition/data/coli_genes.txt");
        String rbs_data = FileUtils.readResourceFile("composition/data/rbs_options.txt");

        //split data into array of lines
        String[] gene_lines = gene_data.split("\\r|\\r?\\n");
        String[] rbs_lines = rbs_data.split("\\r|\\r?\\n");

        //check line by line for match between rbs options native genes and ecoli genes list,
        //create rbsoption to add to rbss
        for (int i = 0; i < rbs_lines.length; i++) {
            String rbs_line = rbs_lines[i]; //single line in rbs_options
            String[] rbs_fields = rbs_line.split("\t");

            for (int j = 0; j < gene_lines.length; j++) {
                String gene_line = gene_lines[j]; //single line in coli_genes
                String[] gene_fields = gene_line.split("\t");

                String src_name = gene_fields[1];
                // match between name of rbs native gene and gene in coli gene list
                if (src_name.equals(rbs_fields[0])) {
                    // define rbsoption attributes
                    String src_desc = gene_fields[0];
                    String src_rbs = rbs_fields[1];
                    String src_cds = gene_fields[6];

                    //translate first 6 aa's from cds
                    String src_first6aas = translator.run(src_cds.substring(0, 18));

                    //add rbsoption to rbss (list of rbs option and corresponding native gene cds)
                    RBSOption rbs_match = new RBSOption(src_name, src_desc, src_rbs, src_cds, src_first6aas);
                    rbss.add(rbs_match);
                }

            }
        }
    }

    /**
     * this private method take a peptide and reverse translates it using the hard coded codon table
     * @param peptide
     * @return cds - reverse translated dna seq
     */

//    private String rev_translate(String peptide) {
//        String cds = "";
//
//        // add codons to cds
//        for(int i=0; i<peptide.length(); i++) {
//            char aa = peptide.charAt(i);
//            String codon = translate.table.get(aa);
//            cds += codon;
//        }
//
//        return cds;
//    }

    /**
     * Provided a protein of sequence 'peptide', this computes the best ribosome
     * binding site to use from a list of options.
     *
     * It also permits a list of options to exclude.
     *
     * @param cds The protein sequence, ie MSKGEE...
     * @param ignores The list of RBS's to exclude
     * @return
     * @throws Exception
     */
    public RBSOption run(String cds, String peptide, Set<RBSOption> ignores) throws Exception {
        //TODO:  Fill in
        double min_score = Integer.MAX_VALUE; //keeps track of lowest score in rbs options
        RBSOption best_rbs = rbss.get(0); //rbsoption corresponding to lowest score, initialized to first rbs in list

        int aa_score;
        double second_struct_score;
        double total_score;

        // iterate and score through rbsoption list to find best rbs
        for (RBSOption rbs_option : rbss) {
            if (!ignores.contains(rbs_option)) {
                //score the similarity between first 6 aas
                String rbs_6_aas = rbs_option.getFirst6aas();
                String pep_6_aas = peptide.substring(0, 6);
                aa_score = six_aa_scorer.run(rbs_6_aas, pep_6_aas);

                //score the secondary structure formation between rbs and cds (only hairpins)
                //String cds = rev_translate(cds);
                String combined_seq = rbs_option.getRbs() + cds;
                second_struct_score = second_struct_scorer.run(combined_seq);

                //define hairpin as #hbonds per nuc, this will acct for diff cds+rbs lengths,
                //and make score similar to aa score (<10)
                second_struct_score = second_struct_score / combined_seq.length();

                total_score = (aa_score) + (second_struct_score);
                // account for a perfect initial 6 aa match possibly by lowering total score by 1?
                // there probably is a better way, but scores are low enough for this to work well
                if (aa_score == 0) {
                    total_score -= 1;
                }

                if (total_score < min_score) {
                    min_score = total_score;
                    best_rbs = rbs_option;
                }

            }
        }
        return best_rbs;
    }


    public static void main(String[] args) throws Exception {
        //Create an example
//        String peptide = "MLSDTIDTKQQQQQLHVLFIDSYDSFTYNVVRLIEQQTDISPGVNAVHVTTVHSDTFQSMDQLLPLLPLFDAIVVGPGPGNPNNGAQDMGIISELFENANGKLDEVPILGICLGFQAMCLAQGADVSELNTIKHGQVYEMHLNDAARACGLFSGYPDTFKSTRYHSLHVNAEGIDTLLPLCTTEDENGILLMSAQTKNKPWFGVQYHPESCCSELGGLLVSNFLKLSFINNVKTGRWEKKKLNGEFSDILSRLDRTIDRDPIYKVKEKYPKGEDTTYVKQFEVSEDPKLTFEICNIIREEKFVMSSSVISENTGEWSIIALPNSASQVFTHYGAMKKTTVHYWQDSEISYTLLKKCLDGQDSDLPGSLEVIHEDKSQFWITLGKFMENKIIDNHREIPFIGGLVGILGYEIGQYIACGRCNDDENSLVPDAKLVFINNSIVINHKQGKLYCISLDNTFPVALEQSLRDSFVRKKNIKQSLSWPKYLPEEIDFIITMPDKLDYAKAFKKCQDYMHKGDSYEMCLTTQTKVVPSAVIEPWRIFQTLVQRNPAPFSSFFEFKDIIPRQDETPPVLCFLSTSPERFLKWDADTCELRPIKGTVKKGPQMNLAKATRILKTPKEFGENLMILDLIRNDLYELVPDVRVEEFMSVQEYATVYQLVSVVKAHGLTSASKKTRYSGIDVLKHSLPPGSMTGAPKKITVQLLQDKIESKLNKHVNGGARGVYSGVTGYWSVNSNGDWSVNIRCMYSYNGGTSWQLGAGGAITVLSTLDGELEEMYNKLESNLQIFM";
//        String peptide = "MSQTVHFQGNPVTVANSIPQ";
//
//        //Initiate the chooser
//        RBSChooser chooser = new RBSChooser();
//        chooser.initiate();
//
//        //Make the first choice with an empty Set of ignores
//        Set<RBSOption> ignores = new HashSet<>();
//        RBSOption selected1 = chooser.run(peptide, ignores);
//
//        //Add the first selection to the list of things to ignore
//        ignores.add(selected1);
//
//        //Choose again with an ignore added
//        RBSOption selected2 = chooser.run(peptide, ignores);
//
//        //Print out the two options, which should be different
//        System.out.println("Protein starts with:");
//        System.out.println(peptide.substring(0, 6));
//        System.out.println();
//        System.out.println("Selected1:\n");
//        System.out.println(selected1.toString());
//        System.out.println();
//        System.out.println("Selected2:\n");
//        System.out.println(selected2.toString());
    }
}