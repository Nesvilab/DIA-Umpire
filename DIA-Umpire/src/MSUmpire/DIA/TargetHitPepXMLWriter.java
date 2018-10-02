/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package MSUmpire.DIA;

import MSUmpire.BaseDataStructure.UmpireInfo;
import MSUmpire.SeqUtility.FastaParser;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.Utility.DateTimeTag;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

/**
 * pepXML writer for targeted re-extraction
 * @author Chih-Chiang Tsou
 */
public class TargetHitPepXMLWriter {

    TargetMatchScoring Tscoring;
    StringBuilder sb = new StringBuilder();
    String Filename;
    String Fasta;
    String Decoytag;
    
    public class ParentProtein {

        String Acc;
        char PreAA = '-';
        char NextAA = '-';
        String Des;
    }

    public TargetHitPepXMLWriter(String Filename, String Fasta, String Decoytab, TargetMatchScoring Tscoring) throws IOException {
        this.Filename = Filename;
        this.Fasta = Fasta;
        this.Decoytag=Decoytab;
        this.Tscoring = Tscoring;
        write();
    }

    public void write() throws IOException {

        Logger.getRootLogger().info("Writing "+Filename);
        int minlength = Integer.MAX_VALUE;
        int maxlength = 0;
        int maxmiss = 0;

        for (UmpireSpecLibMatch match : Tscoring.libTargetMatches) {
            PepIonID pepIonID = match.pepIonID;
            if (pepIonID.Sequence.length() > maxlength) {
                maxlength = pepIonID.Sequence.length();
            }
            if (pepIonID.Sequence.length() < minlength) {
                minlength = pepIonID.Sequence.length();
            }
            if (pepIonID.getNMissedCleavages() > maxmiss) {
                maxmiss = pepIonID.getNMissedCleavages();
            }
        }

        FastaParser fastaparser = new FastaParser(Fasta);
        fastaparser.RemoveDecoy(Decoytag);
        //FastaParser_V2 fastaparser = FastaParser.FasterSerialzationRead(Fasta);        
        fastaparser.digestion(maxmiss, minlength, maxlength,Decoytag);

        Header();
        sb.append("<msms_run_summary base_name=\"" + FilenameUtils.getBaseName(Filename) + "\" search_engine=\"interprophet\" >\n");
        SearchSummary();
        int index = 0;
        for (UmpireSpecLibMatch match : Tscoring.libTargetMatches) {            
            //<editor-fold defaultstate="collapsed" desc="Forward hit">            
            PeakGroupScore bestscore = match.BestHit;
            String Sequence=match.pepIonID.Sequence;
            if (bestscore!=null) {
                ArrayList<ParentProtein> parentprots = new ArrayList<>();
                if (!fastaparser.PeptideList.containsKey(Sequence)) {
                    //Logger.getRootLogger().warn(Sequence + " is not found as tryptic peptides in " + Fasta + ". Searching in protein sequences..");
                    for (String acc : fastaparser.ProteinList.keySet()) {
                        String ProtSeq = fastaparser.ProteinList.get(acc).Seq;
                        if (ProtSeq.contains(Sequence)) {
                            ParentProtein prot = new ParentProtein();
                            prot.Acc = acc;
                            prot.Des = fastaparser.ProteinList.get(acc).Des;
                            int aaindex = ProtSeq.indexOf(Sequence) - 1;
                            if (aaindex >= 0) {
                                prot.PreAA = ProtSeq.charAt(aaindex);
                            }
                            aaindex = ProtSeq.indexOf(Sequence) + Sequence.length();
                            if (aaindex < ProtSeq.length()) {
                                prot.NextAA = ProtSeq.charAt(aaindex);
                            }
                            parentprots.add(prot);
                        }
                    }
                } else {
                    for (String acc : fastaparser.PeptideList.get(Sequence).Proteins) {
                        String ProtSeq = fastaparser.ProteinList.get(acc).Seq;
                        if (ProtSeq.contains(Sequence)) {
                            ParentProtein prot = new ParentProtein();
                            prot.Acc = acc;
                            prot.Des = fastaparser.ProteinList.get(acc).Des;
                            int aaindex = ProtSeq.indexOf(Sequence) - 1;
                            if (aaindex >= 0) {
                                prot.PreAA = ProtSeq.charAt(aaindex);
                            }
                            aaindex = ProtSeq.indexOf(Sequence) + Sequence.length();
                            if (aaindex < ProtSeq.length()) {
                                prot.NextAA = ProtSeq.charAt(aaindex);
                            }
                            parentprots.add(prot);
                        }
                    }
                }

                if (parentprots.isEmpty()) {
                    Logger.getRootLogger().warn(Sequence + " is not found in " + Fasta);
                } else {
                    sb.append("<spectrum_query spectrum=\"" +FilenameUtils.getBaseName(Filename)+"."+ index+"."+ index+"."+match.pepIonID.Charge+ "\" start_scan=\""+index+"\" end_scan=\""+index+"\" precursor_neutral_mass=\"" + bestscore.PrecursorNeutralMass+ "\" assumed_charge=\"" + match.pepIonID.Charge + "\" index=\"" + (index++) + "\" retention_time_sec=\"" +bestscore.PrecursorRT* 60f + "\">\n"
                            + "<search_result>\n"
                            + "<search_hit hit_rank=\"1\" peptide=\"" + Sequence + "\" peptide_prev_aa=\"" + parentprots.get(0).PreAA + "\" peptide_next_aa=\"" + parentprots.get(0).NextAA + "\" protein=\"" + parentprots.get(0).Acc + "\" protein_descr=\"" + parentprots.get(0).Des + "\" num_tot_proteins=\"" + parentprots.size() + "\" num_matched_ions=\"" + (bestscore.NoMatchB+bestscore.NoMatchY) + "\" tot_num_ions=\"" + 2 * (Sequence.length() - 1) + "\" calc_neutral_pep_mass=\"" + match.pepIonID.CalcNeutralPepMass() + "\" massdiff=\"" + (match.pepIonID.CalcNeutralPepMass() - bestscore.PrecursorNeutralMass) + "\"  num_tol_term=\"2\" num_missed_cleavages=\"" + match.pepIonID.getNMissedCleavages() + "\" is_rejected=\"0\">\n");

                    for (int i = 1; i < parentprots.size(); i++) {
                        sb.append("<alternative_protein protein=\"" + parentprots.get(i).Acc + "\" protein_descr=\"" + parentprots.get(i).Des + "\" num_tol_term=\"2\" peptide_prev_aa=\"" + parentprots.get(i).PreAA + "\" peptide_next_aa=\"" + parentprots.get(i).NextAA + "\"/>");
                    }
                    //                        + "<modification_info modified_peptide=\"KITIADCGQLE\">\n"
                    //                        + "<mod_aminoacid_mass position=\"7\" mass=\"160.0306\"/>\n"
                    //                        + "</modification_info>\n"
                    //                        + "<analysis_result analysis=\"peptideprophet\">\n"
                    //                        + "<peptideprophet_result probability=\"0.8297\" all_ntt_prob=\"(0.0059,0.3657,0.8297)\">\n"
                    //                        + "<search_score_summary>\n"
                    //                        + "<parameter name=\"fval\" value=\"1.6872\"/>\n"
                    //                        + "<parameter name=\"ntt\" value=\"2\"/>\n"
                    //                        + "<parameter name=\"nmc\" value=\"0\"/>\n"
                    //                        + "<parameter name=\"massd\" value=\"4.011\"/>\n"
                    //                        + "<parameter name=\"isomassd\" value=\"0\"/>\n"
                    //                        + "</search_score_summary>\n"
                    //                        + "</peptideprophet_result>\n"
                    //                        + "</analysis_result>\n"
                    sb.append("<analysis_result analysis=\"interprophet\">\n"
                            + "<interprophet_result probability=\"" + bestscore.MixtureModelLocalProb + "\" all_ntt_prob=\"("+ bestscore.MixtureModelLocalProb +","+ bestscore.MixtureModelLocalProb +","+ bestscore.MixtureModelLocalProb +")\">\n"
                            //                        + "<search_score_summary>\n"
                            //                        + "<parameter name=\"nrs\" value=\"0\"/>\n"
                            //                        + "<parameter name=\"nsi\" value=\"0\"/>\n"
                            //                        + "<parameter name=\"nsm\" value=\"0\"/>\n"
                            //                        + "</search_score_summary>\n"
                            + "</interprophet_result>\n"
                            + "</analysis_result>\n"
                            + "</search_hit>\n"
                            + "</search_result>\n"
                            + "</spectrum_query>\n");
                }
            }
                //</editor-fold>

            //<editor-fold defaultstate="collapsed" desc="Decoy hit">
            bestscore = match.BestDecoyHit;
            if(!fastaparser.PeptideList.containsKey(Sequence)){
                continue;
            }
            String DecoySeq=fastaparser.PeptideList.get(Sequence).Decoy;
             if (bestscore!=null) {
                ArrayList<ParentProtein> parentprots = new ArrayList<>();
                if (!fastaparser.PeptideList.containsKey(Sequence)) {
                    //Logger.getRootLogger().warn(Sequence + " is not found as tryptic peptides in " + Fasta + ". Searching in protein sequences..");
                    for (String acc : fastaparser.ProteinList.keySet()) {
                        String ProtSeq = fastaparser.ProteinList.get(acc).Seq;
                        if (ProtSeq.contains(Sequence)) {
                            ParentProtein prot = new ParentProtein();
                            prot.Acc = acc;
                            prot.Des = fastaparser.ProteinList.get(acc).Des;
                            int aaindex = ProtSeq.indexOf(Sequence) - 1;
                            if (aaindex >= 0) {
                                prot.PreAA = ProtSeq.charAt(aaindex);
                            }
                            aaindex = ProtSeq.indexOf(Sequence) + Sequence.length();
                            if (aaindex < ProtSeq.length()) {
                                prot.NextAA = ProtSeq.charAt(aaindex);
                            }
                            parentprots.add(prot);
                        }
                    }
                } else {
                    for (String acc : fastaparser.PeptideList.get(Sequence).Proteins) {
                        String ProtSeq = fastaparser.ProteinList.get(acc).Seq;
                        if (ProtSeq.contains(Sequence)) {
                            ParentProtein prot = new ParentProtein();
                            prot.Acc = acc;
                            prot.Des = fastaparser.ProteinList.get(acc).Des;
                            int aaindex = ProtSeq.indexOf(Sequence) - 1;
                            if (aaindex >= 0) {
                                prot.PreAA = ProtSeq.charAt(aaindex);
                            }
                            aaindex = ProtSeq.indexOf(Sequence) + Sequence.length();
                            if (aaindex < ProtSeq.length()) {
                                prot.NextAA = ProtSeq.charAt(aaindex);
                            }
                            parentprots.add(prot);
                        }
                    }
                }

                if (parentprots.isEmpty()) {
                    Logger.getRootLogger().warn(Sequence + " is not found in " + Fasta);
                } else {
                    sb.append("<spectrum_query spectrum=\""+Decoytag+"_" + Sequence+ "\" precursor_neutral_mass=\"" + bestscore.PrecursorNeutralMass + "\" assumed_charge=\"" + match.pepIonID.Charge + "\" index=\"" + (index++) + "\" retention_time_sec=\"" +bestscore.PrecursorRT* 60f + "\">\n"
                            + "<search_result>\n"
                            + "<search_hit hit_rank=\"1\" peptide=\"" + DecoySeq + "\" peptide_prev_aa=\"" + parentprots.get(0).PreAA + "\" peptide_next_aa=\"" + parentprots.get(0).NextAA + "\" protein=\""+Decoytag+"_" + parentprots.get(0).Acc + "\" protein_descr=\""+Decoytag+"_" + parentprots.get(0).Des + "\" num_tot_proteins=\"" + parentprots.size() + "\" num_matched_ions=\"" + (bestscore.NoMatchB+bestscore.NoMatchY) + "\" tot_num_ions=\"" + 2 * (Sequence.length() - 1) + "\" calc_neutral_pep_mass=\"" + match.pepIonID.CalcNeutralPepMass() + "\" massdiff=\"" + (match.pepIonID.CalcNeutralPepMass() - bestscore.PrecursorNeutralMass) + "\" num_missed_cleavages=\"" + match.pepIonID.getNMissedCleavages() + "\" is_rejected=\"0\">\n");

                    for (int i = 1; i < parentprots.size(); i++) {
                        sb.append("<alternative_protein protein=\""+Decoytag+"_" + parentprots.get(i).Acc + "\" protein_descr=\""+Decoytag+"_" + parentprots.get(i).Des + "\" num_tol_term=\"2\" peptide_prev_aa=\"" + parentprots.get(i).PreAA + "\" peptide_next_aa=\"" + parentprots.get(i).NextAA + "\"/>");
                    }
                    //                        + "<modification_info modified_peptide=\"KITIADCGQLE\">\n"
                    //                        + "<mod_aminoacid_mass position=\"7\" mass=\"160.0306\"/>\n"
                    //                        + "</modification_info>\n"
                    //                        + "<analysis_result analysis=\"peptideprophet\">\n"
                    //                        + "<peptideprophet_result probability=\"0.8297\" all_ntt_prob=\"(0.0059,0.3657,0.8297)\">\n"
                    //                        + "<search_score_summary>\n"
                    //                        + "<parameter name=\"fval\" value=\"1.6872\"/>\n"
                    //                        + "<parameter name=\"ntt\" value=\"2\"/>\n"
                    //                        + "<parameter name=\"nmc\" value=\"0\"/>\n"
                    //                        + "<parameter name=\"massd\" value=\"4.011\"/>\n"
                    //                        + "<parameter name=\"isomassd\" value=\"0\"/>\n"
                    //                        + "</search_score_summary>\n"
                    //                        + "</peptideprophet_result>\n"
                    //                        + "</analysis_result>\n"
                    sb.append("<analysis_result analysis=\"interprophet\">\n"
                            + "<interprophet_result probability=\"" + bestscore.MixtureModelLocalProb + "\" all_ntt_prob=\"("+ bestscore.MixtureModelLocalProb +","+ bestscore.MixtureModelLocalProb +","+ bestscore.MixtureModelLocalProb +")\">\n"
                            //                        + "<search_score_summary>\n"
                            //                        + "<parameter name=\"nrs\" value=\"0\"/>\n"
                            //                        + "<parameter name=\"nsi\" value=\"0\"/>\n"
                            //                        + "<parameter name=\"nsm\" value=\"0\"/>\n"
                            //                        + "</search_score_summary>\n"
                            + "</interprophet_result>\n"
                            + "</analysis_result>\n"
                            + "</search_hit>\n"
                            + "</search_result>\n"
                            + "</spectrum_query>\n");
                }
            }
            
            //</editor-fold>
        }
        sb.append("</msms_run_summary>\n");
        sb.append("</msms_pipeline_analysis>\n");
        FileWriter writer = new FileWriter(Filename);
        writer.write(sb.toString());
        writer.close();
    }

    private void SearchSummary() {
        sb.append("<sample_enzyme name=\"trypsin\">\n"
                + "<specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>\n"
                + "</sample_enzyme>\n"
                + "<search_summary base_name=\"" + FilenameUtils.getBaseName(Filename) + "\" search_engine=\"Umpire\" precursor_mass_type=\"monoisotopic\" fragment_mass_type=\"monoisotopic\" search_id=\"1\">\n"
                + "<search_database local_path=\"" + Fasta + "\" type=\"AA\"/>\n"
                + "<enzymatic_search_constraint enzyme=\"trypsin\" max_num_internal_cleavages=\"1\" min_number_termini=\"2\"/>\n"
                + "    </search_summary>\n");
    }

    private void Header() {
        sb.append("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
        sb.append("<msms_pipeline_analysis date=\"" + DateTimeTag.GetTag() + "\" summary_xml=\"" + Filename + "\" xmlns=\"http://regis-web.systemsbiology.net/pepXML\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v18.xsd\">\n");
        sb.append("<analysis_summary analysis=\"Umpire\" version=\"" + UmpireInfo.GetInstance().Version + "\" time=\"" + DateTimeTag.GetTag() + "\"/>\n");
    }
}
