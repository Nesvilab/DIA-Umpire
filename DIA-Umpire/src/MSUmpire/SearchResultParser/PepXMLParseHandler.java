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
package MSUmpire.SearchResultParser;

import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.ModStringConvert;
import MSUmpire.PSMDataStructure.ModificationInfo;
import MSUmpire.PSMDataStructure.PSM;
import MSUmpire.PSMDataStructure.PTMManager;
import com.compomics.util.experiment.biology.AminoAcid;
import com.compomics.util.experiment.biology.AminoAcidPattern;
import com.compomics.util.experiment.biology.AminoAcidSequence;
import com.compomics.util.experiment.biology.PTM;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.vseravno.solna.SolnaHandler;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import org.apache.log4j.Logger;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepXMLParseHandler implements SolnaHandler<Element> {

    public PepXMLParseHandler(LCMSID singleLCMSID, float StartRT, float EndRT, float threshold,boolean CorrectMassDiff) {
        this.StartRT = StartRT;
        this.EndRT = EndRT;
        this.singleLCMSID = singleLCMSID;
        this.threshold = threshold;
        this.CorrectMassDiff=CorrectMassDiff;
    }
    float StartRT;
    float EndRT;
    float threshold;
    boolean CorrectMassDiff=true;
    LCMSID singleLCMSID;

    @Override
    public void handle(Element node) throws Exception {
        switch (node.getNodeName()) {
            case "spectrum_query": {
                ParseSpectrumNode(node);
                break;
            }
            case "search_summary": {
                ParseSearchSummary(node);
                break;
            }
        }
    }

    private void ParseSearchSummary(Element node) throws IOException {
        if (node.getAttributes().getNamedItem("search_engine") != null) {
            singleLCMSID.SearchEngine = node.getAttributes().getNamedItem("search_engine").getNodeValue();
        }
        if (node.getAttributes().getNamedItem("msDetector") != null) {
            singleLCMSID.msDetector = node.getAttributes().getNamedItem("msDetector").getNodeValue();
        }
        if (node.getAttributes().getNamedItem("msIonization") != null) {
            singleLCMSID.msIonization = node.getAttributes().getNamedItem("msIonization").getNodeValue();
        }
        if (node.getAttributes().getNamedItem("msManufacturer") != null) {
            singleLCMSID.msManufacturer = node.getAttributes().getNamedItem("msManufacturer").getNodeValue();
        }
        if (node.getAttributes().getNamedItem("msMassAnalyzer") != null) {
            singleLCMSID.msMassAnalyzer = node.getAttributes().getNamedItem("msMassAnalyzer").getNodeValue();
        }
        if (node.getAttributes().getNamedItem("msModel") != null) {
            singleLCMSID.msModel = node.getAttributes().getNamedItem("msModel").getNodeValue();
        }
        for (int k = 0; k < node.getChildNodes().getLength(); k++) {
            if ("search_database".equals(node.getChildNodes().item(k).getNodeName())) {
                singleLCMSID.DataBase = node.getChildNodes().item(k).getAttributes().getNamedItem("local_path").getNodeValue();
            }
            if ("aminoacid_modification".equals(node.getChildNodes().item(k).getNodeName())) {
                if (node.getChildNodes().item(k).getAttributes().getNamedItem("aminoacid") != null) {
                    String site = node.getChildNodes().item(k).getAttributes().getNamedItem("aminoacid").getNodeValue();
                    float mass = (float) Math.round(Float.parseFloat(node.getChildNodes().item(k).getAttributes().getNamedItem("mass").getNodeValue()) * 1000) / 1000;
                    float massdiff2 = Float.parseFloat(node.getChildNodes().item(k).getAttributes().getNamedItem("massdiff").getNodeValue());
                    AminoAcid aa = AminoAcid.getAminoAcid(site.charAt(0));
                    float massdiff = mass - (float)aa.getMonoisotopicMass();

                    if (massdiff != 0f && Math.abs(massdiff - massdiff2) < 0.1f) {
                        CheckAndAddModification(site, massdiff);

                    } else {
                        if (Math.abs(massdiff2 + 17.0265) < 0.01 && Math.abs(mass - 143.0041f) < 0.001 && "C".equals(site)) {
                            CheckAndAddModification(site, massdiff);
                        } else {
                            Logger.getRootLogger().warn("Warning! modification in pepxml : amino acid " + site + "(mass: " + mass + ", massdiff:" + massdiff2 + ") ignored.");
                        }
                    }
                }
            }
            if ("terminal_modification".equals(node.getChildNodes().item(k).getNodeName())) {

                if (node.getChildNodes().item(k).getAttributes().getNamedItem("terminus") != null) {
                    String site = "";
                    if ("c".equals(node.getChildNodes().item(k).getAttributes().getNamedItem("terminus").getNodeValue().toLowerCase())) {
                        site = "C-term";
                    }
                    if ("n".equals(node.getChildNodes().item(k).getAttributes().getNamedItem("terminus").getNodeValue().toLowerCase())) {
                        site = "N-term";
                    }
                    float massdiff = Float.parseFloat(node.getChildNodes().item(k).getAttributes().getNamedItem("massdiff").getNodeValue());
                    CheckAndAddModification(site, massdiff);
                }
            }
        }
    }

    private void CheckAndAddModification(String site, float massdiff) throws IOException {
        PTM ptm = PTMManager.GetInstance().GetPTM(site, massdiff);
        if (ptm == null) {
            Logger.getRootLogger().warn("Warning! modification in pepxml : amino acid " + site + "(mass diff:" + massdiff + ") cannot be found in the library.");
            massdiff= (float) (Math.floor(massdiff*10)/10);
            Logger.getRootLogger().warn("Creating a custom modification type called \"" + massdiff + "@" + site + "\"");
            ArrayList<String> list = new ArrayList<>();
            list.add(site);
            ptm = new PTM(PTM.MODAA, massdiff + "@" + site, (double) massdiff, list);
            PTMManager.GetInstance().AddPTM(ptm);
        }
        singleLCMSID.AddModification(ptm, site);
    }

    private void ParseSpectrumNode(Element spectrum) throws IOException {
                
        PSM psm = new PSM();
        psm.SpecNumber = spectrum.getAttributes().getNamedItem("spectrum").getNodeValue();
        psm.ObserPrecursorMass = Float.parseFloat(spectrum.getAttributes().getNamedItem("precursor_neutral_mass").getNodeValue());
        psm.Charge = Integer.parseInt(spectrum.getAttributes().getNamedItem("assumed_charge").getNodeValue());
        psm.ScanNo = Integer.parseInt(spectrum.getAttributes().getNamedItem("start_scan").getNodeValue());
        if (spectrum.getAttributes().getNamedItem("retention_time_sec") != null) {
            psm.RetentionTime = Float.parseFloat(spectrum.getAttributes().getNamedItem("retention_time_sec").getNodeValue()) / 60f;
        }
        psm.NeighborMaxRetentionTime = psm.RetentionTime;

        psm.RawDataName = psm.SpecNumber.substring(0, psm.SpecNumber.indexOf("."));

        for (int k = 0; k < spectrum.getChildNodes().getLength(); k++) {
            Node resultNode = spectrum.getChildNodes().item(k);
            if ("search_result".equals(resultNode.getNodeName())) {
                for (int l = 0; l < resultNode.getChildNodes().getLength(); l++) {
                    Node hitNode = resultNode.getChildNodes().item(l);
                    if ("search_hit".equals(hitNode.getNodeName()) && "1".equals(hitNode.getAttributes().getNamedItem("hit_rank").getNodeValue())) {
                        psm.NeutralPepMass = Float.parseFloat(hitNode.getAttributes().getNamedItem("calc_neutral_pep_mass").getNodeValue());
                        float error = Float.parseFloat(hitNode.getAttributes().getNamedItem("massdiff").getNodeValue().replace("+-", ""));
                        float error0 = Math.abs(error);
                        float error1 = Math.abs(error - 1f);
                        float error2 = Math.abs(error - 2f);
                        float error3 = Math.abs(error - 3f);
                        psm.MassError = error;

                        if (CorrectMassDiff && Math.min(Math.min(Math.min(error0, error1), error2), error3) < 0.5f) {
                            if (error1 < error0 && error1 < error2 && error1 < error3) {
                                psm.MassError = error - 1f;
                                psm.ObserPrecursorMass -= 1f;
                            } else if (error2 < error0 && error2 < error1 && error2 < error3) {
                                psm.MassError = error - 2f;
                                psm.ObserPrecursorMass -= 2f;
                            } else if (error3 < error0 && error3 < error1 && error3 < error2) {
                                psm.MassError = error - 3f;
                                psm.ObserPrecursorMass -= 3f;
                            }
                        }

                        if (hitNode.getAttributes().getNamedItem("peptide_prev_aa") != null) {
                            psm.PreAA = hitNode.getAttributes().getNamedItem("peptide_prev_aa").getNodeValue();
                        }
                        if (hitNode.getAttributes().getNamedItem("peptide_next_aa") != null) {
                            psm.NextAA = hitNode.getAttributes().getNamedItem("peptide_next_aa").getNodeValue();
                        }
                        if (hitNode.getAttributes().getNamedItem("num_missed_cleavages") != null) {
                            psm.MissedCleavage = Integer.parseInt(hitNode.getAttributes().getNamedItem("num_missed_cleavages").getNodeValue());
                        }
                        psm.Sequence = hitNode.getAttributes().getNamedItem("peptide").getNodeValue();
                        psm.ModSeq = psm.Sequence;
                        psm.TPPModSeq = psm.Sequence;

                        String ProtACC = hitNode.getAttributes().getNamedItem("protein").getNodeValue();
                        if (!"".equals(ProtACC)) {
                            psm.AddParentProtein(ProtACC);
                        }
                        String altproACC = "";
                        boolean iprophet = false;
                        for (int m = 0; m < hitNode.getChildNodes().getLength(); m++) {
                            Node hitModNode = hitNode.getChildNodes().item(m);

                            switch (hitModNode.getNodeName()) {
                                case ("modification_info"): {
                                    GetModificationInfo(psm, hitModNode);
                                    break;
                                }
                                case ("analysis_result"): {
                                    switch (hitModNode.getAttributes().getNamedItem("analysis").getNodeValue()) {
                                        case "peptideprophet": {
                                            if (!iprophet && hitModNode.getChildNodes().item(1).getAttributes().getNamedItem("probability") != null) {
                                                psm.Probability = Float.parseFloat(hitModNode.getChildNodes().item(1).getAttributes().getNamedItem("probability").getNodeValue());
                                            }
                                            break;
                                        }
                                        case "interprophet": {
                                            iprophet = true;
                                            if (hitModNode.getChildNodes().item(1).getAttributes().getNamedItem("probability") != null) {
                                                psm.Probability = Float.parseFloat(hitModNode.getChildNodes().item(1).getAttributes().getNamedItem("probability").getNodeValue());
                                            }
                                            break;
                                        }
                                        case "percolator": {
                                            if (hitModNode.getChildNodes().item(1).getAttributes().getNamedItem("probability") != null) {
                                                psm.Probability = Float.parseFloat(hitModNode.getChildNodes().item(1).getAttributes().getNamedItem("probability").getNodeValue());
                                            }
                                            break;
                                        }
                                    }
                                    break;
                                }
                                case ("search_score"): {
                                    switch (hitModNode.getAttributes().item(0).getNodeValue()) {
                                        case ("hyperscore"): {
                                            psm.hyperscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("nextscore"): {
                                            psm.nextscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("bscore"): {
                                            psm.bscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("yscore"): {
                                            psm.yscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("zscore"): {
                                            psm.zscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("ascore"): {
                                            psm.ascore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("xscore"): {
                                            psm.xscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("expect"): {
                                            psm.expect = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("XCorr"): {
                                            psm.hyperscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                    }
                                    break;
                                }
                                case ("alternative_protein"): {
                                    altproACC = hitModNode.getAttributes().getNamedItem("protein").getNodeValue();
                                    if (!"".equals(altproACC)) {
                                        psm.AddParentProtein(altproACC);
                                    }
                                    break;
                                }
                            }
                        }
                        if (psm.Probability > threshold) {
                            singleLCMSID.AddPSM(psm);
                        }
                    }
                }
            }
        }
    }

    private void GetModificationInfo(PSM psmid, Node node) throws IOException {
        String PepSeq = psmid.Sequence;
        String modseq = psmid.Sequence;
        String TPPmodseq = modseq;
        if (node.getAttributes().getNamedItem("modified_peptide") != null) {
            TPPmodseq = node.getAttributes().getNamedItem("modified_peptide").getNodeValue();
        }
        if (node.getAttributes().getNamedItem("mod_nterm_mass") != null) {
            float mass = Float.parseFloat(node.getAttributes().getNamedItem("mod_nterm_mass").getNodeValue());
            String site = "N-term";
            CheckAndAddModification(site, mass);
            
            ModificationInfo matchmod = null;
            float massdiff = Float.MAX_VALUE;
            for (ModificationInfo mod : singleLCMSID.ModificationList.values()) {
                if (mod.site.equals(site)) {
                    float diff = Math.abs(mod.mass - mass);
                    if (diff < massdiff) {
                        massdiff = diff;
                        matchmod = mod;
                    }
                }
            }
            if (matchmod != null) {
                psmid.Modifications.add(new ModificationMatch(matchmod.modification.getName(), true, 1));
                modseq = ModStringConvert.AddModIntoSeqBeforeSite(modseq, matchmod.GetKey(), -1);
            } else {
                Logger.getRootLogger().warn("Modification [" + mass + " @ nterm] for spectrum: " + psmid.SpecNumber + " not found in the library:");
            }
        }
        if (node.getAttributes().getNamedItem("mod_cterm_mass") != null) {
            float mass = Float.parseFloat(node.getAttributes().getNamedItem("mod_cterm_mass").getNodeValue());
            String site = "C-term";
            CheckAndAddModification(site, mass);

            ModificationInfo matchmod = null;
            float massdiff = Float.MAX_VALUE;
            for (ModificationInfo mod : singleLCMSID.ModificationList.values()) {
                if (mod.site.equals(site)) {
                    float diff = Math.abs(mod.mass - mass);
                    if (diff < massdiff) {
                        massdiff = diff;
                        matchmod = mod;
                    }
                }
            }            
            if (matchmod != null) {
                psmid.Modifications.add(new ModificationMatch(matchmod.modification.getName(), true, psmid.Sequence.length()));
                modseq = ModStringConvert.AddModIntoSeqBeforeSite(modseq, matchmod.GetKey(), psmid.Sequence.length() - 1);
            } else {
                Logger.getRootLogger().warn("Modification [" + mass + " @ cterm] for spectrum: " + psmid.SpecNumber + " not found in the library:");
            }
        }

        for (int i = 0; i < node.getChildNodes().getLength(); i++) {
            if ("mod_aminoacid_mass".equals(node.getChildNodes().item(i).getNodeName())) {
                int idx = Integer.parseInt(node.getChildNodes().item(i).getAttributes().getNamedItem("position").getNodeValue());
                String site = String.valueOf(PepSeq.charAt(idx - 1));
                float mass = Float.parseFloat(node.getChildNodes().item(i).getAttributes().getNamedItem("mass").getNodeValue());
                AminoAcid aa = AminoAcid.getAminoAcid(site.charAt(0));
                float massdiff = mass - (float) aa.getMonoisotopicMass();
                CheckAndAddModification(site, massdiff);
                                
                ModificationInfo matchmod = null;
                massdiff = Float.MAX_VALUE;
                for (ModificationInfo mod : singleLCMSID.ModificationList.values()) {
                    if (mod.site.equals(site)) {
                        float diff = Math.abs(mod.mass - mass);
                        if (diff < massdiff) {
                            massdiff = diff;
                            matchmod = mod;
                        }
                    }
                }
                if (matchmod != null) {
                    psmid.Modifications.add(new ModificationMatch(matchmod.modification.getName(), true, idx));
                    modseq = ModStringConvert.AddModIntoSeqBeforeSite(modseq, matchmod.GetKey(), idx - 1);
                } else {
                    Logger.getRootLogger().warn("Modification [" + mass + " @ " + site + "] for spectrum: " + psmid.SpecNumber + " not found in the library:");
                }
            }
        }
        psmid.ModSeq = modseq;
        psmid.TPPModSeq = TPPmodseq;
    }
}
