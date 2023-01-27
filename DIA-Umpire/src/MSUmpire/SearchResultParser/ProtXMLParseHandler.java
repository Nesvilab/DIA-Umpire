/*
 * This file is part of DIA-Umpire.
 *
 * DIA-Umpire is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or (at your option) any later 
 * version.
 *
 * DIA-Umpire is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with DIA-Umpire. If not, see <https://www.gnu.org/licenses/>.
 *
 */
package MSUmpire.SearchResultParser;

import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.ModStringConvert;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PSMDataStructure.ProtID;
import com.vseravno.solna.SolnaHandler;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ProtXMLParseHandler implements SolnaHandler<Element> {

    public ProtXMLParseHandler(LCMSID SingleLCMSID, float threshold) {        
        this.SingleLCMSID=SingleLCMSID;
        this.threshold=threshold;
    }
    float threshold;
    LCMSID SingleLCMSID;

    @Override
    public void handle(Element node) throws Exception {
        switch(node.getNodeName()){
            case "protein_group": {
                ParseProteinGroupNode(node);
                break;
            }
        }
    }
    
    private void ParseProteinGroupNode(Element protgroupnode) throws IOException, ClassNotFoundException, InterruptedException {

        int groupindex = Integer.parseInt(protgroupnode.getAttributes().getNamedItem("group_number").getNodeValue());
        
        float groupprob=0f;
        if (protgroupnode.getAttributes().getNamedItem("probability") != null) {
            groupprob = Float.parseFloat(protgroupnode.getAttributes().getNamedItem("probability").getNodeValue());
        }
        else if (protgroupnode.getAttributes().getNamedItem("probability") != null) {
            groupprob = Float.parseFloat(protgroupnode.getAttributes().getNamedItem("group_probability").getNodeValue());
        }
        
        int alphabet = 1;
        
        for (int j = 0; j < protgroupnode.getChildNodes().getLength(); j++) {
            Node protnode = protgroupnode.getChildNodes().item(j);
            if ("protein".equals(protnode.getNodeName())) {
                ProtID proid = new ProtID();
                proid.setAccNo(protnode.getAttributes().getNamedItem("protein_name").getNodeValue());
                proid.IndisProteins.add(proid.getAccNo());
                proid.UniProtID = proid.getAccNo();
                proid.Probability = Float.parseFloat(protnode.getAttributes().getNamedItem("probability").getNodeValue());
                proid.GroupProb=groupprob;              
                
                if (proid.Probability < threshold) {
                    continue;
                }
                proid.ProteinGroup = String.valueOf(groupindex) + "_" + (alphabet++);

                for (int k = 0; k < protnode.getChildNodes().getLength(); k++) {
                    Node child = protnode.getChildNodes().item(k);

                    switch (child.getNodeName()) {
                        case ("indistinguishable_protein"): {
                            proid.IndisProteins.add(child.getAttributes().getNamedItem("protein_name").getNodeValue());
                            proid.IndisProtDes.add(child.getChildNodes().item(1).getAttributes().getNamedItem("protein_description").getNodeValue().replace(",", "_"));
                            break;
                        }
                        case ("annotation"): {
                            proid.SetDescription(child.getAttributes().getNamedItem("protein_description").getNodeValue().replace(",", "_"));
                            proid.IndisProtDes.add(child.getAttributes().getNamedItem("protein_description").getNodeValue().replace(",", "_"));
                            break;
                        }
                        case ("parameter"): {
                            if ("prot_length".equals(child.getAttributes().getNamedItem("name").getNodeValue())) {
                                proid.ProteinLength = Integer.parseInt(child.getAttributes().getNamedItem("value").getNodeValue());
                            }
                            break;
                        }
                        case ("peptide"): {
                            //iProphet
                            if (child.getAttributes().getNamedItem("charge") == null || "0".equals(child.getAttributes().getNamedItem("charge").getNodeValue())) {
                                ArrayList<String> prots=new ArrayList<>();
                                for (int l = 0; l < child.getChildNodes().getLength(); l++) {
                                    Node chidpep = child.getChildNodes().item(l);                                    
                                    if ("peptide_parent_protein".equals(chidpep.getNodeName())) {
                                        prots.add(chidpep.getAttributes().getNamedItem("protein_name").getNodeValue());                                              
                                    }
                                    if ("indistinguishable_peptide".equals(chidpep.getNodeName())) {
                                        PepIonID pepIonID = new PepIonID();
                                        pepIonID.Sequence = child.getAttributes().getNamedItem("peptide_sequence").getNodeValue();
                                        pepIonID.MaxProbability = Float.parseFloat(child.getAttributes().getNamedItem("initial_probability").getNodeValue());
                                        pepIonID.Is_NonDegenerate = "Y".equals(child.getAttributes().getNamedItem("is_nondegenerate_evidence").getNodeValue());
                                        pepIonID.Weight = Float.parseFloat(child.getAttributes().getNamedItem("weight").getNodeValue());
                                        pepIonID.GroupWeight = Float.parseFloat(child.getAttributes().getNamedItem("group_weight").getNodeValue());
                                        pepIonID.Charge = Integer.parseInt(chidpep.getAttributes().getNamedItem("charge").getNodeValue());
                                        pepIonID.ParentProtString_ProtXML=prots;
                                        String modseq = pepIonID.Sequence;
                                        for (int n = 0; n < chidpep.getChildNodes().getLength(); n++) {
                                            Node mod = chidpep.getChildNodes().item(n);
                                            if ("modification_info".equals(mod.getNodeName())) {
                                                if (mod.getAttributes().getNamedItem("mod_nterm_mass") != null) {
                                                    float n_termmass = Float.parseFloat(mod.getAttributes().getNamedItem("mod_nterm_mass").getNodeValue());
                                                    String n_term = new DecimalFormat("#.###").format(n_termmass);
                                                    modseq = ModStringConvert.AddModIntoSeqBeforeSite(modseq, "[" + n_term + "(N-term)]", -1);
                                                }
                                                if (mod.getAttributes().getNamedItem("mod_cterm_mass") != null) {
                                                    float c_termmass = Float.parseFloat(mod.getAttributes().getNamedItem("mod_cterm_mass").getNodeValue());
                                                    String c_term = new DecimalFormat("#.###").format(c_termmass);
                                                    modseq = ModStringConvert.AddModIntoSeqBeforeSite(modseq, "[" + c_term + "(C-term)]", pepIonID.Sequence.length() - 1);
                                                }

                                                for (int m = 0; m < mod.getChildNodes().getLength(); m++) {
                                                    Node aa = mod.getChildNodes().item(m);
                                                    if ("mod_aminoacid_mass".equals(aa.getNodeName())) {
                                                        float tmpmass = Float.parseFloat(aa.getAttributes().getNamedItem("mass").getNodeValue());
                                                        String tmp = new DecimalFormat("#.###").format(tmpmass);
                                                        int idx = Integer.parseInt(aa.getAttributes().getNamedItem("position").getNodeValue());
                                                        modseq = ModStringConvert.AddModIntoSeqBeforeSite(modseq, "[" + tmp + "(" + pepIonID.Sequence.charAt(idx - 1) + ")]", idx - 1);
                                                    }
                                                }
                                            }
                                        }
                                        pepIonID.ModSequence = modseq;
                                        if (proid.MaxIniProb < pepIonID.MaxProbability) {
                                            proid.MaxIniProb = pepIonID.MaxProbability;
                                        }
                                        SingleLCMSID.ProtXMLPepIonList.put(pepIonID.GetKey(), pepIonID);
                                        proid.ProtPeptideID.put(pepIonID.GetKey(), pepIonID);
                                        if (!proid.ProtPepSeq.contains(pepIonID.Sequence)) {
                                            proid.ProtPepSeq.add(pepIonID.Sequence);
                                        }
                                        break;
                                    }
                                }
                            } else {                               
                                PepIonID pepIonID = new PepIonID();
                                pepIonID.Sequence = child.getAttributes().getNamedItem("peptide_sequence").getNodeValue();
                                
                                if (child.getAttributes().getNamedItem("initial_probability") != null) {
                                    pepIonID.MaxProbability = Float.parseFloat(child.getAttributes().getNamedItem("initial_probability").getNodeValue());
                                }
                                if (child.getAttributes().getNamedItem("is_nondegenerate_evidence") != null) {
                                    pepIonID.Is_NonDegenerate = "Y".equals(child.getAttributes().getNamedItem("is_nondegenerate_evidence").getNodeValue());
                                }
                                if (child.getAttributes().getNamedItem("weight") != null) {
                                    pepIonID.Weight = Float.parseFloat(child.getAttributes().getNamedItem("weight").getNodeValue());
                                }
                                if (child.getAttributes().getNamedItem("group_weight") != null) {
                                    pepIonID.GroupWeight = Float.parseFloat(child.getAttributes().getNamedItem("group_weight").getNodeValue());
                                }
                                pepIonID.Charge = Integer.parseInt(child.getAttributes().getNamedItem("charge").getNodeValue());

                                //pepIonID.NeutralPepMass = Float.parseFloat(child.getAttributes().getNamedItem("calc_neutral_pep_mass").getNodeValue());
                                String modseq = pepIonID.Sequence;

                                for (int l = 0; l < child.getChildNodes().getLength(); l++) {
                                    Node mod = child.getChildNodes().item(l);

                                    if ("modification_info".equals(mod.getNodeName())) {
                                        if (mod.getAttributes().getNamedItem("mod_nterm_mass") != null) {
                                            float n_termmass = Float.parseFloat(mod.getAttributes().getNamedItem("mod_nterm_mass").getNodeValue());
                                            String n_term = new DecimalFormat("#.###").format(n_termmass);
                                            modseq = ModStringConvert.AddModIntoSeqBeforeSite(modseq, "[" + n_term + "(N-term)]", -1);
                                        }
                                        if (mod.getAttributes().getNamedItem("mod_cterm_mass") != null) {
                                            float c_termmass = Float.parseFloat(mod.getAttributes().getNamedItem("mod_cterm_mass").getNodeValue());
                                            String c_term = new DecimalFormat("#.###").format(c_termmass);
                                            modseq = ModStringConvert.AddModIntoSeqBeforeSite(modseq, "[" + c_term + "(C-term)]", pepIonID.Sequence.length() - 1);
                                        }

                                        for (int m = 0; m < mod.getChildNodes().getLength(); m++) {
                                            Node aa = mod.getChildNodes().item(m);
                                            if ("mod_aminoacid_mass".equals(aa.getNodeName())) {
                                                float tmpmass = Float.parseFloat(aa.getAttributes().getNamedItem("mass").getNodeValue());
                                                String tmp = new DecimalFormat("#.###").format(tmpmass);
                                                int idx = Integer.parseInt(aa.getAttributes().getNamedItem("position").getNodeValue());

                                                modseq = ModStringConvert.AddModIntoSeqBeforeSite(modseq, "[" + tmp + "(" + pepIonID.Sequence.charAt(idx - 1) + ")]", idx - 1);
                                            }
                                        }
                                    }
                                }
                                pepIonID.ModSequence = modseq;
                                if (proid.MaxIniProb < pepIonID.MaxProbability) {
                                    proid.MaxIniProb = pepIonID.MaxProbability;
                                }
                                SingleLCMSID.ProtXMLPepIonList.put(pepIonID.GetKey(), pepIonID);
                                proid.ProtPeptideID.put(pepIonID.GetKey(), pepIonID);
                                if (!proid.ProtPepSeq.contains(pepIonID.Sequence)) {
                                    proid.ProtPepSeq.add(pepIonID.Sequence);
                                }
                                break;
                            }
                        }
                    }
                }
                if (proid.ProtPeptideID.size() > 0) {
                    SingleLCMSID.AddProtID(proid);
                    if (!SingleLCMSID.ProteinGroups.containsKey(groupindex)) {
                        SingleLCMSID.ProteinGroups.put(groupindex, new ArrayList<ProtID>());
                    }
                    SingleLCMSID.ProteinGroups.get(groupindex).add(proid);
                }
            }
        }        
    }
    
}
