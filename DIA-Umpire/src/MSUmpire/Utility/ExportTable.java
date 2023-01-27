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
package MSUmpire.Utility;

import MSUmpire.PSMDataStructure.FragmentPeak;
import MSUmpire.PSMDataStructure.FragmentSelection;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PSMDataStructure.ProtID;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import org.apache.commons.io.FilenameUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ExportTable {

    HashMap<String, HashMap<String, FragmentPeak>> IDSummaryFragments = new HashMap<>();
    HashMap<String, Float> ProteinFragMap = new HashMap<>();
    HashSet<String> IdentifiedPepMap = new HashSet<>();
    ArrayList<LCMSID> FileList;
    LCMSID CombineProtID;
    String WorkFolder;
    FragmentSelection fragselection;

    public ExportTable(String folder, ArrayList<LCMSID> FileList, HashMap<String, HashMap<String, FragmentPeak>> IDSummaryFragments, LCMSID CombineProtID, FragmentSelection fragselection) {
        this.WorkFolder = folder;
        this.FileList = FileList;
        this.CombineProtID = CombineProtID;
        this.IDSummaryFragments = IDSummaryFragments;
        this.fragselection = fragselection;
    }

    public void Export(int TopNPep, int TopNFrag, float Freq) throws IOException {
        LogManager.getRootLogger().info("Writing reports...............");

        if (CombineProtID != null) {
            ProteinLevelExport(TopNPep, TopNFrag, Freq);
        } else {
            PeptideLevelExport(TopNPep, TopNFrag, Freq);
        }
        LogManager.getRootLogger().info("Job done");
        LogManager.getRootLogger().info("=================================================================================================");
    }

    private void ProteinLevelExport(int TopNPep, int TopNFrag, float Freq) throws IOException {
        FileWriter proWriter = new FileWriter(WorkFolder + "ProtSummary_" + DateTimeTag.GetTag() + ".xls");
        FileWriter pepWriter = new FileWriter(WorkFolder + "PeptideSummary_" + DateTimeTag.GetTag() + ".xls");
        FileWriter NOWriter = new FileWriter(WorkFolder + "IDNoSummary_" + DateTimeTag.GetTag() + ".xls");
        FileWriter fragWriter = new FileWriter(WorkFolder + "FragSummary_" + DateTimeTag.GetTag() + ".xls");
        //Fragment Summary//////////////////////////////////////////
        for (LCMSID IDsummary : FileList) {
            HashMap<String, FragmentPeak> FragMap = IDSummaryFragments.get(FilenameUtils.getBaseName(IDsummary.mzXMLFileName));
            if(FragMap==null){
                LogManager.getRootLogger().error("Cannot find fragment map for "+FilenameUtils.getBaseName(IDsummary.mzXMLFileName));
                 LogManager.getRootLogger().debug("Printing all fragment maps");
                for(String key : FragMap.keySet()){
                    LogManager.getRootLogger().debug(key);
                }
            }
            for (String key : CombineProtID.ProteinList.keySet()) {
                if (IDsummary.ProteinList.containsKey(key)) {
                    ProtID protein = IDsummary.ProteinList.get(key);
                    for (PepIonID pep : protein.PeptideID.values()) {
                        for (FragmentPeak frag : pep.FragmentPeaks) {
                            final String mapkey =
                                    // key + ";" + pep.GetKey() + ";" + frag.IonType;
                                    key + ";" + pep.GetKey() + ";" + frag.IonType + ";+" + frag.Charge;
//                            ProteinFragMap.putIfAbsent(mapkey, frag.FragMZ);// needs Java 8
                            if(!ProteinFragMap.containsKey(mapkey))
                                ProteinFragMap.put(mapkey, frag.FragMZ);
                            frag.Prob1 = pep.MaxProbability;
                            frag.Prob2 = pep.TargetedProbability();
                            frag.RT = pep.PeakRT;
                            FragMap.put(mapkey, frag);
                        }
                    }
                }
            }
        }

        NOWriter.write("File\tNo. Proteins\tNo. peptide ions (Spec-centric)\tNo. peptide ions (Pep-centric)\tNo. proein assoc. ions\n");

        //ProteinSummary/////////////
        proWriter.write("Protein Key\t");
        proWriter.write("Selected_peptides\t");
        for (LCMSID IDSummary : FileList) {
            NOWriter.write(FilenameUtils.getBaseName(IDSummary.mzXMLFileName) + "\t" + IDSummary.ProteinList.size() + "\t" + IDSummary.GetPepIonList().size() + "\t" + IDSummary.GetMappedPepIonList().size() + "\t" + IDSummary.AssignedPepIonList.size() + "\n");
            String file = FilenameUtils.getBaseName(IDSummary.mzXMLFileName);
            proWriter.write(file + "_Prob\t" + file + "_Peptides\t" + file + "_PSMs\t" + file + "_MS1_iBAQ\t" + file + "_Top" + TopNPep + "pep/Top" + TopNFrag + "fra ,Freq>" + Freq + "\t");
        }
        proWriter.write("\n");

        NOWriter.close();
        if (CombineProtID != null) {
            for (String key : CombineProtID.ProteinList.keySet()) {
                proWriter.write(key + "\t");
                for (final LCMSID IDsummary : FileList) {
                    if (IDsummary.ProteinList.containsKey(key)) {
                        final ProtID protein = IDsummary.ProteinList.get(key);
                        proWriter.write(String.join("|", fragselection.TopPeps.get(protein.getAccNo())));
                        break;
                    }
                }
                proWriter.write("\t");
                for (LCMSID IDsummary : FileList) {
                    if (IDsummary.ProteinList.containsKey(key)) {
                        ProtID protein = IDsummary.ProteinList.get(key);
                        proWriter.write(protein.Probability + "\t" + protein.PeptideID.size() + "\t" + protein.GetSpectralCount() + "\t" + protein.GetAbundanceByMS1_IBAQ() + "\t" + protein.GetAbundanceByTopCorrFragAcrossSample(fragselection.TopPeps.get(protein.getAccNo()), fragselection.TopFrags) + "\t");
                    } else {
                        proWriter.write("\t\t\t\t\t");
                    }
                }
                proWriter.write("\n");
            }
            proWriter.close();
        }

        fragWriter.write("Fragment Key\tProtein\tPeptide\tFragment\tFragMz\t");
        for (LCMSID IDSummary : FileList) {
            String file = FilenameUtils.getBaseName(IDSummary.mzXMLFileName);
            fragWriter.write(file + "_RT\t" + file + "_Spec_Centric_Prob\t" + file + "_Pep_Centric_Prob\t" + file + "_Intensity\t" + file + "_Corr\t" + file + "_PPM\t");
        }
        fragWriter.write("\n");
          for (final Map.Entry<String, Float> ent : new java.util.TreeMap<>(ProteinFragMap).entrySet()) {
            final String key = ent.getKey();
            fragWriter.write(key + "\t" + key.split(";")[0] + "\t" + key.split(";")[1] + "\t" + key.split(";")[2] + key.split(";")[3] + "\t" + ent.getValue() + "\t");
            for (LCMSID IDSummary : FileList) {
                if (IDSummaryFragments.get(FilenameUtils.getBaseName(IDSummary.mzXMLFileName)).containsKey(key)) {
                    FragmentPeak fragmentPeak = IDSummaryFragments.get(FilenameUtils.getBaseName(IDSummary.mzXMLFileName)).get(key);
                    fragWriter.write(fragmentPeak.RT + "\t" + fragmentPeak.Prob1 + "\t" + fragmentPeak.Prob2 + "\t" + fragmentPeak.intensity + "\t" + fragmentPeak.corr + "\t" + fragmentPeak.ppm + "\t");
                } else {
                    fragWriter.write("\t\t\t\t\t\t");
                }
            }
            fragWriter.write("\n");
        }
        fragWriter.close();

        ////PepSummary///////////////////////////////////
        for (LCMSID IDsummary : FileList) {
            HashMap<String, FragmentPeak> FragMap = new HashMap<>();
            IDSummaryFragments.put(FilenameUtils.getBaseName(IDsummary.mzXMLFileName), FragMap);

            for (String key : IDsummary.GetPepIonList().keySet()) {
                if (!IdentifiedPepMap.contains(key)) {
                    IdentifiedPepMap.add(key);
                }
            }
            for (String key : IDsummary.GetMappedPepIonList().keySet()) {
                if (!IdentifiedPepMap.contains(key)) {
                    IdentifiedPepMap.add(key);
                }
            }
        }

        pepWriter.write("Peptide Key\tSequence\tModSeq\tProteins\tmz\tCharge\tMaxProb\t");
        pepWriter.write("Selected_fragments\t");
        for (LCMSID IDSummary : FileList) {
            String file = FilenameUtils.getBaseName(IDSummary.mzXMLFileName);
            pepWriter.write(file + "_Spec_Centric_Prob\t" + file + "_Pep_Centric_Prob\t" + file + "_PSMs\t" + file + "_RT\t" + file + "_MS1\t" + file + "_Top" + TopNFrag + "fra\t");
        }
        pepWriter.write("\n");

        for (String key : IdentifiedPepMap) {
            pepWriter.write(key + "\t");
            float maxprob = 0f;
            boolean output = false;
            for (LCMSID IDSummary : FileList) {
                if (IDSummary.GetPepIonList().containsKey(key)) {
                    PepIonID peptide = IDSummary.GetPepIonList().get(key);                    
                    if (!output) {
                        pepWriter.write(peptide.Sequence + "\t" + peptide.ModSequence + "\t" + peptide.ParentProteins() + "\t" + peptide.ObservedMz+ "\t" + peptide.Charge + "\t");
                        output = true;
                    }
                    if (peptide.MaxProbability > maxprob) {
                        maxprob = peptide.MaxProbability;
                    }
                }
                if (IDSummary.GetMappedPepIonList().containsKey(key)) {
                    PepIonID peptide = IDSummary.GetMappedPepIonList().get(key);
                    if (!output) {
                        pepWriter.write(peptide.Sequence + "\t" + peptide.ModSequence + "\t" + peptide.ParentProteins() + "\t" + peptide.ObservedMz + "\t" + peptide.Charge + "\t");
                        output = true;
                    }
                    if (peptide.TargetedProbability() > maxprob) {
                        maxprob = peptide.TargetedProbability();
                    }
                }
            }
            pepWriter.write(maxprob + "\t");

            for (final LCMSID IDSummary : FileList) {
                if (IDSummary.GetPepIonList().containsKey(key)) {
                    final PepIonID peptide = IDSummary.GetPepIonList().get(key);
                    pepWriter.write(String.join("|", fragselection.TopFrags.get(peptide.GetKey())));
                    break;
                } else if (IDSummary.GetMappedPepIonList().containsKey(key)) {
                    final PepIonID peptide = IDSummary.GetMappedPepIonList().get(key);
                    pepWriter.write(String.join("|", fragselection.TopFrags.get(peptide.GetKey())));
                    break;
                } else {
                }
            }
            pepWriter.write("\t");

            for (LCMSID IDSummary : FileList) {
                if (IDSummary.GetPepIonList().containsKey(key)) {
                    PepIonID peptide = IDSummary.GetPepIonList().get(key);
                    pepWriter.write(peptide.MaxProbability + "\t-1\t" + peptide.GetSpectralCount() + "\t" + peptide.PeakRT + "\t" + peptide.PeakHeight[0] + "\t" + peptide.GetPepAbundanceByTopCorrFragAcrossSample(fragselection.TopFrags.get(peptide.GetKey())) + "\t");
                } else if (IDSummary.GetMappedPepIonList().containsKey(key)) {
                    PepIonID peptide = IDSummary.GetMappedPepIonList().get(key);
                    pepWriter.write("-1\t" + peptide.TargetedProbability() + "\t" + peptide.GetSpectralCount() + "\t" + peptide.PeakRT + "\t" + peptide.PeakHeight[0] + "\t" + peptide.GetPepAbundanceByTopCorrFragAcrossSample(fragselection.TopFrags.get(peptide.GetKey())) + "\t");
                } else {
                    pepWriter.write("\t\t\t\t\t\t");
                }
            }
            pepWriter.write("\n");
        }
        pepWriter.close();
    }

    public void PeptideLevelExport(int TopNPep, int TopNFrag, float Freq) throws IOException {
        FileWriter pepWriter = new FileWriter(WorkFolder + "PeptideSummary_" + DateTimeTag.GetTag() + ".xls");
        FileWriter NOWriter = new FileWriter(WorkFolder + "IDNoSummary_" + DateTimeTag.GetTag() + ".xls");
        FileWriter fragWriter = new FileWriter(WorkFolder + "FragSummary_" + DateTimeTag.GetTag() + ".xls");
        //Fragment Summary//////////////////////////////////////////
        for (LCMSID IDsummary : FileList) {
            if (IDSummaryFragments!=null && IDSummaryFragments.containsKey(FilenameUtils.getBaseName(IDsummary.mzXMLFileName))) {
                HashMap<String, FragmentPeak> FragMap = IDSummaryFragments.get(FilenameUtils.getBaseName(IDsummary.mzXMLFileName));
                for (PepIonID pep : IDsummary.GetPepIonList().values()) {
                    GetFragments(pep, FragMap);
                }
                for (PepIonID pep : IDsummary.GetMappedPepIonList().values()) {
                    GetFragments(pep, FragMap);
                }
            }
        }

        fragWriter.write("Fragment Key\tPeptide\tFragment\tFragMz\t");
        for (LCMSID IDSummary : FileList) {
            String file = FilenameUtils.getBaseName(IDSummary.mzXMLFileName);
            fragWriter.write(file + "_RT\t" + file + "_Spec_Centric_Prob\t" + file + "_Pep_Centric_Prob\t" + file + "_Intensity\t" + file + "_Corr\t" + file + "_PPM\t");
        }
        fragWriter.write("\n");
        for (String key : ProteinFragMap.keySet()) {
            fragWriter.write(key + "\t" + key.split(";")[0] + "\t" + key.split(";")[1] + key.split(";")[2] + "\t" + ProteinFragMap.get(key) + "\t");
            for (LCMSID IDSummary : FileList) {
                if (IDSummaryFragments!=null && IDSummaryFragments.get(FilenameUtils.getBaseName(IDSummary.mzXMLFileName)).containsKey(key)) {
                    FragmentPeak fragmentPeak = IDSummaryFragments.get(FilenameUtils.getBaseName(IDSummary.mzXMLFileName)).get(key);
                    fragWriter.write(fragmentPeak.RT + "\t" + fragmentPeak.Prob1 + "\t" + fragmentPeak.Prob2 + "\t" + fragmentPeak.intensity + "\t" + fragmentPeak.corr + "\t" + fragmentPeak.ppm + "\t");
                } else {
                    fragWriter.write("\t\t\t\t\t\t");
                }
            }
            fragWriter.write("\n");
        }
        fragWriter.close();

        //ProteinSummary/////////////
        NOWriter.write("File\tNo. peptide ions (Spec-centric)\tNo. peptide ions (Pep-centric)\n");

        for (LCMSID IDSummary : FileList) {
            NOWriter.write(FilenameUtils.getBaseName(IDSummary.mzXMLFileName) + "\t" + IDSummary.GetPepIonList().size() + "\t" + IDSummary.GetMappedPepIonList().size() + "\n");
        }

        NOWriter.close();

        ////PepSummary///////////////////////////////////
        for (LCMSID IDsummary : FileList) {
            HashMap<String, FragmentPeak> FragMap = new HashMap<>();
            if (IDSummaryFragments != null && IDSummaryFragments.containsKey(FilenameUtils.getBaseName(IDsummary.mzXMLFileName))) {
                IDSummaryFragments.put(FilenameUtils.getBaseName(IDsummary.mzXMLFileName), FragMap);
            }
            for (String key : IDsummary.GetPepIonList().keySet()) {
                if (!IdentifiedPepMap.contains(key)) {
                    IdentifiedPepMap.add(key);
                }
            }
            for (String key : IDsummary.GetMappedPepIonList().keySet()) {
                if (!IdentifiedPepMap.contains(key)) {
                    IdentifiedPepMap.add(key);
                }
            }
        }

        pepWriter.write("Peptide Key\tSequence\tModSeq\tmz\tCharge\tMaxProb\t");

        for (LCMSID IDSummary : FileList) {
            String file = FilenameUtils.getBaseName(IDSummary.mzXMLFileName);
            pepWriter.write(file + "_Spec_Centric_Prob\t" + file + "_Pep_Centric_Prob\t" + file + "_PSMs\t" + file + "_RT\t" + file + "_MS1\t" + file + "_Top" + TopNFrag + "fra\t");
        }
        pepWriter.write("\n");

        for (String key : IdentifiedPepMap) {
            pepWriter.write(key + "\t");
            float maxprob = 0f;
            boolean output = false;
            for (LCMSID IDSummary : FileList) {
                if (IDSummary.GetPepIonList().containsKey(key)) {
                    PepIonID peptide = IDSummary.GetPepIonList().get(key);
                    if (!output) {
                        pepWriter.write(peptide.Sequence + "\t" + peptide.ModSequence + "\t" + peptide.ObservedMz + "\t" + peptide.Charge + "\t");
                        output = true;
                    }
                    if (peptide.MaxProbability > maxprob) {
                        maxprob = peptide.MaxProbability;
                    }
                }
                if (IDSummary.GetMappedPepIonList().containsKey(key)) {
                    PepIonID peptide = IDSummary.GetMappedPepIonList().get(key);
                    if (!output) {
                        pepWriter.write(peptide.Sequence + "\t" + peptide.ModSequence + "\t" + peptide.ParentProteins() + "\t" + peptide.ObservedMz + "\t" + peptide.Charge + "\t");
                        output = true;
                    }
                    if (peptide.TargetedProbability() > maxprob) {
                        maxprob = peptide.TargetedProbability();
                    }
                }
            }
            pepWriter.write(maxprob + "\t");
            for (LCMSID IDSummary : FileList) {
                if (IDSummary.GetPepIonList().containsKey(key)) {
                    PepIonID peptide = IDSummary.GetPepIonList().get(key);
                    float pepabun=0f;
                    if (fragselection != null) {
                        pepabun=peptide.GetPepAbundanceByTopCorrFragAcrossSample(fragselection.TopFrags.get(peptide.GetKey()));
                    }
                    pepWriter.write(peptide.MaxProbability + "\t-1\t" + peptide.GetSpectralCount() + "\t" + peptide.PeakRT + "\t" + peptide.PeakHeight[0] + "\t" + pepabun + "\t");
                } else if (IDSummary.GetMappedPepIonList().containsKey(key)) {
                    PepIonID peptide = IDSummary.GetMappedPepIonList().get(key);
                    float pepabun=0f;
                    if (fragselection != null) {
                        pepabun=peptide.GetPepAbundanceByTopCorrFragAcrossSample(fragselection.TopFrags.get(peptide.GetKey()));
                    }
                    pepWriter.write("-1\t" + peptide.TargetedProbability() + "\t" + peptide.GetSpectralCount() + "\t" + peptide.PeakRT + "\t" + peptide.PeakHeight[0] + "\t" + pepabun+ "\t");
                } else {
                    pepWriter.write("\t\t\t\t\t\t");
                }
            }
            pepWriter.write("\n");
        }
        pepWriter.close();

    }

    private void GetFragments(PepIonID pep, HashMap<String, FragmentPeak> FragMap) {
        LogManager.getRootLogger().info("MSUmpire.Utility.ExportTable.GetFragments()");
        System.exit(123);// never executed
        for (FragmentPeak frag : pep.FragmentPeaks) {
            ProteinFragMap.putIfAbsent(pep.GetKey() + ";" + frag.IonType, frag.FragMZ);
            frag.Prob1 = pep.MaxProbability;
            frag.Prob2 = pep.TargetedProbability();
            frag.RT = pep.PeakRT;
            FragMap.put(pep.GetKey() + ";" + frag.IonType, frag);
        }
    }

    public void ExportPep() throws IOException {
        LogManager.getRootLogger().info("Writing reports...............");
        FileWriter pepWriter = new FileWriter(WorkFolder + "PeptideSummary_" + DateTimeTag.GetTag() + ".xls");

        ////PepSummary///////////////////////////////////
        for (LCMSID IDsummary : FileList) {
            HashMap<String, FragmentPeak> FragMap = new HashMap<>();

            for (String key : IDsummary.GetPepIonList().keySet()) {
                if (!IdentifiedPepMap.contains(key)) {
                    IdentifiedPepMap.add(key);
                }
            }
            for (String key : IDsummary.GetMappedPepIonList().keySet()) {
                if (!IdentifiedPepMap.contains(key)) {
                    IdentifiedPepMap.add(key);
                }
            }
        }

        pepWriter.write("Peptide Key\tSequence\tModSeq\tProteins\tmz\tCharge\tMaxProb\t");

        for (LCMSID IDSummary : FileList) {
            String file = FilenameUtils.getBaseName(IDSummary.mzXMLFileName);
            pepWriter.write(file + "_RT\t" + file + "_MS1\t");
        }
        pepWriter.write("\n");

        for (String key : IdentifiedPepMap) {
            pepWriter.write(key + "\t");
            float maxprob = 0f;
            boolean output = false;
            for (LCMSID IDSummary : FileList) {
                if (IDSummary.GetPepIonList().containsKey(key)) {
                    PepIonID peptide = IDSummary.GetPepIonList().get(key);
                    if (!output) {
                        pepWriter.write(peptide.Sequence + "\t" + peptide.ModSequence + "\t" + peptide.ParentProteins() + "\t" + peptide.ObservedMz + "\t" + peptide.Charge + "\t");
                        output = true;
                    }
                    if (peptide.MaxProbability > maxprob) {
                        maxprob = peptide.MaxProbability;
                    }
                }
            }
            pepWriter.write(maxprob + "\t");
            for (LCMSID IDSummary : FileList) {
                if (IDSummary.GetPepIonList().containsKey(key)) {
                    PepIonID peptide = IDSummary.GetPepIonList().get(key);
                    pepWriter.write(peptide.PeakRT + "\t" + peptide.GetTotalPeakHeight() + "\t");
                } else {
                    pepWriter.write("\t\t");
                }
            }
            pepWriter.write("\n");
        }
        pepWriter.close();

        LogManager.getRootLogger().info("Job done");
        LogManager.getRootLogger().info("=================================================================================================");
    }
}
