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
package MSUmpire.PSMDataStructure;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.SeqUtility.FastaParser;
import MSUmpire.PeakDataStructure.PeakCluster;
import com.compomics.util.experiment.biology.AminoAcid;
import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.PTM;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import com.compomics.util.experiment.biology.ions.PeptideFragmentIon;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.log4j.Logger;

/**
 * Identification data structure for a LC-MS run
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class LCMSID implements Serializable {

    private static final long serialVersionUID = 25494643464616L;

    public HashMap<String, PSM> PSMList;
    public HashMap<String, PSM> LowScorePSMByPepKey;
    public transient HashMap<String, PepIonID> LowScorePep;
    public transient HashMap<String, PSM> LowScorePSM;
    public transient HashMap<Integer, ArrayList<ProtID>> ProteinGroups;
    private HashMap<String, PepIonID> PepIonList;
    private HashMap<Integer, PepIonID> PepIonIndexList;
    public HashMap<String, HashMap<String, PepIonID>> PeptideList;

    private HashMap<String, PepIonID> MappedPepIonList;
    private HashMap<Integer, PepIonID> MappedPepIonIndexList;
    public HashMap<String, HashMap<String, PepIonID>> MappedPeptideList;
    
    public HashMap<String, PepIonID> AssignedPepIonList;
    public HashMap<String, PepIonID> ProtXMLPepIonList;
    
    public HashMap<String, ProtID> ProteinList;
    public HashMap<String, ProtID> IndisProteinIDList;
    
    public HashMap<String, ModificationInfo> ModificationList;
    public HashMap<String, ProtID> PepXMLProteinList;
    
    
    public String DataBase;
    public String SearchEngine;
    public String msModel;
    public String msManufacturer;
    public String msIonization;
    public String msMassAnalyzer;
    public String msDetector;
    public String mzXMLFileName;
    public float FDR = Float.POSITIVE_INFINITY;
    public float ProteinFDR = Float.POSITIVE_INFINITY;
    public float ExpectThreshold = Float.POSITIVE_INFINITY;
    public float SpecProbThreshold = 0f;
    public float PepProbThreshold = 0f;
    public float ProteinProbThreshold = 0f;
    public String DecoyTag = "rev_";
    public String FastaPath;
    private float NorFactor = 1f;
    private transient FastaParser fastaParser;
    public HashMap<String, String> LuciphorResult;
    public String Filename; 

    private FastaParser GetFastaParser() {
        if (fastaParser == null) {
            fastaParser = new FastaParser(FastaPath);
            fastaParser.RemoveDecoy(DecoyTag);
        }
        return fastaParser;
    }
  
    public void WriteLCMSIDSerialization(String filepath) {
        WriteLCMSIDSerialization(filepath, "");
    }

    public void WriteLCMSIDSerialization(String filepath, String tag) {        
        if (!FSWrite(filepath, tag)) {
            Logger.getRootLogger().debug("Writing LCMSID FS failed. writing standard serialization instead");            
        }
    }

    public void RemoveLowWeightPep(float threshold) {
        for (ProtID protein : ProteinList.values()) {
            protein.RemoveLowWeightPepID(threshold);
        }
    }

    public LCMSID CreateEmptyLCMSID() {
        LCMSID newLcmsid = new LCMSID(mzXMLFileName, DecoyTag, FastaPath);
        newLcmsid.SearchEngine = SearchEngine;
        newLcmsid.msDetector = msDetector;
        newLcmsid.msIonization = msIonization;
        newLcmsid.msManufacturer = msManufacturer;
        newLcmsid.msMassAnalyzer = msMassAnalyzer;
        newLcmsid.msModel = msModel;
        newLcmsid.DataBase = DataBase;
        for (ModificationInfo mod : ModificationList.values()) {
            newLcmsid.ModificationList.put(mod.GetKey(), mod);
        }
        return newLcmsid;
    }
    
    public HashMap<String,LCMSID> GetLCMSIDFileMap(){
        
        HashMap<String,LCMSID> lcmsidmap=new HashMap<>();
        
        for(PSM psm : this.PSMList.values()){
            if (!lcmsidmap.containsKey(psm.GetRawNameString())) {
                LCMSID newLcmsid = CreateEmptyLCMSID();
                newLcmsid.mzXMLFileName = psm.GetRawNameString();
                lcmsidmap.put(psm.GetRawNameString(), newLcmsid);
            }            
            lcmsidmap.get(psm.GetRawNameString()).AddPSM(psm);            
        }        
        return lcmsidmap;
    }
  
    private boolean FSWrite(String filepath, String tag) {
        try {
            if (!tag.equals("")) {
                tag = "_" + tag;
            }
            Logger.getRootLogger().info("Writing ID results to file:" + FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS", false);
            FSTObjectOutput out = new FSTObjectOutput(fout);
            ReduceMemoryUsage();
            out.writeObject(this, LCMSID.class);
            out.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    private static LCMSID FS_Read(String filepath, String tag) throws Exception {
        if (!tag.equals("")) {
            tag = "_" + tag;
        }
        if (!new File(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS").exists()) {
            return null;
        }
        try {
            Logger.getRootLogger().info("Reading ID results from file:" + FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            LCMSID lcmsid = (LCMSID) in.readObject(LCMSID.class);
            in.close();
            fileIn.close();
            return lcmsid;

        } catch (Exception ex) {
            Logger.getRootLogger().info("Reading LCMSID FS results failed.");
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
    }

    public static LCMSID ReadLCMSIDSerialization(String filepath) throws Exception {
        return ReadLCMSIDSerialization(filepath, "");
    }

    public static LCMSID ReadLCMSIDSerialization(String filepath, String tag) throws Exception {
        LCMSID lcmsid = FS_Read(filepath, tag);        
        return lcmsid;
    }
   
    public float GetNorFactor() {
        if (NorFactor == 1f) {
            NorFactor = 0f;
            for (PepIonID pepIonID : GetPepIonList().values()) {
                NorFactor += pepIonID.PeakHeight[0];
            }
            NorFactor /= GetPepIonList().size();
        }
        return NorFactor;
    }

    public LCMSID(String mzXMLFileName, String Decoytag, String Fasta) {
        ModificationList = new HashMap<>();
        ProteinList = new HashMap<>();
        PSMList = new HashMap<>();
        PepIonList = new HashMap<>();
        AssignedPepIonList = new HashMap<>();
        ProtXMLPepIonList = new HashMap<>();
        IndisProteinIDList = new HashMap<>();
        ProteinGroups = new HashMap<>();
        PeptideList = new HashMap<>();
        MappedPepIonList = new HashMap<>();
        this.mzXMLFileName = mzXMLFileName;
        this.FastaPath = Fasta;
        this.DecoyTag = Decoytag;
    }

    public HashMap<String, PepIonID> GetPepIonList() {
        return PepIonList;
    }

    public HashMap<String, PepIonID> GetMappedPepIonList() {
        return MappedPepIonList;
    }

    public void ResetMappedPepProb() {
        for (PepIonID pepIonID : MappedPepIonList.values()) {
            pepIonID.UScoreProbability_MS1 = 0f;
            pepIonID.MS1AlignmentProbability = 0f;
            pepIonID.UScoreProbability_MS2 = 0f;
            pepIonID.MS2AlignmentProbability = 0f;
        }
    }

    public void SetMappedPepIonList(HashMap<String, PepIonID> list) {
        MappedPepIonList = list;
    }

    public void GenerateFragmentPeakForPepIonByMSMS(ScanCollection scanCollection, float fragPPM) {

        double protonMass = ElementaryIon.proton.getTheoreticMass();
        for (PepIonID pepIonID : GetPepIonList().values()) {
            for (PSM psm : pepIonID.GetPSMList()) {
                ScanData scan = scanCollection.GetScan(psm.ScanNo);
                for (Ion frag : pepIonID.GetFragments()) {
                    XYData closetPeak = null;
                    float targetmz = (float) (frag.getTheoreticMass() + protonMass);
                    closetPeak = scan.GetHighestPeakInMzWindow(targetmz, fragPPM);
                    if (closetPeak != null) {
                        FragmentPeak fragmentpeak = new FragmentPeak();
                        fragmentpeak.ObservedMZ = closetPeak.getX();
                        fragmentpeak.FragMZ = targetmz;
                        fragmentpeak.intensity = closetPeak.getY();
                        fragmentpeak.Charge = 1;
                        fragmentpeak.ppm = InstrumentParameter.CalcSignedPPM(closetPeak.getX(), targetmz);
                        fragmentpeak.IonType = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber();
                        pepIonID.FragmentPeaks.add(fragmentpeak);
                    }

                    targetmz = (float) (frag.getTheoreticMass() + protonMass * 2) / 2;
                    closetPeak = scan.GetHighestPeakInMzWindow(targetmz, fragPPM);

                    if (closetPeak != null) {
                        FragmentPeak fragmentpeak = new FragmentPeak();
                        fragmentpeak.ObservedMZ = closetPeak.getX();
                        fragmentpeak.FragMZ = targetmz;
                        fragmentpeak.intensity = closetPeak.getY();
                        fragmentpeak.Charge = 2;
                        fragmentpeak.ppm = InstrumentParameter.CalcSignedPPM(closetPeak.getX(), targetmz);
                        fragmentpeak.IonType = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber();
                        pepIonID.FragmentPeaks.add(fragmentpeak);
                    }

                }
            }
        }
    }


    public void AssignProtForPepIon() {
        if (PeptideList != null) {
            for (String pepseq : PeptideList.keySet()) {                 
                for (ProtID protein : ProteinList.values()) {
                    if (protein.Sequence.contains(pepseq) || (protein.ProtPepSeq != null && protein.ProtPepSeq.contains(pepseq))) {
                        for (PepIonID pepIonID : PeptideList.get(pepseq).values()) {
                            protein.PeptideID.put(pepIonID.GetKey(), pepIonID);
                            if (!pepIonID.ParentProtID_ProtXML.contains(protein)) {
                                pepIonID.ParentProtID_ProtXML.add(protein);
                            }
                        }
                    }
                }
            }
        }
    }

    public void AssignProtForMappedIon() {
        if (MappedPeptideList != null) {
            for (String pepseq : MappedPeptideList.keySet()) {
                for (ProtID protein : ProteinList.values()) {
                    if (protein.Sequence.contains(pepseq) || (protein.ProtPepSeq != null && protein.ProtPepSeq.contains(pepseq))) {
                        for (PepIonID pepIonID : MappedPeptideList.get(pepseq).values()) {
                            protein.PeptideID.put(pepIonID.GetKey(), pepIonID);
                            if (!pepIonID.ParentProtID_ProtXML.contains(protein)) {
                                pepIonID.ParentProtID_ProtXML.add(protein);
                            }
                        }
                    }
                }
            }
        }
    }
    
    public void ExportMappedPepID() throws SQLException, IOException {
        ExportMappedPepIonCSV();
    }

     public void ExportPepID() throws SQLException, IOException {
         ExportPepID(null);
     }
    public void ExportPepID(String folder) throws IOException {
        ExportPepIonCSV(folder);
        ExportPepPSMCSV(folder);
    }

    public void ExportProtID() throws SQLException, IOException {
        ExportProtID(null);
    }
    public void ExportProtID(String folder) throws SQLException, IOException {
        ExportProtIDCSV(folder);
    }

    public void ExportPepFragmentPeak() throws SQLException, IOException {
        ExportPepFragmentCSV();
    }

    public void ExportMappedPepFragmentPeak() throws SQLException, IOException {
        ExportMappedPepFragmentCSV();
    }

    private void ExportPepPSMCSV(String folder) throws IOException {
        if(folder==null | "".equals(folder)){
            folder=FilenameUtils.getFullPath(mzXMLFileName);
        }        
        Logger.getRootLogger().info("Writing PSM result to file:" + folder + FilenameUtils.getBaseName(mzXMLFileName) + "_PSMs.csv...");
        FileWriter writer = new FileWriter(folder + FilenameUtils.getBaseName(mzXMLFileName) + "_PSMs.csv");
        writer.write("SpecID,Sequence,ModSeq,TPPModSeq,Modification,Charge,mz,NeutralPepMass,ObservedMass,RT,AdjustedRT,Rank,ScanNo,PreAA,NextAA,MissedCleavage,ExpectValue,MassError,Prob,Rawname,ParentPepIndex,MS1Quant\n");
        for (PepIonID pepion : PepIonList.values()) {
            for (PSM psm : pepion.GetPSMList()) {
                writer.write(psm.SpecNumber + "," + psm.Sequence + "," + psm.ModSeq + "," + psm.TPPModSeq + "," +  psm.GetModificationString() + "," + psm.Charge + "," + psm.ObserPrecursorMz() + "," + psm.NeutralPepMass + "," + psm.ObserPrecursorMass + "," + psm.RetentionTime + "," + psm.NeighborMaxRetentionTime + "," + psm.Rank + "," + psm.ScanNo + "," + psm.PreAA + "," + psm.NextAA + "," + psm.MissedCleavage + "," + psm.expect + "," + psm.MassError + "," + psm.Probability + "," + psm.RawDataName + "," + pepion.Index + "," + pepion.GetMS1() + "\n");
            }
        }
        writer.close();
    }

    private void ExportMappedPepIonCSV() throws IOException {
        Logger.getRootLogger().info("Writing MappedPepIonIDs result to file:" + FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs.csv");
        writer.write("PepIndex,Sequence,ModSeq,TPPModSeq,ModInfo,Charge,mz,PredictRT, PeakRT,MS1ClusIndex,MS2ClusIndex,PeakScore,PeakHeight1,PeakHeight2,PeakHeight3,PeakArea1,PeakArea2,PeakArea3,MS1AlignmentProb,MS1AlignmentLProb,MS2AlignmentProb,MS2AlignmentLProb\n");
        int index = 1;
        for (PepIonID pepion : MappedPepIonList.values()) {
            writer.write((index++) + "," + pepion.Sequence + "," + pepion.ModSequence + "," + pepion.TPPModSeq + "," + pepion.GetModificationString() + "," + pepion.Charge + "," + pepion.NeutralPrecursorMz() + "," + pepion.PredictRTString() + "," + pepion.PeakRT + "," + pepion.GetMS1ClusIndex() + "," + pepion.GetMS2ClusIndex() + "," + pepion.PeakClusterScore + "," + pepion.PeakHeight[0] + "," + pepion.PeakHeight[1] + "," + pepion.PeakHeight[2] + "," + pepion.PeakArea[0] + "," + pepion.PeakArea[1] + "," + pepion.PeakArea[2] + "," + pepion.MS1AlignmentProbability + "," + pepion.UScoreProbability_MS1 + "," + pepion.MS2AlignmentProbability + "," + pepion.UScoreProbability_MS2 + "\n");
        }
        writer.close();
    }

    private void ExportPepIonCSV(String folder) throws IOException {
        if(folder==null | "".equals(folder)){
            folder=FilenameUtils.getFullPath(mzXMLFileName);
        }
        
        Logger.getRootLogger().info("Writing PepIon result to file:" + folder + FilenameUtils.getBaseName(mzXMLFileName) + "_PepIonIDs.csv...");
        FileWriter writer = new FileWriter(folder + FilenameUtils.getBaseName(mzXMLFileName) + "_PepIonIDs.csv");
        writer.write("PepIndex,Sequence,ModSeq,TPPModSeq,IsNonDegenerate,Charge,mz,IDRT,PeakRT,NoPSMs,MS1ClusIndex,MS2ClusIndex,PeakScore,PeakHeight1,PeakHeight2,PeakHeight3,PeakArea1,PeakArea2,PeakArea3\n");
        for (PepIonID pepion : PepIonList.values()) {
            writer.write(pepion.Index + "," + pepion.Sequence + "," + pepion.ModSequence + "," + pepion.TPPModSeq + "," + (pepion.Is_NonDegenerate ? 1 : 0) + "," + pepion.Charge + "," + pepion.NeutralPrecursorMz() + "," + pepion.GetIDRT() + "," + pepion.PeakRT + "," + pepion.GetSpectralCount() + "," + pepion.GetMS1ClusIndex() + "," + pepion.GetMS2ClusIndex() + "," + pepion.PeakClusterScore + "," + pepion.PeakHeight[0] + "," + pepion.PeakHeight[1] + "," + pepion.PeakHeight[2] + "," + pepion.PeakArea[0] + "," + pepion.PeakArea[1] + "," + pepion.PeakArea[2] + "\n");
        }
        writer.close();
    }

    private void ExportProtIDCSV(String folder) throws IOException {
        if(folder==null | "".equals(folder)){
            folder=FilenameUtils.getFullPath(mzXMLFileName);
        }
        Logger.getRootLogger().info("Writing ProteinID result to file:" + folder + FilenameUtils.getBaseName(mzXMLFileName) + "_ProtIDs.csv...");
        FileWriter writer = new FileWriter(folder + FilenameUtils.getBaseName(mzXMLFileName) + "_ProtIDs.csv");
        writer.write("AccNo,UniProtID,ProteinLength,ProteinGroup,IndisProt,Description,Mass,Score,Peptides,Sequence\n");
        for (ProtID protein : ProteinList.values()) {
            String pepstring = "";
            for (PepIonID pep : protein.PeptideID.values()) {
                pepstring += pep.GetKey() + ";";
            }
            String IndisProt = "";
            for (String indisprot : protein.IndisProteins) {
                IndisProt += indisprot + ";";
            }
            writer.write(protein.getAccNo() + "," + protein.UniProtID + "," + protein.ProteinLength + "," + protein.ProteinGroup + "," + IndisProt + "," + protein.Description + "," + protein.Mass + "," + protein.Probability + "," + pepstring + "," + protein.Sequence + "\n");
        }
        writer.close();
    }

    private void ExportMappedPepFragmentCSV() throws IOException {
        Logger.getRootLogger().info("Writing PepFragment result to file:" + FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepFragments.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepFragments.csv");
        writer.write("PepIndex,IonType,fragMZ,ObservedMZ,Charge,Intensity,Correlation,PPM,ApexDelta,RTOverlapP\n");
        for (PepIonID pepion : MappedPepIonList.values()) {
            for (FragmentPeak frag : pepion.FragmentPeaks) {
                writer.write(pepion.Index + "," + frag.IonType + "," + frag.FragMZ + "," + frag.ObservedMZ + "," + frag.Charge + "," + frag.intensity + "," + frag.corr + "," + frag.ppm + "," + frag.ApexDelta + "," + frag.RTOverlapP + "\n");
            }
        }
        writer.close();
    }

    private void ExportPepFragmentCSV() throws IOException {
        Logger.getRootLogger().info("Writing PepFragment result to file:" + FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_PepFragments.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_PepFragments.csv");
        writer.write("PepIndex,IonType,fragMZ,ObservedMZ,Charge,Intensity,Correlation,PPM,ApexDelta,RTOverlapP\n");
        for (PepIonID pepion : PepIonList.values()) {
            for (FragmentPeak frag : pepion.FragmentPeaks) {
                writer.write(pepion.Index + "," + frag.IonType + "," + frag.FragMZ + "," + frag.ObservedMZ + "," + frag.Charge + "," + frag.intensity + "," + frag.corr + "," + frag.ppm + "," + frag.ApexDelta + "," + frag.RTOverlapP + "\n");
            }
        }
        writer.close();
    }

    public PepIonID GetPepID(PSM psm) {
        return PepIonList.get(psm.GetPepKey());
    }

    public void AddPSM(PSM psm) {
        PSMList.put(psm.SpecNumber, psm);
        PepIonID pepIonID = GetPepID(psm);
        if (pepIonID == null) {
            pepIonID = new PepIonID();
            pepIonID.SetInfobyPSM(psm);
            AddPeptideID(pepIonID);
        }
        pepIonID.AddPSM(psm);
    }

    public void LoadSequence() {
        for (ProtID protID : ProteinList.values()) {
            if (protID.Sequence == null || "".equals(protID.Sequence) || "null".equals(protID.Sequence)) {
                if (!protID.IsDecoy(DecoyTag)) {
                    if (FastaPath != null && !"".equals(FastaPath)) {
                        //String Sequence = GetSequenceFactory().getProtein(protID.UniProtID).getSequence();                    
                        //String Sequence = GetFastaParser().ProteinList.get(protID.getAccNo()).Seq;
                        try {
                            String Sequence = GetFastaParser().GetProtSeq(protID.getAccNo());
                            if (Sequence != null) {
                                protID.SetSequence(Sequence);
                            } else {
                                Logger.getRootLogger().error("Can't find sequence in fasta file for protein:" + protID.getAccNo());
                            }
                        } catch (Exception ex) {
                            Logger.getRootLogger().error("Can't find sequence in fasta file for protein:" + protID.getAccNo());
                            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                        }
                    }
                }
            }
        }
    }

    public void AddProtID(ProtID protID) {
        if (!ProteinList.containsKey(protID.getAccNo())) {
            ProteinList.put(protID.getAccNo(), protID);
        }
    }

    public ArrayList<PSM> FindPsmsBymzRT(float mz, float RT, int charge, float PPM, float RTtol) {
        ArrayList<PSM> psms = new ArrayList<>();
        for (PSM psm : PSMList.values()) {
            if (psm.Charge == charge && InstrumentParameter.CalcPPM(mz, psm.ObserPrecursorMz()) < PPM && Math.abs(psm.RetentionTime - RT) < RTtol) {
                psms.add(psm);
            }
        }
        return psms;
    }

    private void FindMaxIniProbThresholdByFDR() {
        if (ProteinList.isEmpty()) {
            return;
        }
        SortedProteinListMaxIniProb sortedlist = new SortedProteinListMaxIniProb();
        sortedlist.addAll(ProteinList.values());

        //writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName)+"/" + FilenameUtils.getBaseName(mzXMLFileName)+"_Pro.txt");
        int positive = 0;
        int negative = 0;
        ProtID protein = sortedlist.get(0);
        if (protein.IsDecoy(DecoyTag)) {
            negative++;
        } else {
            positive++;
        }
        for (int i = 1; i < sortedlist.size(); i++) {
            protein = sortedlist.get(i);
            if (protein.IsDecoy(DecoyTag)) {
                negative++;
                //System.out.println(protein.getAccNo()+"-"+protein.ProteinGroup+"-Decoy");
            } else {
                positive++;
                //System.out.println(protein.getAccNo()+"-"+ protein.ProteinGroup);
            }
            if (i + 1 == sortedlist.size() ||
                    (protein.MaxIniProb > sortedlist.get(i + 1).MaxIniProb && (float) negative / (float) (positive) >= ProteinFDR)) {
                ProteinProbThreshold = protein.MaxIniProb;
                Logger.getRootLogger().info("Protein maxiniprob threshold=" + ProteinProbThreshold + " Estimated raw protein FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
                return;
            }
        }
        Logger.getRootLogger().info("Protein maxiniprob threshold=" + ProteinProbThreshold + " Estimated raw protein FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
    }

    public void ClearPeakData() {
        for (PepIonID pepIonID : PepIonList.values()) {
            for (PeakCluster peak : pepIonID.MS1PeakClusters) {
                peak.IsoPeaksCurves = null;
                peak.MonoIsotopePeak = null;
            }
            for (PeakCluster peak : pepIonID.MS2UnfragPeakClusters) {
                peak.IsoPeaksCurves = null;
                peak.MonoIsotopePeak = null;
            }
        }
    }

    public void ClearAssignPeakCluster() {
        for (PepIonID pepIonID : PepIonList.values()) {
            pepIonID.MS1PeakClusters.clear();
            pepIonID.MS2UnfragPeakClusters.clear();
        }
    }

    public void ClearAssignPeakClusterMappedion() {
        for (PepIonID pepIonID : MappedPepIonList.values()) {
            pepIonID.MS1PeakClusters.clear();
            pepIonID.MS2UnfragPeakClusters.clear();
        }
    }

    public void ReduceMemoryUsage() {
        fastaParser = null;
        for (PepIonID pepIonID : PepIonList.values()) {
            pepIonID.ClearPepFragFactory();
            for (PeakCluster peakCluster : pepIonID.MS1PeakClusters) {
                peakCluster.GroupedFragmentPeaks.clear();
                peakCluster.MonoIsotopePeak = null;
            }
            for (PeakCluster peakCluster : pepIonID.MS2UnfragPeakClusters) {
                peakCluster.GroupedFragmentPeaks.clear();
                peakCluster.MonoIsotopePeak = null;
            }
        }
        for (PepIonID pepIonID : MappedPepIonList.values()) {
            pepIonID.ClearPepFragFactory();
            for (PeakCluster peakCluster : pepIonID.MS1PeakClusters) {
                peakCluster.GroupedFragmentPeaks.clear();
                peakCluster.MonoIsotopePeak = null;
            }
            for (PeakCluster peakCluster : pepIonID.MS2UnfragPeakClusters) {
                peakCluster.GroupedFragmentPeaks.clear();
                peakCluster.MonoIsotopePeak = null;
            }
        }
    }

    public void ClearPSMs() {
        PSMList.clear();
        for (PepIonID pepIonID : PepIonList.values()) {
            pepIonID.CalcRT();
            pepIonID.GetSpectralCount();
            pepIonID.GetPSMList().clear();
        }
    }

    public void ReMapProPep() {
        ClearProPeplist();
        GeneratePepSeqList();
        AssignProtForPepIon();
        GenerateMappedPepSeqList();
        AssignProtForMappedIon();
        GenearteAssignIonList();
    }

   
    public void SetFilterByGroupWeight() {
        for (PepIonID pep : GetPepIonList().values()) {
            pep.FilteringWeight = pep.GroupWeight;
        }
        for (PepIonID pep : GetMappedPepIonList().values()) {
            pep.FilteringWeight = pep.GroupWeight;
        }
    }

    public void SetFilterByWeight() {
        for (PepIonID pep : GetPepIonList().values()) {
            pep.FilteringWeight = pep.Weight;
        }
        for (PepIonID pep : GetMappedPepIonList().values()) {
            pep.FilteringWeight = pep.Weight;
        }
    }

    public void FindPepProbThresholdByFDR() {
        if (PepIonList.isEmpty()) {
            return;
        }
        SortedPepListProb sortedlist = new SortedPepListProb();
        sortedlist.addAll(PepIonList.values());

        int positive = 0;
        int negative = 0;

        PepIonID pep = sortedlist.get(0);
        if (pep.IsDecoy(DecoyTag)) {
            negative++;
        } else {
            positive++;
        }
        for (int i = 1; i < sortedlist.size(); i++) {
            pep = sortedlist.get(i);
            if (pep.MaxProbability < 0.1f) {
                break;
            }
            if (pep.IsDecoy(DecoyTag)) {
                negative++;
            } else {
                positive++;
            }
            if (pep.MaxProbability < sortedlist.get(i - 1).MaxProbability && ((float) negative / (float) (positive) >= FDR)) {
                PepProbThreshold = pep.MaxProbability;
                Logger.getRootLogger().info("Probability threshold=" + PepProbThreshold + " Estimated FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
                return;
            }
        }
        Logger.getRootLogger().info("Probability threshold=" + PepProbThreshold + " Estimated FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
    }

    public void GenearteAssignIonList() {
        AssignedPepIonList.clear();
        for (ProtID protein : ProteinList.values()) {
            for (PepIonID pepIonID : protein.PeptideID.values()) {
                AssignedPepIonList.put(pepIonID.GetKey(), pepIonID);
            }
        }
    }

    public void FilterByProteinDecoyFDRUsingLocalPW(String DecoyTag, float fdr) {
        this.DecoyTag = DecoyTag;
        this.ProteinFDR = fdr;
        FindLocalPWThresholdByFDR();
        RemoveLowProbProteinDecoy();
    }
    
    private void RemoveLowProbProteinDecoy() {

        ArrayList<ProtID> removelist = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            if (protein.Probability < ProteinProbThreshold || protein.IsDecoy(DecoyTag)) {
                removelist.add(protein);
            }
        }
        for (ProtID protein : removelist) {
            ProteinList.remove(protein.getAccNo());
        }
        GenearteAssignIonList();
    }
    
    private void FindLocalPWThresholdByFDR() {
        //FileWriter writer = null;
        //try {
        if (ProteinList.isEmpty()) {
            return;
        }
        SortedProteinListProb sortedlist = new SortedProteinListProb();
        sortedlist.addAll(ProteinList.values());

        //writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName)+"/" + FilenameUtils.getBaseName(mzXMLFileName)+"_Pro.txt");
        int positive = 0;
        int negative = 0;
        ProtID protein = sortedlist.get(0);
        if (protein.IsDecoy(DecoyTag)) {
            negative++;
        } else {
            positive++;
        }
        for (int i = 1; i < sortedlist.size(); i++) {
            protein = sortedlist.get(i);
            if (protein.IsDecoy(DecoyTag)) {
                negative++;
                //System.out.println(protein.getAccNo()+"-"+protein.ProteinGroup+"-Decoy");
            } else {
                positive++;
                //System.out.println(protein.getAccNo()+"-"+ protein.ProteinGroup);
            }
            if (protein.Probability < sortedlist.get(i - 1).Probability && (float) negative / (float) (positive) >= ProteinFDR) {
                ProteinProbThreshold = protein.Probability;
                Logger.getRootLogger().info("Protein probability threshold=" + ProteinProbThreshold + " Estimated raw protein FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
                return;
            }
        }
    }
    
    public void RemoveLowLocalPWProtein(float LocalPW) {
        ArrayList<ProtID> removelist = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            if (protein.Probability < LocalPW) {
                removelist.add(protein);
            }
        }
        for (ProtID protein : removelist) {
            ProteinList.remove(protein.getAccNo());
        }
        GenearteAssignIonList();
    }

    public void RemoveLowMaxIniProbProteinDecoy() {

        ArrayList<ProtID> removelist = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            if (protein.MaxIniProb < ProteinProbThreshold || protein.IsDecoy(DecoyTag)) {
                removelist.add(protein);
            }
        }
        for (ProtID protein : removelist) {
            ProteinList.remove(protein.getAccNo());
        }
        GenearteAssignIonList();
    }

    public void RemoveLowMaxIniProbProtein(float maxiniprob) {

        ArrayList<ProtID> removelist = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            if (protein.MaxIniProb < maxiniprob) {
                removelist.add(protein);
            }
        }
        for (ProtID protein : removelist) {
            ProteinList.remove(protein.getAccNo());
        }
    }

    public void DetermineAssignIonListByProtPepSeq() {
        for (ProtID proid : ProteinList.values()) {
            for (String seq : proid.ProtPepSeq) {
                if (PeptideList.containsKey(seq)) {
                    for (PepIonID pep : PeptideList.get(seq).values()) {
                        proid.AddPeptideID(pep);
                        if (!pep.ParentProtID_ProtXML.contains(proid)) {
                            pep.ParentProtID_ProtXML.add(proid);
                        }
                        if (!AssignedPepIonList.containsKey(pep.GetKey())) {
                            AssignedPepIonList.put(pep.GetKey(), pep);
                        }
                    }
                }
            }
        }
    }

    public void RemoveLowProbPep() {

        // = new HashMap<>();
        ArrayList<PepIonID> removelist = new ArrayList<>();
        for (PepIonID pep : GetPepIonList().values()) {
            if (pep.MaxProbability < PepProbThreshold) {
                removelist.add(pep);
            }
        }
        for (PepIonID pep : removelist) {
            GetPepIonList().remove(pep.GetKey());
            //LowScorePep.put(pep.GetKey(), pep);
        }
    }

    public void RemoveDecoyPep() {
        ArrayList<PepIonID> removelist = new ArrayList<>();
        for (PepIonID pep : PepIonList.values()) {
            if (pep.IsDecoy(DecoyTag)) {
                removelist.add(pep);
            }
        }
        for (PepIonID pep : removelist) {
            PepIonList.remove(pep.GetKey());
        }
    }

    public void RemoveDecoyProtein() {
        ArrayList<ProtID> removelist = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            if (protein.IsDecoy(DecoyTag)) {
                removelist.add(protein);
            }
        }
        for (ProtID protein : removelist) {
            ProteinList.remove(protein.getAccNo());
        }
    }

    public void FilterByPepDecoyFDR(String DecoyTag, float fdr) {
        this.DecoyTag = DecoyTag;
        this.FDR = fdr;        
        FindPepProbThresholdByFDR();
        RemoveLowProbPep();
        RemoveDecoyPep();
        GeneratePepSeqList();
    }

    private void GeneratePepSeqList() {
        PeptideList = new HashMap();
        for (PepIonID pepID : GetPepIonList().values()) {
            if (!PeptideList.containsKey(pepID.Sequence)) {
                PeptideList.put(pepID.Sequence, new HashMap<String, PepIonID>());
            }
            if (!PeptideList.get(pepID.Sequence).containsKey(pepID.GetKey())) {
                PeptideList.get(pepID.Sequence).put(pepID.GetKey(), pepID);
            }
        }
    }

    public void GenerateMappedPepSeqList() {
        MappedPeptideList = new HashMap<>();
        for (PepIonID pepID : GetMappedPepIonList().values()) {
            if (!MappedPeptideList.containsKey(pepID.Sequence)) {
                MappedPeptideList.put(pepID.Sequence, new HashMap<String, PepIonID>());
            }
            if (!MappedPeptideList.get(pepID.Sequence).containsKey(pepID.GetKey())) {
                MappedPeptideList.get(pepID.Sequence).put(pepID.GetKey(), pepID);
            }
        }
    }

    //Generate protein list according to mapping of peptide ions to a master protein list
    public void GenerateProteinByRefIDByPepSeq(LCMSID RefID, boolean UseMappedIon) {
        ProteinList.clear();
        AssignedPepIonList.clear();
        for (PepIonID pepIonID : GetPepIonList().values()) {
            pepIonID.Weight = 0f;
            pepIonID.GroupWeight = 0f;
        }
        for (ProtID protein : RefID.ProteinList.values()) {            
            if (!protein.IsDecoy(DecoyTag)) {
                ProtID newprot = protein.CloneProtein();
                newprot.Probability = protein.Probability;
                newprot.ProtPepSeq = (ArrayList<String>) protein.ProtPepSeq.clone();

                for (PepIonID pep : protein.ProtPeptideID.values()) {                        
                    if (PeptideList.containsKey(pep.Sequence)) {
                        for (PepIonID pepIonID : PeptideList.get(pep.Sequence).values()) {
                            pepIonID.Weight = pep.Weight;
                            pepIonID.GroupWeight = pep.GroupWeight;
                            newprot.AddPeptideID(pepIonID);
                            newprot.IDByDBSearch = true;
                            if (!pepIonID.ParentProtID_ProtXML.contains(newprot)) {
                                pepIonID.ParentProtID_ProtXML.add(newprot);
                            }
                        }
                    }
                    if (UseMappedIon) {
                        if (MappedPeptideList.containsKey(pep.Sequence)) {
                            for (PepIonID pepIonID : MappedPeptideList.get(pep.Sequence).values()) {
                                pepIonID.Weight = pep.Weight;
                                pepIonID.GroupWeight = pep.GroupWeight;
                                newprot.AddPeptideID(pepIonID);
                                if (!pepIonID.ParentProtID_ProtXML.contains(newprot)) {
                                    pepIonID.ParentProtID_ProtXML.add(newprot);
                                }
                            }
                        }
                    }
                }
                if (!newprot.PeptideID.isEmpty()) {
                    ProteinList.put(newprot.getAccNo(), newprot);
                }
            }
        }
    }

    public void UpdateProteinKey() {
        ArrayList<ProtID> temp = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            temp.add(protein);
        }
        ProteinList.clear();
        for (ProtID protein : temp) {
            AddProtID(protein);
        }
    }

    public void FixProteinWithDecoyHead() {
        for (ProtID protein : ProteinList.values()) {
            if (protein.IsDecoy(DecoyTag)) {
                for (int i = 0; i < protein.IndisProteins.size(); i++) {
                    if (!(protein.IndisProteins.get(i).startsWith(DecoyTag)|protein.IndisProteins.get(i).endsWith(DecoyTag))) {
                        protein.setAccNo(protein.IndisProteins.get(i));
                        protein.SetDescription(protein.IndisProtDes.get(i));
                        break;
                    }
                }
            }
        }
        UpdateProteinKey();
    }

    public void UpdateDecoyMaxIniProb() {
        for (ProtID protein : ProteinList.values()) {
            if (protein.IsDecoy(DecoyTag)) {
                protein.MaxIniProb = 0f;
                for (PepIonID pep : protein.ProtPeptideID.values()) {
                    boolean include = true;
                    for (String prot : pep.ParentProtString_ProtXML) {
                        if (!(prot.startsWith(DecoyTag)|prot.endsWith(DecoyTag))) {
                            include = false;
                        }
                    }
                    if (include && protein.MaxIniProb < pep.MaxProbability) {
                        protein.MaxIniProb = pep.MaxProbability;
                    }
                }
            }
        }
    }

    public void SetGroupProbForNonDecoyGroupHead() {
        for (ArrayList<ProtID> group : ProteinGroups.values()) {
            boolean allzero = true;
            for (ProtID protein : group) {
                if (protein.Probability > 0) {
                    allzero = false;
                    break;
                }
            }
            if (allzero) {
                for (ProtID protein : group) {
                    if (!protein.IsDecoy(DecoyTag)) {
                        protein.Probability = protein.GroupProb;
                        break;
                    }
                }
            }
        }
    }

    public void FilterByProteinDecoyFDRUsingMaxIniProb(String DecoyTag, float fdr) {
        this.DecoyTag = DecoyTag;
        this.ProteinFDR = fdr;
        FindMaxIniProbThresholdByFDR();
        RemoveLowMaxIniProbProteinDecoy();
    }

    
    public void AddPeptideID(PepIonID pepID) {
//        if (!PepIonList.containsKey(pepID.GetKey())) {
//            pepID.Index = PepIonList.size();
//            PepIonList.put(pepID.GetKey(), pepID);
//        }
        {
            final PepIonID pep = PepIonList.get(pepID.GetKey());
            if (pep == null || pep.MaxProbability < pepID.MaxProbability) {
                // if not in map or a higher probability is found
                PepIonList.put(pepID.GetKey(), pepID);
                pepID.Index = PepIonList.size() - 1;
            }
        }
        for (PSM psm : pepID.GetPSMList()) {
            PSMList.put(psm.SpecNumber, psm);
        }
    }

    public PSM GetPSM(int SpecNum) {
        return PSMList.get(SpecNum);
    }

    public void AddModification(PTM ptm, String site) {
        ModificationInfo modification = new ModificationInfo();
        modification.site = site;
        modification.massdiff = (float) ptm.getMass();
        modification.mass = (float) (ptm.getMass() + AminoAcid.getAminoAcid(site).monoisotopicMass);
        modification.modification = ptm;
        if (!ModificationList.containsKey(ptm.getName() + "_" + site)) {
            ModificationList.put(ptm.getName() + "_" + site, modification);
        }
    }

    public void CreateInstanceForAllPepIon() {
        for (PepIonID pepIonID : GetPepIonList().values()) {
            pepIonID.CreateQuantInstance(4);
        }
    }

    public void ClearMappedPep() {
        MappedPepIonList.clear();
        MappedPeptideList.clear();
        if (MappedPepIonIndexList != null) {
            MappedPepIonIndexList.clear();
        }
    }

    public void ClearProPeplist() {
        for (ProtID protein : ProteinList.values()) {
            protein.PeptideID.clear();
        }
        for (PepIonID pep : GetPepIonList().values()) {
            pep.ParentProtID_ProtXML.clear();
        }
        for (PepIonID pep : GetMappedPepIonList().values()) {
            pep.ParentProtID_ProtXML.clear();
        }
    }

    public void RemoveLowProbMappedIon(float ProbThreshold) {
        //LowScorePep = new HashMap<>();
        ArrayList<PepIonID> removelist = new ArrayList<>();
        for (PepIonID pep : GetMappedPepIonList().values()) {
            if (pep.TargetedProbability() < ProbThreshold) {
                removelist.add(pep);
            }
        }
        for (PepIonID pep : removelist) {
            GetMappedPepIonList().remove(pep.GetKey());
            //LowScorePep.put(pep.GetKey(), pep);
        }
        GenerateMappedPepSeqList();
    }

    public void ReleaseIDs() {
        PSMList = null;
        LowScorePSMByPepKey = null;
        LowScorePep = null;
        LowScorePSM = null;
        PepIonList = null;
        PepIonIndexList = null;
        MappedPepIonIndexList = null;
        AssignedPepIonList = null;
        ProtXMLPepIonList = null;
        ProteinList = null;
        IndisProteinIDList = null;
        MappedPepIonList = null;
        PepXMLProteinList = null;
        PeptideList = null;
        MappedPeptideList = null;
    }

    public float GetRFactor(float ProbThreshold) {
        int forward = 0;
        int decoy = 0;
        for (ProtID protID : ProteinList.values()) {
            if (protID.Probability < ProbThreshold) {
                if (protID.IsDecoy(DecoyTag)) {
                    decoy++;
                } else {
                    forward++;
                }
            }
        }
        float rf = (float) forward / decoy;
        Logger.getRootLogger().info("Caculating R factor: probability threshold =" + ProbThreshold);
        Logger.getRootLogger().info("R factor=" + rf + " (forward/decoy=" + forward + "/" + decoy + ")");
        return rf;
    }
}
