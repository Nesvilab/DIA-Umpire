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

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.FragmentLib.FragmentLibManager;
import MSUmpire.LCMSPeakStructure.LCMSPeakDIAMS2;
import MSUmpire.LCMSPeakStructure.LCMSPeakMS1;
import MSUmpire.BaseDataStructure.DBSearchParam;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PSM;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.SearchResultParser.PepXMLParser;
import MSUmpire.SpectrumParser.DIA_Setting;
import MSUmpire.SpectrumParser.SpectrumParserBase;
import MSUmpire.SpectrumParser.mzXMLParser;
import MSUmpire.Utility.MSConvert;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;

/**
 * Main data structure which represents a DIA file
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class DIAPack {

    private InstrumentParameter parameter;
    //private ConnectionManager connectionManager;
    public String Filename;
    private int NoCPUs = 4;
    public LCMSPeakMS1 MS1FeatureMap;
    public ArrayList<LCMSPeakDIAMS2> DIAWindows;
    public DIA_Setting dIA_Setting = new DIA_Setting();
    private SpectrumParserBase SpectrumParser;
    public LCMSID IDsummary;
    public HashMap<Integer, Integer> ScanClusterMap_Q1;
    public HashMap<Integer, Integer> ScanClusterMap_Q2;
    public HashMap<Integer, String> ScanClusterMap_Q3;
    public ArrayList<String> iProphPepXMLs = new ArrayList<String>();
    public boolean ExportPrecursorPeak = false;
    public boolean ExportFragmentPeak = false;
    public HashMap<Integer, Double> FactorialTable;
    TargetMatchScoring TScoring;
    public boolean Resume;

    
    private void RemoveIndexFile() {
        new File(FilenameUtils.removeExtension(Filename) + ".ScanPosFS").delete();
        new File(FilenameUtils.removeExtension(Filename) + ".RTidxFS").delete();
        new File(FilenameUtils.removeExtension(Filename) + ".ScanRTFS").delete();
        new File(FilenameUtils.removeExtension(Filename) + ".ScanidxFS").delete();
        new File(FilenameUtils.removeExtension(Filename) + ".DIAWindowsFS").delete();
    }
    
    public void FixScanidx() {
        RemoveIndexFile();
        GetSpectrumParser();
    }

    public int Q1Scan = 0;
    public int Q2Scan = 0;
    public int Q3Scan = 0;
  
    //Whether to use IDs from targeted re-extraction for quantification analysis
    public boolean UseMappedIon = false;
    //Whether to use probability threshold for targeted re-extraction to filter IDs
    public boolean FilterMappedIonByProb = true;
    //Probability threshold for targeted re-extraction
    public float MappedIonProbThreshold = 0.95f;

    public DIAPack(String Filename, int NoCPUs) throws FileNotFoundException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, DataFormatException {
        super();
        this.Filename = Filename;
        this.NoCPUs = NoCPUs;
    }

    public String GetBaseName() {
        return FilenameUtils.getBaseName(Filename);
    }

    public void SetNoCPUs(int cpu) {
        this.NoCPUs = cpu;
    }

    public void SetParameter(InstrumentParameter parameter) {
        this.parameter = parameter;
    }

    public void SetDataType(SpectralDataType.DataType datatype) {
        this.dIA_Setting.dataType = datatype;
    }

    public void AddVariableWindow(XYData window) {
        if (dIA_Setting.DIAWindows == null) {
            dIA_Setting.DIAWindows = new TreeMap<>();
        }
        dIA_Setting.DIAWindows.put(window, new ArrayList<Integer>());
    }

    //Primarily for WiSIM
    public void AddMS1Window(XYData window) {
        if (dIA_Setting.MS1Windows == null) {
            dIA_Setting.MS1Windows = new TreeMap<>();
        }
        dIA_Setting.MS1Windows.put(window, new ArrayList<Integer>());
    }

    public void SetWindowSize(float size) {
        dIA_Setting.F_DIA_WindowSize = size;
    }

    public SpectrumParserBase GetSpectrumParser() {
        if (SpectrumParser == null) {
            try {
                SpectrumParser = SpectrumParserBase.GetInstance(Filename, parameter, dIA_Setting.dataType, dIA_Setting, NoCPUs);
            } catch (Exception ex) {
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                Logger.getRootLogger().error("Read spectral file:" + Filename + " failed.");
                System.exit(2);
            }
            dIA_Setting = SpectrumParser.dIA_Setting;
            SaveDIASetting();            
        }
        return SpectrumParser;
    }
 
    //Entry of processing of signal extraction and generating pseudo MS/MS spectra
    public void process() throws SQLException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, FileNotFoundException, DataFormatException, Exception {
        BuildDIAWindows();
        MS1PeakDetection();
        DIAMS2PeakDetection();
    }

    //Building empty data structure for MS2 feature maps 
    public void BuildDIAWindows() throws IOException, DataFormatException, IOException, IOException, IOException, InterruptedException {
        if (dIA_Setting.DIAWindows == null || dIA_Setting.DIAWindows.isEmpty()) {
            GetSpectrumParser();
        }
        DIAWindows = new ArrayList<>();
        Object[] WindowRange = dIA_Setting.DIAWindows.keySet().toArray();
        for (int i = 0; i < WindowRange.length; i++) {
            XYData DiaWinMz = (XYData) WindowRange[i];
            XYData LastWinMz = null;
            if (i < WindowRange.length - 1) {
                LastWinMz = (XYData) WindowRange[i + 1];
            }
            
            LCMSPeakDIAMS2 diawindow = new LCMSPeakDIAMS2(Filename, this, parameter, DiaWinMz, LastWinMz, GetSpectrumParser(), NoCPUs);
            diawindow.Resume=Resume;
            //pass the settings through to MS2 feature map
            diawindow.datattype = dIA_Setting.dataType;
            diawindow.ExportPeakCurveTable = ExportFragmentPeak;
            diawindow.ExportPeakClusterTable = ExportPrecursorPeak;
            DIAWindows.add(diawindow);
        }
    }

    //Read scan number mapping
    public void ReadScanNoMapping() throws FileNotFoundException, IOException {
        ScanClusterMap_Q1 = new HashMap<>();
        ScanClusterMap_Q2 = new HashMap<>();
        ScanClusterMap_Q3 = new HashMap<>();
        BufferedReader reader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q1"));
        BufferedReader reader2 = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q2"));
        BufferedReader reader4 = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q3"));

        String line = "";
        int StartNo = 0;
        if (new File(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mzXML").exists()) {
            BufferedReader mzReader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mzXML"));
            while ((line = mzReader.readLine()) != null) {
                if (line.contains("<scan num=")) {
                    String substr = line.substring(line.indexOf("<scan num=") + 11);
                    StartNo = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    break;
                }
            }
        }

        line = reader.readLine();
        int offset = StartNo - Integer.parseInt(line.split("_")[0]);
        Integer ScanNo = Integer.parseInt(line.split("_")[0]) + offset;
        Integer ClusterIndex = Integer.parseInt(line.split("_")[1]);
        ScanClusterMap_Q1.put(ScanNo, ClusterIndex);

        while ((line = reader.readLine()) != null) {
            ScanNo = Integer.parseInt(line.split("_")[0]) + offset;
            ClusterIndex = Integer.parseInt(line.split("_")[1]);
            ScanClusterMap_Q1.put(ScanNo, ClusterIndex);
        }

        line = "";
        StartNo = 0;
        if (new File(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mzXML").exists()) {
            BufferedReader mzReader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mzXML"));
            while ((line = mzReader.readLine()) != null) {
                if (line.contains("<scan num=")) {
                    String substr = line.substring(line.indexOf("<scan num=") + 11);
                    StartNo = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    break;
                }
            }
        }
        line = reader2.readLine();
        offset = StartNo - Integer.parseInt(line.split("_")[0]);
        ScanNo = Integer.parseInt(line.split("_")[0]) + offset;
        ClusterIndex = Integer.parseInt(line.split("_")[1]);
        ScanClusterMap_Q2.put(ScanNo, ClusterIndex);

        while ((line = reader2.readLine()) != null) {
            ScanNo = Integer.parseInt(line.split("_")[0]) + offset;
            ClusterIndex = Integer.parseInt(line.split("_")[1]);
            ScanClusterMap_Q2.put(ScanNo, ClusterIndex);
        }

        line = "";
        StartNo = 0;
        if (new File(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML").exists()) {
            BufferedReader mzReader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML"));
            while ((line = mzReader.readLine()) != null) {
                if (line.contains("<scan num=")) {
                    String substr = line.substring(line.indexOf("<scan num=") + 11);
                    StartNo = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    break;
                }
            }
        }
        line = reader4.readLine();
        if (line != null) {
            if (line.split(";").length == 3) {
                offset = StartNo - Integer.parseInt(line.split(";")[0]);
                ScanNo = Integer.parseInt(line.split(";")[0]) + offset;
                String WindowClusterIndex = line.split(";")[1] + ";" + line.split(";")[2];
                ScanClusterMap_Q3.put(ScanNo, WindowClusterIndex);
            } else {
                String ClusterIndexS = line.split("_")[1];
                ScanNo = Integer.parseInt(line.split("_")[0]);
                ScanClusterMap_Q3.put(ScanNo, ClusterIndexS);
            }
            while ((line = reader4.readLine()) != null) {
                if (line.split(";").length == 3) {
                    ScanNo = Integer.parseInt(line.split(";")[0]) + offset;
                    String WindowClusterIndex = line.split(";")[1] + ";" + line.split(";")[2];
                    ScanClusterMap_Q3.put(ScanNo, WindowClusterIndex);
                } else {
                    ScanNo = Integer.parseInt(line.split("_")[0]);
                    String ClusterIndexS = line.split("_")[1];
                    ScanClusterMap_Q3.put(ScanNo, ClusterIndexS);
                }
            }
        }
        reader.close();
        reader2.close();
        reader4.close();
    }

     public String GetForLibQ1Name() {
        //return FilenameUtils.getBaseName(ScanCollectionName) + ".Q1";
        return FilenameUtils.getBaseName(Filename) + ".ForLibQ1";
    }

    public String GetForLibQ2Name() {
        //return FilenameUtils.getBaseName(ScanCollectionName) + ".Q2";
        return FilenameUtils.getBaseName(Filename) + ".ForLibQ2";
    }

    public String GetForLibQ3Name() {
        //return FilenameUtils.getBaseName(ScanCollectionName) + ".Q3";
        return FilenameUtils.getBaseName(Filename) + ".ForLibQ3";
    }
    
    public String GetQ1Name() {
        //return FilenameUtils.getBaseName(ScanCollectionName) + ".Q1";
        return FilenameUtils.getBaseName(Filename) + "_Q1";
    }

    public String GetQ2Name() {
        //return FilenameUtils.getBaseName(ScanCollectionName) + ".Q2";
        return FilenameUtils.getBaseName(Filename) + "_Q2";
    }

    public String GetQ3Name() {
        //return FilenameUtils.getBaseName(ScanCollectionName) + ".Q3";
        return FilenameUtils.getBaseName(Filename) + "_Q3";
    }

    public String GetiProphExtPepxml(String tag) {
        return FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_"+tag+".iproph.pep.xml");
    }

    public String GetQ1Pepxml() {
        return FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(GetQ1Name()) + ".pep.xml");
    }

    public String GetQ2Pepxml() {
        return FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(GetQ2Name()) + ".pep.xml");
    }

    public String GetQ3Pepxml() {
        return FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(GetQ3Name()) + ".pep.xml");
    }

    public String GetIPROPHETPepXML() {
        return FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(Filename) + ".iproph.pep.xml";
    }

    public String GetIPROPHETProtXML() {
        return FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(Filename) + ".iproph.prot.xml";
    }

    public String GetCombinePepXML() {
        return FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(Filename) + "_Combine.pep.xml";
    }

    public String GetCombineProtXML() {
        return FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(Filename) + "_Combine.prot.xml";
    }

    public String GetCombineiProphProtXML() {
        return FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(Filename) + "_Combine.iproph.prot.xml";
    }

    //Building scan number mapping from MGF files
    public void GenerateClusterScanNomapping() throws IOException {
        if (!new File(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q1").exists()) {
            String mgfname1 = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mgf";
            String mgfname2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mgf";
            String mgfname4 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf";

            BufferedReader reader1 = new BufferedReader(new FileReader(mgfname1));
            BufferedReader reader2 = new BufferedReader(new FileReader(mgfname2));
            BufferedReader reader4 = new BufferedReader(new FileReader(mgfname4));

            FileWriter writer = new FileWriter(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q1");
            FileWriter writer2 = new FileWriter(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q2");
            FileWriter writer4 = new FileWriter(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q3");

            BufferedReader mzReader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mzXML"));
            String line = "";
            int ScanNo = 0;
            while ((line = mzReader.readLine()) != null) {
                if (line.contains("<scan num=")) {
                    String substr = line.substring(line.indexOf("<scan num=") + 11);
                    ScanNo = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    break;
                }
            }
            while ((line = reader1.readLine()) != null) {
                if (line.startsWith("TITLE=")) {
                    int ClusterIndex = Integer.parseInt(line.split("ClusterIndex:")[1].split(",")[0]);
                    writer.write(ScanNo + "_" + ClusterIndex + "\n");
                    ScanNo++;
                }
            }
            mzReader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mzXML"));
            line = "";
            ScanNo = 0;
            while ((line = mzReader.readLine()) != null) {
                if (line.contains("<scan num=")) {
                    String substr = line.substring(line.indexOf("<scan num=") + 11);
                    ScanNo = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    break;
                }
            }
            while ((line = reader2.readLine()) != null) {
                if (line.startsWith("TITLE=")) {
                    int ClusterIndex = Integer.parseInt(line.split("ClusterIndex:")[1].split(",")[0]);
                    writer2.write(ScanNo + "_" + ClusterIndex + "\n");
                    ScanNo++;
                }
            }

            mzReader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML"));
            line = "";
            ScanNo = 0;
            while ((line = mzReader.readLine()) != null) {
                if (line.contains("<scan num=")) {
                    String substr = line.substring(line.indexOf("<scan num=") + 11);
                    ScanNo = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    break;
                }
            }
            while ((line = reader4.readLine()) != null) {
                if (line.startsWith("TITLE=")) {
                    int ClusterIndex = Integer.parseInt(line.split("ClusterIndex:")[1].split(",")[0]);
                    String DIAwindow = line.split("TITLE=")[1].split(";")[0];
                    writer4.write(ScanNo + ";" + DIAwindow + ";" + ClusterIndex + "\n");
                    ScanNo++;
                }
            }
            writer.close();
            writer2.close();
            writer4.close();
            reader1.close();
            reader2.close();
            mzReader.close();
            reader4.close();
        }
        ReadScanNoMapping();
    }

    public void AssignQuant() throws IOException, SQLException {
        AssignQuant(true);
    }

    //Entry of quantification process
    public void AssignQuant(boolean export) throws IOException, SQLException {
        Logger.getRootLogger().info("Assign peak cluster to identified peptides");
        GenerateClusterScanNomapping();
        
        ExecutorService executorPool = null;
        for (PeakCluster cluster : MS1FeatureMap.PeakClusters) {
            cluster.Identified = false;
        }
        
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            pepIonID.MS1PeakClusters = new ArrayList<>();
            pepIonID.MS2UnfragPeakClusters = new ArrayList<>();
        }
        
        //Assign precursor features and grouped fragments for all identified peptide ions for a isolation window
        for (LCMSPeakDIAMS2 DIAWindow : DIAWindows) {            
            DIA_window_Quant dia_w = new DIA_window_Quant(GetQ1Name(), GetQ2Name(), GetQ3Name(), ScanClusterMap_Q1, ScanClusterMap_Q2, ScanClusterMap_Q3, MS1FeatureMap, DIAWindow, IDsummary, NoCPUs);            
            dia_w.run();            
        }

        executorPool = Executors.newFixedThreadPool(NoCPUs);

        //Match fragments and calculate quantification for each peptide ion
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            DIAAssignQuantUnit quantunit = new DIAAssignQuantUnit(pepIonID, MS1FeatureMap, parameter);
            executorPool.execute(quantunit);
        }
        executorPool.shutdown();

        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }

        if (export) {
            ExportID();
        }
    }

    public void AssignMappedPepQuant(boolean export, FragmentLibManager libManager) throws IOException, SQLException {
        TargetedExtractionQuant(export, libManager, 1.1f,-1f);
    }

    //Quantification process for peptide ions from targeted re-extraction
    public void TargetedExtractionQuant(boolean export, FragmentLibManager libManager, float ReSearchProb, float RTWindow) throws IOException, SQLException {
        if (IDsummary.GetMappedPepIonList().isEmpty()) {
            Logger.getRootLogger().error("There is no peptide ion for targeted re-extraction.");
            return;
        }
        parameter.RT_window_Targeted=RTWindow;
        GenerateClusterScanNomapping();
        ExecutorService executorPool = null;
        
        //Targeted re-extraction scoring
        TScoring = new TargetMatchScoring(Filename, libManager.LibID);

        if (parameter.UseOldVersion) {
            TScoring.SetUseOldVersion();
        }
        Logger.getRootLogger().info("No. of identified peptide ions: " + IDsummary.GetPepIonList().size());
        Logger.getRootLogger().info("No. of mapped peptide ions: " + IDsummary.GetMappedPepIonList().size());
        ArrayList<PepIonID> SearchList = new ArrayList<>();
        //For each peptide ions in targeted re-extraction, determine whether to research the peptide ion given a re-search probability threshold
        for (PepIonID pepIonID : IDsummary.GetMappedPepIonList().values()) {
            if (libManager.PeptideFragmentLib.containsKey(pepIonID.GetKey()) && libManager.GetFragmentLib(pepIonID.GetKey()).FragmentGroups.size() >= 3 && pepIonID.TargetedProbability() < ReSearchProb) {
                pepIonID.CreateQuantInstance(parameter.MaxNoPeakCluster);
                pepIonID.MS1PeakClusters = new ArrayList<>();
                pepIonID.MS2UnfragPeakClusters = new ArrayList<>();
                pepIonID.UScoreProbability_MS1 = 0f;
                pepIonID.MS1AlignmentProbability = 0f;
                pepIonID.UScoreProbability_MS2 = 0f;
                pepIonID.MS2AlignmentProbability = 0f;
                pepIonID.TPPModSeq="Ext";
                SearchList.add(pepIonID);
            }
        }
        Logger.getRootLogger().info("No. of searchable peptide ions: " + SearchList.size());

        for (LCMSPeakDIAMS2 DIAWindow : DIAWindows) {
            Logger.getRootLogger().info("Assigning clusters for peak groups in MS2 isolation window:" + FilenameUtils.getBaseName(DIAWindow.ScanCollectionName));

            if (!DIAWindow.ReadPeakCluster()|| !DIAWindow.ReadPrecursorFragmentClu2Cur()) {
                Logger.getRootLogger().warn("Reading results for " + DIAWindow.ScanCollectionName + " failed");
                continue;
            }

            executorPool = Executors.newFixedThreadPool(NoCPUs);
            //For each target peptide  ion
            for (PepIonID pepIonID : SearchList) {
                if (DIAWindow.DIA_MZ_Range.getX() <= pepIonID.NeutralPrecursorMz() && DIAWindow.DIA_MZ_Range.getY() >= pepIonID.NeutralPrecursorMz()) {
                    //If the spectrum of peptide ion in the spectral library has more than three fragment peaks
                    if (libManager.GetFragmentLib(pepIonID.GetKey()).FragmentGroups.size() >= 3) {
                        //U-score spectral matching
                        UmpireSpecLibMatch matchunit = new UmpireSpecLibMatch(MS1FeatureMap, DIAWindow, pepIonID, libManager.GetFragmentLib(pepIonID.GetKey()), libManager.GetDecoyFragmentLib(pepIonID.GetKey()), parameter);
                        executorPool.execute(matchunit);                        
                        TScoring.libTargetMatches.add(matchunit);
                    } else {
                        Logger.getRootLogger().warn("skipping " + pepIonID.GetKey() + ", it has only " + libManager.GetFragmentLib(pepIonID.GetKey()).FragmentGroups.size() + " matched fragments");
                    }
                }
            }

            //For each identified peptide ion, calculate their U-score for LDA training
            for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
                if (libManager.PeptideFragmentLib.containsKey(pepIonID.GetKey()) && DIAWindow.DIA_MZ_Range.getX() <= pepIonID.NeutralPrecursorMz() && DIAWindow.DIA_MZ_Range.getY() >= pepIonID.NeutralPrecursorMz()) {
                   //If the spectrum of peptide ion in the spectral library has more than three fragment peaks
                    if (libManager.GetFragmentLib(pepIonID.GetKey()).FragmentGroups.size() >= 3) {
                        //U-score spectral matching
                        UmpireSpecLibMatch matchunit = new UmpireSpecLibMatch(MS1FeatureMap, DIAWindow, pepIonID, libManager.GetFragmentLib(pepIonID.GetKey()), libManager.GetDecoyFragmentLib(pepIonID.GetKey()), parameter);
                        matchunit.IdentifiedPeptideIon = true;
                        executorPool.execute(matchunit);
                        TScoring.libIDMatches.add(matchunit);
                    } else {
                        Logger.getRootLogger().warn("skipping " + pepIonID.GetKey() + ", it has only " + libManager.GetFragmentLib(pepIonID.GetKey()).FragmentGroups.size() + " matched fragments");
                    }
                }
            }
            executorPool.shutdown();

            try {
                executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                Logger.getRootLogger().info("interrupted..");
            }
            DIAWindow.ClearAllPeaks();
        }
        
        Logger.getRootLogger().info("Removing entries with no precursor signal hits: total target entries: "+TScoring.libTargetMatches.size());
        ArrayList<UmpireSpecLibMatch> newlist=new ArrayList<>();
        for(UmpireSpecLibMatch match : TScoring.libTargetMatches){
            if (!match.DecoyHits.isEmpty() || !match.TargetHits.isEmpty()) {
                newlist.add(match);
            }
        }
        TScoring.libTargetMatches=newlist;
        Logger.getRootLogger().info("Remaining entries: "+TScoring.libTargetMatches.size());
        
        //U-score and probablilty calculatation  
        TScoring.Process();    
        TargetHitPepXMLWriter pepxml=new TargetHitPepXMLWriter(GetiProphExtPepxml(libManager.LibID), IDsummary.FastaPath, IDsummary.DecoyTag, TScoring);
        TScoring = null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        
        //Assign precursor peak cluster, extract fragments and do quantification
        for (PepIonID pepIonID : IDsummary.GetMappedPepIonList().values()) {
            DIAAssignQuantUnit quantunit = new DIAAssignQuantUnit(pepIonID, MS1FeatureMap, parameter);
            executorPool.execute(quantunit);
        }
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            DIAAssignQuantUnit quantunit = new DIAAssignQuantUnit(pepIonID, MS1FeatureMap, parameter);
            executorPool.execute(quantunit);
        }
        executorPool.shutdown();

        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }
        
        if (export) {
            ExportID();
        }
    }

    public boolean ReadSerializedLCMSID() throws Exception {
        return ReadSerializedLCMSID("");
    }

    public boolean ReadSerializedLCMSID(String tag) throws Exception {
        this.IDsummary = LCMSID.ReadLCMSIDSerialization(Filename, tag);
        if (this.IDsummary == null) {
            return false;
        }
        if (MS1FeatureMap != null) {
            MS1FeatureMap.IDsummary = IDsummary;
        }
        this.IDsummary.Filename = Filename;
        this.IDsummary.mzXMLFileName=Filename;
        return true;
    }
    
    public void SetPepXMLPath(){
        iProphPepXMLs = new ArrayList<String>();
        String PepXMLPath1 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ1Name() + ".pep.xml");
        String PepXMLPath2 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ2Name() + ".pep.xml");
        String PepXMLPath3 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ3Name() + ".pep.xml");
        iProphPepXMLs.add(PepXMLPath1);
        iProphPepXMLs.add(PepXMLPath2);
        iProphPepXMLs.add(PepXMLPath3);
    }
    
    public void ParsePepXML(DBSearchParam searchPara, LCMSID refID) throws ParserConfigurationException, SAXException, IOException {

        SetPepXMLPath();
        IDsummary = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename), searchPara.DecoyPrefix, searchPara.FastaPath);
        for (String pepxml : iProphPepXMLs) {
            LCMSID pepxmlid = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename), searchPara.DecoyPrefix, searchPara.FastaPath);
            PepXMLParser pepxmlparser = new PepXMLParser(pepxmlid, pepxml, 0f);
            if (refID == null) {
                pepxmlid.FilterByPepDecoyFDR(searchPara.DecoyPrefix, searchPara.PepFDR);
            }
            Logger.getRootLogger().info("No. of peptide ions:" + pepxmlid.GetPepIonList().size() + "; Peptide level threshold: " + pepxmlid.PepProbThreshold);
            for (PepIonID pepID : pepxmlid.GetPepIonList().values()) {
                if (refID != null) {
                    if (refID.GetPepIonList().containsKey(pepID.GetKey())) {
                        IDsummary.AddPeptideID(pepID);
                    }
                } else {
                    IDsummary.AddPeptideID(pepID);
                }
            }
        }
        IDsummary.ReMapProPep();
        Logger.getRootLogger().info("Total number of peptide ions:" + IDsummary.GetPepIonList().size());
        CheckPSMRT();
        if (MS1FeatureMap != null) {
            this.MS1FeatureMap.IDsummary = IDsummary;
        }
    }

    private String GetSkylineFolder(){
        return FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_Skyline/";
    }
    
    public void CreateSkylingImportFolder() {
        new File(GetSkylineFolder()).mkdir();
    }
    
    public boolean RawMGFExist(){
        return (new File(GetSkylineFolder() + GetForLibQ1Name() + ".mgf").exists() && new File(GetSkylineFolder() + GetForLibQ2Name() + ".mgf").exists() && new File(GetSkylineFolder() + GetForLibQ3Name() + ".mgf").exists());            
    }
    
    public void ConvertRawMGF(String msconvertpath){        
        MSConvert mSConvert = new MSConvert(GetSkylineFolder() + GetForLibQ1Name() + ".mgf");
        mSConvert.msconvertpath=msconvertpath;
        mSConvert.Convert();
        mSConvert.SpectrumPath = GetSkylineFolder() + GetForLibQ2Name() + ".mgf";
        mSConvert.Convert();
        mSConvert.SpectrumPath = GetSkylineFolder() + GetForLibQ3Name() + ".mgf";
        mSConvert.Convert();        
    }
    
    //Generate pseudo MS/MS spectra with raw fragment intensities
    public void GenerateRawMGF() throws IOException, Exception {
        
        if(RawMGFExist()){
            return;
        }        
        Logger.getRootLogger().info("Extracting pseudo MS/MS spectra with raw intensity");
        HashMap<Integer, ArrayList<PseudoMSMSProcessing>> ScanList = new HashMap<>();
        HashMap<String, PseudoMSMSProcessing> UnfragScanList = new HashMap<>();
        parameter.BoostComplementaryIon = false;
        ExecutorService executorPool = Executors.newFixedThreadPool(NoCPUs);
        for (LCMSPeakDIAMS2 DIAwindow : DIAWindows) {
            DIAwindow.ReadPeakCluster();
            DIAwindow.ReadPrecursorFragmentClu2Cur();
            DIAwindow.BuildFragmentMS1ranking();
            DIAwindow.FilterByCriteria();
            DIAwindow.BuildFragmentUnfragranking();
            DIAwindow.FilterByCriteriaUnfrag();
            for (PeakCluster ms1cluster : MS1FeatureMap.PeakClusters) {
                if (DIAwindow.DIA_MZ_Range.getX() <= ms1cluster.GetMaxMz() && DIAwindow.DIA_MZ_Range.getY() >= ms1cluster.TargetMz() && DIAwindow.FragmentsClu2Cur.containsKey(ms1cluster.Index)) {
                    DIAwindow.ExtractFragmentForPeakCluser(ms1cluster);   
                    if (DIAwindow.Last_MZ_Range == null || DIAwindow.Last_MZ_Range.getY() < ms1cluster.TargetMz()) {
                        PseudoMSMSProcessing mSMSProcessing = new PseudoMSMSProcessing(ms1cluster, parameter);
                        executorPool.execute(mSMSProcessing);
                        if (!ScanList.containsKey(ms1cluster.Index)) {
                            ScanList.put(ms1cluster.Index, new ArrayList<PseudoMSMSProcessing>());
                        }
                        ScanList.get(ms1cluster.Index).add(mSMSProcessing);
                    }
                }
            }
            
            for (PeakCluster ms2cluster : DIAwindow.PeakClusters) {
                if (DIAwindow.DIA_MZ_Range.getX() <= ms2cluster.TargetMz() && DIAwindow.DIA_MZ_Range.getY() >= ms2cluster.TargetMz() && DIAwindow.UnFragIonClu2Cur.containsKey(ms2cluster.Index)) {
                    DIAwindow.ExtractFragmentForUnfragPeakCluser(ms2cluster);
                    PseudoMSMSProcessing mSMSProcessing = new PseudoMSMSProcessing(ms2cluster, parameter);
                    executorPool.execute(mSMSProcessing);
                    UnfragScanList.put(DIAwindow.WindowID + ";" + ms2cluster.Index, mSMSProcessing);
                }
            }
            DIAwindow.ClearAllPeaks();
            System.gc();
            Logger.getRootLogger().info("(Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)");
        }
        executorPool.shutdown();

        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }
        ReadScanNoMapping();
        String mgffile = GetSkylineFolder() + GetForLibQ1Name() + ".mgf";
        FileWriter mgfWriter = new FileWriter(mgffile, false);

        for (final int ScanNo : new java.util.TreeSet<>(ScanClusterMap_Q1.keySet())) {
            int ClusterIndex = ScanClusterMap_Q1.get(ScanNo);
            XYPointCollection Scan = new XYPointCollection();
            PseudoMSMSProcessing mSMSProcessing = null;
            for (PseudoMSMSProcessing MS2Processing : ScanList.get(ClusterIndex)) {
                mSMSProcessing = MS2Processing;
                for (PrecursorFragmentPairEdge fragmentClusterUnit : MS2Processing.fragments) {
                    Scan.AddPointKeepMaxIfValueExisted(fragmentClusterUnit.FragmentMz, fragmentClusterUnit.Intensity);
                }
            }
            StringBuilder mgfString = new StringBuilder();
            mgfString.append("BEGIN IONS\n");
            mgfString.append("PEPMASS=" + mSMSProcessing.Precursorcluster.TargetMz() + "\n");
            mgfString.append("CHARGE=" + mSMSProcessing.Precursorcluster.Charge + "+\n");
            mgfString.append("RTINSECONDS=" + mSMSProcessing.Precursorcluster.PeakHeightRT[0] * 60f + "\n");
            mgfString.append("TITLE=ClusterIndex:" + mSMSProcessing.Precursorcluster.Index + "\n");
            for (int i = 0; i < Scan.PointCount(); i++) {
                mgfString.append(Scan.Data.get(i).getX()).append(" ").append(Scan.Data.get(i).getY()).append("\n");
            }
            mgfString.append("END IONS\n\n");
            mgfWriter.write(mgfString.toString());
        }
        mgfWriter.close();

        ////////////////////////////////////////////////////////////////////////////////
        String mgffile2 = GetSkylineFolder() + GetForLibQ2Name() + ".mgf";
        FileWriter mgfWriter2 = new FileWriter(mgffile2, false);

        for (final int ScanNo: new java.util.TreeSet<>(ScanClusterMap_Q2.keySet())) {
            int ClusterIndex = ScanClusterMap_Q2.get(ScanNo);
            XYPointCollection Scan = new XYPointCollection();
            PseudoMSMSProcessing mSMSProcessing = null;
            for (PseudoMSMSProcessing MS2Processing : ScanList.get(ClusterIndex)) {
                mSMSProcessing = MS2Processing;
                for (PrecursorFragmentPairEdge fragmentClusterUnit : MS2Processing.fragments) {
                    Scan.AddPointKeepMaxIfValueExisted(fragmentClusterUnit.FragmentMz, fragmentClusterUnit.Intensity);
                }
            }
            StringBuilder mgfString = new StringBuilder();
            mgfString.append("BEGIN IONS\n");
            mgfString.append("PEPMASS=" + mSMSProcessing.Precursorcluster.TargetMz() + "\n");
            mgfString.append("CHARGE=" + mSMSProcessing.Precursorcluster.Charge + "+\n");
            mgfString.append("RTINSECONDS=" + mSMSProcessing.Precursorcluster.PeakHeightRT[0] * 60f + "\n");
            mgfString.append("TITLE=ClusterIndex:" + mSMSProcessing.Precursorcluster.Index + "\n");
            for (int i = 0; i < Scan.PointCount(); i++) {
                mgfString.append(Scan.Data.get(i).getX()).append(" ").append(Scan.Data.get(i).getY()).append("\n");
            }
            mgfString.append("END IONS\n\n");
            mgfWriter2.write(mgfString.toString());
        }

        mgfWriter2.close();

        ////////////////////////////////
        String mgffile3 = GetSkylineFolder() + GetForLibQ3Name() + ".mgf";
        FileWriter mgfWriter3 = new FileWriter(mgffile3, false);
        mzXMLParser Q3mzxml = new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
        Q3mzxml.GetAllScanCollectionByMSLabel(false, true, false, false);
        for (final int ScanNo: new java.util.TreeSet<>(ScanClusterMap_Q3.keySet())) {
            String key = ScanClusterMap_Q3.get(ScanNo);
            XYPointCollection Scan = new XYPointCollection();
            PseudoMSMSProcessing mSMSProcessing = UnfragScanList.get(key);

            for (PrecursorFragmentPairEdge fragmentClusterUnit : mSMSProcessing.fragments) {
                Scan.AddPointKeepMaxIfValueExisted(fragmentClusterUnit.FragmentMz, fragmentClusterUnit.Intensity);
            }

            StringBuilder mgfString = new StringBuilder();
            mgfString.append("BEGIN IONS\n");
            mgfString.append("PEPMASS=" + mSMSProcessing.Precursorcluster.TargetMz() + "\n");
            mgfString.append("CHARGE=" + mSMSProcessing.Precursorcluster.Charge + "+\n");
            mgfString.append("RTINSECONDS=" + mSMSProcessing.Precursorcluster.PeakHeightRT[0] * 60f + "\n");
            mgfString.append("TITLE=ClusterIndex:" + mSMSProcessing.Precursorcluster.Index + "\n");
            for (int i = 0; i < Scan.PointCount(); i++) {
                mgfString.append(Scan.Data.get(i).getX()).append(" ").append(Scan.Data.get(i).getY()).append("\n");
            }
            mgfString.append("END IONS\n\n");
            mgfWriter3.write(mgfString.toString());
        }
        mgfWriter3.close();
    }

    //Check each PSM, if they don't have retention time from pepXML, fill the RT from SpectrumParser file.
    private void FindPSMRT(){        
        try {
            if(!new File(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mzXML").exists()){
                Logger.getRootLogger().warn(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mzXML doesn't exsit." );
                return;
            }
            if(!new File(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mzXML").exists()){
                Logger.getRootLogger().warn(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mzXML doesn't exsit." );
                return;
            }
            if(!new File(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML").exists()){
                Logger.getRootLogger().warn(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML doesn't exsit." );
                return;
            }
            
            mzXMLParser mgfname1 = new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
            mzXMLParser mgfname2 =  new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
            mzXMLParser mgfname3 = new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
            HashMap<String, mzXMLParser> mgfs=new HashMap<>();
            mgfs.put(GetQ1Name(), mgfname1);
            mgfs.put(GetQ2Name(), mgfname2);
            mgfs.put(GetQ3Name(), mgfname3);
            for(PSM psm : IDsummary.PSMList.values()){
                if(psm.RetentionTime==-1f){                    
                    psm.RetentionTime=mgfs.get(psm.GetRawNameString()).GetScanElutionTimeMap().get(psm.ScanNo);
                    psm.NeighborMaxRetentionTime=psm.RetentionTime;
                }
            }
        } catch (Exception ex) {
            Logger.getRootLogger().warn("Exception trying to fill PSM RTs");
        }
    }
   
    public void CheckPSMRT() {
        boolean PSMRT0=false;
        for(PSM psm : IDsummary.PSMList.values()){
            if(psm.RetentionTime==-1f){
                PSMRT0=true;
                break;
            }
        }
        if(PSMRT0){
            FindPSMRT();
        }
    }
        
    public void ClearStructure() {
        MS1FeatureMap = null;
        DIAWindows = null;
        dIA_Setting = null;
        SpectrumParser = null;
    }

    public void BuildStructure() throws SQLException, FileNotFoundException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, DataFormatException {
        
        LoadDIASetting();
        MS1FeatureMap = new LCMSPeakMS1(Filename, NoCPUs);
        MS1FeatureMap.datattype = dIA_Setting.dataType;
        MS1FeatureMap.SetParameter(parameter);
        
        if (IDsummary != null) {
            MS1FeatureMap.IDsummary = IDsummary;
        }
        BuildDIAWindows();
    }

    public boolean MGFgenerated() {
        if (new File(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf").exists()) {
            return true;
        }
        return false;
    }

    //Perform MS1 feature detection
    private void MS1PeakDetection() throws SQLException, InterruptedException, ExecutionException, IOException, ParserConfigurationException, SAXException, FileNotFoundException, Exception {
        //Remove existing pseudo MS/MS MGF files
        RemoveMGF();
                
        MS1FeatureMap = new LCMSPeakMS1(Filename, NoCPUs);
        MS1FeatureMap.datattype = dIA_Setting.dataType;
        MS1FeatureMap.SetParameter(parameter);
        MS1FeatureMap.Resume=Resume;
        //Assign MS1 feature maps
        MS1FeatureMap.SetMS1Windows(dIA_Setting.MS1Windows);
        MS1FeatureMap.CreatePeakFolder();
        MS1FeatureMap.ExportPeakCurveTable = false;
        
        MS1FeatureMap.SetSpectrumParser(GetSpectrumParser());
        Logger.getRootLogger().info("Processing MS1 peak detection");        
        MS1FeatureMap.ExportPeakClusterTable = ExportPrecursorPeak;
        //Start MS1 feature detection
        MS1FeatureMap.PeakClusterDetection();

        Logger.getRootLogger().info("==================================================================================");
    }

    //Remove pseudo MS/MS spectra MGF files
    private void RemoveMGF() {
        String mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mgf";
        String mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mgf";
        String mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf";
        File file = new File(mgffile);
        if (file.exists()) {
            file.delete();
        }
        file = new File(mgffile2);
        if (file.exists()) {
            file.delete();
        }
        file = new File(mgffile3);
        if (file.exists()) {
            file.delete();
        }
        mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mgf.temp";
        mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mgf.temp";
        mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf.temp";
        file = new File(mgffile);
        if (file.exists()) {
            file.delete();
        }
        file = new File(mgffile2);
        if (file.exists()) {
            file.delete();
        }
        file = new File(mgffile3);
        if (file.exists()) {
            file.delete();
        }
    }

    
    public void ExportID() throws SQLException, IOException {
        ExportID("");
    }

    //Export LCMSID instrance to serialization file
    public void ExportID(String tag) throws SQLException, IOException {
        if (IDsummary == null) {
            return;
        }
        IDsummary.WriteLCMSIDSerialization(Filename, tag);
    }

//    private static final TreeSortedMap<String, BufferedWriter> file_handlers = new TreeSortedMap<>();
//
//    private static final org.eclipse.collections.api.block.function.Function<String, BufferedWriter> FO = new org.eclipse.collections.api.block.function.Function<String, BufferedWriter>() {
//        @Override
//        public BufferedWriter valueOf(final String filename) {
//            try {
//                return Files.newBufferedWriter(Paths.get(filename), StandardCharsets.ISO_8859_1);
//            } catch (IOException ex) {
//                throw new RuntimeException(ex);
//            }
//        }
//    };
//
//    public static BufferedWriter get_file(final String filename) {
//        return DIAPack.file_handlers.getIfAbsentPutWithKey(filename, DIAPack.FO);
//    }
    public static enum OutputFile{
        ScanClusterMapping_Q1,
        Mgf_Q1,
        ScanClusterMapping_Q2,
        Mgf_Q2,
        ScanClusterMapping_Q3,
        Mgf_Q3;
    }
    private static final EnumMap<OutputFile, BufferedWriter> file_handlers = new EnumMap<>(OutputFile.class);

    public static BufferedWriter get_file(final OutputFile f,final String filename) {
        final BufferedWriter bw = DIAPack.file_handlers.get(f);
        if (bw == null) {
            final BufferedWriter bw2;
            try {
                bw2 = Files.newBufferedWriter(Paths.get(filename), StandardCharsets.ISO_8859_1);
            } catch (final IOException ex) {
                throw new UncheckedIOException(ex);
            }
            DIAPack.file_handlers.put(f, bw2);
            return bw2;
        }
        return bw;
    }

    //Perform MS2 fragment feature detection
    public void DIAMS2PeakDetection() throws SQLException, IOException, InterruptedException, ExecutionException, FileNotFoundException, Exception {
        int count = 1;
        //CreateSWATHTables();
        for (LCMSPeakDIAMS2 DIAwindow : DIAWindows) {
            Logger.getRootLogger().info("Processing DIA MS2 (mz range):" + DIAwindow.DIA_MZ_Range.getX() + "_" + DIAwindow.DIA_MZ_Range.getY() + "( " + (count++) + "/" + GetSpectrumParser().dIA_Setting.DIAWindows.size() + " )");
            DIAwindow.ExportPeakCurveTable = ExportFragmentPeak;
            DIAwindow.ExportPeakClusterTable = ExportPrecursorPeak;
            DIAwindow.PeakDetectionPFGrouping(MS1FeatureMap);
            DIAwindow.ClearAllPeaks();
            Logger.getRootLogger().info("==================================================================================");
        }
        for(final BufferedWriter bw : DIAPack.file_handlers.values())
            bw.close();
        RenameMGF("");
        convertMGF("");
    }

    private void RenameMGF(String tag) {
        String mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mgf.temp";
        String mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mgf.temp";
        String mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf.temp";
        File file = new File(mgffile);

        if (file.exists()) {
            file.renameTo(new File(file.getAbsolutePath().replace(".mgf.temp", tag + ".mgf")));
        }
        file = new File(mgffile2);
        if (file.exists()) {
            file.renameTo(new File(file.getAbsolutePath().replace(".mgf.temp", tag + ".mgf")));
        }
        file = new File(mgffile3);
        if (file.exists()) {
            file.renameTo(new File(file.getAbsolutePath().replace(".mgf.temp", tag + ".mgf")));
        }
    }

    private void convertMGF(String tag) throws java.security.NoSuchAlgorithmException, IOException {
        String mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + tag + ".mgf";
        String mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + tag + ".mgf";
        String mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + tag + ".mgf";
        {
            final Path file = Paths.get(mgffile);
            if (Files.exists(file)) {
                final String fn = file.getFileName().toString();
                try (final BufferedWriter bw = Files.newBufferedWriter(file.getParent().resolve(fn.substring(0, fn.length() - 4) + ".mzML"), StandardCharsets.US_ASCII)) {
                    MSUmpire.SpectrumParser.mzXMLParser.to_mzML(file, bw);
                }
            }
        }
        {
            final Path file = Paths.get(mgffile2);
            if (Files.exists(file)) {
                final String fn = file.getFileName().toString();
                try (final BufferedWriter bw = Files.newBufferedWriter(file.getParent().resolve(fn.substring(0, fn.length() - 4) + ".mzML"), StandardCharsets.US_ASCII)) {
                    MSUmpire.SpectrumParser.mzXMLParser.to_mzML(file, bw);
                }
            }
        }
        {
            final Path file = Paths.get(mgffile3);
            if (Files.exists(file)) {
                final String fn = file.getFileName().toString();
                try (final BufferedWriter bw = Files.newBufferedWriter(file.getParent().resolve(fn.substring(0, fn.length() - 4) + ".mzML"), StandardCharsets.US_ASCII)) {
                    MSUmpire.SpectrumParser.mzXMLParser.to_mzML(file, bw);
                }
            }
        }


    }

    //Generate mass calibration model
    public void GenerateMassCalibrationRTMap() {
        try {
            MS1FeatureMap.GenerateMassCalibrationRTMap();
        } catch (IOException ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }

        for (LCMSPeakDIAMS2 DIAwindow : DIAWindows) {
            DIAwindow.Masscalibrationfunction = MS1FeatureMap.Masscalibrationfunction;
        }
    }

    public void SaveParams() {
        parameter.WriteParamSerialization(Filename);
    }

    public void SaveDIASetting() {
        dIA_Setting.WriteDIASettingSerialization(Filename);
    }

    public boolean LoadParams() {
        parameter = InstrumentParameter.ReadParametersSerialization(Filename);
        return parameter != null;
    }

    public boolean LoadDIASetting() {
        dIA_Setting = DIA_Setting.ReadDIASettingSerialization(Filename);
        return dIA_Setting != null;
    }

    public InstrumentParameter GetParameter() {
        return parameter;
    }
}
