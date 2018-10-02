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
package MSUmpire.LCMSPeakStructure;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.PSMDataStructure.PSM;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PeakCurve;
import MSUmpire.PeakDataStructure.SortedCurveCollectionApexRT;
import MSUmpire.SpectrumParser.SpectrumParserBase;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.*;
import java.util.logging.Level;
import net.sf.javaml.core.kdtree.KDTree;
import net.sf.javaml.core.kdtree.KeyDuplicateException;
import net.sf.javaml.core.kdtree.KeySizeException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.log4j.Logger;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;


/**
 * Parent peak data structure related to an LCMS map
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class LCMSPeakBase {

    public ArrayList<PeakCluster> PeakClusters = new ArrayList<>(1000);            
    private KDTree PeakClusterMassSearchTree=null;    
    private KDTree PeakCurveSearchTree=null;
    public String ScanCollectionName;
    public String ParentmzXMLName;
    protected SpectrumParserBase SpectrumParser;
    public int MaxNoPeakCluster;
    public int MinNoPeakCluster;
    public int StartCharge;
    public int EndCharge;
    public float MiniIntensity;
    public float SNR;
    private SortedCurveCollectionApexRT PeakCurveListRT = null;
    
    public ArrayList<PeakCurve> UnSortedPeakCurves;        
    public InstrumentParameter parameter;
    public SpectralDataType.DataType datattype;
    public boolean Resume = true;
    public boolean ExportPeakClusterTable=true;
    public boolean SaveSerializationFile=true;
    public boolean ExportPeakCurveTable=false;
    int NoCPUs = Runtime.getRuntime().availableProcessors() - 2;
    public PolynomialSplineFunction Masscalibrationfunction;

    public void ClearAllPeaks(){
        BaseClearAllPeaks();
    }
    public void BaseClearAllPeaks() {        
        PeakClusters = null;
        PeakCurveListRT = null;
        UnSortedPeakCurves=null;
        PeakCurveSearchTree=null;
        PeakClusterMassSearchTree=null;
    }

    public void ClearRawPeaks() {
        for (PeakCurve peakCurve : UnSortedPeakCurves) {
//            peakCurve.GetPeakList().clear();
            peakCurve.nullifyPeakList();
        }
    }
     
    //Get all peak clusters given PSMs related to a peptide ion
    public ArrayList<PeakCluster> FindAllPeakClustersForPepByPSM(PepIonID pep) {
        ArrayList<PeakCluster> allclusterList = new ArrayList<>();
        for (PSM psm : pep.GetPSMList()) {
            ArrayList<PeakCluster> clusterList = FindPeakClustersByMassRTTol(psm.ObserPrecursorMass, psm.Charge, psm.RetentionTime, parameter.MaxCurveRTRange / 2);
            for (PeakCluster cluster : clusterList) {
                if (!allclusterList.contains(cluster)) {
                    allclusterList.add(cluster);
                }
            }
        }
        return allclusterList;
    }

    //Get mass error given a RT according to mass calibration model
    private float GetMassError(float RT) {
        if (RT > Masscalibrationfunction.getKnots()[Masscalibrationfunction.getN()]) {
            RT = (float) Masscalibrationfunction.getKnots()[Masscalibrationfunction.getN()];
        }
        if (RT < Masscalibrationfunction.getKnots()[Masscalibrationfunction.getN()]) {
            RT = (float) Masscalibrationfunction.getKnots()[0];
        }
        return (float) Masscalibrationfunction.value(RT);
    }

    public ArrayList<PeakCluster> FindAllPeakClustersForMappedPep(PepIonID pep) {
        ArrayList<PeakCluster> allclusterList = new ArrayList<>();        
        float idrt=pep.GetRT();
        if (idrt != -1) {
            float calibratedmass = InstrumentParameter.GetMzByPPM(pep.CalcNeutralPepMass(), 1, -GetMassError(idrt));
            ArrayList<PeakCluster> clusterList = FindPeakClustersByMassIDTime(calibratedmass, pep.Charge, idrt, parameter.RTtol);
            return clusterList;
        }
        for (float rt : pep.PredictRT) {
            float calibratedmass = InstrumentParameter.GetMzByPPM(pep.CalcNeutralPepMass(), 1, -GetMassError(rt));
            float rtrange=parameter.RT_window_Targeted;            
            if (rtrange == -1) {
                rtrange = Math.min(5, Math.max(pep.GetRTSD() * 2, parameter.RTtol));
            }
            ArrayList<PeakCluster> clusterList = FindPeakClustersByMassRTTol(calibratedmass, pep.Charge, rt, rtrange);
            for (PeakCluster cluster : clusterList) {
                if (!cluster.Identified && !allclusterList.contains(cluster)) {
                    allclusterList.add(cluster);
                }
            }
        }
        return allclusterList;
    }
    
    public ArrayList<PeakCluster> FindPeakClustersByMassIDTime(float mass, int charge, float RT, float RTtol) {
        ArrayList<PeakCluster> ReturnList = new ArrayList<>();
        float lowrt = RT - RTtol;
        float highrt = RT + RTtol;
        float lowmass = InstrumentParameter.GetMzByPPM(mass, 1, parameter.MS1PPM);
        float highmass = InstrumentParameter.GetMzByPPM(mass, 1, -parameter.MS1PPM);
       
        Object[] found=null;
        try {
            found = GetPeakClusterMassSearchTree().range(new double[]{lowrt,lowmass}, new double[]{highrt,highmass});
        } catch (KeySizeException ex) {
            
        }
        if(found==null || found.length==0){
            return ReturnList;
        }
        for (Object obj : found) {
            PeakCluster candidateCluster = (PeakCluster) obj;            
            if (candidateCluster.Charge == charge && Math.abs(candidateCluster.PeakHeightRT[0] - RT) <= RTtol) {                
                ReturnList.add(candidateCluster);
            }
        }
        return ReturnList;
    }

    public ArrayList<PeakCluster> FindPeakClustersByMassRTTol(float mass, int charge, float RT, float RTtol) {
        
        ArrayList<PeakCluster> Clusters = new ArrayList<>();
        ArrayList<PeakCluster> TightRangeClusters = new ArrayList<>();
        
        float lowrt = RT - RTtol;
        float highrt = RT + RTtol;
        float lowmass = InstrumentParameter.GetMzByPPM(mass, 1, parameter.MS1PPM);
        float highmass = InstrumentParameter.GetMzByPPM(mass, 1, -parameter.MS1PPM);
       
        Object[] found=null;
        try {
            found = GetPeakClusterMassSearchTree().range(new double[]{lowrt,lowmass}, new double[]{highrt,highmass});
        } catch (KeySizeException ex) {
            
        }
        if(found==null || found.length==0){
            return Clusters;
        }
        for (Object obj : found) {
            PeakCluster candidateCluster = (PeakCluster) obj;
            if (candidateCluster.Charge == charge && candidateCluster.startRT - RTtol <= RT && candidateCluster.endRT + RTtol >= RT) {
                Clusters.add(candidateCluster);
            }
            if (candidateCluster.Charge == charge && candidateCluster.startRT <= RT && candidateCluster.endRT >= RT) {
                TightRangeClusters.add(candidateCluster);
            }
        }
        if (!TightRangeClusters.isEmpty()) {
            return TightRangeClusters;
        }
        return Clusters;
    }

    public ArrayList<PeakCluster> FindPeakClustersByMassRTRange(float mass, int charge, float startrt, float endrt) {
        ArrayList<PeakCluster> Clusters = new ArrayList<>();
        float lowrt = startrt;
        float highrt = endrt;
        float lowmass = InstrumentParameter.GetMzByPPM(mass, 1, parameter.MS1PPM);
        float highmass = InstrumentParameter.GetMzByPPM(mass, 1, -parameter.MS1PPM);

        Object[] found = null;
        try {
            found = GetPeakClusterMassSearchTree().range(new double[]{lowrt, lowmass}, new double[]{highrt, highmass});
        } catch (KeySizeException ex) {

        }
        if (found == null || found.length == 0) {
            return Clusters;
        }
        for (Object obj : found) {
            PeakCluster candidateCluster = (PeakCluster) obj;

            if (candidateCluster.Charge == charge && candidateCluster.PeakHeightRT[0] >= startrt && candidateCluster.PeakHeightRT[0] <= endrt) {
                Clusters.add(candidateCluster);
            }

        }
        return Clusters;
    }

    public void ExportPeakClusterResultCSV() throws IOException {

        Logger.getRootLogger().info("Writing PeakCluster CSV:" + FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.csv");

        String mzstring = "";
        String Idxstring = "";
        String CorrString = "";
        String SNRString = "";
        String PeakheightString = "";
        String PeakheightRTString = "";
        String PeakAreaString = "";
        String IdentifiedString = "0";

        for (int i = 0; i < MaxNoPeakCluster; i++) {
            mzstring += ",mz" + (i + 1);
            Idxstring += ",PeakIdx" + (i + 1);
            if (i > 0) {
                CorrString += ",Corr" + (i + 1);
            }
            SNRString += ",SNR" + (i + 1);
            PeakheightString += ",PeakHeight" + (i + 1);
            PeakheightRTString += ",PeakHeightRT" + (i + 1);
            PeakAreaString += ",PeakArea" + (i + 1);
        }

        writer.write("Cluster_Index,StartRT,EndRT,StartScan,EndScan,Identified,Charge" + mzstring + Idxstring + CorrString + SNRString + PeakheightString + PeakheightRTString + PeakAreaString  + ",IsoMapProb,ConflictCorr,LeftInt,RightInt,NoRidges,MS1Score,MS1Prob,MS1LProb\n");

        for (PeakCluster cluster : PeakClusters) {
            IdentifiedString = "0";
            if (cluster.Identified) {
                IdentifiedString = "1";
            }

            String statementString = cluster.Index + "," + cluster.startRT + "," + cluster.endRT + "," + cluster.StartScan + "," + cluster.EndScan + "," + IdentifiedString + "," + cluster.Charge + ",";

            mzstring = "";
            Idxstring = "";
            CorrString = "";
            SNRString = "";
            PeakheightString = "";
            PeakheightRTString = "";
            PeakAreaString = "";

            for (int i = 0; i < MaxNoPeakCluster; i++) {
                mzstring += cluster.mz[i] + ",";
                Idxstring += cluster.IsoPeakIndex[i] + ",";
                if (i > 0) {
                    CorrString += cluster.Corrs[i - 1] + ",";
                }
                SNRString += cluster.GetSNR(i) + ",";
                PeakheightString += cluster.PeakHeight[i] + ",";
                PeakheightRTString += cluster.PeakHeightRT[i] + ",";
                PeakAreaString += cluster.PeakArea[i] + ",";
            }
            statementString += mzstring + Idxstring + CorrString + SNRString + PeakheightString + PeakheightRTString + PeakAreaString  + cluster.IsoMapProb + "," + cluster.GetConflictCorr() + "," + cluster.LeftInt + "," + cluster.RightInt + "," + cluster.NoRidges + "," + cluster.MS1Score + "," + cluster.MS1ScoreProbability + "," + cluster.MS1ScoreLocalProb + "\n";
            writer.write(statementString);
        }
        writer.close();
        //System.out.print("Finished multithreading\n");
    }
    
    public void ExportPeakCluster() throws IOException {
        WritePeakClusterSerialization();      
    }

    //Default path to store peak data serialization files  
     public void CreatePeakFolder(){
        new File(FilenameUtils.getFullPath(ParentmzXMLName)+FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/").mkdir();        
    }
    
    public void WritePeakClusterSerialization() {        
        FS_PeakClusterWrite();
    }

    private void FS_PeakClusterWrite() {
        try {
            Logger.getRootLogger().info("Writing PeakCluster serialization to file:" +  FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS", false);
            FSTObjectOutput out = new FSTObjectOutput(fout);
            out.writeObject(PeakClusters);
            out.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            JavaSerializationPeakClusterWrite();
        }
    }
    
    private void JavaSerializationPeakClusterWrite() {
        try {
            Logger.getRootLogger().info("Writing PeakCluster serialization to file:" +  FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.ser...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.ser", false);
            ObjectOutputStream oos = new ObjectOutputStream(fout);
            oos.writeObject(PeakClusters);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

     public boolean ReadPeakCluster() {
        return ReadPeakClusterSerialization();
    }
    private boolean ReadPeakClusterSerialization() {
        if(!FS_PeakClusterRead()){
            if (JavaSerializationPeakClusterRead()) {
                return true;
            }
            return false;
        }
        for(PeakCluster cluster : PeakClusters){
            cluster.CreateLock();
        }
        return true;        
    }

    private boolean JavaSerializationPeakClusterRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.ser").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().info("Reading PeakCluster serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.ser...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            PeakClusters = (ArrayList<PeakCluster>) in.readObject();
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    private boolean FS_PeakClusterRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().info("Reading PeakCluster serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName)+ "_PeakCluster.serFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            PeakClusters = (ArrayList<PeakCluster>) in.readObject();
            in.close();
            fileIn.close();            
        } catch (Exception ex) {            
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            if(FS_PeakClusterRead_Old()){
                WritePeakClusterSerialization();                
                return true;
            }
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        } 
        return true;
    }
    
    private boolean FS_PeakClusterRead_Old() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().info("Old PeakCluster serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS");
            org.nustaq.serialization.FSTObjectInput in = new org.nustaq.serialization.FSTObjectInput(fileIn);
            PeakClusters = (ArrayList<PeakCluster>) in.readObject();
            in.close();
            fileIn.close();            
        } catch (Exception ex) {
            Logger.getRootLogger().error("Old version reader still failed.");
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        } 
        return true;
    }
    
    public void ExportPeakClusterRegionTXT() throws IOException {
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakClusterRange.txt");
        for (PeakCluster cluster : PeakClusters) {
            writer.write(cluster.startRT + "\t" + cluster.endRT + "\t" + (cluster.mz[0] + 1f) + "\t" + (cluster.mz[0]) + "\n");
        }
        writer.close();
    }

    public void ClearMonoisotopicPeakOfCluster() {
        for (PeakCluster peak : PeakClusters) {
            peak.IsoPeaksCurves = null;
            peak.MonoIsotopePeak = null;
        }
    }
    
    public KDTree GetPeakClusterMassSearchTree() {
        if (PeakClusterMassSearchTree == null) {
            Logger.getRootLogger().info("Building PeakCluster Mass-RT KD tree");
            PeakClusterMassSearchTree = new KDTree(2);
            for (int i = 0; i < PeakClusters.size(); i++) {
                PeakCluster peakCluster = PeakClusters.get(i);
                try {
                    PeakClusterMassSearchTree.insert(new double[]{peakCluster.PeakHeightRT[0], peakCluster.NeutralMass()}, peakCluster);
                } catch (KeyDuplicateException | KeySizeException ex) {
                    Logger.getRootLogger().error(ex.getMessage());
                }
            }
        }
        return PeakClusterMassSearchTree;
    }
       
    public KDTree GetPeakCurveSearchTree(){
        if(PeakCurveSearchTree==null){
            Logger.getRootLogger().info("Building PeakCurve Mass-RT KD tree");
            PeakCurveSearchTree=new KDTree(2);
            for(PeakCurve peakCurve : UnSortedPeakCurves){
                try {
                    PeakCurveSearchTree.insert(new double[]{peakCurve.ApexRT,peakCurve.TargetMz}, peakCurve);
                } catch (KeySizeException ex) {
                    java.util.logging.Logger.getLogger(LCMSPeakBase.class.getName()).log(Level.SEVERE, null, ex);
                } catch (KeyDuplicateException ex) {
                    java.util.logging.Logger.getLogger(LCMSPeakBase.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        return PeakCurveSearchTree;
    }

    public SortedCurveCollectionApexRT GetPeakCurveListRT() {
        if (PeakCurveListRT == null) {
            PeakCurveListRT=new SortedCurveCollectionApexRT();
            PeakCurveListRT.addAll(UnSortedPeakCurves);
            PeakCurveListRT.Finalize();
        }
        return PeakCurveListRT;
    }
}
