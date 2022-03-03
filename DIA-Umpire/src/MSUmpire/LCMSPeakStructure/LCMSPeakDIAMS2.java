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
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.DIA.DIAPack;
import MSUmpire.DIA.PseudoMSMSProcessing;
import MSUmpire.MathPackage.MassDefect;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PeakCurve;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.PeptidePeakClusterDetection.PDHandlerDIAMS2;
import MSUmpire.SpectrumParser.SpectrumParserBase;
import MSUmpire.SpectrumParser.mzXMLParser;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.eclipse.collections.impl.list.mutable.primitive.FloatArrayList;
import org.eclipse.collections.impl.map.mutable.primitive.IntObjectHashMap;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;

/**
 * MS2 peak feature map related to a DIA MS2 isolation window 
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class LCMSPeakDIAMS2 extends LCMSPeakBase {

    public XYData DIA_MZ_Range;
    public XYData Last_MZ_Range;
    public String WindowID;
    public HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> FragmentsClu2Cur;
    public HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> UnFragIonClu2Cur;
    public HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> MatchedFragmentMap = new HashMap<>();
    private DIAPack parentDIA;

    public boolean FragmentGroupByCluster = false;

    public LCMSPeakDIAMS2(String Filename, DIAPack parentDIA, InstrumentParameter parameter, XYData WindowMZ, XYData LastWindowMZ, SpectrumParserBase spectrumparser, int NoCPUs) {
        this.DIA_MZ_Range = WindowMZ;
        this.Last_MZ_Range = LastWindowMZ;
        this.WindowID = (int) Math.floor(WindowMZ.getX()) + "_" + (int) Math.floor(WindowMZ.getY());
        this.SpectrumParser = spectrumparser;
        this.ScanCollectionName = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + (int) Math.floor(WindowMZ.getX()) + "_" + (int) Math.floor(WindowMZ.getY());
        this.ParentmzXMLName = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename);
        this.parentDIA = parentDIA;
        this.parameter = parameter;
        this.MaxNoPeakCluster = parameter.MaxMS2NoPeakCluster;
        this.MinNoPeakCluster = parameter.MinMS2NoPeakCluster;
        this.StartCharge = parameter.MS2StartCharge;
        this.EndCharge = parameter.MS2EndCharge;
        this.MiniIntensity = parameter.MinMSMSIntensity;
        this.SNR = parameter.MS2SNThreshold;
        this.NoCPUs = NoCPUs;
    }

    public void ClearAllPeaks() {        
        BaseClearAllPeaks();
        FragmentsClu2Cur = null;
        UnFragIonClu2Cur = null;
//        FragmentMS1Ranking=null;
//        FragmentUnfragRanking=null;
        MatchedFragmentMap=null;
        System.gc();        
        //System.out.print("Peak data is released (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)\n");
    }

    public void RemoveFragmentPeakByMassDefect(){
        MassDefect MD=new MassDefect();
     LogManager.getRootLogger().info("Performing mass defect filter on fragment peaks");
     LogManager.getRootLogger().info("No. of fragment peaks: "+UnSortedPeakCurves.size());
        ArrayList<PeakCurve> newlist=new ArrayList<>();
        for (PeakCurve peakCurve : UnSortedPeakCurves) {
            for (int charge = 1; charge <= 2; charge++) {
                float mass = charge * (float)((peakCurve.TargetMz - ElementaryIon.proton.getTheoreticMass()));
                if (MD.InMassDefectRange(mass, parameter.MassDefectOffset)) {
                    newlist.add(peakCurve);
                    break;
                }
            }
        }
        UnSortedPeakCurves=newlist;
        LogManager.getRootLogger().info("No. of remaining fragment peaks: "+UnSortedPeakCurves.size());
    }

    public void PeakDetectionPFGrouping(LCMSPeakMS1 ms1lcms) throws InterruptedException, ExecutionException, IOException, FileNotFoundException, Exception {
        if (!(Resume && ReadIfProcessed())) {
            PDHandlerDIAMS2 swathdetection = new PDHandlerDIAMS2(this, NoCPUs, ms1lcms, parameter.MS2PPM);
            swathdetection.MSlevel = 2;            
            if (datattype != SpectralDataType.DataType.pSMART) {
                swathdetection.DetectPeakCurves(GetScanCollection());
                if (UnSortedPeakCurves.isEmpty()) {
                    LogManager.getRootLogger().info("No peak detected...................");
                    return;
                }
                
                ExportPeakCluster();
                if(parameter.MassDefectFilter){
                    RemoveFragmentPeakByMassDefect();
                }
                swathdetection.FragmentGrouping();
            } else {
                //////pSMART////////////////
                swathdetection.pSMARTGrouping(GetScanCollection());
            }
        }
        GenerateMGF(ms1lcms);
    }

    public String GetQ1Name() {
        return FilenameUtils.getBaseName(ParentmzXMLName) + "_Q1";
    }

    public String GetQ2Name() {
        return FilenameUtils.getBaseName(ParentmzXMLName) + "_Q2";
    }

    public String GetQ3Name() {
        return FilenameUtils.getBaseName(ParentmzXMLName) + "_Q3";
    }

    private void PrepareMGF_MS1Cluster(LCMSPeakMS1 ms1lcms) throws IOException {

        ArrayList<PseudoMSMSProcessing> ScanList = new ArrayList<>();
        ExecutorService executorPool = Executors.newFixedThreadPool(NoCPUs);
        for (PeakCluster ms1cluster : ms1lcms.PeakClusters) {
            final ArrayList<PrecursorFragmentPairEdge> frags = FragmentsClu2Cur.get(ms1cluster.Index);
            if (frags!=null && DIA_MZ_Range.getX() <= ms1cluster.GetMaxMz() && DIA_MZ_Range.getY() >= ms1cluster.TargetMz()) {
//            if (DIA_MZ_Range.getX() <= ms1cluster.GetMaxMz() && DIA_MZ_Range.getY() >= ms1cluster.TargetMz() && FragmentsClu2Cur.containsKey(ms1cluster.Index)) {
//                ArrayList<PrecursorFragmentPairEdge> frags = FragmentsClu2Cur.get(ms1cluster.Index);
                ms1cluster.GroupedFragmentPeaks.addAll(frags);
                if (Last_MZ_Range == null || Last_MZ_Range.getY() < ms1cluster.TargetMz()) {
                    PseudoMSMSProcessing mSMSProcessing = new PseudoMSMSProcessing(ms1cluster, parameter);
                    ScanList.add(mSMSProcessing);
                }
            }
        }

        for (PseudoMSMSProcessing proc : ScanList) {
            executorPool.execute(proc);
        }
        executorPool.shutdown();

        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            LogManager.getRootLogger().info("interrupted..");
        }

        String mgffile = FilenameUtils.getFullPath(ParentmzXMLName) + GetQ1Name() + ".mgf.temp";
        String mgffile2 = FilenameUtils.getFullPath(ParentmzXMLName) + GetQ2Name() + ".mgf.temp";
//        FileWriter mapwriter = new FileWriter(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + ".ScanClusterMapping_Q1", true);
//        FileWriter mapwriter2 = new FileWriter(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + ".ScanClusterMapping_Q2", true);

//        FileWriter mgfWriter = new FileWriter(mgffile, true);
//        FileWriter mgfWriter2 = new FileWriter(mgffile2, true);
        final BufferedWriter
                mapwriter = parameter.Q1 ? DIAPack.get_file(DIAPack.OutputFile.ScanClusterMapping_Q1,
                FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + ".ScanClusterMapping_Q1") : null,
                mapwriter2 = parameter.Q2 ? DIAPack.get_file(DIAPack.OutputFile.ScanClusterMapping_Q2,
                                 FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + ".ScanClusterMapping_Q2") : null,
                mgfWriter = parameter.Q1 ? DIAPack.get_file(DIAPack.OutputFile.Mgf_Q1, mgffile) : null,
                mgfWriter2 = parameter.Q2 ? DIAPack.get_file(DIAPack.OutputFile.Mgf_Q2, mgffile2) : null;
        
        for (PseudoMSMSProcessing mSMSProcessing : ScanList) {
            if (MatchedFragmentMap.size() > 0) {
                mSMSProcessing.RemoveMatchedFrag(MatchedFragmentMap);
            }

            XYPointCollection Scan = mSMSProcessing.GetScan();

            if (Scan != null && Scan.PointCount() > parameter.MinFrag) {
//                StringBuilder mgfString = new StringBuilder();
                
                if (mSMSProcessing.Precursorcluster.IsotopeComplete(3)) {
                    if (parameter.Q1) {
                        final BufferedWriter mgfString = mgfWriter;
                        parentDIA.Q1Scan++;
                        mgfString.append("BEGIN IONS\n");
                        mgfString.append("PEPMASS=").append(String.valueOf(mSMSProcessing.Precursorcluster.TargetMz())).append("\n");
                        mgfString.append("CHARGE=").append(String.valueOf(mSMSProcessing.Precursorcluster.Charge)).append("+\n");
                        mgfString.append("RTINSECONDS=").append(String.valueOf(mSMSProcessing.Precursorcluster.PeakHeightRT[0] * 60f)).append("\n");
                        mgfString.append("TITLE=").append(GetQ1Name()).append(".").append(String.valueOf(parentDIA.Q1Scan)).append(".").append(String.valueOf(parentDIA.Q1Scan)).append(".").append(String.valueOf(mSMSProcessing.Precursorcluster.Charge)).append("\n");
                        for (int i = 0; i < Scan.PointCount(); i++) {
                            mgfString.append(String.valueOf(Scan.Data.get(i).getX())).append(" ").append(String.valueOf(Scan.Data.get(i).getY())).append("\n");
                        }
                        mgfString.append("END IONS\n\n");
                        mapwriter.write(parentDIA.Q1Scan + "_" + mSMSProcessing.Precursorcluster.Index + "\n");
                    }
//                    mgfWriter.write(mgfString.toString());
                    //} else if (mSMSProcessing.Precursorcluster.IsotopeComplete(2)) {
                } else {
                    if (parameter.Q2) {
                        final BufferedWriter mgfString = mgfWriter2;
                        parentDIA.Q2Scan++;
                        mgfString.append("BEGIN IONS\n");
                        mgfString.append("PEPMASS=").append(String.valueOf(mSMSProcessing.Precursorcluster.TargetMz())).append("\n");
                        mgfString.append("CHARGE=").append(String.valueOf(mSMSProcessing.Precursorcluster.Charge)).append("+\n");
                        mgfString.append("RTINSECONDS=").append(String.valueOf(mSMSProcessing.Precursorcluster.PeakHeightRT[0] * 60f)).append("\n");
                        mgfString.append("TITLE=").append(GetQ2Name()).append(".").append(String.valueOf(parentDIA.Q2Scan)).append(".").append(String.valueOf(parentDIA.Q2Scan)).append(".").append(String.valueOf(mSMSProcessing.Precursorcluster.Charge)).append("\n");
                        for (int i = 0; i < Scan.PointCount(); i++) {
                            mgfString.append(String.valueOf(Scan.Data.get(i).getX())).append(" ").append(String.valueOf(Scan.Data.get(i).getY())).append("\n");
                        }
                        mgfString.append("END IONS\n\n");
                        mapwriter2.write(parentDIA.Q2Scan + "_" + mSMSProcessing.Precursorcluster.Index + "\n");
//                    mgfWriter2.write(mgfString.toString());
                    }
                }
            }
            mSMSProcessing.Precursorcluster.GroupedFragmentPeaks.clear();
        }
//        mgfWriter2.close();
//        mgfWriter.close();
//        mapwriter.close();
//        mapwriter2.close();
    }

    private void PrepareMGF_UnfragmentIon() throws IOException {
        String mgffile4 = FilenameUtils.getFullPath(ParentmzXMLName) + GetQ3Name() + ".mgf.temp";
//        FileWriter mgfWriter4 = new FileWriter(mgffile4, true);
        final BufferedWriter mgfWriter4 = parameter.Q3 ? DIAPack.get_file(DIAPack.OutputFile.Mgf_Q3, mgffile4) : null;

//        FileWriter mapwriter3 = new FileWriter(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + ".ScanClusterMapping_Q3", true);
        final BufferedWriter mapwriter3 = parameter.Q3 ? DIAPack.get_file(DIAPack.OutputFile.ScanClusterMapping_Q3,
                FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + ".ScanClusterMapping_Q3") : null;
        
        ArrayList<PseudoMSMSProcessing> ScanList = new ArrayList<>();
        ExecutorService executorPool = Executors.newFixedThreadPool(NoCPUs);
        for (PeakCluster ms2cluster : PeakClusters) {
            ArrayList<PrecursorFragmentPairEdge> frags = UnFragIonClu2Cur.get(ms2cluster.Index);
            if (frags!=null && DIA_MZ_Range.getX() <= ms2cluster.TargetMz() && DIA_MZ_Range.getY() >= ms2cluster.TargetMz()) {
//            if (DIA_MZ_Range.getX() <= ms2cluster.TargetMz() && DIA_MZ_Range.getY() >= ms2cluster.TargetMz() && UnFragIonClu2Cur.containsKey(ms2cluster.Index)) {
//                ArrayList<PrecursorFragmentPairEdge> frags = UnFragIonClu2Cur.get(ms2cluster.Index);
                ms2cluster.GroupedFragmentPeaks.addAll(frags);
                PseudoMSMSProcessing mSMSProcessing = new PseudoMSMSProcessing(ms2cluster, parameter);
                executorPool.execute(mSMSProcessing);
                ScanList.add(mSMSProcessing);
            }
        }
        executorPool.shutdown();
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            LogManager.getRootLogger().info("interrupted..");
        }

        for (PseudoMSMSProcessing mSMSProcessing : ScanList) {
            if (MatchedFragmentMap.size() > 0) {
                mSMSProcessing.RemoveMatchedFrag(MatchedFragmentMap);
            }
            XYPointCollection Scan = mSMSProcessing.GetScan();

            if (Scan != null && Scan.PointCount() > parameter.MinFrag) {
                parentDIA.Q3Scan++;
//                StringBuilder mgfString = new StringBuilder();
//                mgfString.append("BEGIN IONS\n");
//                mgfString.append("PEPMASS=" + mSMSProcessing.Precursorcluster.TargetMz() + "\n");
//                mgfString.append("CHARGE=" + mSMSProcessing.Precursorcluster.Charge + "+\n");
//                mgfString.append("RTINSECONDS=" + mSMSProcessing.Precursorcluster.PeakHeightRT[0] * 60f + "\n");
//                mgfString.append("TITLE=").append(GetQ3Name()).append(".").append(parentDIA.Q3Scan).append(".").append(parentDIA.Q3Scan).append(".").append(mSMSProcessing.Precursorcluster.Charge).append("\n");
//                //mgfString.append("TITLE=" + WindowID + ";ClusterIndex:" + mSMSProcessing.ms2cluster.Index + "\n");
//                //mgfString.append("TITLE=" GetQ3Name() + WindowID + ";ClusterIndex:" + mSMSProcessing.ms2cluster.Index + "\n");
//
//                for (int i = 0; i < Scan.PointCount(); i++) {
//                    mgfString.append(Scan.Data.get(i).getX()).append(" ").append(Scan.Data.get(i).getY()).append("\n");
//                }
//                mgfString.append("END IONS\n\n");
//                mgfWriter4.write(mgfString.toString());
                if (parameter.Q3) {
                    mgfWriter4.append("BEGIN IONS\n")
                            .append("PEPMASS=" + mSMSProcessing.Precursorcluster.TargetMz() + "\n")
                            .append("CHARGE=" + mSMSProcessing.Precursorcluster.Charge + "+\n")
                            .append("RTINSECONDS=" + mSMSProcessing.Precursorcluster.PeakHeightRT[0] * 60f + "\n")
                            .append("TITLE=").append(GetQ3Name()).append(".").append(Integer.toString(parentDIA.Q3Scan)).append(".").append(Integer.toString(parentDIA.Q3Scan)).append(".").append(Integer.toString(mSMSProcessing.Precursorcluster.Charge)).append("\n");
                    //mgfWriter4.append("TITLE=" + WindowID + ";ClusterIndex:" + mSMSProcessing.ms2cluster.Index + "\n");
                    //mgfWriter4.append("TITLE=" GetQ3Name() + WindowID + ";ClusterIndex:" + mSMSProcessing.ms2cluster.Index + "\n");

                    for (int i = 0; i < Scan.PointCount(); i++) {
                        mgfWriter4.append(Float.toString(Scan.Data.get(i).getX())).append(" ").append(Float.toString(Scan.Data.get(i).getY())).append("\n");
                    }
                    mgfWriter4.append("END IONS\n\n");

                    mapwriter3.write(parentDIA.Q3Scan + ";" + WindowID + ";" + mSMSProcessing.Precursorcluster.Index + "\n");
                }
            }
            mSMSProcessing.Precursorcluster.GroupedFragmentPeaks.clear();
        }
//        mgfWriter4.close();
//        mapwriter3.close();
    }
    
    public void GenerateMGF(LCMSPeakMS1 ms1lcms) throws IOException, InterruptedException {
        PrepareMGF_MS1Cluster(ms1lcms);
        PrepareMGF_UnfragmentIon();
    }

    public ScanCollection GetScanCollection() throws InterruptedException, ExecutionException, IOException {
        return SpectrumParser.GetScanCollectionDIAMS2(DIA_MZ_Range, true,parameter.startRT, parameter.endRT);
    }

    public boolean ReadIfProcessed() {

        return (ReadPrecursorFragmentClu2Cur() & ReadPeakCluster());
    }

//    HashMap<Integer, ArrayList<Float>> FragmentMS1Ranking;
//    HashMap<Integer, ArrayList<Float>> FragmentUnfragRanking;
    
    public void FilterByCriteriaUnfrag() {
        
        HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> templist=new HashMap<>();
        
        for (int clusterindex : UnFragIonClu2Cur.keySet()) {
            ArrayList<PrecursorFragmentPairEdge> newlist = new ArrayList<>();
            final FloatArrayList CorrArrayList = new FloatArrayList();
           HashMap<PrecursorFragmentPairEdge,Float> ScoreList=new HashMap<>();
            for (PrecursorFragmentPairEdge fragmentClusterUnit : UnFragIonClu2Cur.get(clusterindex)) {                
                float score=fragmentClusterUnit.Correlation*fragmentClusterUnit.Correlation*(float)Math.log(fragmentClusterUnit.Intensity);
                ScoreList.put(fragmentClusterUnit,score);
                CorrArrayList.add(score);
            }
//            Collections.sort(CorrArrayList);
//            Collections.reverse(CorrArrayList);
            CorrArrayList.trimToSize();
            CorrArrayList.sortThis().reverseThis();

            for (PrecursorFragmentPairEdge fragmentClusterUnit : UnFragIonClu2Cur.get(clusterindex)) {
                int CorrRank = 0;                
                for (int intidx = 0; intidx < CorrArrayList.size(); intidx++) {
                    if (CorrArrayList.get(intidx) <= ScoreList.get(fragmentClusterUnit)) {
                        CorrRank = intidx + 1;
                        break;
                    }
                }
                if (fragmentClusterUnit.Correlation >= parameter.CorrThreshold && CorrRank <= parameter.FragmentRank && fragmentClusterUnit.FragmentMS1Rank <= parameter.PrecursorRank && fragmentClusterUnit.ApexDelta <= parameter.ApexDelta) {
                    newlist.add(fragmentClusterUnit);
                }
            }
            templist.put(clusterindex, newlist);
            
        }
        UnFragIonClu2Cur=templist;
    }
        
    public void FilterByCriteria() {        
        HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> templist=new HashMap<>();        
        for (int clusterindex : FragmentsClu2Cur.keySet()) {
            ArrayList<PrecursorFragmentPairEdge> newlist = new ArrayList<>();
//            ArrayList<Float> CorrArrayList = new ArrayList<>();
            final FloatArrayList CorrArrayList = new FloatArrayList();
            HashMap<PrecursorFragmentPairEdge, Float> ScoreList = new HashMap<>();
            for (PrecursorFragmentPairEdge fragmentClusterUnit : FragmentsClu2Cur.get(clusterindex)) {
                float score = fragmentClusterUnit.Correlation * fragmentClusterUnit.Correlation * (float) Math.log(fragmentClusterUnit.Intensity);
                ScoreList.put(fragmentClusterUnit, score);
                CorrArrayList.add(score);
            }
//            Collections.sort(CorrArrayList);
//            Collections.reverse(CorrArrayList);
            CorrArrayList.trimToSize();
            CorrArrayList.sortThis().reverseThis();

            for (PrecursorFragmentPairEdge fragmentClusterUnit : FragmentsClu2Cur.get(clusterindex)) {
                int CorrRank = 0;
                for (int intidx = 0; intidx < CorrArrayList.size(); intidx++) {
                    if (CorrArrayList.get(intidx) <= ScoreList.get(fragmentClusterUnit)) {
                        CorrRank = intidx + 1;
                        break;
                    }
                }
                if (fragmentClusterUnit.Correlation >= parameter.CorrThreshold && CorrRank <= parameter.FragmentRank && fragmentClusterUnit.FragmentMS1Rank <= parameter.PrecursorRank && fragmentClusterUnit.ApexDelta <= parameter.ApexDelta) {
                    newlist.add(fragmentClusterUnit);
                }
            }
            templist.put(clusterindex, newlist);            
        }
        FragmentsClu2Cur=templist;
    }
    
    //Calculate precursor-fragment MS1 ranking (described in DIA-Umpire paper)
    public void BuildFragmentMS1ranking() {
//        FragmentMS1Ranking = new HashMap<>();
        final IntObjectHashMap<FloatArrayList> FragmentMS1Ranking = new IntObjectHashMap<>();
        for (final ArrayList<PrecursorFragmentPairEdge> e: FragmentsClu2Cur.values()) {
            for (final PrecursorFragmentPairEdge framentClusterUnit : e) {
                FragmentMS1Ranking
                        .getIfAbsentPut(framentClusterUnit.PeakCurveIndexB, new FloatArrayList())
                        .add(framentClusterUnit.Correlation);
            }
        }

        for (final FloatArrayList scorelist : FragmentMS1Ranking.values()){
            scorelist.trimToSize();
            scorelist.sortThis().reverseThis();
        }

        for (final ArrayList<PrecursorFragmentPairEdge> e: FragmentsClu2Cur.values()) {
            for (final PrecursorFragmentPairEdge framentClusterUnit : e) {
                final FloatArrayList scorelist = FragmentMS1Ranking.get(framentClusterUnit.PeakCurveIndexB);
                for (int intidx = 0; intidx < scorelist.size(); intidx++) {
                    if (scorelist.get(intidx) <= framentClusterUnit.Correlation) {
                        framentClusterUnit.FragmentMS1Rank = intidx + 1;
                        framentClusterUnit.FragmentMS1RankScore = (float) framentClusterUnit.FragmentMS1Rank / (float) scorelist.size();
                        break;
                    }
                }
            }
        }
    }
    
    //Calculate precursor-fragment MS2 unfragmented ion ranking (described in DIA-Umpire paper)
    public void BuildFragmentUnfragranking() {
//        FragmentUnfragRanking = new HashMap<>();
        final IntObjectHashMap<FloatArrayList> FragmentUnfragRanking = new IntObjectHashMap();
        for (final ArrayList<PrecursorFragmentPairEdge> e: UnFragIonClu2Cur.values()) {
            for (final PrecursorFragmentPairEdge framentClusterUnit : e) {
                FragmentUnfragRanking
                        .getIfAbsentPut(framentClusterUnit.PeakCurveIndexB, new FloatArrayList())
                        .add(framentClusterUnit.Correlation);
            }
        }
        for (final FloatArrayList scorelist : FragmentUnfragRanking.values()){
            scorelist.trimToSize();
            scorelist.sortThis().reverseThis();
        }

        for (final ArrayList<PrecursorFragmentPairEdge> e: UnFragIonClu2Cur.values()) {
            for (final PrecursorFragmentPairEdge framentClusterUnit : e) {
                final FloatArrayList scorelist = FragmentUnfragRanking.get(framentClusterUnit.PeakCurveIndexB);
                for (int intidx = 0; intidx < scorelist.size(); intidx++) {
                    if (scorelist.get(intidx) <= framentClusterUnit.Correlation) {
                        framentClusterUnit.FragmentMS1Rank = intidx + 1;
                        framentClusterUnit.FragmentMS1RankScore = (float) framentClusterUnit.FragmentMS1Rank / (float) scorelist.size();
                        break;
                    }
                }
            }
        }
    }
    
    public boolean ReadPrecursorFragmentClu2Cur() {
        return ReadCluster2CurveCorrSerialization() & ReadUnfragmentedCluster2CurveCorrSerialization();
    }
  
    public void ExportUnfragmentedClusterCurve() throws IOException, SQLException {
        WriteUnfragmentedCluster2CurveCorrSerialization();
    }
    
    public void WritePrecursorFragmentGrouping() {
        WriteCluster2CurveCorrSerialization();
        WriteUnfragmentedCluster2CurveCorrSerialization();
    }

    private void WriteCluster2CurveCorrSerialization() {
        FSCluster2CurveWrite();
    }

    private void FSCluster2CurveWrite() {
        try {
            LogManager.getRootLogger().debug("Writing PrecursorFragmentCorr serialization to file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.serFS...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.serFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(FragmentsClu2Cur);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    private boolean ReadCluster2CurveCorrSerialization() {        
        if (!FSCluster2CurveRead()) {
            if (JavaSerializationCluster2CurveRead()) {
                FSCluster2CurveWrite();
                return true;
            }
            return false;
        }
        return true;
    }

    private boolean FSCluster2CurveRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.serFS").exists()) {
            return false;
        }
        try {
            LogManager.getRootLogger().debug("Reading PrecursorFragmentCorr serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.serFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.serFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            FragmentsClu2Cur = (HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>>) in.readObject();
            in.close();
            fileIn.close();

        } catch (Exception ex) {            
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }
    private boolean JavaSerializationCluster2CurveRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.ser").exists()) {
            return false;
        }
        try {
            LogManager.getRootLogger().debug("Reading PrecursorFragmentCorr serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.ser...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            FragmentsClu2Cur = (HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>>) in.readObject();
            in.close();
            fileIn.close();

        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    private void WriteUnfragmentedCluster2CurveCorrSerialization() {        
        FSCluster2CurveUnfragWrite();
    }

    private void FSCluster2CurveUnfragWrite() {
        try {
            LogManager.getRootLogger().debug("Writing UnfragPrecursorFragCorr serialization to file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.serFS...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.serFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(UnFragIonClu2Cur);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    private boolean ReadUnfragmentedCluster2CurveCorrSerialization() {
        if (!FSCluster2CurveUnfragRead()) {
            if (JavaSerializationCluster2CurveUnfragRead()) {
                FSCluster2CurveUnfragWrite();
                return true;
            }
            return false;
        }
        return true;
    }

    private boolean FSCluster2CurveUnfragRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.serFS").exists()) {
            return false;
        }
        try {
            LogManager.getRootLogger().debug("Reading UnfragPrecursorFragCorr serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.serFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.serFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            UnFragIonClu2Cur = (HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>>) in.readObject();
            in.close();
            fileIn.close();

        } catch (Exception ex) {            
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    private boolean JavaSerializationCluster2CurveUnfragRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.ser").exists()) {
            return false;
        }
        try {
            LogManager.getRootLogger().debug("Reading UnfragPrecursorFragCorr serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.ser...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            UnFragIonClu2Cur = (HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>>) in.readObject();
            in.close();
            fileIn.close();

        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    public void ExportCluster2CurveCorr() throws IOException {
        WriteCluster2CurveCorrSerialization();
    }


    public void ExtractFragmentForPeakCluser(PeakCluster Cluster) {
        if (FragmentsClu2Cur!=null && FragmentsClu2Cur.containsKey(Cluster.Index)) {
            Cluster.fraglock.writeLock().lock();
            try {
                for (PrecursorFragmentPairEdge fragmentClusterUnit : FragmentsClu2Cur.get(Cluster.Index)) {
                    if (!Cluster.GroupedFragmentPeaks.contains(fragmentClusterUnit)) {
                        Cluster.GroupedFragmentPeaks.add(fragmentClusterUnit);
                    }
                }
            } finally {
                Cluster.fraglock.writeLock().unlock();
            }
        }
    }

    public void ExtractFragmentForUnfragPeakCluser(PeakCluster Cluster) {
        if (UnFragIonClu2Cur!=null && UnFragIonClu2Cur.containsKey(Cluster.Index)) {
            Cluster.fraglock.writeLock().lock();
            try {
                for (PrecursorFragmentPairEdge fragmentClusterUnit : UnFragIonClu2Cur.get(Cluster.Index)) {
                    if (!Cluster.GroupedFragmentPeaks.contains(fragmentClusterUnit)) {
                        Cluster.GroupedFragmentPeaks.add(fragmentClusterUnit);
                    }
                }
            } finally {
                Cluster.fraglock.writeLock().unlock();
            }
        }
    }
}
