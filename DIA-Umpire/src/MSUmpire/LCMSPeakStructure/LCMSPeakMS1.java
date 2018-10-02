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
import MSUmpire.BaseDataStructure.DBSearchParam;
import MSUmpire.MathPackage.ChiSquareGOF;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PSM;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeptidePeakClusterDetection.PDHandlerMS1;
import MSUmpire.SearchResultParser.PepXMLParser;
import MSUmpire.SpectrumParser.SpectrumParserBase;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.geom.Ellipse2D;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.log4j.Logger;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.xml.sax.SAXException;

/**
 * MS1 peak feature data related to a LCMS1 feature map
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class LCMSPeakMS1 extends LCMSPeakBase {

    public LCMSID IDsummary;
    public boolean UseMappedIon = false;
    public TreeMap<XYData, ArrayList<Integer>> MS1Windows;

     public LCMSPeakMS1(String Filename, InstrumentParameter parameter, int NoCPUs) {
        this.ScanCollectionName = Filename;
        this.ParentmzXMLName = Filename;
        this.datattype = SpectralDataType.DataType.DDA;
        this.NoCPUs = NoCPUs;
        SetParameter(parameter);
    }
    
    public LCMSPeakMS1(String Filename, int NoCPUs) {
        this.ScanCollectionName = Filename;
        this.ParentmzXMLName=Filename;
        this.datattype=SpectralDataType.DataType.DDA;        
        this.NoCPUs = NoCPUs;
    }
        
    public void SetMS1Windows(TreeMap<XYData, ArrayList<Integer>> MS1Windows){
        this.MS1Windows=MS1Windows;
    }

    public void SetParameter(InstrumentParameter parameter){
        this.parameter = parameter;
        this.MaxNoPeakCluster = parameter.MaxNoPeakCluster;
        this.MinNoPeakCluster = parameter.MinNoPeakCluster;
        ChiSquareGOF.GetInstance(this.MaxNoPeakCluster);
        this.StartCharge = parameter.StartCharge;
        this.EndCharge = parameter.EndCharge;
        this.MiniIntensity = parameter.MinMSIntensity;
        this.SNR = parameter.SNThreshold;
    }

    public void GeneratePeptideFragmentByMSMS(DBSearchParam para) throws InterruptedException, ExecutionException, IOException, SQLException {
        IDsummary.GenerateFragmentPeakForPepIonByMSMS(GetSpectrumParser().GetAllScanCollectionByMSLabel(false, true, false, true), para.FragPPM);
}

    public void WriteLCMSIDSerialization(){
        WriteLCMSIDSerialization("");
    }
    public void WriteLCMSIDSerialization(String tag){
        this.IDsummary.WriteLCMSIDSerialization(ParentmzXMLName,tag);
    }
    public boolean ReadSerializedLCMSID() throws Exception{
        return ReadSerializedLCMSID("");
    }
    public boolean ReadSerializedLCMSID(String tag) throws Exception{
        this.IDsummary=LCMSID.ReadLCMSIDSerialization(ParentmzXMLName,tag);            
        if(this.IDsummary==null){
            return false;
        }
        this.IDsummary.mzXMLFileName=ParentmzXMLName;
        this.IDsummary.Filename=ParentmzXMLName;
        return true;
    }
       
    public String GetDefaultPepXML() {
        return FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + "interact-" + FilenameUtils.getBaseName(ScanCollectionName) + ".pep.xml");
    }

    public String GetDefaultProtXML() {
        return GetDefaultPepXML().replace(".pep.xml",".prot.xml");
    }

    public void SetDefaultPath(DBSearchParam para) {
        para.PepXMLPath = GetDefaultPepXML();
        para.ProtXMLPath = GetDefaultProtXML();
    }
    
     public void AssignMS1Cluster() {
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            for (String indexint : pepIonID.MS1ClusIndex.split(";")) {
                if (!"".equals(indexint)) {
                    pepIonID.MS1PeakClusters.add(PeakClusters.get(Integer.parseInt(indexint) - 1));
                }
            }
        }
        if (IDsummary.GetMappedPepIonList() != null) {
            for (PepIonID pepIonID : IDsummary.GetMappedPepIonList().values()) {
                for (String indexint : pepIonID.MS1ClusIndex.split(";")) {
                    if (!"".equals(indexint)) {
                        pepIonID.MS1PeakClusters.add(PeakClusters.get(Integer.parseInt(indexint) - 1));
                    }
                }
            }
        }
    }
  
    public void SetSpectrumParser(SpectrumParserBase mzxml) {
        this.SpectrumParser = mzxml;
    }

    public SpectrumParserBase GetSpectrumParser() {
        if (SpectrumParser == null) {
            try {
                //SpectrmParser = new mzXMLParser(ScanCollectionName, parameter, datattype, null, NoCPUs);
                SpectrumParser = SpectrumParserBase.GetInstance(ScanCollectionName, parameter, datattype, null, NoCPUs);
            } catch (Exception ex) {
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            }
        }
        return SpectrumParser;
    }

    //MS1 feature detection
    public void PeakClusterDetection() throws IOException, InterruptedException, ExecutionException {

        if (Resume && ReadIfProcessed()) {
            return;
        }
        CreatePeakFolder();
        ArrayList<ScanCollection> scanCollections = new ArrayList<>();
        //Calculate how many point per minute when we do B-spline peak smoothing
        parameter.NoPeakPerMin = (int) (parameter.SmoothFactor / GetSpectrumParser().GetMS1CycleTime());
        Logger.getRootLogger().info("MS1 average cycle time : "+GetSpectrumParser().GetMS1CycleTime()*60+ " seconds");
        
        if (MS1Windows == null || MS1Windows.isEmpty()) {
            //The data has only one MS1 scan set
            ScanCollection scanCollection=GetSpectrumParser().GetAllScanCollectionByMSLabel(true, true, true, false,parameter.startRT, parameter.endRT);
            scanCollections.add(scanCollection);            
        }
        else{
            //Get MS1 ScanCollection for each MS1 scan set
            for( XYData window : MS1Windows.keySet()){
                scanCollections.add(GetSpectrumParser().GetScanCollectionMS1Window(window, true));
            }
        }
        
        PDHandlerMS1 detection = new PDHandlerMS1(this, NoCPUs, parameter.MS1PPM);        
        //Set the flag to tell the program to do targeted peak detection.
        if(parameter.TargetIDOnly){
            detection.SetTargetedDetectionOnly();
        }
        
        if (IDsummary != null) {
            //If identifications are present, assign mz-RT for each PSM into targeted peak detection list.
            for (PSM psm : IDsummary.PSMList.values()) {
                for (int i = 0; i < parameter.MinNoPeakCluster; i++) {
                    detection.AddToInclusionList(psm.GetObsrIsotopicMz(i), psm.RetentionTime);
                }
            }
            Logger.getRootLogger().info("No. of targeted PSM IDs:"+IDsummary.PSMList.size());
        }
        //Start detect MS1 peak curve and isotope clusters
        detection.DetectPeakClusters(scanCollections);
        
        //Release scan collection data structure
        scanCollections.clear();
        scanCollections=null;
        
        //Find closet scan number for each detected isotope peak cluster
        MapScanNoForPeakClusters();
        
        //Export peak cluster CSV table
        if (ExportPeakClusterTable) {
            ExportPeakClusterResultCSV();
        }
        
        //Export peak cluster serialization file
        if (SaveSerializationFile) {
            ExportPeakCluster();
        }
    }

    public void MapScanNoForPeakClusters() {
        for (PeakCluster peakCluster : PeakClusters) {
            peakCluster.StartScan = GetSpectrumParser().GetScanNoByRT(peakCluster.startRT);
            peakCluster.EndScan = GetSpectrumParser().GetScanNoByRT(peakCluster.endRT);
        }
    }
  
    public void AssignQuant() throws IOException{
        AssignQuant(true);
    }

    public void AssignQuant(boolean export) throws IOException {
        if (IDsummary == null) {
            return;
        }
        if (!FilenameUtils.getBaseName(IDsummary.mzXMLFileName).equals(FilenameUtils.getBaseName(ScanCollectionName))) {
            return;
        }
        Logger.getRootLogger().info("Assigning peak cluster to peptide IDs......");
        Logger.getRootLogger().info("Total peptide ions to be quantified..:"+IDsummary.GetPepIonList().size());
        int count=0;
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            pepIonID.CreateQuantInstance(MaxNoPeakCluster);
            pepIonID.MS1PeakClusters.clear();
            ArrayList<PeakCluster> clusterList = FindAllPeakClustersForPepByPSM(pepIonID);
            if (!clusterList.isEmpty()) {
                count++;
                PeakCluster targetCluster = clusterList.get(0);
                float Intensity = targetCluster.PeakHeight[0];

                for (int i = 1; i < clusterList.size(); i++) {
                    if (clusterList.get(i).PeakHeight[0] > Intensity) {
                        targetCluster = clusterList.get(i);
                        Intensity = clusterList.get(i).PeakHeight[0];
                    }
                }
                pepIonID.PeakArea = targetCluster.PeakArea;
                pepIonID.PeakHeight = targetCluster.PeakHeight;
                if (!pepIonID.MS1PeakClusters.contains(targetCluster)) {
                    pepIonID.MS1PeakClusters.add(targetCluster);
                }
                pepIonID.PeakClusterScore = targetCluster.MS1Score;
                pepIonID.PeakRT = targetCluster.PeakHeightRT[0];
                pepIonID.ObservedMz = targetCluster.mz[0];
                targetCluster.AssignedPepIon = pepIonID.GetKey();
                targetCluster.Identified = true;
            }
            else {
                Logger.getRootLogger().warn("Feature for "+pepIonID.ModSequence+ " not found: Charge="+pepIonID.Charge+", mz="+pepIonID.ObservedMz+", RT="+pepIonID.GetRT());
            }
        }
        DecimalFormat df = new DecimalFormat("###.##");
              
        Logger.getRootLogger().info("No. of peptide ions quantified..:"+count+"("+ df.format((float)count/IDsummary.GetPepIonList().size()*100)+"%)");

        if (export) {
            ExportID();
        }
    }

    //Perform MS1 quant for mapped peptides
    public void AssignMappedPepQuant() throws SQLException, IOException {
        if (IDsummary == null || IDsummary.GetMappedPepIonList().isEmpty()) {
            return;
        }
        if (!FilenameUtils.getBaseName(IDsummary.mzXMLFileName).equals(FilenameUtils.getBaseName(ScanCollectionName))) {
            return;
        }
        Logger.getRootLogger().info("Assigning peak cluster to mapped peptide IDs......");
        for (PepIonID pepIonID : IDsummary.GetMappedPepIonList().values()) {
            pepIonID.CreateQuantInstance(MaxNoPeakCluster);

            ArrayList<PeakCluster> clusterList = FindAllPeakClustersForMappedPep(pepIonID);

            if (!clusterList.isEmpty()) {
                PeakCluster targetCluster = null;
                float Score = 0f;

                for (int i = 0; i < clusterList.size(); i++) {
                    PeakCluster cluster = clusterList.get(i);
                    if ("".equals(cluster.AssignedPepIon)) {
                        Score = cluster.PeakHeight[0] * cluster.MS1Score;
                        if (targetCluster == null || clusterList.get(i).MS1Score > Score) {
                            targetCluster = cluster;
                            Score = cluster.MS1Score;
                        }
                    }
                }
                if (targetCluster != null) {
                    pepIonID.PeakArea = targetCluster.PeakArea;
                    pepIonID.PeakHeight = targetCluster.PeakHeight;
                    pepIonID.MS1PeakClusters.add(targetCluster);
                    pepIonID.PeakClusterScore = targetCluster.MS1Score;
                    pepIonID.PeakRT = targetCluster.PeakHeightRT[0];
                    pepIonID.ObservedMz = targetCluster.mz[0];
                    targetCluster.AssignedPepIon = pepIonID.GetKey();
                }
            }
        }
    }

    public boolean ReadIfProcessed(){
        return ReadPeakCluster();    
    }
    
    public void ExportMappedIDQuant() throws SQLException, IOException {
        if (IDsummary == null) {
            return;
        }
        if (!FilenameUtils.getBaseName(IDsummary.mzXMLFileName).equals(FilenameUtils.getBaseName(ScanCollectionName))) {
            return;
        }
        IDsummary.ExportMappedPepID();
    }

    public void ExportID()throws IOException{
        ExportID("");
    }
            
    public void ExportID(String tag) throws IOException {
        if (IDsummary == null) {
            return;
        }
        if (!FilenameUtils.getBaseName(IDsummary.mzXMLFileName).equals(FilenameUtils.getBaseName(ScanCollectionName))) {
            return;
        }
        IDsummary.WriteLCMSIDSerialization(ParentmzXMLName,tag);
    }
    
    public void GenerateMassCalibrationRTMap() throws IOException {
        String pngfile = FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_masscaliRT.png";
        XYSeries series = new XYSeries("PSM");
        XYSeriesCollection xySeriesCollection = new XYSeriesCollection();
        LoessInterpolator loessInterpolator = new LoessInterpolator(
                0.75,//bandwidth,
                2//robustnessIters
        );

        for (PSM psm : this.IDsummary.PSMList.values()) {
            float ppm = InstrumentParameter.CalcSignedPPM(psm.ObserPrecursorMass, psm.NeutralPepMass);            
            series.add(new XYDataItem(psm.RetentionTime, ppm));
        }
        double x[] = new double[IDsummary.PSMList.size()];
        double y[] = new double[x.length];
        double currentmin = 0f;
        for (int i = 0; i < series.getItemCount(); i++) {
            x[i] = (double) series.getX(i);
            if (x[i] <= currentmin) {
                x[i] = currentmin + 0.0001f;
            }
            currentmin = x[i];
            y[i] = (double) series.getY(i);
        }

        Masscalibrationfunction = loessInterpolator.interpolate(x, y);
        XYSeries smoothline = new XYSeries("Loess Regression");

        double xvalue = x[0];

        while (xvalue < x[x.length-1]) {
            smoothline.add(xvalue, Masscalibrationfunction.value(xvalue));
            xvalue += 0.05d;
        }
        xySeriesCollection.addSeries(smoothline);
        xySeriesCollection.addSeries(series);

        JFreeChart chart = ChartFactory.createScatterPlot("Mass calibration", "RT", "Mass error (ppm)", xySeriesCollection,
                PlotOrientation.VERTICAL, true, true, false);
        XYPlot xyPlot = (XYPlot) chart.getPlot();
        xyPlot.setDomainCrosshairVisible(true);
        xyPlot.setRangeCrosshairVisible(true);

        XYItemRenderer renderer = xyPlot.getRenderer();
        renderer.setSeriesPaint(1, Color.blue);
        renderer.setSeriesPaint(0, Color.BLACK);
        renderer.setSeriesShape(1, new Ellipse2D.Double(0, 0, 3, 3));
        renderer.setSeriesStroke(1, new BasicStroke(1.0f));
        xyPlot.setBackgroundPaint(Color.white);
        ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 1000, 600);
    }
       
    public String GetIPROPHETPepXML(){
        return FilenameUtils.getFullPath(ParentmzXMLName) + "interact-"+FilenameUtils.getBaseName(ParentmzXMLName) + ".iproph.pep.xml";
    }
    
    public String GetIPROPHETProtXML(){
        return FilenameUtils.getFullPath(ParentmzXMLName) + "interact-"+FilenameUtils.getBaseName(ParentmzXMLName) + ".iproph.prot.xml";
    }
    
    public void ParsePepXML(DBSearchParam param, float prob) throws ParserConfigurationException, IOException, SAXException {
        IDsummary = new LCMSID(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName),param.DecoyPrefix,param.FastaPath);                
        PepXMLParser pepxmlparser = new PepXMLParser(IDsummary, param.InteractPepXMLPath, prob);
        CheckPSMRT();
    }
  
    private void CheckPSMRT() {
        for(PSM psm : IDsummary.PSMList.values()){
            if(psm.RetentionTime==-1f){
                psm.RetentionTime=GetSpectrumParser().GetScanElutionTimeMap().get(psm.ScanNo);
                psm.NeighborMaxRetentionTime=psm.RetentionTime;
            }
        }
    }

    public void RemoveContaminantPeaks(float proportion) {
        Logger.getRootLogger().info("Removing peak clusters whose m/z appear more than " +proportion*100+ "% chromatography. No. of peak clusters : "+PeakClusters.size());
        float minmz=Float.MAX_VALUE;
        float maxmz=Float.MIN_VALUE;
        float minrt=Float.MAX_VALUE;
        float maxrt=Float.MIN_VALUE;
        for(PeakCluster peak : PeakClusters){
            if(peak.TargetMz()>maxmz){
                maxmz=peak.TargetMz();
            }
            if(peak.TargetMz()<minmz){
                minmz=peak.TargetMz();
            }
            if(peak.endRT>maxrt){
                maxrt=peak.endRT;
            }
            if(peak.startRT<minrt){
                minrt=peak.startRT;
            }
        }
        HashMap<Integer,ArrayList<PeakCluster>> map=new HashMap<>();
        
        float[] MzBin=new float[(int)Math.ceil((maxmz-minmz)*10)+1];       
        for(PeakCluster peak : PeakClusters){
            int binkey=(int)Math.ceil((peak.TargetMz()-minmz)*10);
            MzBin[binkey]+=peak.endRT-peak.startRT;
            if(!map.containsKey(binkey)){
                map.put(binkey, new ArrayList<PeakCluster>());
            }
            map.get(binkey).add(peak);
        }
        float threshold=proportion*(maxrt-minrt);
        for(int i=0;i<MzBin.length;i++){
            if(MzBin[i]>threshold){
                for (PeakCluster peakCluster : map.get(i)){ 
                    PeakClusters.remove(peakCluster);
                    //Logger.getRootLogger().debug("Removing the cluster m/z: "+ peakCluster.TargetMz()+", StartRT: "+ peakCluster.startRT +", EndRT: "+peakCluster.endRT);
                }
            }
        }
        Logger.getRootLogger().info("Remaining peak clusters : "+PeakClusters.size());
    }
}
