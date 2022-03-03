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
import MSUmpire.MathPackage.PiecewiseRegression;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.BaseDataStructure.XYZData;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.geom.Ellipse2D;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import org.apache.commons.io.FilenameUtils;
import org.apache.logging.log4j.LogManager;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Cross assign peptide ions from two LCMS runs
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class RTAlignedPepIonMapping implements Runnable{

    private PiecewiseRegression regression;
    private LCMSID LCMSA;
    private LCMSID LCMSB;
    InstrumentParameter parameter;
    String Workfolder;

    public RTAlignedPepIonMapping(String Workfolder, InstrumentParameter parameter, LCMSID LCMSA, LCMSID LCMSB) {
        this.parameter=parameter;
        this.LCMSA = LCMSA;
        this.LCMSB = LCMSB;
        this.Workfolder=Workfolder;
    }

    public void GenerateModel() throws IOException {
        
        XYPointCollection points = new XYPointCollection();
        XYSeries series = new XYSeries("Peptide ions");
        XYSeriesCollection xySeriesCollection = new XYSeriesCollection();

        for (PepIonID pepA : LCMSA.GetPepIonList().values()) {
            if (LCMSB.GetPepIonList().containsKey(pepA.GetKey())) {
                PepIonID pepB = LCMSB.GetPepIonList().get(pepA.GetKey());
                points.AddPoint(pepA.GetRT(), pepB.GetRT());
                series.add(new XYDataItem(pepA.GetRT(), pepB.GetRT()));
            }
        }
        regression = new PiecewiseRegression(parameter.MaxCurveRTRange,parameter.MaxCurveRTRange);
        regression.SetData(points);
        float R2 = regression.GetR2();
        LogManager.getRootLogger().info("Retention time prediction model:(" + FilenameUtils.getBaseName(LCMSA.mzXMLFileName) + "-" + FilenameUtils.getBaseName(LCMSB.mzXMLFileName) + ")..R2=" + R2 + "(No. of commonly identified peptide ions="+points.PointCount()+")");
        
        GenerateRTMapPNG(xySeriesCollection, series, R2);
    }

    private void GenerateRTMapPNG(XYSeriesCollection xySeriesCollection, XYSeries series, float R2) throws IOException {
        new File(Workfolder+ "/RT_Mapping/").mkdir();
        String pngfile = Workfolder+ "/RT_Mapping/" + FilenameUtils.getBaseName(LCMSA.mzXMLFileName).substring(0,Math.min(120, FilenameUtils.getBaseName(LCMSA.mzXMLFileName).length()-1)) + "_" + FilenameUtils.getBaseName(LCMSB.mzXMLFileName).substring(0,Math.min(120, FilenameUtils.getBaseName(LCMSB.mzXMLFileName).length()-1)) + "_RT.png";
        
        XYSeries smoothline = new XYSeries("RT fitting curve");
        for (XYZData data : regression.PredictYList) {
            smoothline.add(data.getX(), data.getY());
        }
        xySeriesCollection.addSeries(smoothline);
        xySeriesCollection.addSeries(series);
        JFreeChart chart = ChartFactory.createScatterPlot("Retention time mapping: R2=" + R2, "RT:" + FilenameUtils.getBaseName(LCMSA.mzXMLFileName), "RT:" + FilenameUtils.getBaseName(LCMSB.mzXMLFileName), xySeriesCollection,
                PlotOrientation.VERTICAL, true, true, false);
        XYPlot xyPlot = (XYPlot) chart.getPlot();
        xyPlot.setDomainCrosshairVisible(true);
        xyPlot.setRangeCrosshairVisible(true);

        XYItemRenderer renderer = xyPlot.getRenderer();
        renderer.setSeriesPaint(1, Color.blue);
        renderer.setSeriesPaint(0, Color.BLACK);
        renderer.setSeriesShape(1, new Ellipse2D.Double(0, 0, 3, 3));
        renderer.setSeriesStroke(1, new BasicStroke(3.0f));
        renderer.setSeriesStroke(0, new BasicStroke(3.0f));
        xyPlot.setBackgroundPaint(Color.white);        
        ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 1000, 600);
    }

    public void GenerateMappedPepIon() {
        LogManager.getRootLogger().info("Mapping predicted peptide ions for " + FilenameUtils.getBaseName(LCMSB.mzXMLFileName) + "...");

        if(!regression.valid()){
            return;
        }
        
        for (PepIonID pepion : LCMSA.GetPepIonList().values()) {
            PepIonID predictedPepIon = null;
            if (!LCMSB.GetPepIonList().containsKey(pepion.GetKey())) {
                if (LCMSB.GetMappedPepIonList().containsKey(pepion.GetKey())) {
                    predictedPepIon = LCMSB.GetMappedPepIonList().get(pepion.GetKey());
                } else {
                    predictedPepIon = pepion.ClonePepIonID();
                    LCMSB.GetMappedPepIonList().put(pepion.GetKey(), predictedPepIon);
                }
            } else {
                predictedPepIon = LCMSB.GetPepIonList().get(pepion.GetKey());
            }
            XYZData predict=regression.GetPredictTimeSDYByTimelist(pepion.GetIDRT());            
            float PRT = predict.getY();
            boolean added = true;
            for (float rt : predictedPepIon.PredictRT) {
                if (Math.abs(PRT - rt) < 0.1f) {
                    added = false;
                }
            }
            if (added) {
                predictedPepIon.PredictRT.add(PRT);
            }
            predictedPepIon.SetRTSD(predict.getZ());
        }
        
        LogManager.getRootLogger().info("Mapping predicted peptide ions for " + FilenameUtils.getBaseName(LCMSA.mzXMLFileName) + "...");

        for (PepIonID pepion : LCMSB.GetPepIonList().values()) {
            PepIonID predictedPepIon = null;
            if (!LCMSA.GetPepIonList().containsKey(pepion.GetKey())) {                
                if (LCMSA.GetMappedPepIonList().containsKey(pepion.GetKey())) {
                    predictedPepIon = LCMSA.GetMappedPepIonList().get(pepion.GetKey());
                } else {
                    predictedPepIon = pepion.ClonePepIonID();
                    LCMSA.GetMappedPepIonList().put(pepion.GetKey(), predictedPepIon);
                }            
            } else {
                predictedPepIon = LCMSA.GetPepIonList().get(pepion.GetKey());
            }
            XYZData predict=regression.GetPredictTimeSDXByTimelist(pepion.GetIDRT());
            float PRT = predict.getY();
            boolean added = true;
            for (float rt : predictedPepIon.PredictRT) {
                if (Math.abs(PRT - rt) < 0.1f) {
                    added = false;
                }
            }
            if (added) {
                predictedPepIon.PredictRT.add(PRT);
            }
            predictedPepIon.SetRTSD(predict.getZ());
        }        
    }

    public void GenerateMappedPepIonToList(HashMap<String, PepIonID>ListA, HashMap<String, PepIonID> ListB) {
        LogManager.getRootLogger().info("Mapping predicted peptide ions for " + FilenameUtils.getBaseName(LCMSB.mzXMLFileName) + "...");

        if(!regression.valid()){
            return;
        }
        
        for (PepIonID pepion : LCMSA.GetPepIonList().values()) {
            if (!LCMSB.GetPepIonList().containsKey(pepion.GetKey())) {
                PepIonID predictedPepIon = null;
                if (ListB.containsKey(pepion.GetKey())) {
                    predictedPepIon = ListB.get(pepion.GetKey());
                } else {
                    predictedPepIon = pepion.ClonePepIonID();
                    ListB.put(pepion.GetKey(), predictedPepIon);
                }
                //System.out.println(pepion.GetIDRT());
                XYZData predict=regression.GetPredictTimeSDYByTimelist(pepion.GetIDRT());
                float PRT = predict.getY();
                boolean added = true;
                for (float rt : predictedPepIon.PredictRT) {
                    if (Math.abs(PRT - rt) < 0.1f) {
                        added = false;
                    }
                }
                if (added) {
                    predictedPepIon.PredictRT.add(PRT);
                }              
                predictedPepIon.SetRTSD(predict.getZ());
            }
        }
        
        LogManager.getRootLogger().info("Mapping predicted peptide ions for " + FilenameUtils.getBaseName(LCMSA.mzXMLFileName) + "...");

        for (PepIonID pepion : LCMSB.GetPepIonList().values()) {
            if (!LCMSA.GetPepIonList().containsKey(pepion.GetKey())) {
                PepIonID predictedPepIon = null;
                if (ListA.containsKey(pepion.GetKey())) {
                    predictedPepIon = ListA.get(pepion.GetKey());
                } else {
                    predictedPepIon = pepion.ClonePepIonID();
                    ListA.put(pepion.GetKey(), predictedPepIon);
                }
               XYZData predict=regression.GetPredictTimeSDXByTimelist(pepion.GetIDRT());
               float PRT = predict.getY();
                boolean added = true;
                for (float rt : predictedPepIon.PredictRT) {
                    if (Math.abs(PRT - rt) < 0.1f) {
                        added = false;
                    }
                }
                if (added) {
                    predictedPepIon.PredictRT.add(PRT);
                }          
                predictedPepIon.SetRTSD(predict.getZ());
            }
        }        
    }
    
    public void ExportMappedPepIon() throws SQLException, IOException{
        LCMSB.ExportMappedPepID();
        LCMSA.ExportMappedPepID();
    }
    
    @Override
    public void run() {
        try {
            GenerateModel();
            GenerateMappedPepIon();
        } catch (Exception ex) {
            System.out.println(ex.toString());
        }
    }
}
