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
package MSUmpire.DIA;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.MathPackage.PiecewiseRegression;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.BaseDataStructure.XYZData;
import MSUmpire.FragmentLib.FragmentLibManager;
import MSUmpire.PSMDataStructure.PepFragmentLib;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.geom.Ellipse2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
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
 * Generate peptide ions from external library
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class RTMappingExtLib implements Runnable{

    private PiecewiseRegression regression;
    private LCMSID TargetLCMS;
    private FragmentLibManager libManager;
    InstrumentParameter parameter;

    public RTMappingExtLib(LCMSID TargetLCMS, FragmentLibManager libManager, InstrumentParameter parameter) {
       this.TargetLCMS=TargetLCMS;
       this.parameter=parameter;
       this.libManager=libManager;
    }

    public void GenerateModel() throws IOException {
        
        XYPointCollection points = new XYPointCollection();
        XYSeries series = new XYSeries("Peptide ions");
        XYSeriesCollection xySeriesCollection = new XYSeriesCollection();

        FileWriter writer= new FileWriter(FilenameUtils.getFullPath(TargetLCMS.mzXMLFileName) + "/"+FilenameUtils.getBaseName(TargetLCMS.mzXMLFileName)+"_"+libManager.LibID+"_RTMapPoints.txt");        
        
        for (String pepkey : libManager.PeptideFragmentLib.keySet()) {
            if (TargetLCMS.GetPepIonList().containsKey(pepkey)) {
                PepFragmentLib peplib = libManager.GetFragmentLib(pepkey);
                for (float rt : peplib.RetentionTime) {
                    float y = TargetLCMS.GetPepIonList().get(pepkey).GetRT();
                    points.AddPoint(rt, y);
                    series.add(new XYDataItem(rt, y));
                    writer.write(rt + "\t" + y + "\n");
                }
            }
        }        
        writer.close();
        regression = new PiecewiseRegression(parameter.MaxCurveRTRange,parameter.MaxCurveRTRange);
        regression.SetData(points);
        float R2 = regression.GetR2();
        LogManager.getRootLogger().info("Retention time prediction model:(" + FilenameUtils.getBaseName(TargetLCMS.mzXMLFileName) + "..R2=" + R2 + "(No. of commonly identified peptide ions="+points.PointCount()+")");
        
        GenerateRTMapPNG(xySeriesCollection, series, R2);
    }

    private void GenerateRTMapPNG(XYSeriesCollection xySeriesCollection, XYSeries series, float R2) throws IOException {        
        String pngfile = FilenameUtils.getFullPath(TargetLCMS.mzXMLFileName) + "/"+FilenameUtils.getBaseName(TargetLCMS.mzXMLFileName)+"_"+libManager.LibID+"_RTMap.png";        
        FileWriter writer= new FileWriter(FilenameUtils.getFullPath(TargetLCMS.mzXMLFileName) + "/"+FilenameUtils.getBaseName(TargetLCMS.mzXMLFileName)+"_"+libManager.LibID+"_RTMap.txt");        
                
        XYSeries smoothline = new XYSeries("RT fitting curve");
        for (XYZData data : regression.PredictYList) {
            smoothline.add(data.getX(), data.getY());
            writer.write(data.getX()+"\t"+data.getY()+"\n");
        }
        writer.close();
        xySeriesCollection.addSeries(smoothline);
        xySeriesCollection.addSeries(series);
        JFreeChart chart = ChartFactory.createScatterPlot("Retention time mapping: R2=" + R2, "Normalized RT ("+libManager.LibID+")", "RT:" + FilenameUtils.getBaseName(TargetLCMS.mzXMLFileName), xySeriesCollection,
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
        LogManager.getRootLogger().info("Mapping peptide ions for " + FilenameUtils.getBaseName(TargetLCMS.mzXMLFileName) + "...");

        if(!regression.valid()){
            return;
        }
        int cont=0;
        for (String pepkey : libManager.PeptideFragmentLib.keySet()) {
            PepFragmentLib peplib = libManager.GetFragmentLib(pepkey);
            PepIonID predictedPepIon = null;
            if (!TargetLCMS.GetPepIonList().containsKey(pepkey)) {
                if (TargetLCMS.GetMappedPepIonList().containsKey(pepkey)) {
                    predictedPepIon = TargetLCMS.GetMappedPepIonList().get(pepkey);
                } else {
                    predictedPepIon = new PepIonID();
                    predictedPepIon.ModSequence = peplib.ModSequence;
                    predictedPepIon.Sequence = peplib.Sequence;
                    predictedPepIon.Modifications = (ArrayList<ModificationMatch>) peplib.Modifications.clone();
                    predictedPepIon.Charge = peplib.Charge;
                    TargetLCMS.GetMappedPepIonList().put(peplib.GetKey(), predictedPepIon);
                    cont++;
                }
            } else {
                predictedPepIon = TargetLCMS.GetPepIonList().get(pepkey);
            }
            for (float librt : peplib.RetentionTime) {
                XYZData predict=regression.GetPredictTimeSDYByTimelist(librt);
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
        LogManager.getRootLogger().info("No. of peptide ions added:"+cont);
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
