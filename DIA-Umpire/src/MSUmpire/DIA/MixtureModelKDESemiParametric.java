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

import MSUmpire.MathPackage.KMeans;
import ExternalPackages.jMEF.ExpectationMaximization1D;
import ExternalPackages.jMEF.MixtureModel;
import ExternalPackages.jMEF.PVector;
import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Vector;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import umontreal.ssj.probdist.*;

/**
 * This class provides semi-parametric mixture modeling
 */
public class MixtureModelKDESemiParametric {

    private EmpiricalDist TargetEmpiricalDist;
    private EmpiricalDist DecoyEmpiricalDist;
    private EmpiricalDist IDEmpiricalDist;
    double bandWidth;
    NormalDist kern = new NormalDist();
    public float[][] MixtureModelProb;
    MixtureModel mmc;
    public int NoBinPoints = 1000;
    float max;
    float min;
    double[] p;
    double[] f1;
    double[] f0;
    double weight_incorrect;

    double weight_correct;
    double pisum = 0d;
    int MAX_ITERATIONS = 50;
    double[] model_kde_x;
    double[] model_kde_y;
    double[] decoy_kde_y;
    double[] correct_kde_y;
    double[] inicorrect_kde_y;

    public void GeneratePlot(String pngfile) throws IOException {
        String modelfile = FilenameUtils.getFullPath(pngfile) + "/" + FilenameUtils.getBaseName(pngfile)+ "_ModelPoints.txt";
        FileWriter writer=new FileWriter(modelfile);
        
        double[] IDObs = new double[IDEmpiricalDist.getN()];
        double[] DecoyObs = new double[DecoyEmpiricalDist.getN()];

        for (int i = 0; i < IDEmpiricalDist.getN(); i++) {
            IDObs[i] = IDEmpiricalDist.getObs(i);
        }
        for (int i = 0; i < DecoyEmpiricalDist.getN(); i++) {
            DecoyObs[i] = DecoyEmpiricalDist.getObs(i);
        }

        XYSeries model1 = new XYSeries("Incorrect matches");
        XYSeries model2 = new XYSeries("Correct matches");
        XYSeries model3 = new XYSeries("All target hits");

        writer.write("UScore\tModel\tCorrect\tDecoy\n");
        for (int i = 0; i < NoBinPoints; i++) {
            model1.add(model_kde_x[i], decoy_kde_y[i]);
            model2.add(model_kde_x[i], correct_kde_y[i]);
            model3.add(model_kde_x[i], model_kde_y[i]);
            writer.write(model_kde_x[i]+"\t"+model_kde_y[i]+"\t"+correct_kde_y[i]+"\t"+decoy_kde_y[i]+"\n");
        }
        writer.close();
        
        MixtureModelProb = new float[NoBinPoints + 1][3];
        float positiveaccu = 0f;
        float negativeaccu = 0f;

        MixtureModelProb[0][0] = (float) model2.getMaxX() + Float.MIN_VALUE;
        MixtureModelProb[0][1] = 1f;
        MixtureModelProb[0][2] = 1f;

        for (int i = 1; i < NoBinPoints + 1; i++) {
            double positiveNumber = correct_kde_y[NoBinPoints-i];
            double negativeNumber = decoy_kde_y[NoBinPoints-i];
            MixtureModelProb[i][0] = (float) model_kde_x[NoBinPoints-i];
            positiveaccu += positiveNumber;
            negativeaccu += negativeNumber;
            MixtureModelProb[i][2] = 0.999999f * (float) (positiveNumber / (negativeNumber + positiveNumber));
            MixtureModelProb[i][1] = 0.999999f * (float) (positiveaccu / (negativeaccu + positiveaccu));
            }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(model1);
        dataset.addSeries(model2);
        dataset.addSeries(model3);

        HistogramDataset histogramDataset = new HistogramDataset();
        histogramDataset.setType(HistogramType.SCALE_AREA_TO_1);
        histogramDataset.addSeries("ID hits", IDObs, 100);
        histogramDataset.addSeries("Decoy hits", DecoyObs, 100);
        //histogramDataset.addSeries("Model hits", ModelObs, 100);

        JFreeChart chart = ChartFactory.createHistogram(FilenameUtils.getBaseName(pngfile), "Score", "Hits", histogramDataset, PlotOrientation.VERTICAL, true, false, false);
        XYPlot plot = chart.getXYPlot();

        NumberAxis domain = (NumberAxis) plot.getDomainAxis();
        domain.setRange(min, max);
        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);
        plot.setForegroundAlpha(0.8f);
        chart.setBackgroundPaint(Color.white);

        XYLineAndShapeRenderer render = new XYLineAndShapeRenderer();

        plot.setDataset(1, dataset);
        plot.setRenderer(1, render);
        plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);
        try {
            ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 1000, 600);
        } catch (IOException e) {
        }
    }

    public void Modeling() {
        GenerateInitialModel();
        p = new double[TargetEmpiricalDist.getN()];
        f1 = new double[TargetEmpiricalDist.getN()];
        f0 = new double[TargetEmpiricalDist.getN()];
        //double[] plast = new double[TargetEmpiricalDist.getN()];
        double miniIDscore=IDEmpiricalDist.getMean()-IDEmpiricalDist.getStandardDeviation()*3;
        double meandecoy=DecoyEmpiricalDist.getSampleMean();
        if (meandecoy > miniIDscore) {
            miniIDscore = meandecoy;
        }
        
        double targetlowidcdf=TargetEmpiricalDist.cdf(miniIDscore);
        double decoylowidcdf=DecoyEmpiricalDist.cdf(miniIDscore);
        weight_incorrect = targetlowidcdf/ decoylowidcdf;
        weight_correct = 1-weight_incorrect;
        
        //initialization
        GenerateIniCorrectDensity();
        GenerateDecoyDensity();
        pisum = 0d;
        for (int i = 0; i < TargetEmpiricalDist.getN(); i++) {
            f1[i] = IniCorrectKDELookUp(TargetEmpiricalDist.getObs(i));
            f0[i] = DecoyKDELookUp(TargetEmpiricalDist.getObs(i));
            //p[i] = f1[i] / TargetKDELookUp(TargetEmpiricalDist.getObs(i));
            p[i] = f1[i] / (f1[i]+f0[i]);
            pisum += p[i];
        }

        ///EM
        GenerateCorrectDensity();
        double likelihood = 0d;
        double correctdensity = 0f;
        double decoydensity = 0f;
        for (int i = 0; i < TargetEmpiricalDist.getN(); i++) {
            decoydensity = DecoyKDELookUp(TargetEmpiricalDist.getObs(i));
            correctdensity = CorrectKDELookUp(TargetEmpiricalDist.getObs(i));
            likelihood += Math.log(decoydensity + correctdensity);            
        }

        int iterations = 0;
        double logLikelihoodNew = likelihood;
        double logLikelihoodThreshold = Math.abs(logLikelihoodNew) * 0.00001;
        double logLikelihoodOld;

        Logger.getRootLogger().debug("EM mixture modeling");
        do {

            logLikelihoodOld = logLikelihoodNew;
            
            // Update of iterations and log likelihood value
            Update();
            iterations++;
            logLikelihoodNew=0;
            for (int i = 0; i < TargetEmpiricalDist.getN(); i++) {
                decoydensity = DecoyKDELookUp(TargetEmpiricalDist.getObs(i));
                correctdensity = CorrectKDELookUp(TargetEmpiricalDist.getObs(i));
                logLikelihoodNew += Math.log(decoydensity + correctdensity);    
            }
            // Display
            Logger.getRootLogger().debug("Iterations:"+iterations + " LogLikelihood :"+logLikelihoodNew+ " total prob="+pisum);
            //System.out.printf("%2d : %12.6f\n", iterations, logLikelihoodNew);
        } while (Math.abs(logLikelihoodNew - logLikelihoodOld) > logLikelihoodThreshold && iterations < MAX_ITERATIONS);
        weight_correct = pisum / TargetEmpiricalDist.getN();
        weight_incorrect = 1 - weight_correct;
        Update();
        GenerateDecoyDensity();
    }

    private void Update() {
        pisum=0d;
        for (int i = 0; i < TargetEmpiricalDist.getN(); i++) {
            f1[i] = CorrectKDELookUp(TargetEmpiricalDist.getObs(i));
            f0[i] = DecoyKDELookUp(TargetEmpiricalDist.getObs(i));
            p[i] = f1[i] / (f1[i]+f0[i]);
            pisum += p[i];
        }
        GenerateCorrectDensity(); 
    }

    public void SetData(double[] targetdata, double[] decoydata, double[] iddata) {
        Arrays.sort(targetdata);
        Arrays.sort(decoydata);
        Arrays.sort(iddata);
        min = (float) Math.min(targetdata[0], decoydata[0]);
        max = (float) Math.max(targetdata[targetdata.length - 1], iddata[iddata.length - 1]);
        TargetEmpiricalDist = new EmpiricalDist(targetdata);
        DecoyEmpiricalDist = new EmpiricalDist(decoydata);
        IDEmpiricalDist = new EmpiricalDist(iddata);
        //Silverman's ‘rule of thumb’ (Scott Variation uses factor = 1.06)
        bandWidth = 0.99 * Math.min(TargetEmpiricalDist.getSampleStandardDeviation(), (TargetEmpiricalDist.getInterQuartileRange() / 1.34)) / Math.pow(TargetEmpiricalDist.getN(), 0.2);
        
        float intv = (max - min) / NoBinPoints;
        //bandWidth=intv*5;
        model_kde_x = new double[NoBinPoints];
        for (int i = 0; i < NoBinPoints; i++) {
            model_kde_x[i] = min + i * intv;
        }
        GenerateTargetDensity();        
    }

    private void GenerateIniCorrectDensity() {
        inicorrect_kde_y = new double[NoBinPoints];
        for (int i = 0; i < NoBinPoints; i++) {
            PVector point = new PVector(2);
            point.array[0] = model_kde_x[i];
            inicorrect_kde_y[i] = mmc.EF.density(point, mmc.param[1]) * weight_correct;
        }
    }

    private void GenerateCorrectDensity() {
        correct_kde_y = new double[NoBinPoints];
        for (int i = 0; i < NoBinPoints; i++) {
            correct_kde_y[i] = weight_correct * ProbBasedKDE(TargetEmpiricalDist, model_kde_x[i]);
        }
    }

    private void GenerateDecoyDensity() {
        decoy_kde_y = new double[NoBinPoints];
        for (int i = 0; i < NoBinPoints; i++) {
            decoy_kde_y[i] = weight_incorrect * DecoyKDE(model_kde_x[i]);
        }
    }

    private void GenerateTargetDensity() {
        model_kde_y = new double[NoBinPoints];
        for (int i = 0; i < NoBinPoints; i++) {
            model_kde_y[i] = TargetKDE(model_kde_x[i]);
        }
    }

    private void GenerateInitialModel() {
        PVector[] points = new PVector[TargetEmpiricalDist.getN()];
        PVector[] centroids = new PVector[2];

        for (int i = 0; i < TargetEmpiricalDist.getN(); i++) {
            points[i] = new PVector(1);
            points[i].array[0] = TargetEmpiricalDist.getObs(i);
        }

        centroids[0] = new PVector(1);
        centroids[0].array[0] = DecoyEmpiricalDist.getMedian();
        centroids[1] = new PVector(1);
        centroids[1].array[0] = IDEmpiricalDist.getMedian();
        Vector<PVector>[] clusters = KMeans.run(points, 2, centroids);
        MixtureModel mm = ExpectationMaximization1D.initialize(clusters);
        mmc = ExpectationMaximization1D.run(points, mm);
        DecimalFormat df = new DecimalFormat("#.####");
        Logger.getRootLogger().debug("----------------------------------------------------------------------------------------");
        Logger.getRootLogger().debug("No. of modeling points=" + TargetEmpiricalDist.getN());
        Logger.getRootLogger().debug("ID hits mean=" + df.format(IDEmpiricalDist.getMean()));
        Logger.getRootLogger().debug("Decoy hits mean=" + df.format(DecoyEmpiricalDist.getMean()));
        Logger.getRootLogger().debug("Initial model of incorrect dist., mean=" + df.format(((PVector) mmc.param[0]).array[0]) + " variance=" + df.format(((PVector) mmc.param[0]).array[1]) + " weight=" + df.format(mmc.weight[0]));
        Logger.getRootLogger().debug("Initial model of correct dist., mean=" + df.format(((PVector) mmc.param[1]).array[0]) + " variance=" + df.format(((PVector) mmc.param[1]).array[1]) + " weight=" + df.format(mmc.weight[1]));
    }

    private double IniCorrectKDELookUp(double x) {
        if (x <= model_kde_x[0]) {
            return inicorrect_kde_y[0];
        }
        for (int i = 0; i < NoBinPoints - 1; i++) {
            if (x >= model_kde_x[i] && x < model_kde_x[i + 1]) {
                return inicorrect_kde_y[i];
            }
        }
        return inicorrect_kde_y[NoBinPoints - 1];
    }

    private double CorrectKDELookUp(double x) {
        if (x <= model_kde_x[0]) {
            return correct_kde_y[0];
        }
        for (int i = 0; i < NoBinPoints - 1; i++) {
            if (x >= model_kde_x[i] && x < model_kde_x[i + 1]) {
                return correct_kde_y[i];
            }
        }
        return correct_kde_y[NoBinPoints - 1];
    }

    private double DecoyKDELookUp(double x) {
        if (x <= model_kde_x[0]) {
            return decoy_kde_y[0];
        }
        for (int i = 0; i < NoBinPoints - 1; i++) {
            if (x >= model_kde_x[i] && x < model_kde_x[i + 1]) {
                return decoy_kde_y[i];
            }
        }
        return decoy_kde_y[NoBinPoints - 1];
    }

    private double ProbBasedKDE(EmpiricalDist dist, double x) {
        // Computes and returns the kernel density estimate at $y$, where the 
        // kernel is the density kern.density(x), and the bandwidth is $h$.
        double z;
        double a = kern.getXinf();       // lower limit of density
        double b = kern.getXsup();       // upper limit of density
        double sum = 0;
        int n = dist.getN();
        for (int i = 0; i < n; i++) {
            z = (x - dist.getObs(i)) / bandWidth;
            if ((z >= a) && (z <= b)) {
                sum += kern.density(z) * p[i];
            }
        }
        sum /= (bandWidth * pisum);
        return sum;
    }

    private double KDE(EmpiricalDist dist, double y) {
        // Computes and returns the kernel density estimate at $y$, where the 
        // kernel is the density kern.density(x), and the bandwidth is $h$.
        double z;
        double a = kern.getXinf();       // lower limit of density
        double b = kern.getXsup();       // upper limit of density
        double sum = 0;
        int n = dist.getN();
        for (int i = 0; i < n; i++) {
            z = (y - dist.getObs(i)) / bandWidth;
            if ((z >= a) && (z <= b)) {
                sum += kern.density(z);
            }
        }
        sum /= (bandWidth * n);
        return sum;
    }

    private double TargetKDE(double x) {
        return KDE(TargetEmpiricalDist, x);
    }

    private double DecoyKDE(double x) {
        return KDE(DecoyEmpiricalDist, x);
    }    
}
