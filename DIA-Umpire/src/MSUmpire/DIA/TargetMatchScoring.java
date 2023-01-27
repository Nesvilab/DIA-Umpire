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

import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.MathPackage.KMeans;
import MSUmpire.MathPackage.KernelDensityEstimator;
import MSUmpire.MathPackage.Regression;
import ExternalPackages.jMEF.ExpectationMaximization1D;
import ExternalPackages.jMEF.MixtureModel;
import ExternalPackages.jMEF.PVector;
import java.awt.Color;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Vector;
import javastat.multivariate.DiscriminantAnalysis;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;

/**
 * Targeted re-extraction scoring
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class TargetMatchScoring implements Serializable {

    private static final long serialVersionUID = 3845269874265L;
    public ArrayList<UmpireSpecLibMatch> libTargetMatches;
    public ArrayList<UmpireSpecLibMatch> libIDMatches;
    private boolean Terminate = false;
    String Filename;
    String LibID;
    public transient int NoBinPoints = 1000;
    MatchSubscore matchSubscore = new MatchSubscore(0);
    float[][] MixtureModelProb;
    ArrayList<PeakGroupScore> decoyModelingList = new ArrayList<>();    
    private transient boolean UseOldVersion=false;

    public void SetUseOldVersion(){
        UseOldVersion=true;
        matchSubscore=new MatchSubscore(1);
    }

    public TargetMatchScoring(String Filename, String LibID) {
        libTargetMatches = new ArrayList<>();
        libIDMatches = new ArrayList<>();
        this.Filename = Filename;
        this.LibID = LibID;
    }

    public void RemoveLowScoreEntry() {
        if (Terminate) {
            return;
        }
        LogManager.getRootLogger().info("Removing low probability (<0.2) hits: total target entries: " + libTargetMatches.size());
        ArrayList<UmpireSpecLibMatch> newlist = new ArrayList<>();
        for (UmpireSpecLibMatch match : libTargetMatches) {
            if ((match.BestDecoyHit != null && match.BestDecoyHit.MixtureModelLocalProb > 0.2) || (match.BestHit != null && match.BestHit.MixtureModelLocalProb > 0.2)) {
                newlist.add(match);
            }
        }
        libTargetMatches = newlist;
        LogManager.getRootLogger().info("Remaining entries: " + libTargetMatches.size());
        decoyModelingList.clear();
    }

    public void Process() throws IOException {
        //ClearGroupFragments();
        //LibraryMatchWrite();

        if(UseOldVersion){
            EM_LDATraining();
            MixtureModeling();
            AssignMixtureModelProb();
        }
        else {
            EM_LDATraining();
            //MixtureModeling();
            MixtureModelingSemiSupervised();
            //GeneratePrecursorCentrolRanking();
            AssignMixtureModelProb();
        }
        ExportAlignmentResult();
        RemoveLowScoreEntry();
        ClearPepFrag();
        LibraryMatchWrite();
    }

    public void ClearPepFrag() {
        if (Terminate) {
            return;
        }
        for (UmpireSpecLibMatch match : libTargetMatches) {
            match.pepIonID.ClearPepFragFactory();
        }
    }

    public void GeneratePrecursorCentrolRanking() {
        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peakGroupScore : match.TargetHits) {
                peakGroupScore.AddScore(peakGroupScore.UmpireScore);
            }
            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
                peakGroupScore.AddScore(peakGroupScore.UmpireScore);
            }
        }
        for (UmpireSpecLibMatch match : libIDMatches) {
            for (PeakGroupScore peakGroupScore : match.TargetHits) {
                peakGroupScore.AddScore(peakGroupScore.UmpireScore);
            }
            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
                peakGroupScore.AddScore(peakGroupScore.UmpireScore);
            }
        }
        /////////////////////////////////calculate delta
        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peakGroupScore : match.TargetHits) {
                peakGroupScore.PrecursorCentralRank = peakGroupScore.GetScoreRank(peakGroupScore.UmpireScore);
            }
            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
                peakGroupScore.PrecursorCentralRank = peakGroupScore.GetScoreRank(peakGroupScore.UmpireScore);
            }
        }
        for (UmpireSpecLibMatch match : libIDMatches) {
            for (PeakGroupScore peakGroupScore : match.TargetHits) {
                peakGroupScore.PrecursorCentralRank = peakGroupScore.GetScoreRank(peakGroupScore.UmpireScore);
            }
            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
                peakGroupScore.PrecursorCentralRank = peakGroupScore.GetScoreRank(peakGroupScore.UmpireScore);
            }
        }
    }

    public void AssignModelingDecoy() {
        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
                decoyModelingList.add(peakGroupScore);
            }
        }
        for (UmpireSpecLibMatch match : libIDMatches) {
            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
                decoyModelingList.add(peakGroupScore);
            }
        }
    }

    private void ExportAlignmentResult() throws IOException {
        try {
            FileWriter writer = new FileWriter(FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_TargetReExtraction.xls");
            writer.write("Modseq\tPredictRT\tPrecursorMz\tType\tClusterRT\tMS level\tApexDeltaScore\tCorrScore\tSpecDotProduct\tSpecCorrelation\tSpecConstrastAngle\tPPMScore\tRTOverlap\tIntScore\tMS1Corr\tNoB\tNoY\tRTdiff\tPrcursorPPM\tSumCorrScore\tSumCorrPPMScore\tPrecursorIsoPattern\tPrecursorCentralRank\tRank\tU-Score\tProbability\n");
            for (UmpireSpecLibMatch match : libTargetMatches) {
                int Rank = 1;
                for (PeakGroupScore target : match.TargetHits) {
                    writer.write(match.pepIonID.ModSequence + "\t" + match.pepIonID.PredictRTString() + "\t" + match.pepIonID.NeutralPrecursorMz() + "\tTarget\t" + target.PrecursorRT + "\t" + target.MSlevel + "\t" + target.ApexDeltaScore + "\t" + target.AveCorrScore + "\t" + target.SpecDotProduct + "\t" + target.SpecCorrelation + "\t" + target.ContrastAngle + "\t" + target.PPMScore + "\t" + target.RTOverlapScore + "\t" + target.FragIntAvgScore + "\t" + target.PrecursorCorr + "\t" + target.NoMatchB + "\t" + target.NoMatchY + "\t" + target.RTDiff + "\t" + target.PrecursorPPM + "\t" + target.SumCorrScore + "\t" + target.SumCorrPPMScore + "\t" + target.PrecursorIsoPattern + "\t" + target.PrecursorCentralRank + "\t" + (Rank++) + "\t" + target.UmpireScore + "\t" + target.MixtureModelLocalProb + "\n");
                }
                Rank = 1;
                for (PeakGroupScore target : match.DecoyHits) {
                    writer.write(match.pepIonID.ModSequence + "\t" + match.pepIonID.PredictRTString() + "\t" + match.pepIonID.NeutralPrecursorMz() + "\tTarget_Decoy\t" + target.PrecursorRT + "\t" + target.MSlevel + "\t" + target.ApexDeltaScore + "\t" + target.AveCorrScore + "\t" + target.SpecDotProduct + "\t" + target.SpecCorrelation + "\t" + target.ContrastAngle + "\t" + target.PPMScore + "\t" + target.RTOverlapScore + "\t" + target.FragIntAvgScore + "\t" + target.PrecursorCorr + "\t" + target.NoMatchB + "\t" + target.NoMatchY + "\t" + target.RTDiff + "\t" + target.PrecursorPPM + "\t" + target.SumCorrScore + "\t" + target.SumCorrPPMScore + "\t" + target.PrecursorIsoPattern + "\t" + target.PrecursorCentralRank + "\t" + (Rank++) + "\t" + target.UmpireScore + "\t" + target.MixtureModelLocalProb + "\n");
                }
            }
            for (UmpireSpecLibMatch match : libIDMatches) {
                int Rank = 1;
                for (PeakGroupScore target : match.TargetHits) {
                    writer.write(match.pepIonID.ModSequence + "\t" + match.pepIonID.PredictRTString() + "\t" + match.pepIonID.NeutralPrecursorMz() + "\tID\t" + target.PrecursorRT + "\t" + target.MSlevel + "\t" + target.ApexDeltaScore + "\t" + target.AveCorrScore + "\t" + target.SpecDotProduct + "\t" + target.SpecCorrelation + "\t" + target.ContrastAngle + "\t" + target.PPMScore + "\t" + target.RTOverlapScore + "\t" + target.FragIntAvgScore + "\t" + target.PrecursorCorr + "\t" + target.NoMatchB + "\t" + target.NoMatchY + "\t" + target.RTDiff + "\t" + target.PrecursorPPM + "\t" + target.SumCorrScore + "\t" + target.SumCorrPPMScore + "\t" + target.PrecursorIsoPattern + "\t" + target.PrecursorCentralRank + "\t" + (Rank++) + "\t" + target.UmpireScore + "\t" + target.MixtureModelLocalProb + "\n");
                }
                for (PeakGroupScore target : match.DecoyHits) {
                    writer.write(match.pepIonID.ModSequence + "\t" + match.pepIonID.PredictRTString() + "\t" + match.pepIonID.NeutralPrecursorMz() + "\tID_Decoy\t" + target.PrecursorRT + "\t" + target.MSlevel + "\t" + target.ApexDeltaScore + "\t" + target.AveCorrScore + "\t" + target.SpecDotProduct + "\t" + target.SpecCorrelation + "\t" + target.ContrastAngle + "\t" + target.PPMScore + "\t" + target.RTOverlapScore + "\t" + target.FragIntAvgScore + "\t" + target.PrecursorCorr + "\t" + target.NoMatchB + "\t" + target.NoMatchY + "\t" + target.RTDiff + "\t" + target.PrecursorPPM + "\t" + target.SumCorrScore + "\t" + target.SumCorrPPMScore + "\t" + target.PrecursorIsoPattern + "\t" + target.PrecursorCentralRank + "\t" + (Rank++) + "\t" + target.UmpireScore + "\t" + target.MixtureModelLocalProb + "\n");
                }
            }
            writer.close();
        } catch (Exception e) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(e));
        }
    }

    private void CalcUmpireScore() {
        LogManager.getRootLogger().debug("Calculating Umpire Scores");
        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore target : match.TargetHits) {
                AssignUmpireScore(target);
            }
            for (PeakGroupScore target : match.DecoyHits) {
                AssignUmpireScore(target);
            }
            match.Ranking();
        }
        for (UmpireSpecLibMatch match : libIDMatches) {
            for (PeakGroupScore target : match.TargetHits) {
                AssignUmpireScore(target);
            }
            for (PeakGroupScore target : match.DecoyHits) {
                AssignUmpireScore(target);
            }
            match.Ranking();
        }
    }

    public void AssignUmpireScore(PeakGroupScore peakgroup) {
        float[] Subs = matchSubscore.GetSubScoreArray(peakgroup);
        peakgroup.UmpireScore = 0f;
        int idx = 0;
        for (int i = 0; i < Subs.length; i++) {
            if (matchSubscore.Subscores[i].enable) {
                peakgroup.UmpireScore += (float) (matchSubscore.SubSCoeff[idx++] * Subs[i]);
            }
        }
    }

    private void EM_LDATraining() throws IOException {

        LogManager.getRootLogger().info("Linear discriminant analysis using identified peptides and decoys as training set.");
        Regression regression = new Regression();
        XYPointCollection points = new XYPointCollection();
        matchSubscore.InitializeLDACoeff();
        int NoLDAComp = matchSubscore.NoEnableSubscores;

        CalcUmpireScore();

        float LDASimialrity = 0f;
        float StopThreshold = 0.95f;
        int MaxIterations = 50;
        int iteration = 0;
        float samplingratio = 0.5f;

        ArrayList<PeakGroupScore> IDList = new ArrayList<>();
        ArrayList<PeakGroupScore> decoyList = new ArrayList<>();

        for (UmpireSpecLibMatch match : libTargetMatches) {
            if (match.BestDecoyHit != null) {
                decoyList.add(match.BestDecoyHit);
            }
            if (match.BestMS2DecoyHit != null) {
                decoyList.add(match.BestMS2DecoyHit);
            }
        }
        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null) {
                IDList.add(match.BestHit);
            }
        }

        Collections.shuffle(decoyList);
        ArrayList<PeakGroupScore> decoyTList = new ArrayList<>();
        for (int i = 0; i < decoyList.size() / 2; i++) {
            decoyModelingList.add(decoyList.get(i));
        }
        for (int i = decoyList.size() / 2; i < decoyList.size(); i++) {
            decoyTList.add(decoyList.get(i));
        }
        decoyList = decoyTList;

        int targetNo = (int) (IDList.size() * samplingratio);
        int TrainNo = Math.min(targetNo, decoyList.size());

        LogManager.getRootLogger().info("No. of identified peptide ions:"+IDList.size());
        LogManager.getRootLogger().info("No. of decoys:"+decoyList.size());
        if (TrainNo < 5) {
            Terminate = true;
            LogManager.getRootLogger().warn("No. of training data is less than 5, the training process will exit.");
            return;
        }
      
        while (LDASimialrity < StopThreshold && iteration < MaxIterations) {
            Collections.shuffle(decoyList);
            Collections.shuffle(IDList);

            double[][] testdata = new double[NoLDAComp][2 * TrainNo];
            double[] testgroup = new double[2 * TrainNo];

            int idx = 0;
            for (int i = 0; i < TrainNo; i++) {
                PeakGroupScore peakgroup = IDList.get(i);
                double[] enablesubscore = matchSubscore.GetEnableSubScoreArray(peakgroup);
                for (int j = 0; j < enablesubscore.length; j++) {
                    testdata[j][idx] = enablesubscore[j];
                }
                testgroup[idx] = 1;
                idx++;
            }

            for (int i = 0; i < TrainNo; i++) {
                PeakGroupScore peakgroup = decoyList.get(i);
                double[] enablesubscore = matchSubscore.GetEnableSubScoreArray(peakgroup);
                for (int j = 0; j < enablesubscore.length; j++) {
                    testdata[j][idx] = enablesubscore[j];
                }
                testgroup[idx] = 0;
                idx++;
            }
            DiscriminantAnalysis LDA = new DiscriminantAnalysis();
            int[] group = LDA.predictedGroup(testgroup, testdata, testdata);
            boolean modelvalid = true;
            for (int i = 0; i < NoLDAComp; i++) {
                if (Float.isNaN((float) LDA.linearDiscriminants[0][i])) {
                    modelvalid = false;
                    break;
                }
                points.AddPoint((float) matchSubscore.SubSCoeff[i], (float) -LDA.linearDiscriminants[0][i]);
                matchSubscore.SubSCoeff[i] = -LDA.linearDiscriminants[0][i];
            }

            if (!modelvalid) {
                LogManager.getRootLogger().debug("LDA failed at iteration:" + iteration);
                break;
            }
            regression.SetData(points);

            LDASimialrity = regression.GetR2();
            if (Float.isNaN(LDASimialrity)) {
                LDASimialrity = -1f;
            }

            DecimalFormat df = new DecimalFormat("#.####");
            LogManager.getRootLogger().debug("----------------------------------------------------------------------------------------");
            LogManager.getRootLogger().debug("Iteration:" + (iteration++));
            LogManager.getRootLogger().debug("No of target hits:" + targetNo);
            LogManager.getRootLogger().debug("No of decoy hits:" + decoyList.size());
            LogManager.getRootLogger().debug("Training set size:" + TrainNo * 2);
            LogManager.getRootLogger().debug("Trained weights:");
            for (int i = 0; i < NoLDAComp; i++) {
                LogManager.getRootLogger().debug(matchSubscore.SubSName[i] + ":" + df.format(matchSubscore.SubSCoeff[i]));
            }
            LogManager.getRootLogger().debug("LDA weight similarity to previous iteration:" + df.format(LDASimialrity));

            CalcUmpireScore();
        }
    }
    
    public void MixtureModelingSemiSupervised() throws IOException {

        if (libTargetMatches.isEmpty() || Terminate) {
            return;
        }
        LogManager.getRootLogger().info("Semi-parametric mixture modeling");
        int IDNo = 0;
        int decoyNo = 0;
        int modelNo = 0;

        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null) {
                IDNo++;
            }
        }
        decoyNo = decoyModelingList.size();

        for (UmpireSpecLibMatch match : libTargetMatches) {
            //modelNo+= match.TargetHits.size();
            if (match.BestMS1Hit != null) {
                modelNo++;
            }
            if (match.BestMS2Hit != null) {
                modelNo++;
            }
        }

        double[] IDObs = new double[IDNo];
        double[] DecoyObs = new double[decoyNo];
        double[] ModelObs = new double[modelNo];
        int idx = 0;
        int didx = 0;
        int midx = 0;
        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null) {
                IDObs[idx++] = match.BestHit.UmpireScore;
            }
        }
        for (PeakGroupScore peakGroupScore : decoyModelingList) {
            DecoyObs[didx++] = peakGroupScore.UmpireScore;
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
            if (match.BestMS1Hit != null) {
                ModelObs[midx++] = match.BestMS1Hit.UmpireScore;
            }
            if (match.BestMS2Hit != null) {
                ModelObs[midx++] = match.BestMS2Hit.UmpireScore;
            }
        }

        MixtureModelKDESemiParametric mixkde = new MixtureModelKDESemiParametric();
        mixkde.NoBinPoints = NoBinPoints;
        mixkde.SetData(ModelObs, DecoyObs, IDObs);
        mixkde.Modeling();
        mixkde.GeneratePlot(FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatchModel.png");
        MixtureModelProb = mixkde.MixtureModelProb;
    }

    public void MixtureModeling() throws IOException {
        if (libTargetMatches.isEmpty()) {
            return;
        }
        int IDNo = 0;
        int decoyNo = 0;
        int modelNo = 0;
        double IDmean = 0d;
        double Decoymean = 0d;

        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null) {
                IDNo++;
                IDmean += match.BestHit.UmpireScore;
            }
        }

        decoyNo = decoyModelingList.size();
        for (PeakGroupScore peakGroupScore : decoyModelingList) {
            Decoymean += peakGroupScore.UmpireScore;
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
            //modelNo+= match.TargetHits.size();
            if (match.BestMS1Hit != null) {
                modelNo++;
            }
            if (match.BestMS2Hit != null) {
                modelNo++;
            }
        }

        Decoymean /= decoyNo;
        IDmean /= IDNo;

        PVector[] points = new PVector[modelNo];
        PVector[] centroids = new PVector[2];

        int idx = 0;
        for (UmpireSpecLibMatch match : libTargetMatches) {
            if (match.BestMS1Hit != null) {
                points[idx] = new PVector(1);
                points[idx].array[0] = match.BestMS1Hit.UmpireScore;
                idx++;
            }
            if (match.BestMS2Hit != null) {
                points[idx] = new PVector(1);
                points[idx].array[0] = match.BestMS2Hit.UmpireScore;
                idx++;
            }
//            for(PeakGroupScore peakGroupScore : match.TargetHits){
//                points[idx] = new PVector(1);
//                points[idx].array[0] = match.BestMS2Hit.UmpireScore;
//                idx++;
//            }
        }

        MixtureModel mmc;
        centroids[0] = new PVector(1);
        centroids[0].array[0] = Decoymean;
        centroids[1] = new PVector(1);
        centroids[1].array[0] = IDmean;
        Vector<PVector>[] clusters = KMeans.run(points, 2, centroids);
        MixtureModel mm = ExpectationMaximization1D.initialize(clusters);
        mmc = ExpectationMaximization1D.run(points, mm);
        DecimalFormat df = new DecimalFormat("#.####");
        LogManager.getRootLogger().debug("----------------------------------------------------------------------------------------");
        LogManager.getRootLogger().debug("No. of modeling points=" + modelNo);
        LogManager.getRootLogger().debug("ID hits mean=" + df.format(IDmean));
        LogManager.getRootLogger().debug("Decoy hits mean=" + df.format(Decoymean));
        //System.out.print("T-test: p-value=" + df.format(model.ttest.pValue).toString() + "\n");
        LogManager.getRootLogger().debug("Incorrect hits model mean=" + df.format(((PVector) mmc.param[0]).array[0]) + " variance=" + df.format(((PVector) mmc.param[0]).array[1]) + " weight=" + df.format(mmc.weight[0]));
        LogManager.getRootLogger().debug("Correct hits model mean=" + df.format(((PVector) mmc.param[1]).array[0]) + " variance=" + df.format(((PVector) mmc.param[1]).array[1]) + " weight=" + df.format(mmc.weight[1]));

        if (((PVector) mmc.param[0]).array[0] > ((PVector) mmc.param[1]).array[0]) {
            return;
        }

        float max = (float) (((PVector) mmc.param[1]).array[0] + 4 * Math.sqrt(((PVector) mmc.param[1]).array[1]));
        float min = (float) (((PVector) mmc.param[0]).array[0] - 4 * Math.sqrt(((PVector) mmc.param[0]).array[1]));

        IDNo = 0;
        decoyNo = 0;
        modelNo = 0;

        for (PeakGroupScore peakGroupScore : decoyModelingList) {
            if (peakGroupScore.UmpireScore > min && peakGroupScore.UmpireScore < max) {
                decoyNo++;
            }
        }

        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null && match.BestHit.UmpireScore > min && match.BestHit.UmpireScore < max) {
                IDNo++;
            }
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
            //targetNo += match.TargetHits.size();
            //decoyNo += match.DecoyHits.size();
            if (match.BestMS1Hit != null && match.BestMS1Hit.UmpireScore > min && match.BestMS1Hit.UmpireScore < max) {
                modelNo++;
            }
            if (match.BestMS2Hit != null && match.BestMS2Hit.UmpireScore > min && match.BestMS2Hit.UmpireScore < max) {
                modelNo++;
            }
            //modelNo += match.TargetHits.size();            
        }

        double[] IDObs = new double[IDNo];
        double[] DecoyObs = new double[decoyNo];
        double[] ModelObs = new double[modelNo];
        idx = 0;
        int didx = 0;
        int midx = 0;
        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null && match.BestHit.UmpireScore > min && match.BestHit.UmpireScore < max) {
                IDObs[idx++] = match.BestHit.UmpireScore;
            }
        }
        for (PeakGroupScore peakGroupScore : decoyModelingList) {
            if (peakGroupScore.UmpireScore > min && peakGroupScore.UmpireScore < max) {
                DecoyObs[didx++] = peakGroupScore.UmpireScore;
            }
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
//            for(PeakGroupScore peak : match.TargetHits){
//                ModelObs[midx++]=peak.UmpireScore;
//            }
            if (match.BestMS1Hit != null && match.BestMS1Hit.UmpireScore > min && match.BestMS1Hit.UmpireScore < max) {
                ModelObs[midx++] = match.BestMS1Hit.UmpireScore;
            }
            if (match.BestMS2Hit != null && match.BestMS2Hit.UmpireScore > min && match.BestMS2Hit.UmpireScore < max) {
                ModelObs[midx++] = match.BestMS2Hit.UmpireScore;
            }
        }

        String pngfile = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatchModel.png";
        XYSeries model1 = new XYSeries("Incorrect matches");
        XYSeries model2 = new XYSeries("Correct matches");
        XYSeries model3 = new XYSeries("All target hits");

        String modelfile = FilenameUtils.getFullPath(pngfile) + "/" + FilenameUtils.getBaseName(pngfile)+ "_ModelPoints.txt";
        FileWriter writer=new FileWriter(modelfile);
        writer.write("UScore\tModel\tCorrect\tDecoy\n");
        
        int NoPoints = 1000;
        double[] model_kde_x = new double[NoPoints];
        float intv = (max - min) / NoPoints;
        PVector point = new PVector(2);
        for (int i = 0; i < NoPoints; i++) {
            point.array[0] = max - i * intv;
            model_kde_x[i] = point.array[0];
            point.array[1] = mmc.EF.density(point, mmc.param[0]) * mmc.weight[0];
            model1.add(point.array[0], point.array[1]);
            point.array[1] = mmc.EF.density(point, mmc.param[1]) * mmc.weight[1];
            model2.add(point.array[0], point.array[1]);   
            
        }

        KernelDensityEstimator kde = new KernelDensityEstimator();
        kde.SetData(ModelObs);
        double[] model_kde_y = kde.Density(model_kde_x);

        for (int i = 0; i < NoPoints; i++) {
            if (model_kde_x[i] > min && model_kde_x[i] < max) {
                point.array[0] = max - i * intv;
                model_kde_x[i] = point.array[0];
                model3.add(model_kde_x[i], model_kde_y[i]);
                writer.write(point.array[0] + "\t" + mmc.EF.density(point, mmc.param[0]) * mmc.weight[0]
                        + "\t" + mmc.EF.density(point, mmc.param[1]) * mmc.weight[1]
                        + "\t" + model_kde_y[i] + "\n");
            }
        }
        writer.close();
        MixtureModelProb = new float[NoPoints + 1][3];
        float positiveaccu = 0f;
        float negativeaccu = 0f;

        MixtureModelProb[0][0] = (float) model2.getMaxX() + Float.MIN_VALUE;
        MixtureModelProb[0][1] = 1f;
        MixtureModelProb[0][2] = 1f;

        for (int i = 1; i < NoPoints + 1; i++) {
            float positiveNumber = model2.getY(NoPoints - i).floatValue();
            float negativeNumber = model1.getY(NoPoints - i).floatValue();
            MixtureModelProb[i][0] = model2.getX(NoPoints - i).floatValue();
            positiveaccu += positiveNumber;
            negativeaccu += negativeNumber;
            MixtureModelProb[i][2] = positiveNumber / (negativeNumber + positiveNumber);
            MixtureModelProb[i][1] = positiveaccu / (negativeaccu + positiveaccu);
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
//        render.setSeriesPaint(0, Color.DARK_GRAY);
//        render.setSeriesPaint(1, Color.DARK_GRAY); 
//        render.setSeriesPaint(2, Color.GREEN); 
//        render.setSeriesShape(0, new Ellipse2D.Double(0, 0, 2, 2));
//        render.setSeriesShape(1, new Ellipse2D.Double(0, 0, 2, 2));
//        render.setSeriesShape(2, new Ellipse2D.Double(0, 0, 2.5f, 2.5f));
//        render.setSeriesStroke(1, new BasicStroke(1.0f));
//        render.setSeriesStroke(0, new BasicStroke(1.0f));
//        render.setSeriesStroke(2, new BasicStroke(2.0f));
        plot.setDataset(1, dataset);
        plot.setRenderer(1, render);
        plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);
        try {
            ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 1000, 600);
        } catch (IOException e) {
        }
    }

    public void AssignMixtureModelProb() {
        if (Terminate) {
            return;
        }
        for (UmpireSpecLibMatch match : libTargetMatches) {
            AssignProb(match);
        }
        for (UmpireSpecLibMatch match : libIDMatches) {
            AssignProb(match);
        }
    }

    private void AssignProb(UmpireSpecLibMatch match) {
        for (PeakGroupScore peakscore : match.TargetHits) {
            for (int i = 0; i < NoBinPoints - 1; i++) {
                if (peakscore.UmpireScore >= MixtureModelProb[i][0]) {
                    peakscore.MixtureModelProb = MixtureModelProb[i][1];
                    peakscore.MixtureModelLocalProb = MixtureModelProb[i][2];
                    break;
                }
            }
        }
        for (PeakGroupScore peakscore : match.DecoyHits) {
            for (int i = 0; i < NoBinPoints - 1; i++) {
                if (peakscore.UmpireScore >= MixtureModelProb[i][0]) {
                    peakscore.MixtureModelProb = MixtureModelProb[i][1];
                    peakscore.MixtureModelLocalProb = MixtureModelProb[i][2];
                    break;
                }
            }
        }
        match.AssignProbToPepIon();
    }

    public static TargetMatchScoring LibraryMatchRead(String Filename, String LibID) throws FileNotFoundException {

        if (!new File(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.serFS").exists()) {
            return null;
        }
        TargetMatchScoring match = null;
        try {
            LogManager.getRootLogger().info("Loading Target library match results to file:" + FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.serFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.serFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            match = (TargetMatchScoring) in.readObject();
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
        return match;
    }

    public static TargetMatchScoring LibraryMatchReadJS(String Filename, String LibID) throws FileNotFoundException {

        if (!new File(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.serFS").exists()) {
            return null;
        }
        TargetMatchScoring match = null;
        try {
            LogManager.getRootLogger().info("Loading Target library match results to file:" + FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.ser...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            match = (TargetMatchScoring) in.readObject();
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
        return match;
    }

    public boolean LibraryMatchWriteJS() throws FileNotFoundException {
        try {
            LogManager.getRootLogger().info("Writing Target match results to file:" + FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.ser...");
            FileOutputStream fout2 = new FileOutputStream(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.ser", false);
            ObjectOutputStream oos = new ObjectOutputStream(fout2);
            oos.writeObject(this);
            oos.close();
            fout2.close();
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    public boolean LibraryMatchWrite() throws FileNotFoundException {
        try {
            LogManager.getRootLogger().info("Writing Target match results to file:" + FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.serFS...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.serFS", false);
            FSTObjectOutput out = new FSTObjectOutput(fout);
            out.writeObject(this);
            out.close();
            fout.close();
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

}
