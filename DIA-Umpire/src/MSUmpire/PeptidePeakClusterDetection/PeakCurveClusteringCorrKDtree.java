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
package MSUmpire.PeptidePeakClusterDetection;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.MathPackage.MassDefect;
import MSUmpire.PeakDataStructure.IonChargeHashSet;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PeakCurve;
import MSUmpire.PeakDataStructure.SortedCurveCollectionMZ;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import net.sf.javaml.core.kdtree.KDTree;
import net.sf.javaml.core.kdtree.KeySizeException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Isotope peak curve clustering based on peak profile correlation. The searching is using KD tree
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
//public class PeakCurveClusteringCorrKDtree implements Runnable {
public class PeakCurveClusteringCorrKDtree implements Callable<ArrayList<PeakCluster>> {

    PeakCurve peakA;
    InstrumentParameter parameter;
    MassDefect MD=new MassDefect();
    private final KDTree PeakCurveSearchTree;
    private final TreeMap<Float, XYData>[] IsotopePatternMap;
    public ArrayList<PeakCluster> ResultClusters = new ArrayList<>();

    private final int MaxNoOfClusters;
    private final int MinNoOfClusters;
    private final int StartCharge;
    private final int EndCharge;
    
    public PeakCurveClusteringCorrKDtree(PeakCurve targeCurve, KDTree PeakCurveSearchTree, InstrumentParameter parameter, TreeMap<Float, XYData>[] IsotopePatternMap, int StartCharge, int EndCharge, int MaxNoClusters, int MinNoClusters) {
        this.peakA = targeCurve;
        this.PeakCurveSearchTree = PeakCurveSearchTree;
        this.parameter = parameter;
        this.IsotopePatternMap = IsotopePatternMap;
        this.MaxNoOfClusters = MaxNoClusters;
        this.MinNoOfClusters = MinNoClusters;
        this.StartCharge = StartCharge;
        this.EndCharge = EndCharge;
    }

    @Override
//    public void run() {
    public ArrayList<PeakCluster> call() {

        float lowrt = peakA.ApexRT - parameter.ApexDelta;
        float highrt = peakA.ApexRT + parameter.ApexDelta;
        float lowmz = InstrumentParameter.GetMzByPPM(peakA.TargetMz, 1, parameter.MS1PPM);
        float highmz = InstrumentParameter.GetMzByPPM((peakA.TargetMz+ ((float)this.MaxNoOfClusters/this.StartCharge)), 1, -parameter.MS1PPM);
       
        Object[] found=null;
        try {
            found = PeakCurveSearchTree.range(new double[]{lowrt,lowmz}, new double[]{highrt,highmz});
        } catch (KeySizeException ex) {
            StringWriter w = new StringWriter();
            PrintWriter pw = new PrintWriter(w);
            ex.printStackTrace(pw);
            LogManager.getRootLogger().error(w.toString());
        }
        if(found==null || found.length==0){
            return ResultClusters;
        }
        SortedCurveCollectionMZ PeakCurveListMZ=new SortedCurveCollectionMZ();
        for(Object peakCurve : found){
            PeakCurveListMZ.add((PeakCurve)peakCurve);
        }
        
        float Arange = peakA.EndRT() - peakA.StartRT();
        for (int charge = EndCharge; charge >= StartCharge; charge--) {
            float mass=charge * (float)((peakA.TargetMz - ElementaryIon.proton.getTheoreticMass()));
            if(mass<parameter.MinPrecursorMass || mass>parameter.MaxPrecursorMass || (parameter.MassDefectFilter && !MD.InMassDefectRange(mass, parameter.MassDefectOffset))){
                continue;
            }
            PeakCluster peakCluster = new PeakCluster(MaxNoOfClusters, charge);
            peakCluster.IsoPeaksCurves[0] = peakA;
            peakCluster.MonoIsotopePeak=peakA;
            XYData[] Ranges = new XYData[MaxNoOfClusters - 1];
            for (int i = 0; i < MaxNoOfClusters - 1; i++) {
                Entry range = IsotopePatternMap[i].ceilingEntry(peakCluster.NeutralMass());
                if (range == null) {
                    range = IsotopePatternMap[i].lastEntry();
                }
                Ranges[i] = (XYData) range.getValue();
            }

            for (int pkidx = 1; pkidx < MaxNoOfClusters; pkidx++) {
                float ppmthreshold = parameter.MS1PPM + (parameter.MS1PPM * pkidx * 0.5f);
                float lowtheomz = InstrumentParameter.GetMzByPPM(peakA.TargetMz + (pkidx * ((float)ElementaryIon.proton.getTheoreticMass() / charge)), charge, ppmthreshold);
                float uptheomz = InstrumentParameter.GetMzByPPM(peakA.TargetMz + (pkidx * ((float) ElementaryIon.proton.getTheoreticMass() / charge)), charge, -ppmthreshold);
                int startmzidx = PeakCurveListMZ.BinarySearchLower(lowtheomz);

                float theomz = peakA.TargetMz + (pkidx * ((float)ElementaryIon.proton.getTheoreticMass() / charge));
                float maxscore = 0f;
                float maxcorr = 0f;
                float maxoverlap = 0f;
                PeakCurve closestPeak = null;

                for (int mzidx = startmzidx; mzidx < PeakCurveListMZ.size(); mzidx++) {
                    PeakCurve peakB = PeakCurveListMZ.get(mzidx);

                    if (peakB.TargetMz <= peakA.TargetMz) {
                        continue;
                    }
                    if (peakB.TargetMz > uptheomz) {
                        break;
                    }

                    float Brange = peakB.EndRT() - peakB.StartRT();
                    float OverlapP = 0f;
                    if (peakA.StartRT() >= peakB.StartRT() && peakA.StartRT() <= peakB.EndRT() && peakA.EndRT() >= peakB.EndRT()) {
                        OverlapP = (peakB.EndRT() - peakA.StartRT()) / Brange;

                    } else if (peakA.EndRT() >= peakB.StartRT() && peakA.EndRT() <= peakB.EndRT() && peakA.StartRT() <= peakB.StartRT()) {
                        OverlapP = (peakA.EndRT() - peakB.StartRT()) / Brange;

                    } else if (peakA.StartRT() <= peakB.StartRT() && peakA.EndRT() >= peakB.EndRT()) {
                        OverlapP = 1;

                    } else if (peakA.StartRT() >= peakB.StartRT() && peakA.EndRT() <= peakB.EndRT()) {
                        OverlapP = Arange / Brange;
                    }
                    if (parameter.TargetIDOnly || (OverlapP > parameter.MiniOverlapP && (!parameter.CheckMonoIsotopicApex || (peakA.ApexRT >= peakB.StartRT() && peakA.ApexRT <= peakB.EndRT() && peakB.ApexRT >= peakA.StartRT() && peakB.ApexRT <= peakA.EndRT())))) {
                        float ppm = InstrumentParameter.CalcPPM(theomz, peakB.TargetMz);
                        if (ppm < ppmthreshold) {
                            float corr = 0f;
                            try {
                                corr = PeakCurveCorrCalc.CalPeakCorr(peakA, peakB, parameter.NoPeakPerMin);
                            } catch (IOException ex) {
                                StringWriter w = new StringWriter();
                                PrintWriter pw = new PrintWriter(w);
                                ex.printStackTrace(pw);
                                LogManager.getRootLogger().error(w.toString());
                            }
                            if (Float.isNaN(corr)) {
                                corr = 0f;
                                //System.out.print("Corr=NAN\n");
                            }

                            float PeakIntA = peakA.ApexInt;
                            float PeakIntB = peakB.GetMaxIntensityByRegionRange(Math.max(peakA.StartRT(), peakB.StartRT()), Math.min(peakB.EndRT(), peakA.EndRT()));

                            if ((parameter.TargetIDOnly && corr>0.2f) || corr > parameter.IsoCorrThreshold ) {
                            //if (corr > parameter.IsoCorrThreshold || (PeakIntA > PeakIntB*1.5f && PeakIntB>0.1f * PeakIntA && (peakA.EndScan-peakA.StartScan)>4 && (peakB.EndScan-peakB.StartScan)>4)) {
                                float intscore = 0f;
                                float IntRatio = PeakIntB / PeakIntA;

                                if (IntRatio > Ranges[pkidx - 1].getY() && IntRatio <= Ranges[pkidx - 1].getX()) {
                                    intscore = 1f;
                                } else {
                                    if (Math.abs(IntRatio - Ranges[pkidx - 1].getY()) > Math.abs(IntRatio - Ranges[pkidx - 1].getX())) {
                                        intscore = 1 - Math.abs(IntRatio - Ranges[pkidx - 1].getX());
                                    } else {
                                        intscore = 1 - Math.abs(IntRatio - Ranges[pkidx - 1].getY());
                                    }
                                }
                                if (intscore < 0f) {
                                    intscore = 0f;
                                }
                                float score = ((ppmthreshold - ppm) / ppmthreshold) + corr + intscore;

                                if (maxscore < score) {
                                    maxscore = score;
                                    closestPeak = peakB;
                                    maxcorr = corr;
                                    maxoverlap = OverlapP;
                                }
                            }
                        }
                    }
                }

                if (closestPeak != null) {
                    peakCluster.Corrs[pkidx - 1] = maxcorr;
                    peakCluster.IsoPeaksCurves[pkidx] = closestPeak;
                    peakCluster.OverlapRT[pkidx - 1] = maxoverlap;
                    peakCluster.GetSNR(pkidx - 1);

                    if (pkidx == 1) {
                        peakCluster.OverlapP = maxoverlap;
                    }
                }
                if (closestPeak == null) {
                    break;
                }
            }
            if (peakCluster.IsotopeComplete(MinNoOfClusters)) {
                peakCluster.CalcPeakArea_V2();
                peakCluster.UpdateIsoMapProb(IsotopePatternMap);
                peakCluster.AssignConfilictCorr();
                peakCluster.LeftInt = peakA.GetSmoothedList().Data.get(0).getY();
                peakCluster.RightInt = peakA.GetSmoothedList().Data.get(peakA.GetSmoothedList().PointCount() - 1).getY();
                if (parameter.TargetIDOnly ||peakCluster.IsoMapProb > parameter.IsoPattern) {
                    ResultClusters.add(peakCluster);
                    if (!parameter.TargetIDOnly ||(parameter.RemoveGroupedPeaks && peakCluster.Corrs[0] > parameter.RemoveGroupedPeaksCorr && peakCluster.OverlapP > parameter.RemoveGroupedPeaksRTOverlap)) {
                        for (int i = 1; i < peakCluster.IsoPeaksCurves.length; i++) {
                            PeakCurve peak = peakCluster.IsoPeaksCurves[i];
                            if (peak != null && peakCluster.Corrs[i - 1] > parameter.RemoveGroupedPeaksCorr && peakCluster.OverlapRT[i - 1] > parameter.RemoveGroupedPeaksRTOverlap) {
//                                peak.ChargeGrouped.add(charge);
                                peak.ChargeGrouped = IonChargeHashSet.add(peak.ChargeGrouped,charge);
                            }
                        }
                    }
                }
            }
        }

    //System.out.print("....done\n");
        return ResultClusters;
    }
}
