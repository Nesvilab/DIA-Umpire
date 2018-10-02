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
package MSUmpire.PeakDataStructure;

import MSUmpire.BaseDataStructure.SortedListFloat;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.MathPackage.ChiSquareGOF;
import MSUmpire.SpectralProcessingModule.Binning;
import MSUmpire.SpectralProcessingModule.ScoreFunction;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

/**
 * Peak isotope cluster data structure 
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakCluster implements Serializable {
    private static final long serialVersionUID = 3545854121L;

    public int Index;
    private transient SortedListFloat MatchScores;
    public transient PeakCurve[] IsoPeaksCurves;
    public PeakCurve MonoIsotopePeak;
    public int[] IsoPeakIndex;
    public float[] Corrs;
    private float[] SNR;
    public float[] PeakHeight;
    public float[] PeakHeightRT;
    public float[] PeakArea;
    public float[] mz;
    public float startRT = Float.MAX_VALUE;
    public float endRT = Float.MIN_VALUE;
    public transient int StartScan;
    public transient int EndScan;
    public int Charge;
    public float IsoMapProb = -1f;
    private float conflictCorr = -1f;
    public float[] PeakDis;    
    public int NoRidges;
    public float OverlapP;
    public transient float[] OverlapRT;
    public float LeftInt;
    public float RightInt;
    public boolean Identified;
    public String AssignedPepIon = "";
    public ArrayList<PrecursorFragmentPairEdge> GroupedFragmentPeaks = new ArrayList<>();
    public float MS1Score;
    public float MS1ScoreLocalProb;
    public float MS1ScoreProbability;
    private transient XYPointCollection FragmentScan;
    public transient String SpectrumKey;

    public PeakCluster(int IsotopicNum, int Charge) {
        IsoPeaksCurves = new PeakCurve[IsotopicNum];
        Corrs = new float[IsotopicNum - 1];
        SNR = new float[IsotopicNum];
        OverlapRT = new float[IsotopicNum-1];
        PeakHeight = new float[IsotopicNum];
        PeakHeightRT = new float[IsotopicNum];
        PeakArea = new float[IsotopicNum];
        IsoPeakIndex = new int[IsotopicNum];
        //PeakDis=new float[IsotopicNum];
        mz = new float[IsotopicNum];
        for (int i = 0; i < IsotopicNum; i++) {
            SNR[i] = -1f;
        }
        this.Charge = Charge;
    }

    public transient ReadWriteLock fraglock = new ReentrantReadWriteLock();          
    transient ReadWriteLock lock = new ReentrantReadWriteLock();
    
    boolean locked=false;
    public XYPointCollection GetNormalizedFragmentScan() throws InterruptedException {
        if (FragmentScan != null && !locked) {
            return FragmentScan;
        }        
        else {
            lock.writeLock().lock();
            try {
                if (FragmentScan == null) {
                    locked=true;
                    FragmentScan = new XYPointCollection();
                    for (PrecursorFragmentPairEdge fragment : GroupedFragmentPeaks) {
                        FragmentScan.AddPoint(fragment.FragmentMz, fragment.Intensity);
                    }
                    FragmentScan.Data.Finalize();
                    Binning bining = new Binning();
                    if (FragmentScan.PointCount() > 2) {
                        FragmentScan=ScoreFunction.SpectralNormalizationForScan(bining.Binning(FragmentScan, 0f, null));
                    }
                    locked=false;
                }
            }
            finally {
                lock.writeLock().unlock();
            }
        }
        return FragmentScan;
    }
    
    public void AddScore(float score) {
        if (MatchScores == null) {
            MatchScores = new SortedListFloat();
        }
        MatchScores.add(score);
    }
    public int GetScoreRank(float score) {
        if (MatchScores != null) {
            return MatchScores.size() - MatchScores.BinarySearchHigher(score) + 1;
        }
        return -1;
    }
    
    public int GetQualityCategory() {
        if ((IsoPeaksCurves==null ||IsoPeaksCurves[2] == null) && mz[2] == 0.0f) {
            return 2;
        }
        return 1;
    }

    public void SetMz(int pkidx, float value) {
        mz[pkidx] = value;
    }
    float RTVar = 0f;

    public float GetConflictCorr() {
        if (conflictCorr == -1f) {
            conflictCorr = IsoPeaksCurves[0].ConflictCorr;
        }
        return conflictCorr;
    }

    public void SetConflictCorr(float ConflictCorr) {
        conflictCorr = ConflictCorr;
    }

    public float TargetMz() {
        if (mz[0] == 0f) {
            mz[0] = IsoPeaksCurves[0].TargetMz;
        }
        return mz[0];
    }

    public void SetSNR(int pkidx, float _snr) {
        SNR[pkidx] = _snr;
    }

    public float GetSNR(int pkidx) {
        if (SNR[pkidx] == -1) {
            if (IsoPeaksCurves!=null && IsoPeaksCurves[pkidx] != null) {
                SNR[pkidx] = IsoPeaksCurves[pkidx].GetRawSNR();
            }
            else if (pkidx==1){
                //Logger.getRootLogger().error("Failed to get SNR");
            }
        }
        return SNR[pkidx];
    }

    private transient float mass=0f;
    public float NeutralMass() {
        if (mass == 0f) {
            if (MonoIsotopePeak != null) {
                mass = Charge * (float)((MonoIsotopePeak.TargetMz - ElementaryIon.proton.getTheoreticMass()));
            }
            else {
                mass = Charge * (float)((mz[0] - ElementaryIon.proton.getTheoreticMass()));
            }
        }
        return mass;
    }

    public void UpdateIsoMapProb(TreeMap<Float, XYData>[] IsotopePatternMap) {
        if (IsoMapProb == -1) {
            IsoMapProb = GetChiSquareProbByIsoMap(IsotopePatternMap);
        }
    }

    public void AssignConfilictCorr() {
        for (int i = 1; i < IsoPeaksCurves.length; i++) {
            if (IsoPeaksCurves[i] != null) {
                if (Corrs[i - 1] > 0.6f) {
                    IsoPeaksCurves[i].AddConflictScore(Corrs[i - 1]);
                }
            }
        }
    }

    public void CalcPeakArea_V2() {
        int NoOfIsotopic = IsoPeaksCurves.length;

        PeakCurve peakA = MonoIsotopePeak;
        startRT = MonoIsotopePeak.StartRT();
        endRT = MonoIsotopePeak.EndRT();
        
        if (IsoPeaksCurves.length>1 && IsoPeaksCurves[1]!=null) {
            startRT = Math.min(MonoIsotopePeak.StartRT(), IsoPeaksCurves[1].StartRT());
            endRT = Math.max(MonoIsotopePeak.EndRT(), IsoPeaksCurves[1].EndRT());
        }
        
        if(endRT==startRT){
            startRT=MonoIsotopePeak.GetSmoothedList().Data.get(0).getX();
            endRT=MonoIsotopePeak.GetSmoothedList().Data.get(MonoIsotopePeak.GetSmoothedList().PointCount()-1).getX();
        }

        NoRidges = 0;
        if (peakA.RegionRidge != null) {
//            for (Float ridge : peakA.RegionRidge) {
            for (int i=0; i<peakA.RegionRidge.size(); ++i) {
                final float ridge = peakA.RegionRidge.get(i);
                if (ridge >= startRT && ridge <= endRT) {
                    NoRidges++;
                }
            }
        }

        for (int i = 0; i < NoOfIsotopic; i++) {
            PeakCurve peak = IsoPeaksCurves[i];
            if (peak == null) {
                break;
            }
            for (int j = 0; j < peak.GetSmoothedList().PointCount(); j++) {
                XYData pt = peak.GetSmoothedList().Data.get(j);
                if (pt.getX() >= startRT && pt.getX() <= endRT) {
                    PeakArea[i] += pt.getY();
                    if (pt.getY() > PeakHeight[i]) {
                        PeakHeight[i] = pt.getY();
                        PeakHeightRT[i] = pt.getX();
                    }
                }
            }
            mz[i]=peak.TargetMz;
            IsoPeakIndex[i]=peak.Index;
        }        
    }

    //Generate isotope peak distribution
    private void GeneratePeakDis() {
        if (PeakDis != null) {
            return;
        }
        PeakDis = new float[PeakHeight.length];
        float firstPeak = PeakHeight[0];
        for (int i = 0; i < PeakDis.length; i++) {
            if (PeakHeight[i] > 0) {
                PeakDis[i] = PeakHeight[i] / firstPeak;
            }
        }
    }

    //Get isotope pattern range according the mass of this peak cluster
    public XYData[] GetPatternRange(TreeMap<Float, XYData>[] IsotopePatternMap) {
        XYData[] PatternRange = new XYData[IsotopePatternMap.length];
        for (int i = 0; i < IsotopePatternMap.length; i++) {
            Map.Entry range = IsotopePatternMap[i].ceilingEntry(NeutralMass());
            if (range == null) {
                range = IsotopePatternMap[i].lastEntry();
            }
            PatternRange[i] = (XYData) range.getValue();
        }
        return PatternRange;
    }


    private float GetChiSquareProbByIsoMap(TreeMap<Float, XYData>[] IsotopePatternMap) {

        GeneratePeakDis();
        XYData[] PatternRange = new XYData[IsotopePatternMap.length];

        for (int i = 0; i < IsotopePatternMap.length; i++) {
            Map.Entry range = IsotopePatternMap[i].ceilingEntry(NeutralMass());
            if (range == null) {
                range = IsotopePatternMap[i].lastEntry();
            }
            PatternRange[i] = (XYData) range.getValue();
        }
        float[] TheoIso = new float[IsotopePatternMap.length];

        TheoIso[0] = 1f;

        for (int i = 1; i < IsotopePatternMap.length; i++) {
            if (PeakDis[i] >= PatternRange[i - 1].getY() && PeakDis[i] <= PatternRange[i - 1].getX()) {
                TheoIso[i] = PeakDis[i];
            } else {
                if (Math.abs(PeakDis[1] - PatternRange[i - 1].getY()) > Math.abs(PeakDis[i] - PatternRange[i - 1].getX())) {
                    TheoIso[i] = PatternRange[i - 1].getX();
                } else {
                    TheoIso[i] = PatternRange[i - 1].getY();
                }
            }
        }
        float prob = ChiSquareGOF.GetInstance(IsoPeaksCurves.length).GetGoodNessOfFitProb(TheoIso, PeakDis);

        return prob;
    }

    //Check is the number of detected isotope peaks passes the criterion
    public boolean IsotopeComplete(int minIsonum) {
        for (int i = 0; i < minIsonum; i++) {
            if ((IsoPeaksCurves==null || IsoPeaksCurves[i] == null) && mz[i] == 0.0f) {
                return false;
            }
        }
        return true;
    }

    //Create locks for multithreading 
    public void CreateLock() {
        lock=new ReentrantReadWriteLock();
        locked=false;
        fraglock=new ReentrantReadWriteLock();
    }

    public float GetMaxMz() {
        for (int i = mz.length - 1; i > 0; i--) {
            if (mz[i] > 0.0f) {
                return mz[i];
            }
        }
        return mz[0];
    }
}
