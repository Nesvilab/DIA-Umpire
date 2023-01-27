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
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.LCMSPeakStructure.LCMSPeakDIAMS2;
import MSUmpire.LCMSPeakStructure.LCMSPeakMS1;
import MSUmpire.MathPackage.Regression;
import MSUmpire.PSMDataStructure.FragmentPeakGroup;
import MSUmpire.PSMDataStructure.PepFragmentLib;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.SpectralProcessingModule.ScoreFunction;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Spectral matching to calculate subscores given a peptide ion from spectral library and a precursor-fragment group
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class UmpireSpecLibMatch implements Runnable, Serializable{
    private static final long serialVersionUID = 3164978164L;

    transient PepFragmentLib fragmentLib;
    transient PepFragmentLib decoyfragmentLib;
    transient InstrumentParameter parameter;
    transient LCMSPeakMS1 ms1lcms;
    transient LCMSPeakDIAMS2 DIAWindow;
    PepIonID pepIonID;
    public PeakGroupScore BestMS1Hit;
    public PeakGroupScore BestMS2Hit;
    public PeakGroupScore BestHit;
    public PeakGroupScore BestMS1DecoyHit;
    public PeakGroupScore BestMS2DecoyHit;
    public PeakGroupScore BestDecoyHit;
    public ArrayList<PeakGroupScore> TargetHits=new ArrayList<>();
    public ArrayList<PeakGroupScore> DecoyHits=new ArrayList<>();
    public boolean IdentifiedPeptideIon = false;

    public UmpireSpecLibMatch(LCMSPeakMS1 ms1lcms, LCMSPeakDIAMS2 DIAWindow, PepIonID pepIonID, PepFragmentLib fragmentLib, PepFragmentLib decoyfragmentLib, InstrumentParameter parameter) {
        this.ms1lcms = ms1lcms;
        this.DIAWindow = DIAWindow;
        this.fragmentLib = fragmentLib;
        this.decoyfragmentLib = decoyfragmentLib;
        this.parameter = parameter;
        this.pepIonID = pepIonID;
    }

    public class MatchFragment {

        public FragmentPeakGroup libfrag;
        public PrecursorFragmentPairEdge peakfrag;

        public MatchFragment(FragmentPeakGroup libfrag, PrecursorFragmentPairEdge peakfrag) {
            this.libfrag = libfrag;
            this.peakfrag = peakfrag;
        }
    }
    
    public String GetPepKey(){
        return pepIonID.GetKey();
    }

    private void CalMatchScoreOld(PeakCluster cluster, PepFragmentLib fragmentLib, PeakGroupScore peakscore) {
        Regression regression = new Regression();
        XYPointCollection pointset = new XYPointCollection();
        ArrayList<MatchFragment> matchFragments = new ArrayList<>();
        for (FragmentPeakGroup frag : fragmentLib.FragmentGroups.values()) {
            PrecursorFragmentPairEdge bestfragment = null;
            for (PrecursorFragmentPairEdge fragmentClusterUnit : cluster.GroupedFragmentPeaks) {
                if (InstrumentParameter.CalcPPM(frag.FragMZ, fragmentClusterUnit.FragmentMz) <= parameter.MS2PPM) {
                    if (bestfragment == null || fragmentClusterUnit.Correlation > bestfragment.Correlation) {
                        bestfragment = fragmentClusterUnit;
                    }
                }
            }
            matchFragments.add(new MatchFragment(frag, bestfragment));
        }

        for (MatchFragment matchFragment : matchFragments) {
            if (matchFragment.peakfrag != null) {
                pointset.AddPoint(matchFragment.libfrag.GetAvgInt(), matchFragment.peakfrag.Intensity);
                peakscore.FragIntAvgScore += matchFragment.peakfrag.Intensity;
                peakscore.SumCorrScore += matchFragment.peakfrag.Correlation;
                peakscore.ApexDeltaScore += 1 - matchFragment.peakfrag.ApexDelta;
                peakscore.PPMScore += 1 - InstrumentParameter.CalcPPM(matchFragment.peakfrag.FragmentMz, matchFragment.libfrag.FragMZ) / parameter.MS2PPM;
                peakscore.RTOverlapScore += matchFragment.peakfrag.RTOverlapP;
                if (matchFragment.libfrag.IonType.startsWith("b")) {
                    peakscore.NoMatchB++;
                }
                if (matchFragment.libfrag.IonType.startsWith("y")) {
                    peakscore.NoMatchY++;
                }
            } else {
                pointset.AddPoint(matchFragment.libfrag.GetAvgInt(), 0f);
            }
        }
        peakscore.NoFragmentLib = matchFragments.size();

        int NoMatchFrag=peakscore.NoMatchB+peakscore.NoMatchY;
        peakscore.FragIntAvgScore /= NoMatchFrag;
        peakscore.SumCorrScore /= NoMatchFrag;
        peakscore.ApexDeltaScore /= NoMatchFrag;
        peakscore.PPMScore /= NoMatchFrag;
        peakscore.RTOverlapScore /= NoMatchFrag;
                
        float R2 = 0f;
        regression.SetData(pointset);
        if (regression.equation.Mvalue > 0) {
            R2 = regression.GetR2();
        }
        peakscore.SpecCorrelation = R2;
    }

    private void CalMatchScore(PeakCluster cluster, PepFragmentLib fragmentLib, PeakGroupScore peakscore) {
        if (parameter.UseOldVersion) {
            CalMatchScoreOld(cluster, fragmentLib, peakscore);
        } else {
            CalMatchScoreNew(cluster, fragmentLib, peakscore);
        }
    }

    private void CalMatchScoreNew(PeakCluster cluster, PepFragmentLib fragmentLib, PeakGroupScore peakscore) {                
        XYPointCollection pointset = new XYPointCollection();
        ArrayList<MatchFragment> matchFragments = new ArrayList<>();        
        for (FragmentPeakGroup frag : fragmentLib.FragmentGroups.values()) {
            PrecursorFragmentPairEdge bestfragment = null;
            for (PrecursorFragmentPairEdge fragmentClusterUnit : cluster.GroupedFragmentPeaks) {                
                if (InstrumentParameter.CalcPPM(frag.FragMZ, fragmentClusterUnit.FragmentMz) <= parameter.MS2PPM) {
                    if (bestfragment == null || fragmentClusterUnit.Correlation > bestfragment.Correlation) {
                        bestfragment = fragmentClusterUnit;
                    }
                }
            }
            matchFragments.add(new MatchFragment(frag, bestfragment));
        }

        for (MatchFragment matchFragment : matchFragments) {
            if (matchFragment.peakfrag != null) {
                pointset.AddPoint(matchFragment.libfrag.GetAvgInt(), matchFragment.peakfrag.Intensity);
                float ppmscore=1 - InstrumentParameter.CalcPPM(matchFragment.peakfrag.FragmentMz, matchFragment.libfrag.FragMZ) / parameter.MS2PPM;
                peakscore.FragIntAvgScore += matchFragment.peakfrag.Intensity;         
                peakscore.SumCorrScore += matchFragment.peakfrag.Correlation; 
                peakscore.SumCorrPPMScore += matchFragment.peakfrag.Correlation*ppmscore;
                if(peakscore.MaxMatchCorr<matchFragment.peakfrag.Correlation){
                    peakscore.MaxMatchCorr=matchFragment.peakfrag.Correlation;
                }
                peakscore.ApexDeltaScore += 1 - matchFragment.peakfrag.ApexDelta;
                peakscore.PPMScore += ppmscore;
                peakscore.RTOverlapScore += matchFragment.peakfrag.RTOverlapP;
                if (matchFragment.libfrag.IonType.startsWith("b")) {
                    peakscore.NoMatchB++;
                }
                if (matchFragment.libfrag.IonType.startsWith("y")) {
                    peakscore.NoMatchY++;
                }
            } else {
                pointset.AddPoint(matchFragment.libfrag.GetAvgInt(), 0.000000000001f);
            }
        }
        int totmamatch=peakscore.NoMatchB+peakscore.NoMatchY;
        peakscore.NoFragmentLib = matchFragments.size();
        peakscore.FragIntAvgScore = (float) Math.log(peakscore.FragIntAvgScore/totmamatch);
        
        peakscore.AveCorrScore =peakscore.SumCorrScore/ totmamatch;
        peakscore.ApexDeltaScore /= totmamatch;
        peakscore.PPMScore /= totmamatch;
        peakscore.RTOverlapScore /= totmamatch;
        peakscore.PrecursorScore = cluster.MS1Score;
        peakscore.PrecursorCorr = cluster.Corrs[0];
        peakscore.PrecursorIsoPattern=cluster.IsoMapProb;
        peakscore.Peplength=pepIonID.Sequence.length();
                
        XYPointCollection norP=ScoreFunction.SpectralNormalizationForPairCollection(pointset);
        peakscore.SpecCorrelation=ScoreFunction.CalcSpecCorrForPairPointCollection(norP);
        peakscore.SpecDotProduct=ScoreFunction.CalcDotProductForPairPointCollection(norP);
        peakscore.ContrastAngle=ScoreFunction.CalcSpecContrastAngleForPairPointCollection(norP);        
    }

    @Override
    public void run() {
        if (fragmentLib != null) {
            if (IdentifiedPeptideIon) {
                IDPeak();
            } else {
                Target();
                if (decoyfragmentLib != null) {
                    Decoy();
                }
                else{
                    LogManager.getRootLogger().error("decoy spectrum is null : "+pepIonID.GetKey());
                }
            }
        }
        else{
            LogManager.getRootLogger().warn("lib spectrum is null : "+pepIonID.GetKey());
        }
    }

    public void IDPeak() {
        ArrayList<PeakCluster> clusterList = pepIonID.MS1PeakClusters;
        TargetHits = new ArrayList<>();
        for (PeakCluster Assigncluster : clusterList) {
            PeakCluster cluster=ms1lcms.PeakClusters.get(Assigncluster.Index-1);
            DIAWindow.ExtractFragmentForPeakCluser(cluster);            
            PeakGroupScore peakGroupScore = new PeakGroupScore(cluster);
            peakGroupScore.MSlevel = 1;       
             if(!pepIonID.PredictRT.isEmpty()){
                peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetAvgPredictRT());                        
            }
            else {
                peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetIDRT());
            }
            peakGroupScore.PrecursorPPM = InstrumentParameter.CalcPPM(cluster.TargetMz(), pepIonID.NeutralPrecursorMz());            
            CalMatchScore(cluster, fragmentLib, peakGroupScore);
            if (peakGroupScore.NoMatchB+peakGroupScore.NoMatchY > 0) {
                TargetHits.add(peakGroupScore);
            }
        }
        clusterList = pepIonID.MS2UnfragPeakClusters;
        for (PeakCluster Assigncluster : clusterList) {
            if (DIAWindow.PeakClusters.size() >= Assigncluster.Index) {
                PeakCluster cluster = DIAWindow.PeakClusters.get(Assigncluster.Index - 1);
                if (cluster.TargetMz() == Assigncluster.TargetMz() || cluster.Charge == Assigncluster.Charge) {
                    DIAWindow.ExtractFragmentForUnfragPeakCluser(cluster);
                    PeakGroupScore peakGroupScore = new PeakGroupScore(cluster);                    
                    peakGroupScore.MSlevel = 2;
                    if (!pepIonID.PredictRT.isEmpty()) {
                        peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetAvgPredictRT());
                    } else {
                        peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetIDRT());
                    }
                    peakGroupScore.PrecursorPPM = InstrumentParameter.CalcPPM(cluster.TargetMz(), pepIonID.NeutralPrecursorMz());
                    CalMatchScore(cluster, fragmentLib, peakGroupScore);
                    if (peakGroupScore.NoMatchB + peakGroupScore.NoMatchY > 0) {
                        TargetHits.add(peakGroupScore);
                    }
                }
            }
        }
    }

    public void Target() {        
        ArrayList<PeakCluster> clusterList = ms1lcms.FindAllPeakClustersForMappedPep(pepIonID);
        TargetHits = new ArrayList<>();
        for (PeakCluster cluster : clusterList) {
            DIAWindow.ExtractFragmentForPeakCluser(cluster);
            PeakGroupScore peakGroupScore = new PeakGroupScore(cluster);
            peakGroupScore.MSlevel = 1;
            peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetAvgPredictRT());            
            peakGroupScore.PrecursorPPM = InstrumentParameter.CalcPPM(cluster.TargetMz(), pepIonID.NeutralPrecursorMz());            
            CalMatchScore(cluster, fragmentLib, peakGroupScore);
            if (peakGroupScore.NoMatchB+peakGroupScore.NoMatchY > 0) {
                TargetHits.add(peakGroupScore);
            }
        }
        clusterList = DIAWindow.FindAllPeakClustersForMappedPep(pepIonID);
        for (PeakCluster cluster : clusterList) {
            DIAWindow.ExtractFragmentForUnfragPeakCluser(cluster);
            PeakGroupScore peakGroupScore = new PeakGroupScore(cluster);
            peakGroupScore.MSlevel = 2;
            peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetAvgPredictRT());            
            peakGroupScore.PrecursorPPM = InstrumentParameter.CalcPPM(cluster.TargetMz(), pepIonID.NeutralPrecursorMz());
            CalMatchScore(cluster, fragmentLib, peakGroupScore);
            if (peakGroupScore.NoMatchB+peakGroupScore.NoMatchY > 0) {
                TargetHits.add(peakGroupScore);
            }
        }
    }

    public void Decoy() {
        ArrayList<PeakCluster> clusterList = ms1lcms.FindAllPeakClustersForMappedPep(pepIonID);
        DecoyHits = new ArrayList<>();
        for (PeakCluster cluster : clusterList) {
            PeakGroupScore peakGroupScore = new PeakGroupScore(cluster);
            peakGroupScore.MSlevel = 1;      
            if(!pepIonID.PredictRT.isEmpty()){
                peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetAvgPredictRT());                        
            }
            else {
                peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetIDRT());
            }
            peakGroupScore.PrecursorPPM = InstrumentParameter.CalcPPM(cluster.TargetMz(), pepIonID.NeutralPrecursorMz());
            CalMatchScore(cluster, decoyfragmentLib, peakGroupScore);
            if (peakGroupScore.NoMatchB+peakGroupScore.NoMatchY > 0) {
                DecoyHits.add(peakGroupScore);
            }
        }
        clusterList = DIAWindow.FindAllPeakClustersForMappedPep(pepIonID);
        for (PeakCluster cluster : clusterList) {
            PeakGroupScore peakGroupScore = new PeakGroupScore(cluster);
            peakGroupScore.MSlevel = 2;
             if(!pepIonID.PredictRT.isEmpty()){
                peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetAvgPredictRT());                        
            }
            else {
                peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetIDRT());
            }           
            peakGroupScore.PrecursorPPM = InstrumentParameter.CalcPPM(cluster.TargetMz(), pepIonID.NeutralPrecursorMz());            
            CalMatchScore(cluster, decoyfragmentLib, peakGroupScore);
            if (peakGroupScore.NoMatchB+peakGroupScore.NoMatchY> 0) {
                DecoyHits.add(peakGroupScore);
            }
        }
    }

    public void AssignProbToPepIon() {
        pepIonID.MS1PeakClusters.clear();
        pepIonID.MS2UnfragPeakClusters.clear();
        pepIonID.UScoreProbability_MS1 = 0f;
        pepIonID.MS1AlignmentProbability = 0f;
        pepIonID.UScoreProbability_MS2 = 0f;
        pepIonID.MS2AlignmentProbability = 0f;
        if (BestMS1Hit != null) {
            pepIonID.MS1PeakClusters.add(BestMS1Hit.cluster);
            pepIonID.MS1AlignmentProbability = BestMS1Hit.MixtureModelProb;
            pepIonID.UScoreProbability_MS1 = BestMS1Hit.MixtureModelLocalProb;
        }
        if (BestMS2Hit != null) {
            pepIonID.MS2UnfragPeakClusters.add(BestMS2Hit.cluster);
            pepIonID.MS2AlignmentProbability = BestMS2Hit.MixtureModelProb;
            pepIonID.UScoreProbability_MS2 = BestMS2Hit.MixtureModelLocalProb;
        }                
    }

    
    public void Ranking() {
        Collections.sort(TargetHits, new Comparator<PeakGroupScore>() {
            @Override
            public int compare(PeakGroupScore o1, PeakGroupScore o2) {
                return -Float.compare(o1.UmpireScore, o2.UmpireScore);
            }
        });

        Collections.sort(DecoyHits, new Comparator<PeakGroupScore>() {
            @Override
            public int compare(PeakGroupScore o1, PeakGroupScore o2) {
                return -Float.compare(o1.UmpireScore, o2.UmpireScore);
            }
        });

        if (!TargetHits.isEmpty()) {
            BestHit = TargetHits.get(0);
        }
        for (int i = 0; i < TargetHits.size(); i++) {
            if (BestMS1Hit != null && BestMS2Hit != null) {
                break;
            }
            if (TargetHits.get(i).MSlevel == 1 && BestMS1Hit == null) {
                BestMS1Hit = TargetHits.get(i);
            }
            if (TargetHits.get(i).MSlevel == 2 && BestMS2Hit == null) {
                BestMS2Hit = TargetHits.get(i);
            }
        }

         if (!DecoyHits.isEmpty()) {
            BestDecoyHit = DecoyHits.get(0);
        }
        for (int i = 0; i < DecoyHits.size(); i++) {
            if (BestMS1DecoyHit != null && BestMS2DecoyHit != null) {
                break;
            }
            if (DecoyHits.get(i).MSlevel == 1 && BestMS1DecoyHit == null) {
                BestMS1DecoyHit = DecoyHits.get(i);
            }
            if (DecoyHits.get(i).MSlevel == 2 && BestMS2DecoyHit == null) {
                BestMS2DecoyHit = DecoyHits.get(i);
            }
        }
    }
}
