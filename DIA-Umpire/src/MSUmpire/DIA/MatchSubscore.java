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

import java.io.Serializable;

/**
 * Data structure to define which match subscores to use
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MatchSubscore implements Serializable{
    private static final long serialVersionUID = -6347236539892953039L;    
    public subscore[] Subscores;    
    public int NoEnableSubscores;
    public double[] SubSCoeff;
    public String[] SubSName;
    int Version=2;
    
    public class subscore implements Serializable{
        private static final long serialVersionUID = 2801998695766081693L;
        public String name;
        public boolean enable;
        public double initialweight;

        public subscore(String name, boolean enable,double iniweight) {
            this.name = name;
            this.enable = enable;
            this.initialweight=iniweight;
        }
    }
    private void V1_Array(){
         Subscores = new subscore[7];
        Subscores[0] = new subscore("SpecCorrelation", true, 25d);
        Subscores[1] = new subscore("CorrScore", true, -2d);
        Subscores[2] = new subscore("IntScore", true, -0.2d);
        Subscores[3] = new subscore("PPMScore", true, 4d);
        Subscores[4] = new subscore("ApexDeltaScore", true, -2d);
        Subscores[5] = new subscore("RTOverlapScore", true, 1d);
        Subscores[6] = new subscore("NoMatchFragment", true, 1d);  
    }
     private void V2_Array(){
         Subscores = new subscore[18];
        Subscores[0] = new subscore("SpecDotProduct", true, 25d);
        Subscores[1] = new subscore("SpecCorrelation", false, 25d);
        Subscores[2] = new subscore("SpecContrastAngle", false, 25d);
        Subscores[3] = new subscore("CorrScore", true, -2d);
        Subscores[4] = new subscore("IntScore", false, -0.2d);
        Subscores[5] = new subscore("PPMScore", true, 4d);
        Subscores[6] = new subscore("ApexDeltaScore", true, -2d);
        Subscores[7] = new subscore("RTOverlapScore", true, 1d);
        Subscores[8] = new subscore("NoMatchB", false, 1d);
        Subscores[9] = new subscore("NoMatchY", false, 1d);
        Subscores[10] = new subscore("PrecursorCorr", true, 1d);
        Subscores[11] = new subscore("RTDiff", true, 0.1d);
        Subscores[12] = new subscore("PrecursorPPM", true, 1d);
        Subscores[13] = new subscore("WeightTotalMatch", false, 1d);
        Subscores[14] = new subscore("NorTotalMatch", false, 1d);
        Subscores[15] = new subscore("PrecursorIsotope", true, 1d);
        Subscores[16] = new subscore("SumCorr", false, 1d);
        Subscores[17] = new subscore("SumPPMCorr", true, 5d);
        //Subscores[18] = new subscore("DeltaPPMCorr", false, 1d);
    }
    public MatchSubscore(int Version) {
        this.Version=Version;
        switch (Version) {
            case 1: {
                V1_Array();
                break;
            }
            default:{
                V2_Array();
                break;
            }
        }
    }
    
    public void InitializeLDACoeff() {
        NoEnableSubscores = 0;
        for (int i = 0; i < Subscores.length; i++) {
            if (Subscores[i].enable) {
                NoEnableSubscores++;
            }
        }
        SubSCoeff=new double[NoEnableSubscores];
        SubSName=new String[NoEnableSubscores];
        int idx=0;
        for (int i = 0; i < Subscores.length; i++) {
            if (Subscores[i].enable) {
                SubSCoeff[idx] = Subscores[i].initialweight;
                SubSName[idx]=Subscores[i].name;
                idx++;
            }
        }
    }
  
    private float [] GetScoreArray_V1(PeakGroupScore peakgroup){
        float[] Subs = new float[7];
        int idx=0;
        Subs[idx++] = peakgroup.SpecCorrelation;
        Subs[idx++] = peakgroup.AveCorrScore;
        Subs[idx++] = peakgroup.FragIntAvgScore;
        Subs[idx++] = peakgroup.PPMScore;
        Subs[idx++] = peakgroup.ApexDeltaScore;
        Subs[idx++] = peakgroup.RTOverlapScore;
        Subs[idx++] = peakgroup.NoMatchB+peakgroup.NoMatchY;
        return Subs;
    }
    
    private float [] GetScoreArray_V2(PeakGroupScore peakgroup){
        float[] Subs = new float[18];
        int idx=0;
        Subs[idx++] = peakgroup.SpecDotProduct;
        Subs[idx++] = peakgroup.SpecCorrelation;
        Subs[idx++] = peakgroup.ContrastAngle;
        Subs[idx++] = peakgroup.AveCorrScore;
        Subs[idx++] = peakgroup.FragIntAvgScore;
        Subs[idx++] = peakgroup.PPMScore;
        Subs[idx++] = peakgroup.ApexDeltaScore;
        Subs[idx++] = peakgroup.RTOverlapScore;
        Subs[idx++] = peakgroup.NoMatchB;
        Subs[idx++] = peakgroup.NoMatchY;
        Subs[idx++] = peakgroup.PrecursorCorr;
        Subs[idx++] = peakgroup.RTDiff;
        Subs[idx++] = peakgroup.PrecursorPPM;
        Subs[idx++] = peakgroup.MaxMatchCorr*(peakgroup.NoMatchY+peakgroup.NoMatchB);
        Subs[idx++] = (peakgroup.MaxMatchCorr*(peakgroup.NoMatchY+peakgroup.NoMatchB))/peakgroup.Peplength;
        Subs[idx++] = peakgroup.PrecursorIsoPattern;
        Subs[idx++] = peakgroup.SumCorrScore;
        Subs[idx++] = peakgroup.SumCorrPPMScore;
        //Subs[idx++] = peakgroup.PrecursorCentralRank;
        return Subs;
    }
    public float [] GetSubScoreArray(PeakGroupScore peakgroup){
        switch (Version) {
            case 1: {
                return GetScoreArray_V1(peakgroup);
            }
            default: {
                return GetScoreArray_V2(peakgroup);
            }
        }
    }
    
    public double[] GetEnableSubScoreArray(PeakGroupScore peakgroup){
        float[] Subs=GetSubScoreArray(peakgroup);
        double[] enablesubsore=new double[NoEnableSubscores];
        int idx=0;
        for (int i = 0; i < Subscores.length; i++) {
            if (Subscores[i].enable) {
                enablesubsore[idx++] = Subs[i];
            }
        }    
        return enablesubsore;
    }
}
