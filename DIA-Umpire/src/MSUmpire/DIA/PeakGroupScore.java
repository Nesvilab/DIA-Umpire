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

import MSUmpire.PeakDataStructure.PeakCluster;
import java.io.Serializable;

/**
 * Class to calculate U-score
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakGroupScore implements Serializable{
    private static final long serialVersionUID = 164894161684L;

    public transient PeakCluster cluster;
    public int MSlevel = 0;
    public float PrecursorRT;

    public float SpecDotProduct=0f;
    public float SpecCorrelation=0f;
    public float ContrastAngle=0f;
    public float AveCorrScore=0f;
    public float SumCorrScore=0f;
    public float SumCorrPPMScore=0f;
    public float PrecursorCentralRank=0f;
    public float MaxMatchCorr=0f;
    public float FragIntAvgScore=0f;
    public float PPMScore=0f;
    public float ApexDeltaScore=0f;
    public float RTOverlapScore=0f;
    public int NoMatchB = 0;
    public int NoMatchY = 0;
    public float PrecursorCorr=0f;
    public float RTDiff=0f;
    public float PrecursorPPM=0f;
    public int Peplength=0;
    public float PrecursorScore=0f;    
    public float PrecursorIsoPattern=0f;    
    public int NoFragmentLib = 0;
    public float MixtureModelProb;
    public float MixtureModelLocalProb;
    public float UmpireScore;
    public float PrecursorNeutralMass;
    
    public PeakGroupScore(PeakCluster cluster) {
        this.cluster = cluster;
        this.PrecursorRT=cluster.PeakHeightRT[0];
        this.PrecursorNeutralMass=cluster.NeutralMass();
    }

    void AddScore(float UmpireScore) {
        cluster.AddScore(UmpireScore);
    }

    float GetScoreRank(float UmpireScore) {
        return cluster.GetScoreRank(UmpireScore);
    }
    
}
