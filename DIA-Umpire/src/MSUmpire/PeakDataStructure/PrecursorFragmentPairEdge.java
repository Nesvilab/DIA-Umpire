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
package MSUmpire.PeakDataStructure;

import java.io.Serializable;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PrecursorFragmentPairEdge extends PeakOverlapRegion implements Serializable{
    private static final long serialVersionUID = 2890272386178844781L;
    public float FragmentMz;
    public float Intensity;
    public float ApexDelta;
    public float RTOverlapP;    
    public int FragmentMS1Rank = 0;
    public float FragmentMS1RankScore = 1f;
    public float AdjustedFragInt;
    public boolean ComplementaryFragment = false;
    public float MatchedFragMz;
}
