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

import ExternalPackages.SortedListLib.SortedList;
import java.util.Comparator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SortedCorrFrag extends SortedList<PrecursorFragmentPairEdge> {
    private static final long serialVersionUID = -2283085827334564318L;
    
    public SortedCorrFrag() {
        super(new Comparator<PrecursorFragmentPairEdge>() {
            @Override
            public int compare(PrecursorFragmentPairEdge x, PrecursorFragmentPairEdge y) {
                if (x.Correlation == y.Correlation) {
                    return Float.compare(x.FragmentMz, y.FragmentMz);
                }
                return -Float.compare(x.Correlation, y.Correlation);
            }
        });
    }
}
