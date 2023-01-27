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
package MSUmpire.PSMDataStructure;

import ExternalPackages.SortedListLib.SortedList;
import java.util.Comparator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SortedPepListProb extends SortedList<PepIonID> {

    public SortedPepListProb() {
        super(new Comparator<PepIonID>() {
            @Override
            public int compare(PepIonID x, PepIonID y) {
                if (x.MaxProbability == y.MaxProbability) {
                    return -1;
                }
                return Float.compare(y.MaxProbability, x.MaxProbability);
            }
        });
    }

}
