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
package MSUmpire.BaseDataStructure;

import ExternalPackages.SortedListLib.SortedList;
import java.io.Serializable;
import java.util.Comparator;

/**
 * Sorted collection of XYZData, sorted by x value
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SortedXYZCollection extends SortedList<XYZData> implements Serializable {
    private static final long serialVersionUID = 836465863L;

    public SortedXYZCollection() {
        super(new Comparator<XYZData>() {
            @Override
            public int compare(XYZData x, XYZData y) {
                if (x.getX() == y.getX()) {
                    return 1;
                }
                return -Float.compare(y.getX(), x.getX());
            }
        });
    }

    public XYZData GetCloset(float value) {
        return get(BinarySearchClosest(value));
    }
    public int BinarySearchClosest(float value) {

        if (isEmpty()) {
            return 0;
        }
        int lower = 0;
        int upper = size() - 1;

        if (value - get(upper).getX() >= 0) {
            return upper;
        }
        if (value - get(0).getX() <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle).getX();
            if (comparisonResult == 0) {
                return middle;
            } else if (comparisonResult < 0) {
                upper = middle - 1;
            } else {
                lower = middle + 1;
            }
        }

        if (Math.abs(value - get(lower).getX()) > Math.abs(value - get(upper).getX())) {
            return upper;
        } else {
            return lower;
        }
    }

}
