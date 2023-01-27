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
import java.io.Serializable;
import java.util.Comparator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SortedCurveCollectionApexRT extends SortedList<PeakCurve> implements Serializable{

    public SortedCurveCollectionApexRT() {
        super(new Comparator<PeakCurve>() {
            @Override
            public int compare(PeakCurve x, PeakCurve y) {
                if (x.ApexRT == y.ApexRT) {
                    return Float.compare(x.TargetMz, y.TargetMz);
                }
                return Float.compare(x.ApexRT, y.ApexRT);
            }
        });
    }

    public int BinarySearchLower(float value) {
        if (isEmpty()) {
            return 0;
        }
        int lower = 0;
        int upper = size() - 1;

        if (value - get(upper).ApexRT >= 0) {
            return upper;
        }
        if (value - get(0).ApexRT <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle).ApexRT;
            if (comparisonResult == 0) {
                return middle;
            } else if (comparisonResult < 0) {
                upper = middle - 1;
            } else {
                lower = middle + 1;
            }
        }
        if (upper < 0) {
            return 0;
        }
        while (upper > 0 && get(upper).ApexRT >= value) {
            upper--;
        }
        return upper;
    }

    public int BinarySearchHigher(float value) {
        if (isEmpty()) {
            return 0;
        }
        int lower = 0;
        int upper = size() - 1;

        if (value - get(upper).ApexRT >= 0) {
            return upper;
        }
        if (value - get(0).ApexRT <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle).ApexRT;
            if (comparisonResult == 0) {
                return middle;
            } else if (comparisonResult < 0) {
                upper = middle - 1;
            } else {
                lower = middle + 1;
            }
        }
        if (lower > size() - 1) {
            return size() - 1;
        }
        while (upper < size() && get(upper).ApexRT <= value) {
            upper++;
        }
        return upper;
    }

    public int BinarySearchClosest(float value) {
        if (isEmpty()) {
            return 0;
        }
        int lower = 0;
        int upper = size() - 1;

        if (value - get(upper).ApexRT >= 0) {
            return upper;
        }
        if (value - get(0).ApexRT <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle).ApexRT;
            if (comparisonResult == 0) {
                return middle;
            } else if (comparisonResult < 0) {
                upper = middle - 1;
            } else {
                lower = middle + 1;
            }
        }

        if (Math.abs(value - get(lower).ApexRT) > Math.abs(value - get(upper).ApexRT)) {
            return upper;
        } else {
            return lower;
        }
    }
}
