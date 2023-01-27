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
/*
 * 
 */

import ExternalPackages.SortedListLib.SortedList;
import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;

/** 
 * Sorted collection of XYData, sorted by x value
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SortedXYCollectionClass extends SortedList<XYData> implements Serializable{
    private static final long serialVersionUID = 65464643184541L;

//    private float[][] SortedArray;
    private float[] SortedArray;

    public int size;

    public SortedXYCollectionClass() {
        super(new Comparator<XYData>() {
            @Override
            public int compare(XYData x, XYData y) {
                return -x.compareTo(y);
            }
        });
    }

    @Override
    public boolean isEmpty() {
        if (!Finalized) {
            return getRoot() == null;
        }
        return size == 0;
    }

    /**
     * Removes all elements from the list, leaving it empty.
     */
    @Override
    public void clear() {
        FinalizedSortedArray = null;
        SortedArray = null;
    }

    public Object[] toArray(){
        if (Finalized) {
            XYData[] arr = new XYData[size];
            for (int i = 0; i < size; i++) {
                arr[i] = get(i);
            }
            return arr;
        }        
        return super.toArray();
    }
    
    public Iterator<XYData> iterator() {
        if (Finalized) {                        
            return (Iterator<XYData>) Arrays.asList((XYData[])toArray()).iterator();
        }
        return super.iterator();
    }
    private synchronized void writeObject(java.io.ObjectOutputStream stream) throws java.io.IOException {
        if (!Finalized) {
            Finalize();
        }
        stream.defaultWriteObject();
        stream.writeInt(size);
        for (int i = 0; i < size; i++) {
//            stream.writeFloat(SortedArray[0][i]);
//            stream.writeFloat(SortedArray[1][i]);
            stream.writeFloat(SortedArray[0+i*2]);
            stream.writeFloat(SortedArray[1+i*2]);
        }
    }
    private void readObject(java.io.ObjectInputStream in) throws ClassNotFoundException, IOException {
        in.defaultReadObject();
        Finalized = true;
        size=in.readInt();
//        SortedArray = new float[2][size];
        SortedArray = new float[2*size];
        for (int i = 0; i < size; i++) {            
            //XYData xy = (XYData) in.readObject();
//            SortedArray[0][i] = xy.getX();
//            SortedArray[1][i] = xy.getY();
            SortedArray[0+i*2] = in.readFloat();
            SortedArray[1+i*2] = in.readFloat();
        }
    }
    
    @Override
    public void Finalize() {
        FinalizedSortedArray = toArray();
        ClearTree();
        Finalized = true;
        size = FinalizedSortedArray.length;
//        SortedArray = new float[2][size];
        SortedArray = new float[2*size];
        for (int i = 0; i < size; i++) {
            XYData xy = (XYData) FinalizedSortedArray[i];
            SortedArray[0+i*2] = xy.getX();
            SortedArray[1+i*2] = xy.getY();
        }
        FinalizedSortedArray = null;
    }

    @Override
    public XYData get(int index) {
        if (Finalized) {
//            return new XYData(SortedArray[0][index], SortedArray[1][index]);
            return new XYData(SortedArray[0+2*index], SortedArray[1+2*index]);
        }
        return (XYData) findNodeAtIndex(index).getValue();
    }

    @Override
    public int size() {
        if (!Finalized) {
            return (getRoot() == null) ? 0 : 1 + getRoot().GetNumOfChildren();
        }
        return size;
    }

    public int BinarySearchLower(XYData value) {
        return BinarySearchLower(value.getX());
    }

    public int BinarySearchHigher(float value) {
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
                while (middle - 1 >= 0 && get(middle - 1).getX() == value) {
                    middle--;
                }
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
        while (lower < size() - 1 && get(lower).getX() <= value) {
            lower++;
        }
        return lower;
    }

    public int BinarySearchLower(float value) {
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
                while (middle - 1 >= 0 && get(middle - 1).getX() == value) {
                    middle--;
                }
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
        while (upper > 0 && get(upper).getX() >= value) {
            upper--;
        }
        return upper;
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
