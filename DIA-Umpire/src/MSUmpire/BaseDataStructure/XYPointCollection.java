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

import java.io.Serializable;

/*
 * Collection of two dimensional
 */
/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class XYPointCollection implements Serializable {
    private static final long serialVersionUID = 654648165471L;

    public float MaxY;
    public SortedXYCollectionClass Data;
    //XYComparator comparator=new XYComparator();

    public XYPointCollection() {
        Data = new SortedXYCollectionClass();
    }

    public float GetSumX(){
        float sum=0f;
        for(XYData point : Data){
            sum+=point.getX();
        }
        return sum;
    }
    
    public float GetSumY(){
        float sum=0f;
        for(XYData point : Data){
            sum+=point.getY();
        }
        return sum;
    }
    
        
    public void AddPoint(float x, float y) {
        this.Data.add(new XYData(x, y));
        if (MaxY < y) {
            MaxY = y;
        }
    }
    
    public void AddPointKeepMaxIfCloseValueExisted(float x, float y, float ppm) {
        boolean insert = true;
        if (this.Data.size() > 0) {
            int idx = GetClosetIndexOfX(x);
            XYData pt = Data.get(idx);
            if (InstrumentParameter.CalcPPM(pt.getX(),x)<ppm) {
                insert = false;
                if (y < pt.getY()) {
                    pt.setY(y);
                    pt.setX(x);
                }
                if (MaxY < y) {
                    MaxY = y;
                }
            }
        }
        if (insert) {
            AddPoint(x, y);
        }
    }

    public void AddPointKeepMaxIfValueExisted(float x, float y) {
        boolean insert = true;
        if (this.Data.size() > 0) {
            int idx = GetClosetIndexOfX(x);
            XYData pt = Data.get(idx);
            if (pt.getX() == x) {
                insert = false;
                if (y < pt.getY()) {
                    pt.setY(y);
                }
                if (MaxY < y) {
                    MaxY = y;
                }
            }
        }
        if (insert) {
            AddPoint(x, y);
        }
    }

    public void AddPoint(XYData point) {
        this.Data.add(point);
        if (MaxY < point.getY()) {
            MaxY = point.getY();
        }
    }

    public int PointCount() {
        return Data.size();
    }

    public void CentroidingbyLocalMaximum(int Resolution, float MinMZ) {
        if (Data.size() == 0) {
            return;
        }
        int oldcount = Data.size();
        SortedXYCollectionClass DataTemp = Data;
        int startindex=DataTemp.BinarySearchHigher(MinMZ);
        XYData pt = DataTemp.get(startindex);
        float maxintensity = pt.getY();
        float maxmz = pt.getX();
        float gap = pt.getX() / Resolution;
        Data = new SortedXYCollectionClass();
        for (int i = startindex+1; i < oldcount; i++) {
            XYData pti = DataTemp.get(i);
            if (pti.getX() - maxmz < gap) {
                if (pti.getY() > maxintensity) {
                    maxintensity = pti.getY();
                    maxmz = pti.getX();
                    gap = pti.getX() / Resolution;
                }
            } else {
                AddPoint(maxmz, maxintensity);
                maxintensity = pti.getY();
                maxmz = pti.getX();
                gap = pti.getX() / Resolution;
            }
        }
        DataTemp.clear();
        DataTemp = null;
        Data.Finalize();
    }

    public int GetLowerIndexOfX(float x) {
        return Data.BinarySearchLower(x);
    }

    public int GetHigherIndexOfX(float x) {
        return Data.BinarySearchHigher(x);
    }

    public int GetClosetIndexOfX(float x) {
        return Data.BinarySearchClosest(x);
    }

    public XYData GetPoinByXLower(float x) {
        return Data.get(GetLowerIndexOfX(x));
    }

    public XYData GetPoinByXCloset(float x) {
        return Data.get(GetClosetIndexOfX(x));
    }

    public XYData GetPoinByXHigher(float x) {
        return Data.get(GetHigherIndexOfX(x));

    }

    public XYPointCollection GetSubSetByXRange(float xlower, float xupper) {
        if (PointCount() == 0) {
            return null;
        }
        XYPointCollection NewXYCollection = new XYPointCollection();
        int start = GetLowerIndexOfX(xlower);

        if (start < 0) {
            start = 0;
        }

        for (int i = start; i < PointCount(); i++) {
            float x = Data.get(i).getX();
            if (x >= xlower && x <= xupper) {
                NewXYCollection.AddPoint(x, Data.get(i).getY());
            } else if (x > xupper) {
                break;
            }
        }
        return NewXYCollection;
    }
    
}
