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

import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import java.util.Arrays;

/**
 * B-spline smoothing
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class Bspline {
	static long ff_to_long(final float x,final float y){
		final long l=((long)Float.floatToRawIntBits(x)) << 32;
		return l|Float.floatToRawIntBits(y);
	}
        static float get_Y(final long l){
		return Float.intBitsToFloat((int)(l & 0xffffffff));
	}
	static float get_X(final long l){
		return Float.intBitsToFloat((int) (l>>>32));
	}
    static public float[] Run__(final float[] data, final int PtNum, final int smoothDegree) {
        int p = smoothDegree;
        final int data_size = data.length / 2;
        int n = data_size - 1;
        int m = data_size + p;


        if (data_size <= p)
            return data.clone();


        final float[] bspline_T_ = new float[m + p];
        for (int i = 0; i <= n; i++) {
            bspline_T_[i] = 0;
            bspline_T_[m - i] = 1;
        }
        float intv = 1.0f / (m - 2 * p);
        for (int i = 1; i <= (m - 1); i++) {
            bspline_T_[p + i] = bspline_T_[p + i - 1] + intv;
        }
        final float[] ret = new float[2*(PtNum + 1)];

        for (int i = 0; i <= PtNum; i++) {
            final float t = (float) i / PtNum;
//            getbspline
            float x = 0, y = 0;//point
            for (int ii = 0; ii <= n; ii++) {
                final float a = bspline_base__(ii, p, t, bspline_T_);
                x += data[2*ii]*a;
                y += data[2*ii+1]*a;
            }
            ret[2 * i] = x;
            ret[2 * i + 1] = y;
        }
//        final boolean addHead = data[0]<ret[0];
        final boolean addTail = ret[2*(PtNum+1)-2]<data[data.length-2];
//        System.out.println("HT\t"+addHead+"\t"+addTail+"\t"+data[0]+"\t"+ret[0]+"\t"+ret[2*(PtNum+1)-2]+"\t"+data[data.length-2]);
        if(addTail){
//            final float[] tmp=new float[ret.length+2];
//            System.arraycopy(ret, 0, tmp, 0, ret.length);
            final float[] tmp = Arrays.copyOf(ret, ret.length+2);
            tmp[tmp.length-2]=data[data.length-2];
            tmp[tmp.length-1]=data[data.length-1];
            return tmp;
        }
//        if (bsplineCollection.Data.get(bsplineCollection.PointCount() - 1).getX() < data.Data.get(data.PointCount() - 1).getX()) {
//            bsplineCollection.AddPoint(data.Data.get(data.PointCount() - 1));
//        }
//        if (bsplineCollection.Data.get(0).getX() > data.Data.get(0).getX()) {
//            bsplineCollection.AddPoint(data.Data.get(0));
//        }
        return ret;
    }

    static float bspline_base__(int i, int p, float t, float[] bspline_T_) {
        float n, c1, c2;
        float tn1 = 0;
        float tn2 = 0;
        if (p == 0) {
            if (bspline_T_[i] <= t && t < bspline_T_[i + 1] && bspline_T_[i] < bspline_T_[i + 1]) {
                n = 1;
            } else {
                n = 0;
            }
        } else {
            if ((bspline_T_[i + p] - bspline_T_[i]) == 0) {
                c1 = 0;
            } else {
                tn1 = bspline_base__(i, (p - 1), t, bspline_T_);
                c1 = (t - bspline_T_[i]) / (bspline_T_[i + p] - bspline_T_[i]);
            }
            if ((bspline_T_[i + p + 1] - bspline_T_[i + 1]) == 0) {
                c2 = 0;
            } else {
                tn2 = bspline_base__((i + 1), (p - 1), t, bspline_T_);
                c2 = (bspline_T_[i + p + 1] - t) / (bspline_T_[i + p + 1] - bspline_T_[i + 1]);
            }
            n = (c1 * tn1) + (c2 * tn2);
        }
        return n;
    }

    private float[] bspline_T_ = null;

    public XYPointCollection Run(XYPointCollection data, int PtNum, int smoothDegree) {
        XYPointCollection bsplineCollection = new XYPointCollection();
        int p = smoothDegree;
        int n = data.Data.size() - 1;
        int m = data.Data.size() + p;
        bspline_T_ = new float[m + p];

        if (data.Data.size() <= p) {
            return data;
        }

        for (int i = 0; i <= n; i++) {
            bspline_T_[i] = 0;
            bspline_T_[m - i] = 1;
        }
        float intv = 1.0f / (m - 2 * p);
        for (int i = 1; i <= (m - 1); i++) {
            bspline_T_[p + i] = bspline_T_[p + i - 1] + intv;
        }


        for (int i = 0; i <= PtNum; i++) {
            final float t = ((float) i / PtNum);
            XYData pt = getbspline(data, t, n, p);
            bsplineCollection.AddPoint(pt);
        }
        if (bsplineCollection.Data.get(bsplineCollection.PointCount() - 1).getX() < data.Data.get(data.PointCount() - 1).getX()) {
            bsplineCollection.AddPoint(data.Data.get(data.PointCount() - 1));
        }
        if (bsplineCollection.Data.get(0).getX() > data.Data.get(0).getX()) {
            bsplineCollection.AddPoint(data.Data.get(0));
        }
        bsplineCollection.Data.Finalize();
        return bsplineCollection;
    }

    XYData getbspline(XYPointCollection data, float t, int n, int p) {
        float x=0,y=0;
        for (int i = 0; i <= n; i++) {
            final float a=bspline_base(i, p, t);
            final XYData pt=data.Data.get(i);
            x += pt.getX() * a;
            y += pt.getY() * a;
        }
        return new XYData(x, y);
    }

    float bspline_base(int i, int p, float t) {
        float n, c1, c2;
        float tn1 = 0;
        float tn2 = 0;
        if (p == 0) {
            if (bspline_T_[i] <= t && t < bspline_T_[i + 1] && bspline_T_[i] < bspline_T_[i + 1]) {
                n = 1;
            } else {
                n = 0;
            }
        } else {
            if ((bspline_T_[i + p] - bspline_T_[i]) == 0) {
                c1 = 0;
            } else {
                tn1 = bspline_base(i, (p - 1), t);
                c1 = (t - bspline_T_[i]) / (bspline_T_[i + p] - bspline_T_[i]);
            }
            if ((bspline_T_[i + p + 1] - bspline_T_[i + 1]) == 0) {
                c2 = 0;
            } else {
                tn2 = bspline_base((i + 1), (p - 1), t);
                c2 = (bspline_T_[i + p + 1] - t) / (bspline_T_[i + p + 1] - bspline_T_[i + 1]);
            }
            n = (c1 * tn1) + (c2 * tn2);
        }
        return n;
    }
}
