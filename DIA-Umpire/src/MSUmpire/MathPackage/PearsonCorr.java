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
package MSUmpire.MathPackage;

import MSUmpire.BaseDataStructure.XYPointCollection;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PearsonCorr {

    public float CalcCorrNeighborBin(XYPointCollection CollectionA, XYPointCollection CollectionB) {
        Regression regression = new Regression();

        int num = (int) ((Math.min(CollectionA.Data.get(CollectionA.PointCount() - 1).getX(), CollectionB.Data.get(CollectionB.PointCount() - 1).getX()) - Math.max(CollectionA.Data.get(0).getX(), CollectionB.Data.get(0).getX())) * 100);
        float timeinterval = 1 / 100f;

        float[] arrayA = new float[num];
        float[] arrayB = new float[num];

        float start = Math.max(CollectionA.Data.get(0).getX(), CollectionB.Data.get(0).getX());

        for (int i = 0; i < num - 1; i++) {
            float low = start + i * timeinterval;
            float up = start + (i + 1) * timeinterval;

            for (int j = 0; j < CollectionA.PointCount(); j++) {
                if (CollectionA.Data.get(j).getX() >= low && CollectionA.Data.get(j).getX() < up) {
                    float intenlow = CollectionA.Data.get(j).getY() * (1 - (CollectionA.Data.get(j).getX() - low) / timeinterval);
                    float intenup = CollectionA.Data.get(j).getY() * (1 - (up - CollectionA.Data.get(j).getX()) / timeinterval);
                    if (intenlow > arrayA[i]) {
                        arrayA[i] = intenlow;
                    }
                    if (intenup > arrayA[i + 1]) {
                        arrayA[i + 1] = intenup;
                    }
                } else if (CollectionA.Data.get(j).getX() > up) {
                    break;
                }
            }

            for (int j = 0; j < CollectionB.PointCount(); j++) {
                if (CollectionB.Data.get(j).getX() >= low && CollectionB.Data.get(j).getX() < up) {
                    float intenlow = CollectionB.Data.get(j).getY() * (1 - (CollectionB.Data.get(j).getX() - low) / timeinterval);
                    float intenup = CollectionB.Data.get(j).getY() * (1 - (up - CollectionB.Data.get(j).getX()) / timeinterval);
                    if (intenlow > arrayB[i]) {
                        arrayB[i] = intenlow;
                    }
                    if (intenup > arrayB[i + 1]) {
                        arrayB[i + 1] = intenup;
                    }
                } else if (CollectionB.Data.get(j).getX() > up) {
                    break;
                }
            }
        }

        XYPointCollection pointset = new XYPointCollection();
        for (int i = 0; i < num; i++) {
            if (arrayA[i] > 0 && arrayB[i] > 0) {
                pointset.AddPoint(arrayA[i], arrayB[i]);
            }
        }

        float R2 = 0f;

        if (pointset.PointCount() > 5) {
            regression.SetData(pointset);
            if (regression.equation.Mvalue > 0) {
                R2 = regression.GetR2();
            }
        }
        return R2;
    }
    
    public double CalcCorrV2(XYPointCollection CollectionA, XYPointCollection CollectionB, int NoPointPerInterval) {
        SpearmansCorrelation pearsonsCorrelation = new SpearmansCorrelation();
        
        int num = Math.max(CollectionA.PointCount(), CollectionB.PointCount()) / 2;
        float timeinterval = 2f / (float) NoPointPerInterval;
        if (num < 6) {
            return 0f;
        }

        double[] arrayA = new double[num];
        double[] arrayB = new double[num];

        float start = Math.max(CollectionA.Data.get(0).getX(), CollectionB.Data.get(0).getX());

        int i = 0;
        float low = start;
        float up = start + timeinterval;

        for (int j = 0; j < CollectionA.PointCount(); j++) {
            while (CollectionA.Data.get(j).getX() > up) {
                i++;
                low = up;
                up = low + timeinterval;
            }
            if (i >= num) {
                break;
            }
            if (CollectionA.Data.get(j).getX() >= low && CollectionA.Data.get(j).getX() < up) {
                if (CollectionA.Data.get(j).getY() > arrayA[i]) {
                    arrayA[i] = CollectionA.Data.get(j).getY();
                }
            }
        }
        i = 0;
        low = start;
        up = start + timeinterval;
        for (int j = 0; j < CollectionB.PointCount(); j++) {
            while (CollectionB.Data.get(j).getX() > up) {
                i++;
                low = up;
                up = low + timeinterval;
            }
            if (i >= num) {
                break;
            }
            if (CollectionB.Data.get(j).getX() >= low && CollectionB.Data.get(j).getX() < up) {
                if (CollectionB.Data.get(j).getY() > arrayB[i]) {
                    arrayB[i] = CollectionB.Data.get(j).getY();
                }
            }
        }

        if(arrayA[0]==0f){
            arrayA[0]=arrayA[1];
        }
        if(arrayB[0]==0f){
            arrayB[0]=arrayB[1];
        }
        for (int idx = 1; idx < num - 1; idx++) {
            if (arrayA[idx] == 0f) {
                arrayA[idx] = (arrayA[idx - 1] + arrayA[idx + 1]) / 2;
            }
            if (arrayB[idx] == 0f) {
                arrayB[idx] = (arrayB[idx - 1] + arrayB[idx + 1]) / 2;
            }
        }
        
        if(arrayA[num - 1]==0f){
            arrayA[num - 1]=arrayA[num - 2];
        }
        if(arrayB[num - 1]==0f){
            arrayB[num - 1]=arrayB[num - 2];
        }
        double R2 =pearsonsCorrelation.correlation(arrayA, arrayB); 
        return R2;
    }
    
    public float CalcCorr(XYPointCollection CollectionA, XYPointCollection CollectionB, int NoPointPerInterval) {
        Regression regression = new Regression();

        int num = Math.max(CollectionA.PointCount(), CollectionB.PointCount()) / 2;
        float timeinterval = 2f / (float) NoPointPerInterval;
        if (num < 6) {
            return 0f;
        }

        float[] arrayA = new float[num];
        float[] arrayB = new float[num];

        float start = Math.max(CollectionA.Data.get(0).getX(), CollectionB.Data.get(0).getX());

        int i = 0;
        float low = start;
        float up = start + timeinterval;

        for (int j = 0; j < CollectionA.PointCount(); j++) {
            while (CollectionA.Data.get(j).getX() > up) {
                i++;
                low = up;
                up = low + timeinterval;
            }
            if (i >= num) {
                break;
            }
            if (CollectionA.Data.get(j).getX() >= low && CollectionA.Data.get(j).getX() < up) {
                if (CollectionA.Data.get(j).getY() > arrayA[i]) {
                    arrayA[i] = CollectionA.Data.get(j).getY();
                }
            }
        }
        i = 0;
        low = start;
        up = start + timeinterval;
        for (int j = 0; j < CollectionB.PointCount(); j++) {
            while (CollectionB.Data.get(j).getX() > up) {
                i++;
                low = up;
                up = low + timeinterval;
            }
            if (i >= num) {
                break;
            }
            if (CollectionB.Data.get(j).getX() >= low && CollectionB.Data.get(j).getX() < up) {
                if (CollectionB.Data.get(j).getY() > arrayB[i]) {
                    arrayB[i] = CollectionB.Data.get(j).getY();
                }
            }
        }

        for (int idx = 1; idx < num - 1; idx++) {
            if (arrayA[idx] == 0f) {
                arrayA[idx] = (arrayA[idx - 1] + arrayA[idx + 1]) / 2;
            }
            if (arrayB[idx] == 0f) {
                arrayB[idx] = (arrayB[idx - 1] + arrayB[idx + 1]) / 2;
            }
        }

        XYPointCollection pointset = new XYPointCollection();
        for (int idx = 0; idx < num; idx++) {
            if (arrayA[idx] > 0 && arrayB[idx] > 0) {
                pointset.AddPoint(arrayA[idx], arrayB[idx]);
            }
        }

        float R2 = 0f;
        if (pointset.PointCount() > 5) {
            regression.SetData(pointset);
            if (regression.equation.Mvalue > 0) {
                R2 = regression.GetR2();
            }
        }
        return R2;
    }
}
