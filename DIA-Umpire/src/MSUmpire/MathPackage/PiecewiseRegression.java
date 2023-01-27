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

import MSUmpire.BaseDataStructure.SortedXYZCollection;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.BaseDataStructure.XYZData;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PiecewiseRegression extends Regression {

    private float ErrorEstimateRange;
    private float SDEstimateRange;

    public PiecewiseRegression(float ErrorEstimateRange, float SDEstimateRange) {
        this.ErrorEstimateRange = ErrorEstimateRange;
        this.SDEstimateRange=SDEstimateRange;               
    }

    @Override
    public float GetR2() {
        ComputeSST();
        ComputeSSR();
        equation.R2value = (SST - SSR) / SST;
        return equation.R2value;
    }

    @Override
    public void SetData(XYPointCollection pointset) {
        this.pointset = pointset;
        equation = new Equation();
        FindEquation();
        GeneratePredictYList();
        GeneratePredictXList();
        CalculateSD();
    }

    private void ComputeSSR() {
        SSR = 0;
        for (int i = 0; i < pointset.PointCount(); i++) {
            SSR += (pointset.Data.get(i).getY() - (GetPredictTimeSDYByTimelist(pointset.Data.get(i).getX()).getY())) * (pointset.Data.get(i).getY() - (GetPredictTimeSDYByTimelist(pointset.Data.get(i).getX()).getY()));
        }
    }
    public SortedXYZCollection PredictYList;
    public SortedXYZCollection PredictXList;

    public XYZData GetPredictTimeSDYByTimelist(float xvalue) {
        return PredictYList.GetCloset(xvalue);
    }

    public XYZData GetPredictTimeSDXByTimelist(float yvalue) {
        return PredictXList.GetCloset(yvalue);
    }

    private float PredictXByLocalRegion(float Yvalue) {
        float error = 0;
        int count = 0;
        
        for (int i = 0; i < pointset.PointCount(); i++) {            
            if (pointset.Data.get(i).getY() > Yvalue - ErrorEstimateRange && pointset.Data.get(i).getY() < Yvalue + ErrorEstimateRange) {
                error += pointset.Data.get(i).getX() - (float) (GetX(pointset.Data.get(i).getY()));                
                count++;
            }
        }
        if (count < 10) {
            error = 0;
            count = 0;
            for (int i = 0; i < pointset.PointCount(); i++) {
                if (pointset.Data.get(i).getY() > Yvalue - 2 * ErrorEstimateRange && pointset.Data.get(i).getY() < Yvalue + 2 * ErrorEstimateRange) {
                    error += pointset.Data.get(i).getX() - GetX(pointset.Data.get(i).getY());                    
                    count++;
                }
            }
        }
        if (count > 10) {
            error /= count;
        }
        float predictvalue=GetX(Yvalue) + error;
        
        return predictvalue;
    }
    
    private float PredictYByLocalRegion(float Xvalue) {
        float error = 0f;        
        int count = 0;
        
        for (int i = 0; i < pointset.PointCount(); i++) {
            if (pointset.Data.get(i).getX() > Xvalue - ErrorEstimateRange && pointset.Data.get(i).getX() < Xvalue + ErrorEstimateRange) {
                error += pointset.Data.get(i).getY() - (float) (GetY(pointset.Data.get(i).getX()));
                count++;
            }
        }
        if (count < 10) {
             error = 0f;        
             count = 0;
            for (int i = 0; i < pointset.PointCount(); i++) {
                if (pointset.Data.get(i).getX() > Xvalue - 2 * ErrorEstimateRange && pointset.Data.get(i).getX() < Xvalue + 2 * ErrorEstimateRange) {
                    error += pointset.Data.get(i).getY() - GetY(pointset.Data.get(i).getX());
                    count++;
                }
            }
        }
        if (count > 10) {
            error /= count;
        }

        float predictvalue=GetY(Xvalue) + error;        
    
        return predictvalue;
    }
    
    public void CalculateSD()
    {
        int[] countX=new int[PredictXList.size()];
        int[] countY=new int[PredictYList.size()];
        float[] SDX=new float[PredictXList.size()];
        float[] SDY=new float[PredictYList.size()];
        float[] PredictX=new float[pointset.PointCount()];
        float[] PredictY=new float[pointset.PointCount()];
        for (int i = 0; i < pointset.PointCount(); i++) {
            XYData data = pointset.Data.get(i);
            PredictX[i] = PredictXList.GetCloset(data.getY()).getX();
            PredictY[i] = PredictYList.GetCloset(data.getX()).getX();
        }
        for (int i = 0; i < pointset.PointCount(); i++) {
            for(int j=0;j<PredictXList.size();j++){
                if(pointset.Data.get(i).getX() > PredictXList.get(j).getX() - SDEstimateRange && pointset.Data.get(i).getX() < PredictXList.get(j).getX() + SDEstimateRange){
                    countX[j]++;
                    float error=PredictX[i]-pointset.Data.get(i).getX();
                    SDX[j]+=error*error;
                }
            }
            for(int j=0;j<PredictYList.size();j++){
                if(pointset.Data.get(i).getY() > PredictYList.get(j).getX() - SDEstimateRange && pointset.Data.get(i).getY() < PredictYList.get(j).getX() + SDEstimateRange){
                    countY[j]++;
                    float error=PredictY[i]-pointset.Data.get(i).getY();
                    SDY[j]+=error*error;
                }
            }
        }
        for (int j = 0; j < PredictXList.size(); j++) {
            SDX[j] /= countX[j];
            SDX[j] = (float) Math.sqrt(SDX[j]);            
            PredictXList.get(j).setZ(SDX[j]);
        }
        for (int j = 0; j < PredictYList.size(); j++) {
            SDY[j] /= countY[j];
            SDY[j] = (float) Math.sqrt(SDY[j]);            
            PredictYList.get(j).setZ(SDY[j]);
        }
    }

    public void GeneratePredictYList() {
        PredictYList = new SortedXYZCollection();
        float gap = 0.2f;
        for (int i = 0; i < (int) (max_x + 2 - min_x) * 5; i++) {
            float x = min_x - 1 + gap * i;
            PredictYList.add(new XYZData(x, PredictYByLocalRegion(x), 0f));
        }
    }

    public void GeneratePredictXList() {
        PredictXList = new SortedXYZCollection();
        float gap = 0.2f;
        for (int i = 0; i < (int) (max_y + 2 - min_y) * 5; i++) {
            float y = min_y - 1 + gap * i;
            PredictXList.add(new XYZData(y,PredictXByLocalRegion(y),0f));
        }
    }
}
