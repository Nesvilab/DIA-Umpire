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

import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class Regression {

    public Equation equation;
    protected XYPointCollection pointset;
    private float SigXY = 0;
    private float SigX = 0;
    private float SigY = 0;
    private float SigX2 = 0;
    private float SigY2 = 0;
    protected float SST;
    protected float SSR;
    private float SXX;
    private float SYY;
    private float SXY;
    private float MeanY;
    private float MeanX;
    protected float max_x = 0f;
    protected float min_x = Float.MAX_VALUE;
    protected float max_y = 0f;
    protected float min_y = Float.MAX_VALUE;
    public int MinPoint=3;


    /// <summary>
    /// Equation class: y=mx+b
    /// </summary>
    public class Equation {

        public float Bvalue;
        public float Mvalue;
        public float SDvalue;
        public float R2value;
        public int NoPoints;
        public float CorrelationCoffe;        

        public String GetEquationText() {
            return "Y=(" + Math.round(Mvalue * 1000) / 1000 + ")X+" + Math.round(Bvalue * 1000) / 1000;
        }
    }
    
    public boolean valid(){
        return pointset.PointCount()>=MinPoint;
    }

    public void SetData(XYPointCollection pointset) {
        this.pointset = pointset;
        equation = new Equation();
        FindEquation();
    }

    public float GetX(float y) {
        return (y - equation.Bvalue) / equation.Mvalue;
    }

    public float GetY(float x) {
        return equation.Mvalue * x + equation.Bvalue;
    }

    protected void FindEquation() {
        for (int i = 0; i < pointset.PointCount(); i++) {
            XYData point = pointset.Data.get(i);
            SigXY += point.getX() * point.getY();
            SigX += point.getX();
            SigY += point.getY();
            SigX2 += point.getX() * point.getX();
            SigY2 += point.getY() * point.getY();
            if (point.getX() > max_x) {
                max_x = point.getX();
            }
            if (point.getX() < min_x) {
                min_x = point.getX();
            }
            if (point.getY() > max_y) {
                max_y = point.getY();
            }
            if (point.getY() < min_y) {
                min_y = point.getY();
            }
        }
        equation.Mvalue = ((pointset.PointCount() * SigXY) - (SigX * SigY)) / ((pointset.PointCount() * SigX2) - (SigX * SigX));
        equation.Bvalue = (SigY - (equation.Mvalue * SigX)) / pointset.PointCount();
        equation.NoPoints = pointset.PointCount();
        MeanY = SigY / pointset.PointCount();
        MeanX = SigX / pointset.PointCount();
        //ComputeSD();
        //ComputeCorrelationCoff();
    }

    private void ComputeCorrelationCoff() {
        ComputeSXY();
        ComputeSXX();
        ComputeSYY();
        equation.CorrelationCoffe = (float) (SXY / Math.pow((double) SXX * SYY, 0.5));
    }

    private void ComputeSXY() {
        SXY = 0;
        for (int i = 0; i < pointset.PointCount(); i++) {
            SXY += (pointset.Data.get(i).getX() - MeanX) * (pointset.Data.get(i).getY() - MeanY);
        }
    }

    private void ComputeSXX() {
        SXX = 0;
        for (int i = 0; i < pointset.PointCount(); i++) {
            SXX += (pointset.Data.get(i).getX() - MeanX) * (pointset.Data.get(i).getX() - MeanX);
        }
    }

    private void ComputeSYY() {
        SYY = 0;
        for (int i = 0; i < pointset.PointCount(); i++) {
            SYY += (pointset.Data.get(i).getY() - MeanY) * (pointset.Data.get(i).getY() - MeanY);
        }
    }

    private void ComputeSD() {
        equation.SDvalue = (float) Math.sqrt((double) ((((pointset.PointCount() * SigY2) - (SigY * SigY)) - equation.Mvalue * ((pointset.PointCount() * SigXY) - (SigX * SigY))) / pointset.PointCount()));
    }

    private void ComputeR2() {
        ComputeSST();
        ComputeSSR();
        equation.R2value = (SST - SSR) / SST;
    }

    public float GetR2() {
        ComputeR2();
        return equation.R2value;
    }

    protected void ComputeSST() {
        SST = 0;
        for (int i = 0; i < pointset.PointCount(); i++) {
            SST += (pointset.Data.get(i).getY() - MeanY) * (pointset.Data.get(i).getY() - MeanY);
        }
    }

    private void ComputeSSR() {
        SSR = 0;
        for (int i = 0; i < pointset.PointCount(); i++) {
            SSR += (pointset.Data.get(i).getY() - (GetY(pointset.Data.get(i).getX()))) * (pointset.Data.get(i).getY() - (GetY(pointset.Data.get(i).getX())));
        }
    }
}
