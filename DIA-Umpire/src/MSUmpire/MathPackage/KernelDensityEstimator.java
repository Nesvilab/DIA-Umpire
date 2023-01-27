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

import java.util.Arrays;
import umontreal.ssj.gof.KernelDensity;
import umontreal.ssj.probdist.EmpiricalDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.randvar.KernelDensityGen;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class KernelDensityEstimator {

    private EmpiricalDist empiricalDist;
    private KernelDensityGen kernelDensityGen;
    private int ObsDataSize=0;
    
    public void SetData(double [] data){
        Arrays.sort(data);
        empiricalDist = new EmpiricalDist(data);
        ObsDataSize=data.length;
    }
    
    public double[] Density(double[] xdata) {
        
        NormalDist kern = new NormalDist();        
        //Silverman's ‘rule of thumb’ (Scott Variation uses factor = 1.06)
        double bandWidth = 0.99 * Math.min(empiricalDist.getSampleStandardDeviation(), (empiricalDist.getInterQuartileRange() / 1.34)) / Math.pow(ObsDataSize, 0.2);        
        double[] DensityValues = KernelDensity.computeDensity(empiricalDist, kern, bandWidth, xdata);        
        return DensityValues;
    }    
    
}
