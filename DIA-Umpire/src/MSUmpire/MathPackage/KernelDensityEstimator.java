/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
