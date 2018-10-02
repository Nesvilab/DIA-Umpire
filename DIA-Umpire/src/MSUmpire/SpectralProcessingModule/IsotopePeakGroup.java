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
package MSUmpire.SpectralProcessingModule;

import MSUmpire.BaseDataStructure.XYData;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import java.util.ArrayList;

/**
 * Isotope peak group class
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class IsotopePeakGroup {

    public int Charge;
    public ArrayList<XYData> PeakGroupList = new ArrayList<>();

    public IsotopePeakGroup(int charge) {
        this.Charge = charge;
    }

    public float NeutralMass(){
        return Charge * (float)((PrecursorMz() - ElementaryIon.proton.getTheoreticMass()));
    }
    
    public float PrecursorMz() {
        return PeakGroupList.get(0).getX();
    }

    public XYData GetPeakXYPointByPeakidx(int pkidx) {
        return PeakGroupList.get(pkidx);
    }

    public void AddPeak(XYData peakPoint) {
        PeakGroupList.add(peakPoint);
    }
    public float[] GetIsoPattern() {
        float[] Pattern = new float[3];
        for (int i = 0; i < 3; i++) {
            if (i == 0) {
                Pattern[i] = 1;
            }
            Pattern[i] = PeakGroupList.get(i).getY() / PeakGroupList.get(0).getY();
        }
        return Pattern;
    }
        
}
